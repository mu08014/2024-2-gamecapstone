using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Rendering.VirtualTexturing;
using UnityEngine.UIElements;
using static UnityEngine.ParticleSystem;

public class StrandCompliantEuler
{
    private TripletXs m_A_nz;
    private TripletXs m_J_nz;
    private TripletXs m_Jv_nz;
    private TripletXs m_Jxv_nz;
    private TripletXs m_invC_nz;
    private TripletXs m_invCv_nz;

    private TripletXs m_J_inter_nz;
    private TripletXs m_Jv_inter_nz;
    private TripletXs m_invC_inter_nz;
    private TripletXs m_invCv_inter_nz;

    private SparseXs m_A;
    private SparseXs m_J;
    private SparseXs m_Jv;
    private SparseXs m_Jxv;
    private SparseXs m_invC;
    private SparseXs m_invCv;

    private SparseXs m_J_inter;
    private SparseXs m_Jv_inter;
    private SparseXs m_invC_inter;
    private SparseXs m_invCv_inter;

    private Vectors m_A_inv_diag;

    SimplicialLDLT m_solver = new SimplicialLDLT();

    private int m_hidx;

    private StrandCompliantManager m_parent;

    private int m_start_global_dof;
    private int m_num_global_dof;

    private int DIM;

    public StrandCompliantEuler(
        StrandCompliantManager parent, int hidx, int dim = 3)
    {
        m_parent = parent;
        m_hidx = hidx;
        DIM = 3;

        HairFlow flow = m_parent.m_scene.getFilmFlows()[hidx];
        int ndof = m_parent.m_scene.isMassSpring() ? flow.size() * DIM
                                               : flow.size() * 4 - 1;
        m_A = new SparseXs(ndof, ndof);

        m_num_global_dof = ndof;
        m_start_global_dof = m_parent.m_scene.getDof(flow.getParticleIndices()[0]);
    }

    public void computeRHS(ref TwoDScene scene, double dt, ref Vectors b)
    {
        Vectors v = scene.getV().segment(m_start_global_dof, m_num_global_dof);
        HairFlow flow = m_parent.m_scene.getFilmFlows()[m_hidx];
        Vectors m = scene.getInterpolatedM().segment(m_start_global_dof, m_num_global_dof);
        Vectors gradu = m_parent.m_gradU.segment(m_start_global_dof, m_num_global_dof);
        VectorXi constraint_idx = flow.getConstraintIdx();
        VectorXi constraint_num = flow.getNumConstraints();

        Vectors phi = m_parent.m_Phi.segment(constraint_idx[0], constraint_num[0]);

        b.SetSegment(m_start_global_dof, m_num_global_dof,
            v.cwiseProduct(m) - gradu * dt - dt * (m_J.transpose() * (m_invC * phi)));

        if (constraint_num[1] > 0)
        {
            Vectors phi_v =
                m_parent.m_Phi_v.segment(constraint_idx[1], constraint_num[1]);
            b.PlusSegment(m_start_global_dof, m_num_global_dof,
                -dt * (m_Jv.transpose() * (m_invCv * (phi_v - m_Jv * v))));
        }
    }

    public bool stepScene(ref TwoDScene scene, double dt, ref Vectors r, in Vectors b)
    {
        Vectors rr = b.segment(m_start_global_dof, m_num_global_dof);
        r.SetSegment(m_start_global_dof, m_num_global_dof,
            m_solver.solve(rr));  // m_A_diag.cwiseProduct(rr);
        return true;
    }

    public void updateLambda(ref TwoDScene scene, in Vectors dx, in Vectors dv, double dt)
    {
        HairFlow flow = scene.getFilmFlows()[m_hidx];
        List<int> particles = flow.getParticleIndices();

        VectorXi constraint_idx = flow.getConstraintIdx();
        VectorXi constraint_num = flow.getNumConstraints();

        Vectors phi =
            m_parent.m_Phi.segment(constraint_idx[0], constraint_num[0]);

        Vectors ddx = dx.segment(m_start_global_dof, m_num_global_dof);
        Vectors ddv = dv.segment(m_start_global_dof, m_num_global_dof);

        m_parent.m_lambda.SetSegment(constraint_idx[0], constraint_num[0],
            (-1) * m_invC * (m_J * ddx + phi));
        if (constraint_num[1] > 0)
        {
            Vectors phi_v =
                m_parent.m_Phi_v.segment(constraint_idx[1], constraint_num[1]);
            m_parent.m_lambda_v.SetSegment(constraint_idx[1], constraint_num[1],
                (-1) * m_invCv * (m_Jxv * ddx + m_Jv * ddv + phi_v));
        }
    }

    public void computeRHSIncremental(ref TwoDScene scene, double dt, ref Vectors b, in Vectors vplus)
    {
        Vectors dv = (vplus.segment(m_start_global_dof, m_num_global_dof) -
                scene.getV().segment(m_start_global_dof, m_num_global_dof));
        HairFlow flow = m_parent.m_scene.getFilmFlows()[m_hidx];
        Vectors m =
            scene.getInterpolatedM().segment(m_start_global_dof, m_num_global_dof);
        Vectors gradu =
            m_parent.m_gradU.segment(m_start_global_dof, m_num_global_dof);
        VectorXi constraint_idx = flow.getConstraintIdx();
        VectorXi constraint_num = flow.getNumConstraints();

        Vectors phi =
            m_parent.m_Phi.segment(constraint_idx[0], constraint_num[0]);

        b.SetSegment(m_start_global_dof, m_num_global_dof,
            (-1) * dv.cwiseProduct(m) - gradu * dt -
            dt * (m_J.transpose() * (m_invC * phi)));
        if (constraint_num[1] > 0)
        {
            Vectors phi_v =
                m_parent.m_Phi_v.segment(constraint_idx[1], constraint_num[1]);
            b.PlusSegment(m_start_global_dof, m_num_global_dof,
                -dt * (m_Jv.transpose() * (m_invCv * phi_v)));
        }
    }

    public bool stepScene(ref TwoDScene scene, double dt)
    {
        return false;
    }

    public void updateNextV(ref TwoDScene scene, in Vectors vplus)
    {
        Vectors v_next = scene.getV();
        HairFlow flow = scene.getFilmFlows()[m_hidx];
        List<int> particles = flow.getParticleIndices();

        foreach (int pidx in particles)
        {
            if (!scene.isFixed(pidx))
            {
                int numdofs = scene.isMassSpring() || scene.isTip(pidx) ? DIM : 4;
                v_next.SetSegment(scene.getDof(pidx), numdofs,
                    vplus.segment(scene.getDof(pidx), numdofs));
            }
        }
    }

    public void computeAp(in Vectors p, ref Vectors b)
    {
        b.SetSegment(m_start_global_dof, m_num_global_dof, 
            m_A * p.segment(m_start_global_dof, m_num_global_dof));
    }

    public bool PreconditionScene(ref TwoDScene scene, double dt, ref Vectors r, in Vectors b)
    {
        if (m_parent.m_use_preconditioner)
        {
            return stepScene(ref scene, dt, ref r, b);
        }
        else
        {
            Vectors rr = b.segment(m_start_global_dof, m_num_global_dof);
            r.SetSegment(m_start_global_dof, m_num_global_dof,
                m_A_inv_diag.cwiseProduct(rr));
            return true;
        }
    }

    public void preIterate(ref TwoDScene scene, double dt)
    {
        // Foreach local forces:
        //	compute an A matrix (localized version of compute Integration vars),
        HairFlow flow = m_parent.m_scene.getFilmFlows()[m_hidx];
        int ndof = scene.isMassSpring() ? flow.size() * DIM : flow.size() * 4 - 1;
        Vectors m = scene.getInterpolatedM();
        VectorXi constraint_idx = flow.getConstraintIdx();
        VectorXi constraint_num = flow.getNumConstraints();
        List<int> global_local =
            m_parent.m_scene.getParticleToHairLocalIndices();
        List<int> particle_hair =
            m_parent.m_scene.getParticleToHairs();

        m_A_nz = new TripletXs(0);
        m_A_nz.reserve(constraint_num[5] + ndof);
        TripletXs m_A_nz_ref = m_parent.m_A_nz;
        int nK = constraint_num[5];
        int base_K = constraint_idx[5];
        for (int i = 0; i < nK; ++i)
        {
            Triplets t = m_A_nz_ref[base_K + i];
            int ip = scene.getVertFromDof(t.getrow());
            int jp = scene.getVertFromDof(t.getcol());
            int ip_local = global_local[ip];
            int jp_local = global_local[jp];
            int idir = scene.getComponent(t.getrow());
            int jdir = scene.getComponent(t.getcol());
            
            double val = 0;
            if (!scene.isFixed(ip) && !scene.isFixed(jp))
            {
                val += t.getvalue() * dt * dt;
            }

            if (val != 0.0)
            {
                if (scene.isMassSpring())
                    m_A_nz.Add(
                        new Triplets(ip_local * DIM + idir, jp_local * DIM + jdir, val));
                else
                    m_A_nz.Add(
                        new Triplets(ip_local * 4 + idir, jp_local * 4 + jdir, val));
            }
        }

        for (int i = 0; i < ndof; ++i)
        {
            if (scene.isFixed(scene.getVertFromDof(m_start_global_dof + i)))
                m_A_nz.Add(new Triplets(i, i, 1.0));
            else
                m_A_nz.Add(new Triplets(i, i, m[m_start_global_dof + i]));
        }

        m_A.setFromTriplets(m_A_nz);

        int nconstraint = constraint_num[0];
        int nconstraint_v = constraint_num[1];
        int base_constraint = constraint_idx[0];
        int base_constraint_v = constraint_idx[1];

        m_J_nz = new TripletXs(0);
        m_J_nz.reserve(constraint_num[2]);
        int nJ = constraint_num[2];
        int base_J = constraint_idx[2];
        TripletXs m_J_nz_ref = m_parent.m_J_nz;
        for (int i = 0; i < nJ; ++i)
        {
            Triplets t = m_J_nz_ref[base_J + i];
            int jp = scene.getVertFromDof(t.getcol());
            if (scene.isFixed(jp))
                continue;
            int jp_local = global_local[jp];
            int jdir = scene.getComponent(t.getcol());
            if (scene.isMassSpring())
                m_J_nz.Add(new Triplets(t.getrow() - base_constraint,
                                          jp_local * DIM + jdir, t.getvalue()));
            else
                m_J_nz.Add(
                    new Triplets(t.getrow() - base_constraint, jp_local * 4 + jdir, t.getvalue()));
        }
        m_J = new SparseXs(nconstraint, ndof);
        m_J.setFromTriplets(m_J_nz);

        m_Jv_nz = new TripletXs(0);
        m_Jv_nz.reserve(constraint_num[3]);
        int nJv = constraint_num[3];
        int base_Jv = constraint_idx[3];
        TripletXs m_Jv_nz_ref = m_parent.m_Jv_nz;
        for (int i = 0; i < nJv; ++i)
        {
            Triplets t = m_Jv_nz_ref[base_Jv + i];
            int jp = scene.getVertFromDof(t.getcol());
            if (scene.isFixed(jp))
                continue;
            int jp_local = global_local[jp];
            int jdir = scene.getComponent(t.getcol());
            if (scene.isMassSpring())
                m_Jv_nz.Add(new Triplets(t.getrow() - base_constraint_v,
                                           jp_local * DIM + jdir, t.getvalue()));
            else
                m_Jv_nz.Add(new Triplets(t.getrow() - base_constraint_v,
                                           jp_local * 4 + jdir, t.getvalue()));
        }
        m_Jv = new SparseXs(nconstraint_v, ndof);
        m_Jv.setFromTriplets(m_Jv_nz);

        m_Jxv_nz = new TripletXs(0);
        m_Jxv_nz.reserve(constraint_num[4]);
        int nJxv = constraint_num[4];
        int base_Jxv = constraint_idx[4];
        TripletXs m_Jxv_nz_ref = m_parent.m_Jxv_nz;
        for (int i = 0; i < nJxv; ++i)
        {
            Triplets t = m_Jxv_nz_ref[base_Jxv + i];
            int jp = scene.getVertFromDof(t.getcol());
            if (scene.isFixed(jp))
                continue;
            int jp_local = global_local[jp];
            int jdir = scene.getComponent(t.getcol());
            if (scene.isMassSpring())
                m_Jxv_nz.Add(new Triplets(t.getrow() - base_constraint_v,
                                            jp_local * DIM + jdir, t.getvalue()));
            else
                m_Jxv_nz.Add(new Triplets(t.getrow() - base_constraint_v,
                                            jp_local * 4 + jdir, t.getvalue()));
        }
        m_Jxv = new SparseXs(nconstraint_v, ndof);
        m_Jxv.setFromTriplets(m_Jxv_nz);

        m_invC_nz = new TripletXs(nconstraint);
        TripletXs m_invC_nz_ref = m_parent.m_invC_nz;
        for (int i = 0; i < nconstraint; ++i)
        {
            Triplets t = m_invC_nz_ref[base_constraint + i];
            m_invC_nz[i] = new Triplets(i, i, t.getvalue());
        }
        m_invC = new SparseXs(nconstraint, nconstraint);
        m_invC.setFromTriplets(m_invC_nz);

        m_invCv_nz = new TripletXs(nconstraint_v);
        TripletXs m_invCv_nz_ref = m_parent.m_invCv_nz;
        for (int i = 0; i < nconstraint_v; ++i)
        {
            Triplets t = m_invCv_nz_ref[base_constraint_v + i];
            m_invCv_nz[i] = new Triplets(i, i, t.getvalue());
        }
        m_invCv = new SparseXs(nconstraint_v, nconstraint_v);
        m_invCv.setFromTriplets(m_invCv_nz);

        // pre-factor it with simplicial LDLT
        m_A += (new SparseXs(m_J.transpose()) * ((m_invC * dt * dt) * m_J));

        if (nconstraint_v > 0)
        {
            m_A +=
                (new SparseXs(m_Jv.transpose()) * ((m_invCv * dt) * (m_Jv + m_Jxv * dt)));
        }

        m_solver.compute(m_A);

        if (!m_parent.m_use_preconditioner)
        {
            m_A_inv_diag = m_A.diagonal().cwiseInverse();
        }
    }
}

public class SimplicialLDLT
{
    private SparseXs L;
    private Dictionary<int, double> D;

    public SimplicialLDLT() { }

    public void compute(SparseXs matrix)
    {
        int n = matrix.Rows;
        L = new SparseXs(n, n);
        D = new Dictionary<int, double>();

        for (int i = 0; i < n; i++)
        {
            double sum = matrix.GetValue(i, i);
            for (int k = 0; k < i; k++)
            {
                double l_ik = L.GetValue(i, k);
                if (D.TryGetValue(k, out double d_k))
                    sum -= l_ik * l_ik * d_k;
            }
            D[i] = sum;

            for (int j = i + 1; j < n; j++)
            {
                if (matrix.GetValue(j, i) != 0)
                {
                    sum = matrix.GetValue(j, i);
                    for (int k = 0; k < i; k++)
                    {
                        sum -= L.GetValue(j, k) * L.GetValue(i, k) * (D.ContainsKey(k) ? D[k] : 0);
                    }
                    if (D[i] != 0)
                        L.insert(j, i, sum / D[i]);
                }
            }
            L.insert(i, i, 1.0);
        }
    }

    public Vectors solve(Vectors v)
    {
        int n = v.Size;
        double[] y = new double[n];
        Vectors result = new Vectors(n);

        // Forward substitution
        for (int i = 0; i < n; i++)
        {
            double sum = v[i];
            for (int k = 0; k < i; k++)
                sum -= L.GetValue(i, k) * y[k];
            y[i] = sum;
        }

        // Diagonal solve
        double[] z = new double[n];
        for (int i = 0; i < n; i++)
        {
            if (D.TryGetValue(i, out double d_i) && d_i != 0)
                z[i] = y[i] / d_i;
            else
                continue;
        }

        // Back substitution
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = z[i];
            for (int k = i + 1; k < n; k++)
                sum -= L.GetValue(k, i) * result[k];
            result[i] = sum;
        }

        return result;
    }
}
