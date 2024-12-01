using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.Threading.Tasks;
using Unity.Burst;
using Unity.Collections;
using Unity.Jobs;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.UIElements;

public class StrandCompliantManager : SceneStepper
{
    List<StrandCompliantEuler> m_integrators;

    protected double m_dt;
    public TwoDScene m_scene;

    public TripletXs m_A_nz;
    public TripletXs m_J_nz;
    public TripletXs m_Jv_nz;
    public TripletXs m_Jxv_nz;
    protected TripletXs m_Fix_nz;
    protected TripletXs m_M_nz;
    public TripletXs m_invC_nz;
    public TripletXs m_invCv_nz;

    protected TripletXs m_J_interhair_nz;
    protected TripletXs m_Jv_interhair_nz;
    protected SparseXs m_J_interhair;
    protected SparseXs m_Jv_interhair;
    protected SparseXs m_JT_interhair;
    protected SparseXs m_JvT_interhair;
    protected Vectors m_invC_interhair;
    protected Vectors m_invCv_interhair;

    public Vectors m_lambda;
    public Vectors m_lambda_v;
    public Vectors m_gradU;
    public Vectors m_Phi;
    protected Vectors m_Phi_interhair;
    public Vectors m_Phi_v;
    protected Vectors m_Phiv_interhair;
    protected Vectors m_vplus;

    protected Vectors m_dvi;
    protected Vectors m_dxi;
    protected Vectors m_dv;
    protected Vectors m_dx;
    protected Vectors m_dx_scripted;

    protected Vectors m_v;
    protected Vectors m_r;
    protected Vectors m_p;
    protected Vectors m_q;
    protected Vectors m_t;
    protected Vectors m_z;
    protected Vectors m_rhs;

    protected Vectors m_cv_buffer;
    protected Vectors m_c_buffer;

    protected VectorXi m_interhair_idx;
    protected VectorXi m_interhair_num;

    protected int m_max_num_newton;
    protected int m_max_num_iters;
    protected double m_criterion;

    protected bool m_compute_interhair;
    public bool m_use_preconditioner;

    public StrandCompliantManager(
        ref TwoDScene scene, int max_newton, int max_iters, double criterion,
        bool compute_interhair, bool use_preconditioner, int dim = 3)
    {
        DIM = dim;

        m_max_num_newton = max_newton;
        m_max_num_iters = max_iters;
        m_criterion = criterion;
        m_scene = scene;
        m_compute_interhair = true;//compute_interhair;

        int numhairs = scene.getNumFlows();

        m_integrators = new List<StrandCompliantEuler>(numhairs);
        for (int i = 0; i < numhairs; ++i)
        {
            m_integrators.Add(new StrandCompliantEuler(this, i, dim));
        }

        m_dt = 0.0;

        int ndof = scene.getNumDofs();
        m_gradU = new Vectors(ndof);
        m_r = new Vectors(ndof);
        m_v = new Vectors(ndof);
        m_p = new Vectors(ndof);
        m_q = new Vectors(ndof);
        m_t = new Vectors(ndof);
        m_z = new Vectors(ndof);
        m_rhs = new Vectors(ndof);
        m_vplus = new Vectors(ndof);

        m_dv = new Vectors(ndof);
        m_dx = new Vectors(ndof);
        m_dvi = new Vectors(ndof);
        m_dxi = new Vectors(ndof);
        m_dx_scripted = new Vectors(ndof);

        m_timing_statistics = new List<double>(6);
        for (int i = 0; i < 6; i++)
        {
            m_timing_statistics.Add(0.0);
        }

    }

    public bool stepScene(ref TwoDScene scene, double dt,
        bool updatePreCompute)
    {
        return stepSceneLinear(ref scene, dt, updatePreCompute);
        //m_max_num_newton == 0 이 false 이면 stepSceneNonlinear도 해야됨
    }

    public bool stepSceneLinear(ref TwoDScene scene, double dt, bool updatePreCompute)
    {
        double t0 = Time.realtimeSinceStartup;
        double t1;

        m_dt = dt;
        m_scene = scene;

        Vectors x = scene.getX();
        Vectors v = scene.getV();

        m_dx = v * dt;
        m_dv.setZero();
        m_dxi = m_dv.Clone();

        m_dx_scripted.setZero();

        int np = m_scene.getNumParticles();

        for (int i = 0; i < np; ++i)
        {
            if (scene.isFixed(i))
            {
                int numdofs = scene.isMassSpring() || scene.isTip(i) ? DIM : 4;
                m_dx_scripted.SetSegment(scene.getDof(i), numdofs,
                    v.segment(scene.getDof(i), numdofs) * dt);
            }
        }

        scene.preComputeLocal(m_dx_scripted, m_dv, dt);
       
        for (int i = 0; i < m_dv.Size; i++)
        {
            if (m_dv[i] > 1)
                m_dv[i] = 1;
            else if (m_dv[i] < -1)
                m_dv[i] = -1;
        }
        
        t1 = Time.realtimeSinceStartup;
        m_timing_statistics[0] += (t1 - t0);
        t0 = t1;

        localUpdateNumConstraints(m_dx_scripted, m_dv, dt);
        
        for (int i = 0; i < m_dv.Size; i++)
        {
            if (m_dv[i] > 1)
                m_dv[i] = 1;
            else if (m_dv[i] < -1)
                m_dv[i] = -1;
        }
        
        t1 = Time.realtimeSinceStartup;
        m_timing_statistics[1] += (t1 - t0);
        t0 = t1;

        localPreIterate(ref scene, dt);
        computeRHSLocal(ref scene, dt, ref m_rhs);
        t1 = Time.realtimeSinceStartup;
        m_timing_statistics[2] += (t1 - t0);  // local Matrix Composition
        t0 = t1;

        localStepScene(ref scene, dt, ref m_vplus, m_rhs);
        
        for (int i = 0; i < m_vplus.Size; i++)
        {
            if (m_vplus[i] > 1)
                m_vplus[i] = 1;
            else if (m_vplus[i] < -1)
                m_vplus[i] = -1;
        }
        
        m_dv = m_vplus - v;
        m_dx = m_vplus * dt;
        
        for (int i = 0; i < m_dv.Size; i++)
        {
            if (m_dv[i] > 1)
                m_dv[i] = 1;
            else if (m_dv[i] < -1)
                m_dv[i] = -1;
        }
        
        if (!m_compute_interhair)
        {
            localUpdateLambda(ref scene, m_dx, m_dv, dt);

            m_scene.storeLambda(m_lambda, m_lambda_v);

            updateNextV(ref scene, m_vplus);

            m_next_x = x + dt * v;

            return true;
        }

        t0 = Time.realtimeSinceStartup;

        scene.preComputeInterhair(m_dx, m_dv, dt);
        
        for (int i = 0; i < m_dv.Size; i++)
        {
            if (m_dv[i] > 1)
                m_dv[i] = 1;
            else if (m_dv[i] < -1)
                m_dv[i] = -1;
        }
        
        interhairUpdateNumConstraints(m_dxi, m_dv, dt);
        
        for (int i = 0; i < m_dv.Size; i++)
        {
            if (m_dv[i] > 1)
                m_dv[i] = 1;
            else if (m_dv[i] < -1)
                m_dv[i] = -1;
        }
        
        computeRHSIncrementalInterhair(scene, dt, m_rhs, m_vplus);
        
        for (int i = 0; i < m_vplus.Size; i++)
        {
            if (m_vplus[i] > 1)
                m_vplus[i] = 1;
            else if (m_vplus[i] < -1)
                m_vplus[i] = -1;
        }
        
        t1 = Time.realtimeSinceStartup;
        m_timing_statistics[3] += (t1 - t0);  // compute Interhair Variables
        t0 = t1;

        // Ap = Ax0
        computeAp(scene, m_vplus, m_r, dt);
        
        for (int i = 0; i < m_vplus.Size; i++)
        {
            if (m_vplus[i] > 1)
                m_vplus[i] = 1;
            else if (m_vplus[i] < -1)
                m_vplus[i] = -1;
        }
        
        // r = b - Ax0
        m_r = m_rhs - m_r;

        // solve Mz = r
        localPreconditionScene(scene, dt, m_z, m_r);

        // p = z
        m_p = m_z;

        double rho = m_r.dot(m_z);

        double res_norm = m_r.norm();

        int iter = 0;
        if (res_norm < m_criterion)
        {
            localUpdateLambda(ref scene, m_dx, m_dv, dt);

            m_scene.storeLambda(m_lambda, m_lambda_v);

            updateNextV(ref scene, m_vplus);
            m_next_x = x + dt * v;

            return true;
        }
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // q = Ap
        computeAp(scene, m_p, m_q, dt);

        // alpha = rho / (p, q)
        double alpha = rho / m_p.dot(m_q);

        // x = x + alpha*p
        m_vplus += m_p * alpha;

        // r = r - alpha*q
        m_r -= m_q * alpha;

        res_norm = m_r.norm();

        double rho_old, beta;
        for (; iter < m_max_num_iters && res_norm > m_criterion; ++iter)
        {
            rho_old = rho;

            // solve Mz = r
            localPreconditionScene(scene, dt, m_z, m_r);

            rho = m_r.dot(m_z);

            beta = rho / rho_old;

            // p = beta * p + z
            m_p = m_z + m_p * beta;

            // q = Ap
            computeAp(scene, m_p, m_q, dt);

            // alpha = rho / (p, q)
            alpha = rho / m_p.dot(m_q);

            // x = x + alpha*p
            m_vplus += m_p * alpha;

            // r = r - alpha*q
            m_r -= m_q * alpha;

            res_norm = m_r.norm();
        }
        
        for (int i = 0; i < m_vplus.Size; i++)
        {
            if (m_vplus[i] > 1)
                m_vplus[i] = 1;
            else if (m_vplus[i] < -1)
                m_vplus[i] = -1;
        }
        
        m_dv = m_vplus - v;
        m_dx = m_vplus * dt;
        
        for (int i = 0; i < m_dv.Size; i++)
        {
            if (m_dv[i] > 1)
                m_dv[i] = 1;
            else if (m_dv[i] < -1)
                m_dv[i] = -1;
        }
        
        localUpdateLambda(ref scene, m_dx, m_dv, dt);
        
        for (int i = 0; i < m_dv.Size; i++)
        {
            if (m_dv[i] > 1)
                m_dv[i] = 1;
            else if (m_dv[i] < -1)
                m_dv[i] = -1;
        }
        
        m_scene.storeLambda(m_lambda, m_lambda_v);

        updateNextV(ref scene, m_vplus);
        
        for (int i = 0; i < m_vplus.Size; i++)
        {
            if (m_vplus[i] > 1)
                m_vplus[i] = 1;
            else if (m_vplus[i] < -1)
                m_vplus[i] = -1;
        }
        
        m_next_x = x + dt * v;

        return true;
    }

    public void localUpdateNumConstraints(in Vectors dx, in Vectors dv, double dt)
    {
        int num_pos = 0;
        int num_vel = 0;
        int num_J = 0;
        int num_Jv = 0;
        int num_Jxv = 0;
        int num_tildeK = 0;
        m_scene.updateNumConstraintsLocal(ref num_pos, ref num_vel, ref num_J, ref num_Jv, ref num_Jxv,
                                           ref num_tildeK);

        m_lambda = new Vectors(num_pos);
        m_lambda_v = new Vectors(num_vel);
        m_A_nz = new TripletXs(num_tildeK);
        m_J_nz = new TripletXs(num_J);
        m_Jv_nz = new TripletXs(num_Jv);
        m_Jxv_nz = new TripletXs(num_Jxv);
        m_invC_nz = new TripletXs(num_pos);
        m_invCv_nz = new TripletXs(num_vel);
        m_Phi = new Vectors(num_pos);
        m_Phi_v = new Vectors(num_vel);

        m_scene.localPostPreprocess(ref m_lambda, ref m_lambda_v, ref m_J_nz, ref m_Jv_nz, ref m_Jxv_nz,
                                     ref m_A_nz, ref m_invC_nz, ref m_invCv_nz, ref m_Phi, ref m_Phi_v,
                                     dx, dv, dt);
    }

    public void localPreIterate(ref TwoDScene scene, double dt)
    {
        TwoDScene local_scene = scene;
        int nhairs = m_integrators.Count;

        Parallel.For(0, nhairs, hidx => {
            m_integrators[hidx].preIterate(ref local_scene, dt);
        });

        scene = local_scene;
    }

    public void computeRHSLocal(ref TwoDScene scene, double dt, ref Vectors b)
    {
        m_gradU.setZero();
        scene.accumulateExternalGradU(ref m_gradU, new Vectors(0), new Vectors(0));
        zeroFixedDoFs(scene, ref m_gradU);

        TwoDScene local_scene = scene;
        Vectors local_b = b;
        // compute local b
        int nhairs = m_integrators.Count;
        Parallel.For(0, nhairs, hidx => {
            m_integrators[hidx].computeRHS(ref local_scene, dt, ref local_b);
        });

        scene = local_scene;
        b = local_b;
    }

    public void zeroFixedDoFs(in TwoDScene scene, ref Vectors vec)
    {
        int nprts = scene.getNumParticles();
        for (int i = 0; i < nprts; ++i)
        {
            if (scene.isFixed(i))
            {
                int numdofs = scene.isMassSpring() || scene.isTip(i) ? DIM : 4;
                vec.segment(scene.getDof(i), numdofs).setZero();
            }
        }
    }

    public void localPreconditionScene(TwoDScene scene, double dt, Vectors r, Vectors b)
    {
        double t0 = Time.realtimeSinceStartup;
        int nhairs = m_integrators.Count;
        Parallel.For(0, nhairs, hidx => {
            m_integrators[hidx].PreconditionScene(ref scene, dt, ref r, b);
        });
        double t1 = Time.realtimeSinceStartup;
        m_timing_statistics[4] += (t1 - t0);  // local Solve
        t0 = t1;
    }

    public void localStepScene(ref TwoDScene scene,
                                                 double dt, ref Vectors r,
                                                 in Vectors b) {
        double t0 = Time.realtimeSinceStartup;
        int nhairs = m_integrators.Count;
        TwoDScene local_scene = scene;
        Vectors local_r = r;
        Vectors local_b = b;

        Parallel.For(0, nhairs, hidx => {
            m_integrators[hidx].stepScene(ref local_scene, dt, ref local_r, local_b);
        });
        double t1 = Time.realtimeSinceStartup;
        m_timing_statistics[4] += (t1 - t0);  // local Solve
        t0 = t1;

        scene = local_scene;
        r = local_r;
    }

    public void localUpdateLambda(ref TwoDScene scene,
                                                    in Vectors dx,
                                                    in Vectors dv,
                                                    double dt)
    {
        TwoDScene local_scene = scene;
        Vectors local_dx = dx;
        Vectors local_dv = dv;
        int nhairs = m_integrators.Count;
        Parallel.For(0, nhairs, hidx => {
            m_integrators[hidx].updateLambda(ref local_scene, local_dx, local_dv, dt);
        });
    }

    public void updateNextV(ref TwoDScene scene, in Vectors vplus)
    {
        int nhairs = m_integrators.Count;
        TwoDScene local_scene = scene;
        Vectors local_vplus = vplus;
        Parallel.For(0, nhairs, hidx => {
            m_integrators[hidx].updateNextV(ref local_scene, local_vplus);
        });
    }

    public void interhairUpdateNumConstraints(in Vectors dx, in Vectors dv, double dt)
    {
        int num_pos = 0;
        int num_vel = 0;
        int num_J = 0;
        int num_Jv = 0;
        int num_Jxv = 0;
        int num_tildeK = 0;
        m_scene.updateNumConstraintsInterHair(ref num_pos, ref num_vel, ref num_J, ref num_Jv,
                                               ref num_Jxv, ref num_tildeK, ref m_interhair_idx,
                                               ref m_interhair_num);

        m_lambda.conservativeResize(num_pos);
        m_lambda_v.conservativeResize(num_vel);
        m_A_nz.resize(num_tildeK);
        m_J_nz.resize(num_J);
        m_Jv_nz.resize(num_Jv);
        m_Jxv_nz.resize(num_Jxv);
        m_invC_nz.resize(num_pos);
        m_invCv_nz.resize(num_vel);
        m_Phi.conservativeResize(num_pos);
        m_Phi_v.conservativeResize(num_vel);

        m_scene.interhairPostPreprocess(m_lambda, m_lambda_v, m_J_nz, m_Jv_nz,
                                         m_Jxv_nz, m_A_nz, m_invC_nz, m_invCv_nz,
                                         m_Phi, m_Phi_v, dx, dv, dt);

        m_c_buffer = new Vectors(m_interhair_num[0]);
        m_cv_buffer = new Vectors(m_interhair_num[1]);

        int ndof = m_scene.getNumDofs();

        int num_J_inter = m_interhair_num[2];
        m_J_interhair_nz = new TripletXs(num_J_inter);
        Parallel.For(0, num_J_inter, i => {
            Triplets t = m_J_nz[i + m_interhair_idx[2]];
            m_J_interhair_nz[i] =
                new Triplets(t.getrow() - m_interhair_idx[0], t.getcol(), t.getvalue());
        });
        m_J_interhair = new SparseXs(m_interhair_num[0], ndof);
        m_J_interhair.setFromTriplets(m_J_interhair_nz);
        m_JT_interhair = m_J_interhair.transpose();
        //m_J_interhair.makeCompressed();
        //m_JT_interhair.makeCompressed();

        int num_Jv_inter = m_interhair_num[3];
        m_Jv_interhair_nz = new TripletXs(num_Jv_inter);
        Parallel.For(0, num_Jv_inter, i => {
            Triplets t = m_Jv_nz[i + m_interhair_idx[3]];
            m_Jv_interhair_nz[i] =
                new Triplets(t.getrow() - m_interhair_idx[1], t.getcol(), t.getvalue());
        });
        m_Jv_interhair = new SparseXs(m_interhair_num[1], ndof);
        m_Jv_interhair.setFromTriplets(m_Jv_interhair_nz);
        m_JvT_interhair = m_Jv_interhair.transpose();
        //m_Jv_interhair.makeCompressed();
        //m_JvT_interhair.makeCompressed();

        m_invC_interhair = new Vectors(m_interhair_num[0]);
        Parallel.For(0, m_interhair_num[0], i => {
            m_invC_interhair[i] = m_invC_nz[i + m_interhair_idx[0]].getvalue();
        });

        m_invCv_interhair = new Vectors(m_interhair_num[1]);
        Parallel.For(0, m_interhair_num[1], i => {
            m_invCv_interhair[i] = m_invCv_nz[i + m_interhair_idx[1]].getvalue();
        });

        m_Phi_interhair = m_Phi.segment(m_interhair_idx[0], m_interhair_num[0]);
        m_Phiv_interhair = m_Phi_v.segment(m_interhair_idx[1], m_interhair_num[1]);
    }

    public void computeRHSIncrementalInterhair(TwoDScene scene, double dt, Vectors b, Vectors vplus)
    {
        if (m_interhair_num[0] > 0)
        {
            CMath.compute_cwiseProduct(ref m_c_buffer, m_invC_interhair,
                                            m_Phi_interhair);
            CMath.accumulateJTPhi_coeff(ref b, -dt, m_c_buffer, m_J_interhair);
        }

        if (m_interhair_num[1] > 0)
        {
            CMath.compute_cwiseProduct(ref m_cv_buffer, m_invCv_interhair,
                                            m_Phiv_interhair);
            CMath.accumulateJTPhi_coeff(ref b, -dt, m_cv_buffer, m_Jv_interhair);
        }
    }

    public void computeAp(TwoDScene scene, Vectors p, Vectors b, double dt)
    {
        double t0 = Time.realtimeSinceStartup;

        // compute local Ap
        int nhairs = m_integrators.Count;
        Parallel.For(
            0, nhairs, hidx => { m_integrators[hidx].computeAp(p, ref b); });

        if (m_interhair_num[2] > 0)
        {
            CMath.computeJTPhi_coeff(ref m_c_buffer, m_invC_interhair, p,
                                          m_JT_interhair);
            CMath.accumulateJTPhi_coeff(ref b, dt * dt, m_c_buffer, m_J_interhair);
        }

        if (m_interhair_num[3] > 0)
        {
            CMath.computeJTPhi_coeff(ref m_cv_buffer, m_invCv_interhair, p,
                                          m_JvT_interhair);
            CMath.accumulateJTPhi_coeff(ref b, dt, m_cv_buffer, m_Jv_interhair);
        }

        double t1 = Time.realtimeSinceStartup;
        m_timing_statistics[5] += (t1 - t0);  // global CG
        t0 = t1;
    }
}
