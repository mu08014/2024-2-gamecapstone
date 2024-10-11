using System;
using System.Collections.Generic;

public class LinearSpringForce : Force
{
    private Tuple<int, int> m_endpoints;

    private double m_k;
    private double m_l0;
    private double m_b;

    private Vectors m_D;
    private Vectors m_lambda_v;
    private Vectors m_lambda;

    private TwoDScene m_scene;

    public int DIM;

    public LinearSpringForce(ref TwoDScene parent, in Tuple<int, int> endpoints, in double k, in double l0, in double b, int dim) : base()
    {
        m_scene = parent;
        dim = DIM;
        m_endpoints = endpoints;
        m_k = k;
        m_l0 = l0;
        m_b = b;
        m_lambda = new Vectors(dim);
        m_lambda_v = new Vectors(dim);
    }

    public override void preCompute(in VectorXs x, in VectorXs v, in VectorXs m, in double dt)
    {
        m_D = (x.segment(m_scene.getDof(m_endpoints.Item1), DIM) - 
            x.segment(m_scene.getDof(m_endpoints.Item2), DIM))
            .normalized() *
            m_l0;
    }

    public override void computeIntegrationVars(in VectorXs x, in VectorXs v, in VectorXs m, ref VectorXs lambda, ref VectorXs lambda_v, ref TripletXs J, ref TripletXs Jv, ref TripletXs Jxv, ref TripletXs tildeK, ref TripletXs stiffness, ref TripletXs damping, ref VectorXs Phi, ref VectorXs Phiv, in double dt)
    {
        Phi.SetSegment(m_internal_index_pos, DIM,
            x.segment(m_scene.getDof(m_endpoints.Item1), DIM) -
            x.segment(m_scene.getDof(m_endpoints.Item2), DIM) - m_D);
        lambda_v.SetSegment(m_internal_index_pos, DIM, m_lambda);

        if (m_b > 0.0)
        {
            Phiv.SetSegment(m_internal_index_vel, DIM,
                v.segment(m_scene.getDof(m_endpoints.Item1), DIM) -
                v.segment(m_scene.getDof(m_endpoints.Item2), DIM));
            lambda_v.SetSegment(m_internal_index_vel, DIM, m_lambda_v);
        }

        for (int r = 0; r < DIM; r++)
        {
            stiffness[m_internal_index_pos + r] =
                new Triplets(m_internal_index_pos + r, m_internal_index_pos + r, m_k);

            J[m_internal_index_pos + r] = 
                new Triplets(m_internal_index_pos + r, m_internal_index_pos + r, m_k);
            J[m_internal_index_J + DIM + r] =
                new Triplets(m_internal_index_pos + r,
                m_scene.getDof(m_endpoints.Item2) + r, -1.0);

            if (m_b > 0.0)
            {
                damping[m_internal_index_vel + r] =
                    new Triplets(m_internal_index_vel + r, m_internal_index_vel + r, m_b);

                Jv[m_internal_index_Jv + r] =
                    new Triplets(m_internal_index_vel + r,
                    m_scene.getDof(m_endpoints.Item1) + r, 1.0);
                Jv[m_internal_index_Jv + DIM + r] =
                    new Triplets(m_internal_index_vel + r,
                    m_scene.getDof(m_endpoints.Item2) + r, -1.0);
            }
        }
    }

    public override void getAffectedVars(int pidx, ref HashSet<int> vars)
    {
        int idir = m_scene.getComponent(pidx);
        if (idir == DIM)
        {
            return;
        }
        int ip = m_scene.getVertFromDof(pidx);

        if (ip == m_endpoints.Item1 || ip == m_endpoints.Item2)
        {
            for (int r = 0; r < DIM; r++)
            {
                vars.Add(m_scene.getDof(m_endpoints.Item1) + r);
                vars.Add(m_scene.getDof(m_endpoints.Item2) + r);
            }
        }
    }

    public override int getAffectedHair(in List<int> particle_to_hairs)
    {
        return particle_to_hairs[m_endpoints.Item1];
    }

    public override bool isContained(int pidx)
    {
        int idir = m_scene.getComponent(pidx);
        if (idir == DIM)
        {
            return false;
        }
        int ip = m_scene.getVertFromDof(pidx);

        if (ip == m_endpoints.Item1 || ip == m_endpoints.Item2)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    public override int numConstraintPos()
    {
        return DIM;
    }

    public override int numConstraintVel()
    {
        return m_b > 0.0 ? DIM : 0;
    }

    public override int numJ()
    {
        return DIM * 2;
    }

    public override int numJv()
    {
        return m_b > 0.0 ? DIM * 2 : 0;
    }

    public override int numJxv()
    {
        return 0;
    }

    public override int numTildeK()
    {
        return 0;
    }


    public override bool isParallelized() { return false; }

    public override bool isPrecomputationParallelized() { return false; }

    public override void storeLambda(in VectorXs lambda, in VectorXs lamda_v)
    {
        m_lambda = lambda.segment(m_internal_index_pos, DIM);
        if (m_b > 0.0)
        {
            m_lambda_v = lamda_v.segment(m_internal_index_vel, DIM);
        }
    }

    public override string name()
    {
        return "linearspringforce";
    }
}
