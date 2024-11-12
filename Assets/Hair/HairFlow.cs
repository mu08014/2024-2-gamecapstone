using System;
using System.Collections.Generic;

enum LIQUID_SIM_TYPE
{
    LST_SHALLOW,
    LST_COUNT
}

struct HairParticleBridge
{
    Vectors vel;
    double volume;
    double alpha;
    int pidx;
    int eidx;
};

struct HairPlaneIntersection
{
    int lidx;
    int pidx_0;
    int pidx_1;
    int flow_idx;
    int hair_idx;
    int prev_inter_idx;
    int next_inter_idx;
    double dist_to_prev;
    double dist_to_next;
    double flow_height;
    double alpha;
    double pos_y;
    double porosity;
    VectorXs jac;
    MatrixXs hess;
};


public abstract class HairFlow
{
    protected static int m_flow_counter;

    protected TwoDScene m_parent;

    protected List<int> m_particle_indices = new List<int>();
    protected List<bool> m_particle_state = new List<bool>();
    protected List<int> m_edge_indices = new List<int>();
    protected List<Tuple<int, int>> m_internal_edges = new List<Tuple<int, int>>();
    protected List<Tuple<int, int>> m_global_edges = new List<Tuple<int, int>>();
    protected List<List<int>> m_particle_to_edges = new List<List<int>>();
    protected Dictionary<int, int> m_global_to_local = new Dictionary<int, int>();

    protected VectorXs m_liquid_phi;
    protected VectorXs m_porosity;
    protected VectorXs m_eta;
    protected VectorXs m_avg_eta;
    protected VectorXs m_edge_eta;
    protected VectorXs m_u;
    protected VectorXs m_stored_friction_coeff;

    protected MatrixXs m_accel_v;
    protected VectorXs m_area_v;
    protected VectorXs m_area_v_hair;
    protected VectorXs m_area_e;
    protected VectorXs m_normal_e;
    protected VectorXs m_normal_v;
    protected VectorXs m_rad_vec;
    protected VectorXs m_edge_rad_vec;
    protected MatrixXs m_dir_v;
    protected MatrixXs m_dir_f;

    protected MatrixXs m_actual_u_f;
    protected MatrixXs m_actual_u_v;

    protected MatrixXs m_c_v;
    protected MatrixXs m_c_star;

    protected double m_avg_area_e;
    protected double m_min_area_e;
    protected double m_max_area_e;

    protected double m_min_area_v;
    protected double m_max_area_v;
    protected double m_max_eta;
    protected double m_min_eta;

    protected VectorXi m_constraint_starts;
    protected VectorXi m_constraint_length;

    protected List<List<int>> m_edge_bridges = new List<List<int>>();

    protected int m_flow_index;

    public int DIM;

    public enum STATE
    {
        NORMAL,
        ENDPOINT,
        STAR
    }

    public HairFlow(ref TwoDScene parent, in List<int> involved_particles, in VectorXs eta, in List<bool> particle_state, int dim)
    {
        m_parent = parent;
        m_particle_indices = involved_particles;
        m_eta = eta;
        m_particle_state = particle_state;
        m_flow_index = m_flow_index + 1;
        m_porosity = new VectorXs(m_eta.size());
        m_liquid_phi = new VectorXs(m_eta.size());

        VectorXs radii = m_parent.getRadii();

        int np = m_porosity.size();
        for (int i = 0; i < np; i++)
        {
            double r = radii[m_particle_indices[i]];
            m_porosity[i] = 1.0 - r * r / ((m_eta[i] + r) * (m_eta[i] + r));
        }

        m_avg_eta = m_eta;

        DIM = dim;
    }

    public VectorXs getVelocity()
    {
        return m_u;
    }

    public virtual double getAvgAreaE()
    {
        return m_avg_area_e;
    }

    public virtual double getMinAreaE()
    {
        return m_min_area_e;
    }

    public virtual double getMaxAreaE()
    {
        return m_max_area_e;
    }

    public virtual double getMinAreaV()
    {
        return m_min_area_v;
    }
    public virtual double getMaxAreaV()
    {
        return m_max_area_v;
    }

    public virtual double getMinEta()
    {
        return m_min_eta;
    }

    public virtual double getMaxEta()
    {
        return m_max_eta;
    }

    public virtual VectorXs getEta()
    {
        return m_eta;
    }

    public virtual VectorXs getAvgEta()
    {
        return m_avg_eta;
    }

    public virtual MatrixXs getActualVertexVelocity()
    {
        return m_actual_u_v;
    }

    public virtual VectorXs getPorosity()
    {
        return m_porosity;
    }

    public virtual VectorXs getAreaV()
    {
        return m_area_v;
    }

    public virtual VectorXs getAreaVHair()
    {
        return m_area_v_hair;
    }

    public virtual MatrixXs getAccelV()
    {
        return m_accel_v;
    }

    public virtual VectorXs getAreaE()
    {
        return m_area_e;
    }

    public virtual VectorXs getNormalV()
    {
        return m_normal_v;
    }

    public virtual VectorXs getRadiiV()
    {
        return m_rad_vec;
    }

    public virtual VectorXs getRadiiE()
    {
        return m_edge_rad_vec;
    }

    public virtual VectorXs getNormalE()
    {
        return m_normal_e;
    }

    public virtual MatrixXs getTangentV()
    {
        return m_dir_v;
    }

    public virtual MatrixXs getTangentE()
    {
        return m_dir_f;
    }

    public abstract double getTotalLength();

    public virtual void adjustVolumeGlobal(in double prop)
    {
        for (int i = 0; i < m_eta.size(); i++)
        {
            m_eta[i] = Math.Sqrt(Math.Max(0.0, prop) * m_eta[i] * (m_eta[i] + m_rad_vec[i] * 2.0) +
                m_rad_vec[i] * m_rad_vec[i]) - 
                m_rad_vec[i];
        }
    }

    public virtual Dictionary<int, int> getGlobalToLocal()
    {
        return m_global_to_local;
    }

    public virtual List<List<int>> getEdgeBridges()
    {
        return m_edge_bridges;
    }

    public virtual List<Tuple<int, int>> getLocalEdges()
    {
        return m_internal_edges;
    }

    public virtual List<Tuple<int, int>> getGlobalEdges()
    {
        return m_global_edges;
    }

    public virtual List<int> getEdgeIndices()
    {
        return m_edge_indices;
    }

    public virtual List<int> getParticleIndices()
    {
        return m_particle_indices;
    }

    public abstract void updateGeometricState(in VectorXs x, in VectorXs v,
                                    ref FluidSim fluidsim);

    public abstract void updateToFilteredGrid(in VectorXs x, in VectorXs v,
                                    ref FluidSim fluidsim, in double dt,
                                    int ibuffer);

    public abstract void updateFromFilteredGrid(in VectorXs x, ref VectorXs v,
                                      ref FluidSim fluidsim, in double dt);

    public abstract void advance(in VectorXs x, in double dt);

    public abstract void add_force(in VectorXs x, in VectorXs accel,
                         ref FluidSim fluidsim, in double dt);

    public abstract void updateHairMass();

    public abstract void updateHairFlowHeight(in double dt);

    public abstract void preUpdateHairFlowHeight(in double dt);

    public abstract void postUpdateHairFlowHeight(in double dt);

    public abstract void updateReservoir(ref FluidSim fluidsim, in VectorXs x,
                               in VectorXs v, in double dt);

    public virtual List<bool> getState()
    {
        return m_particle_state;
    }

    public abstract void resizeSystem();

    public virtual int index()
    {
        return m_flow_index;
    }

    public virtual int find(int idx_global)
    {
        if (m_global_to_local.TryGetValue(idx_global, out var localval))
        {
            return localval;
        }
        else
        {
            return -1;
        }
    }

    public virtual int size()
    {
        return m_global_to_local.Count;
    }

    public virtual void getAffectedVars(int pidx, ref HashSet<int> vars)
    {
        int ip = m_parent.getVertFromDof(pidx);
        if (m_parent.getComponent(pidx) == 0 && find(ip) != -1)
        {
            vars.Add(pidx);
        }
    }

    public virtual Vectors computeHairLiquidMomentum(in VectorXs v)
    {
        int ne = m_global_to_local.Count;

        double rho = m_parent.getLiquidDensity();

        Vectors momentum = new Vectors(DIM); // 인스턴스 생성 시 0으로 초기화

        for (int i = 0; i < ne; i++)
        {
            int i0 = m_internal_edges[i].Item1;
            int i1 = m_internal_edges[i].Item2;
            int global_i0 = m_global_edges[i].Item1;
            int global_i1 = m_global_edges[i].Item2;
            Vectors v0 = v.segment(m_parent.getDof(global_i0), DIM);
            Vectors v1 = v.segment(m_parent.getDof(global_i1), DIM);

            double H0 = m_eta[i0] + m_rad_vec[i0];
            double H1 = m_eta[i1] + m_rad_vec[i1];

            double vol = WetHairParameter.M_PI *
                (H0 * H0 + H1 * H1 - m_rad_vec[i0] * m_rad_vec[i0] -
                m_rad_vec[i1] * m_rad_vec[i1]) *
                0.5 * m_area_e[i];
            double mass = rho * vol;

            Vectors ve =
                (v0 + v1) * 0.5 + m_dir_f.row(i).normalized() * m_u[i];

            momentum += ve * mass;
        }

        return momentum;
    }

    public abstract VectorXs computeHairLiquidAngularMomentum(
      in VectorXs x, in VectorXs v, ref FluidSim fluidsim);

    public virtual Vectors computeHairDragForce(in VectorXs v)
    {
        int ne = m_global_edges.Count;

        double rho = m_parent.getLiquidDensity();

        Vectors momentum = new Vectors(DIM);

        for (int i = 0; i < ne; i++)
        {
            int i0 = m_internal_edges[i].Item1;
            int i1 = m_internal_edges[i].Item2;
            int global_i0 = m_global_edges[i].Item1;
            int global_i1 = m_global_edges[i].Item2;
            Vectors v0 = v.segment(m_parent.getDof(global_i0), DIM);
            Vectors v1 = v.segment(m_parent.getDof(global_i1), DIM);

            double H0 = m_eta[i0] + m_rad_vec[i0];
            double H1 = m_eta[i1] + m_rad_vec[i1];

            double vol = WetHairParameter.M_PI *
                (H0 * H0 + H1 * H1 - m_rad_vec[i0] * m_rad_vec[i0] -
                m_rad_vec[i1] * m_rad_vec[i1]) *
                0.5 * m_area_e[i];
            double mass = rho * vol;

            Vectors ve = (v0 + v1) * 0.5;

            momentum += ve * mass;
        }

        return momentum;
    }

    public virtual double computeHairLiquidEnergy(in VectorXs v)
    {
        int ne = m_global_edges.Count;

        double rho = m_parent.getLiquidDensity();

        double E = 0;

        for (int i = 0; i < ne; i++)
        {
            int i0 = m_internal_edges[i].Item1;
            int i1 = m_internal_edges[i].Item2;
            int global_i0 = m_global_edges[i].Item1;
            int global_i1 = m_global_edges[i].Item2;
            Vectors v0 = v.segment(m_parent.getDof(global_i0), DIM);
            Vectors v1 = v.segment(m_parent.getDof(global_i1), DIM);

            double H0 = m_eta[i0] + m_rad_vec[i0];
            double H1 = m_eta[i1] + m_rad_vec[i1];

            double vol = WetHairParameter.M_PI *
                (H0 * H0 + H1 * H1 - m_rad_vec[i0] * m_rad_vec[i0] -
                m_rad_vec[i1] * m_rad_vec[i1]) *
                0.5 * m_area_e[i];
            double mass = rho * vol;

            Vectors ve =
                (v0 + v1) * 0.5 + m_dir_f.row(i).normalized() * m_u[i];

            E += mass * ve.squaredNorm();
        }

        return E * 0.5;
    }

    public virtual bool isContained(int pidx)
    {
        if (m_parent.getComponent(pidx) == DIM)
            return false;
        int ip = m_parent.getVertFromDof(pidx);
        return find(ip) != -1;

    }

    public abstract double computeTotalLiquidVol();

    public abstract double computeTotalReservoirVol();

    public abstract double getPoolSize();


    public virtual void read(in double[] data)
    {
        int neta = m_eta.size();
        int k = 0;
        for (int i = 0; i < neta; i++)
        {
            m_eta[i] = data[k++];
        }
        int nu = m_u.size();
        for (int i = 0; i < nu; i++)
        {
            m_u[i] = data[k++];
        }
    }

    public virtual void write(ref List<double> data)
    {
        int neta = m_eta.size();
        for (int i = 0; i < neta; i++)
        {
            data.Add(m_eta[i]);
        }
        int nu = m_u.size();
        for (int i = 0; i < nu; i++)
        {
            data.Add(m_u[i]);
        }
    }

    public virtual void writeReadable()
    {
        // 유니티에 필요 x
    }

    public virtual void readReadable()
    {
        // 유니티에 필요 x
    }

    public virtual int serialized_size()
    {
        return (m_eta.size() + m_u.size()) * sizeof(double);
    }

    public virtual void setConstraintParameters(in VectorXi start,
                                       in VectorXi num)
    {
        m_constraint_starts = start;
        m_constraint_length = num;

    }

    public virtual VectorXi getConstraintIdx()
    {
        return m_constraint_starts;
    }

    public virtual VectorXi getNumConstraints()
    {
        return m_constraint_length;
    }

    public abstract void geodesic_to_local(in double geopos, ref int pidx_low,
                                 ref double alpha);


}
