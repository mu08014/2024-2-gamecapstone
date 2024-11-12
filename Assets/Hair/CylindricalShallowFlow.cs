using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

public class CylindricalShallowFlow : HairFlow
{
    protected SparseXs m_Gv;
    protected SparseXs m_iGv;
    protected SparseXs m_Gv_hair;
    protected SparseXs m_iGv_hair;
    protected SparseXs m_Gf;
    protected SparseXs m_W_fv;
    protected SparseXs m_W_fv_hair;
    protected SparseXs m_gradF;
    protected SparseXs m_divV;
    protected SparseXs m_divV_hair;
    protected SparseXs m_L;
    protected SparseXs m_dir_f_expand;
    protected SparseXs m_D;
    protected SparseXs m_A;
    protected SparseXs m_LHS;
    protected SparseXs m_sp0;
    protected SparseXs m_sp1;

    protected SparseXs m_area_v_hair_flat;  // used by igl::active_set for volume computation

    protected MatrixXs m_accel_f;

    protected VectorXs m_delta_sol;  // used by igl::active_set for result
    protected VectorXs m_edge_x;
    protected VectorXs m_area_e_accu;
    protected VectorXs m_area_e_inv_mapping;
    protected VectorXs m_uvert;
    protected VectorXs m_ustar;
    protected VectorXs m_accel_e;
    protected VectorXs m_cap_e;
    protected MatrixXs m_rhs;
    protected MatrixXs m_rhs_plus;
    protected MatrixXs m_sol;
    protected VectorXs m_mass_v;

    protected VectorXs m_porosity_e;
    protected VectorXs m_pressure;

    protected VectorXs m_old_eta;
    protected VectorXs m_old_area_v;
    protected VectorXs m_old_rad_vec;

    protected VectorXs m_liquid_mass;
    //Eigen::SparseLU<SparseXs> m_solver;

    protected double m_sum_area_e;
    protected double m_pool_liquid;

    public CylindricalShallowFlow(
        ref TwoDScene parent, in List<int> involved_particles,
        in VectorXs eta, in List<bool> particle_state) : base(ref parent, involved_particles, eta, particle_state, 3)
    {
        int np = involved_particles.Count;
        for (int i = 0; i < np; ++i)
        {
            m_global_to_local[involved_particles[i]] = i;
        }

        m_particle_to_edges = new List<List<int>>(np);
        for (int i = 0; i < np; ++i)
        {
            m_particle_to_edges.Add(new List<int>());
        }

        var edges = m_parent.getEdges();
        int ne = edges.Count;
        for (int i = 0; i < ne; ++i)
        {
            var e = edges[i];
            bool itr0_success = m_global_to_local.TryGetValue(e.Item1, out var itr0);
            bool itr1_success = m_global_to_local.TryGetValue(e.Item2, out var itr1);
            if (itr0_success && itr1_success)
            {
                m_edge_indices.Add(i);
                m_global_edges.Add(e);
                m_internal_edges.Add(
                    new Tuple<int, int>(itr0, itr1));
                m_particle_to_edges[itr0].Add(i);
                m_particle_to_edges[itr1].Add(i);
            }
        }

        resizeSystem();
        m_u.setZero();
        m_uvert.setZero();
        m_rad_vec.setZero();
        m_area_v.setZero();

        m_porosity.setZero();
        m_c_v.setZero();
        m_c_star.setZero();

        compute_SP0_matrix(m_sp0);
        compute_SP1_matrix(m_sp1);

        copy_rest_mass(m_parent.getHairRestMass(),
                            m_particle_indices, m_mass_v,
                            m_parent);

        compute_edge_val(m_eta, m_internal_edges, ref m_edge_eta);

        compute_rest_porosity(m_eta, m_parent.getRadii(), involved_particles, ref m_porosity);
        m_pool_liquid = 0.0;
    }

    private void compute_rest_porosity(in VectorXs eta, in VectorXs global_radii, in List<int> indices, ref VectorXs porosity)
    {
        int np = eta.size();
        for (int i = 0; i < np; ++i)
        {
            double ra = global_radii[indices[i]];
            porosity[i] = 1.0 - CMath.clamp(ra * ra / (eta[i] * eta[i]), 0.0, 1.0);
        }
    }

    private void compute_edge_val(in VectorXs vert_val, in List<Tuple<int, int>> edges, ref VectorXs edge_val)
    {
        int ne = edges.Count;

        for (int i = 0; i < ne; ++i)
        {
            int i0 = edges[i].Item1;
            int i1 = edges[i].Item2;
            edge_val[i] = (vert_val[i0] + vert_val[i1]) * 0.5;
        }
    }

    private void compute_SP0_matrix(SparseXs sp0)
    {
        TripletXs tri = new TripletXs();
        int np = sp0.Rows;
        tri.reserve(np);
        for (int i = 1; i < np - 1; ++i)
        {
            tri.Add(new Triplets(i, i, 1.0));
        }

        sp0.setFromTriplets(tri);
    }

    private void compute_SP1_matrix(SparseXs sp1)
    {
        TripletXs tri = new TripletXs();
        int np = sp1.Rows;
        tri.reserve(np);
        for (int i = 1; i < np - 1; ++i)
        {
            tri.Add(new Triplets(i, i, 1.0));
        }
        tri.Add(new Triplets(0, 1, 1.0));
        tri.Add(new Triplets(np - 1, np - 2, -1.0));

        sp1.setFromTriplets(tri);
    }

    private void copy_rest_mass(VectorXs rest_mass, List<int> indices, VectorXs mass_v, TwoDScene parent)
    {
        int np = mass_v.size();

        for (int i = 0; i < np; ++i)
        {
            mass_v[i] = rest_mass[parent.getDof(indices[i])];
        }
    }

    public override void add_force(in VectorXs x, in VectorXs accel, ref FluidSim fluidsim, in double dt)
    {
        throw new System.NotImplementedException();
    }

    public override void advance(in VectorXs x, in double dt)
    {
        throw new System.NotImplementedException();
    }

    public override VectorXs computeHairLiquidAngularMomentum(in VectorXs x, in VectorXs v, ref FluidSim fluidsim)
    {
        throw new System.NotImplementedException();
    }

    public override double computeTotalLiquidVol()
    {
        throw new System.NotImplementedException();
    }

    public override double computeTotalReservoirVol()
    {
        throw new System.NotImplementedException();
    }

    public override void geodesic_to_local(in double geopos, ref int pidx_low, ref double alpha)
    {
        throw new System.NotImplementedException();
    }

    public override double getPoolSize()
    {
        throw new System.NotImplementedException();
    }

    public override double getTotalLength()
    {
        throw new System.NotImplementedException();
    }

    public override void postUpdateHairFlowHeight(in double dt)
    {
        throw new System.NotImplementedException();
    }

    public override void preUpdateHairFlowHeight(in double dt)
    {
        throw new System.NotImplementedException();
    }

    public override void resizeSystem()
    {
        int num_particles = m_particle_indices.Count;
        int num_edges = m_internal_edges.Count;

        m_area_v = new VectorXs(num_particles);
        m_area_v_hair = new VectorXs(num_particles);
        m_area_e = new VectorXs(num_edges);

        m_dir_f = new MatrixXs(num_edges, DIM);
        m_dir_v = new MatrixXs(num_particles, DIM);

        m_actual_u_f = new MatrixXs(num_edges, DIM);
        m_actual_u_v = new MatrixXs(num_particles, DIM);

        m_edge_rad_vec = new VectorXs(num_edges);
        m_rad_vec = new VectorXs(num_particles);
        m_stored_friction_coeff = new VectorXs(num_edges);

        m_area_v_hair_flat = new SparseXs(1, num_particles);
        m_Gv = new SparseXs(num_particles, num_particles);
        m_iGv = new SparseXs(num_particles, num_particles);
        m_Gv_hair = new SparseXs(num_particles, num_particles);
        m_iGv_hair = new SparseXs(num_particles, num_particles);
        m_Gf = new SparseXs(num_edges * DIM, num_edges * DIM);
        m_W_fv = new SparseXs(num_particles, num_edges);
        m_W_fv_hair = new SparseXs(num_particles, num_edges);
        m_gradF = new SparseXs(num_edges * DIM, num_particles);
        m_divV = new SparseXs(num_particles, num_edges * DIM);
        m_dir_f_expand = new SparseXs(num_edges, num_edges * DIM);
        m_D = new SparseXs(num_particles, num_particles);
        m_A = new SparseXs(num_particles, num_particles);
        m_LHS = new SparseXs(num_particles, num_particles);
        m_L = new SparseXs(num_particles, num_particles);
        m_sp0 = new SparseXs(num_particles, num_particles);
        m_sp1 = new SparseXs(num_particles, num_particles);

        m_accel_v = new MatrixXs(num_particles, DIM);
        m_accel_f = new MatrixXs(num_edges, DIM);
        m_porosity_e = new VectorXs(num_edges);
        m_pressure = new VectorXs(num_particles);

        m_c_v = new MatrixXs(num_particles, DIM * DIM);
        m_c_star = new MatrixXs(num_edges, DIM * DIM);

        m_edge_x = new VectorXs(num_edges);
        m_u = new VectorXs(num_edges);
        m_ustar = new VectorXs(num_edges);
        m_accel_e = new VectorXs(num_edges);
        m_cap_e = new VectorXs(num_edges);
        m_cap_e.setZero();
        m_edge_eta = new VectorXs(num_edges);

        m_area_e_accu = new VectorXs(num_particles);
        m_area_e_inv_mapping = new VectorXs(num_particles);
        m_uvert = new VectorXs(num_particles);

        m_old_rad_vec = new VectorXs(num_particles);
        m_old_area_v = new VectorXs(num_particles);
        m_delta_sol = new VectorXs(num_particles);

        if (m_parent.isMassSpring())
        {
            m_rhs = new MatrixXs(num_particles, DIM + 1);
            m_rhs_plus = new MatrixXs(num_particles, DIM + 1);
            m_sol = new MatrixXs(num_particles, DIM + 1);
        }
        else
        {
            m_rhs = new MatrixXs(num_particles, DIM + 2);
            m_rhs_plus = new MatrixXs(num_particles, DIM + 2);
            m_sol = new MatrixXs(num_particles, DIM + 2);
        }

        m_liquid_mass = new VectorXs(num_particles);
        m_mass_v = new VectorXs(num_particles);

        m_old_eta = new VectorXs(num_particles);

        m_edge_bridges = new List<List<int>>(num_edges);
    }

    public override void updateFromFilteredGrid(in VectorXs x, ref VectorXs v, ref FluidSim fluidsim, in double dt)
    {
        throw new System.NotImplementedException();
    }

    public override void updateGeometricState(in VectorXs x, in VectorXs v, ref FluidSim fluidsim)
    {
        throw new System.NotImplementedException();
    }

    public override void updateHairFlowHeight(in double dt)
    {
        throw new System.NotImplementedException();
    }

    public override void updateHairMass()
    {
        throw new System.NotImplementedException();
    }

    public override void updateReservoir(ref FluidSim fluidsim, in VectorXs x, in VectorXs v, in double dt)
    {
        throw new System.NotImplementedException();
    }

    public override void updateToFilteredGrid(in VectorXs x, in VectorXs v, ref FluidSim fluidsim, in double dt, int ibuffer)
    {
        throw new System.NotImplementedException();
    }
}
