using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;
using System.Threading;
using System.Threading.Tasks;
using Unity.Burst;
using Unity.Collections;
using static FurMesh;

public struct Script
{
    enum TARGET
    {
        CAMERA,
        ROOT,
        SOLID,
        SOURCE,
        CURLRADIUS,
        CURLDENSITY,
        ALL,
        TARGET_COUNT
    }

    enum TYPE
    {
        ROTATE,
        TRANSLATE,
        SCALE,
        ARSORB,
        TYPE_COUNT
    }

    enum FUNC
    {
        CUBIC,
        COSINE,
        WENO,
        FUNC_COUNT
    }

    TARGET target;
    TYPE type;
    FUNC func;
    int index;
    Vectors v;
    Vectors origin;
    double start;
    double end;
    double ease_start;
    double ease_end;
    double amplitude;
    double frequencey;
    double base_dt;
    double base_pos;
    bool updateSDF;
    bool transform_global;

    List<double> base_vertices;
}

public enum MASS_UPDATE_MODE
{
    MUM_NONE,
    MUM_MASS_ONLY,
    MUM_DIRECT_DIV,
    MUM_MOMENTUM,
    MUM_COUNT
}

[Serializable]
public class WetHairParameter
{
    public double dt;
    public double rho;
    public double sigma;
    public double theta;
    public double viscosity;
    public double quadratic_dragging;
    public double friction;
    public double max_limited_eta_prop;
    public double latitude;
    public double earth_radius;
    public double earth_rotation;
    public double height_smooth;
    public double regularizer_shell;
    public double capillary_accel_multiplier;
    public double dripping_radius_multiplier;
    public double bulk_threshold_multiplier;
    public double absorptionRate;
    public double airviscosity;
    public double drag_radius_multiplier;
    public double hair_hair_cohesion_multiplier;
    public double hair_solid_cohesion_multiplier;
    public double radius_multiplier;
    public double collision_stiffness;
    public double radius_multiplier_planar;
    public double collision_stiffness_planar;
    public double damping_multiplier;
    public double damping_multiplier_planar;
    public double friction_multiplier_planar;

    public int hairsteps;
    public int swesteps;
    public int fluidcorrectionsteps;

    public bool no_fluids;
    public bool no_swe;
    public bool no_absorb;
    public bool no_dripping;
    public bool no_fictitious;

    public bool drippingnear;
    public bool drippingfar;
    public bool drippingmiddle;
    public bool use_ctcd;
    public bool apply_coriolis;
    public bool global_volume_control;
    public bool individual_transfer;
    public bool volume_summary;
    public bool viscous_solve;

    public MASS_UPDATE_MODE mass_update_mode;
    public Vectors gravity;

    double max_velocity_ratio;

    public static double M_PI = 3.14159265358979323846;
    public const int MUM_MOMENTUM = 3;

    public WetHairParameter()
    { 
        dt = 0.004;
        rho = 1.0;
        sigma = 72.0;
        theta = M_PI / 3;
        viscosity = 8.9e-3;
        friction = 0.36;
        max_limited_eta_prop = 6.0;
        latitude = 0.71226395812;
        earth_radius = 6.37662216306e8;
        earth_rotation = 1.160576e-5 * 2.0 * M_PI;
        height_smooth = 1.0;
        regularizer_shell = 0.0;
        capillary_accel_multiplier = 0.0;
        dripping_radius_multiplier = 0.5;
        bulk_threshold_multiplier = 0.5;
        absorptionRate = 0.001;
        airviscosity = 0.0;
        drag_radius_multiplier = 1.0;
        quadratic_dragging = 1.0;
        hair_hair_cohesion_multiplier = 1.0;
        hair_solid_cohesion_multiplier = 1.0;
        radius_multiplier = 1.6;
        collision_stiffness = 10000.0;
        radius_multiplier_planar = 1.1;
        collision_stiffness_planar = 10000.0;
        damping_multiplier = 0.0;
        damping_multiplier_planar = 0.0;
        friction_multiplier_planar = 0.0;
        hairsteps = 1;
        swesteps = 1;
        fluidcorrectionsteps = 8;
        drippingnear = true;
        drippingfar = true;
        drippingmiddle = true;
        no_fluids = false;
        no_swe = false;
        no_absorb = false;
        no_dripping = false;
        no_fictitious = false;
        use_ctcd = false;
        global_volume_control = true;
        individual_transfer = true;
        volume_summary = false;
        viscous_solve = false;
        apply_coriolis = false;
        mass_update_mode = (MASS_UPDATE_MODE)MUM_MOMENTUM;
        gravity = new Vectors(3);
        gravity[0] = 0.0;
        gravity[1] = -981.0;
        gravity[2] = 0.0;
        max_velocity_ratio = 100.0;
    }
};

public class TwoDScene
{
    private Vectors m_base_x;
    private Vectors m_x;
    private Vectors m_v;
    private Vectors m_m;
    private Vectors m_rest_m = new Vectors();
    private Vectors m_radii;
    private Vectors m_interpolated_m;

    private Vectors m_fluid_drag_buffer = new Vectors();

    private List<int> m_particle_to_dofs = new List<int>();
    private VectorXi m_dofs_to_component = new VectorXi(0);
    private List<bool> m_is_strand_tip;
    private VectorXi m_dofs_to_particle;
    private List<StrandForce> m_strands = new List<StrandForce>();
    private int m_num_strands;
    private List<int> m_edge_to_hair;

    private List<double> m_liquid_on_hair;
    private List<double> m_liquid_free;

    private double m_search_radius;
    private double m_collision_keep_proximity;

    private Vectors m_bb_min;
    private Vectors m_bb_max;

    private FluidSim m_fluid_sim;

    private PolygonalCohesion m_polygonal_cohesion;

    private List<bool> m_fixed;

    private List<Tuple<int, int>> m_edges;
    private Vectors m_edge_radii;
    private Vectors m_edge_rest_radii = new Vectors();
    private Vectors m_edge_rest_length = new Vectors();
    private Vectors m_edge_poisson_ratio = new Vectors();

    private List<Force> m_forces;
    private List<Force> m_external_forces = new List<Force>();
    private List<Force> m_internal_forces = new List<Force>();
    private List<List<Force>> m_hair_internal_forces = new List<List<Force>>(); // hair-> forces only affecting that hair
    private List<Force> m_inter_hair_forces = new List<Force>(); // String 'tags' assinged to  particles. Can be used to identify and single out
    private List<String> m_particle_tags;

    private List<int> m_particle_to_hairs;

    private List<int> m_script_group;

    private List<Vectors> m_scripted_translate;

    private List<Vectors> m_scripted_rotation;

    private List<int> m_particle_to_hair_local_indices;

    private List<HairFlow> m_flows = new List<HairFlow>();

    private List<List<int>> m_particle_to_edge = new List<List<int>>();

    private List<HashSet<int>> m_bp_edge_edge_pairs = new List<HashSet<int>>();
    private List<HashSet<int>> m_bp_particle_edge_pairs = new List<HashSet<int>>();
    private List<HashSet<int>> m_bp_particle_particle_pairs = new List<HashSet<int>>();

    private List<StrandParameters> m_strandParameters = new List<StrandParameters>();
    private List<StrandEquilibriumParameters> m_strandEquilibriumParameters = new List<StrandEquilibriumParameters>();

    private bool m_massSpringSim;

    private double m_volume_hair;
    private double m_volume_particle;
    private double m_volume_particle_insterted;
    private double m_volume_particle_removed;
    private double m_volume_reservoir;

    private double m_volume_hair_old;
    private double m_volume_particle_old;
    private double m_volume_reservoir_old;

    private VectorXi m_constraint_idx;

    private WetHairParameter m_parameters;

    public int DIM = 3;

    public TwoDScene(in bool isMassSpring)
    {
        m_x = new Vectors();
        m_base_x = new Vectors();
        m_v = new Vectors();
        m_m = new Vectors();
        m_interpolated_m = new Vectors();
        m_fixed = new List<bool>();
        m_radii = new Vectors();
        m_edges = new List<Tuple<int, int>>();
        m_edge_radii = new Vectors();
        m_forces = new List<Force>();
        m_particle_tags = new List<string>();
        m_massSpringSim = isMassSpring;
        m_fluid_sim = null;
        m_polygonal_cohesion = null;
        m_massSpringSim = isMassSpring;
        m_parameters = new WetHairParameter();
    }

    public TwoDScene(in TwoDScene otherscene, in bool isMassSpring)
    {
        m_x = otherscene.m_x;
        m_base_x = otherscene.m_base_x;
        m_v = otherscene.m_v;
        m_m = otherscene.m_m;
        m_interpolated_m = otherscene.m_interpolated_m;
        m_rest_m = otherscene.m_rest_m;
        m_fixed = otherscene.m_fixed;
        m_radii = new Vectors();
        m_edges = new List<Tuple<int, int>>();
        m_edge_radii = new Vectors();
        m_forces = new List<Force>(otherscene.m_forces.Count);
        m_particle_tags = new List<String>();
        m_fluid_sim = null;
        m_polygonal_cohesion = null;
        m_massSpringSim = isMassSpring;
        for (int i = 0; i < m_forces.Count; i++)
        {
            m_forces[i] = otherscene.m_forces[i].Clone();
        }
    }

    public StrandParameters getStrandParameters(int index)
    {
        return m_strandParameters[index];
    }

    public bool isTip(int particle)
    {
        if (!m_massSpringSim)
        {
            return m_is_strand_tip[particle];
        }
        return false;
    }

    public Vectors getX()
    {
        return m_x;
    }

    public Vectors getV()
    {
        return m_v;
    }

    public int getDof(int particle)
    {
        if (m_massSpringSim)
        {
            return DIM * particle;
        }
        else
        {
            return m_particle_to_dofs[particle];
        }
    }

    public int getComponent(int dof)
    {
        if (m_massSpringSim)
        {
            return dof % DIM;
        }
        else
        {
            return m_dofs_to_component[dof];
        }
    }

    public int getVertFromDof(int dof)
    {
        if (m_massSpringSim)
        {
            return dof / DIM;
        }
        else
        {
            return m_dofs_to_particle[dof];
        }
    }

    public Vectors getRadii()
    {
        return m_radii;
    }

    public double getLiquidDensity()
    {
        return m_parameters.rho;
    }

    public void insertForce(ref Force newForce)
    {
        m_forces.Add(newForce);
    }

    public void insertStrandParameters(ref StrandParameters newparams)
    {
        m_strandParameters.Add(newparams);
    }

    public void insertStrandEquilibriumParameters(ref StrandEquilibriumParameters newparams)
    {
        m_strandEquilibriumParameters.Add(newparams);
    }

    public void resizeSystem(int num_particles, int num_edges, int num_strands)
    {
        int numDofs = (4 * num_particles) - num_strands;
        m_fluid_drag_buffer.resize(numDofs);
        m_x.resize(numDofs);
        m_base_x.resize(numDofs);
        m_v.resize(numDofs);
        m_m.resize(numDofs);
        m_interpolated_m.resize(numDofs);
        m_rest_m.resize(numDofs);

        m_script_group = new List<int>(new int[num_particles]);
        m_particle_to_dofs = new List<int>(new int[num_particles]);
        m_dofs_to_component.resize(numDofs);
        m_is_strand_tip = new List<bool>(new bool[num_particles]);

        m_fixed = new List<bool>(new bool[num_particles]);
        m_radii.resize(num_particles);
        m_particle_tags = new List<string>(new string[num_particles]);
        m_edges = new List<Tuple<int, int>>(new Tuple<int, int>[num_edges]);
        m_edge_radii.resize(num_edges);
        m_edge_rest_radii.resize(num_edges);
        m_edge_rest_length.resize(num_edges);
        m_edge_poisson_ratio.resize(num_edges);
        m_particle_to_edge = new List<List<int>>(num_particles);
        for (int i = 0; i < num_particles; i++)
        {
            m_particle_to_edge.Add(new List<int>());
        }

        m_num_strands = num_strands;
        m_strands = new List<StrandForce>(new StrandForce[num_particles]);
        m_edge_to_hair = new List<int>(new int[num_edges]);
    }

    public void setVertToDoFMap(in List<int> vert_to_dof,
        in VectorXi dofs_to_vars,
        in List<bool> tipVerts,
        in VectorXi dof_to_vert)
    {
        m_particle_to_dofs = new List<int>();
        foreach (int i in vert_to_dof)
        {
            m_particle_to_dofs.Add(i);
        }
        m_dofs_to_component = dofs_to_vars;
        m_is_strand_tip = new List<bool>();
        foreach (bool b in tipVerts)
        {
            m_is_strand_tip.Add(b);
        }
        m_dofs_to_particle = dof_to_vert;
    }

    public List<string> getParticleTags()
    {
        return m_particle_tags;
    }

    public void setPosition(int particle, in Vectors pos)
    {
        m_x.SetSegment(getDof(particle), DIM, pos);
    }

    public Vectors getPosition(int particle)
    {
        return m_x.segment(getDof(particle), DIM);
    }

    public void setVelocity(int particle, in Vectors vel)
    {
        m_v.SetSegment(getDof(particle), DIM, vel);
    }

    public Vectors getVelocity(int particle)
    {
        return m_v.segment(getDof(particle), DIM);
    }

    public void setMass(int particle, in double mass)
    {
        if (m_massSpringSim)
        {
            m_m.SetSegment(particle * DIM, DIM, m_m.segment(particle * DIM, DIM).setConstant(mass));
            m_rest_m.SetSegment(particle * DIM, DIM, m_rest_m.segment(particle * DIM, DIM).setConstant(mass));
            m_interpolated_m.SetSegment(particle * DIM, DIM, m_interpolated_m.segment(particle * DIM, DIM).setConstant(mass));
        }
    }

    public void setFixed(int particle, bool fixe)
    {
        m_fixed[particle] = fixe;
    }

    public void setEdge(int idx, Tuple<int, int> edge,
        double radius)
    {

        m_edges[idx] = edge;
        m_particle_to_edge[edge.Item1].Add(idx);
        m_particle_to_edge[edge.Item2].Add(idx);

        if (m_massSpringSim)
        {
            m_edge_radii[idx] = radius;
            m_edge_rest_radii[idx] = radius;
        }
    }

    public void setEdgeRestLength(int idx, in double l0)
    {
        m_edge_rest_length[idx] = l0;
    }

    public void setScriptedGroup(int particle, int group_idx)
    {
        m_script_group[particle] = group_idx;
    }

    public void computeMassesAndRadiiFromStrands()
    {
        m_strands = new List<StrandForce>(new StrandForce[m_num_strands]);

        int nStrand = 0;
        for (int f = 0; f < m_forces.Count(); ++f)
        {
            StrandForce strand = (StrandForce)m_forces[f];
            if (strand == null)
                continue;
            m_strands[nStrand] = strand;

            if (strand.m_strandParams.m_straightHairs != 1.0)
            {
                List<Vectors> kappas = strand.alterRestKappas();
                for (int k = 0; k < kappas.Count(); ++k)
                {
                    kappas[k] = kappas[k] * strand.m_strandParams.m_straightHairs;
                }
            }

            for (int v = 0; v < strand.getNumVertices(); ++v)
            {
                int globalVtx = strand.m_verts[v];
                int globalEdx = globalVtx - nStrand;
                int globalDof = getDof(globalVtx);
                double r =
                    strand.m_strandParams.getRadius(v, strand.getNumVertices());

                m_radii[globalVtx] = r;
                m_m.SetSegment(globalDof, DIM, m_m.segment(globalDof, DIM).setConstant(strand.m_vertexMasses[v]));

                if (v < strand.getNumEdges())
                {
                    m_edge_to_hair[globalEdx] = nStrand;
                    // Edge radius, edge's should be indexed the same as
                    m_edge_radii[globalEdx] = r;
                    m_edge_rest_radii[globalEdx] = r;

                    // Twist Mass (Second moment of inertia * length)
                    double mass = strand.m_strandParams.m_density * CMath.M_PI * r * r *
                                        strand.m_restLengths[v];
                    double vtm = 0.25 * mass * 2 * r * r;
                    m_m[globalDof + 3] = vtm;
                }
            }
            ++nStrand;
        }

        m_rest_m = m_m.Clone();
        m_interpolated_m = m_m.Clone();
    }

    public List<HairFlow> getFilmFlows()
    {
        return m_flows;
    }


    public List<Tuple<int, int>> getEdges()
    {
        return m_edges;
    }

    public double getDt()
    {
        return m_parameters.dt;
    }

    public double getHairSteps()
    {
        return m_parameters.hairsteps;
    }

    public WetHairParameter getParameter()
    {
        return m_parameters;
    }

    public int getNumParticles()
    {
        return (m_x.size() + m_num_strands) / 4;
    }

    public void updateStrandParamsTimestep(in double dt)
    {
        for (int p = 0; p < m_strandParameters.Count; ++p)
        {
            m_strandParameters[p].computeViscousForceCoefficients(dt);
        }
    }

    public void updateStrandStartStates()
    {
        for (int s = 0; s < m_num_strands; ++s)
        {
            m_strands[s].updateStartDoFs(m_x);
        }
        /*
        if (m_polygonal_cohesion)
        {
            m_polygonal_cohesion->updateViscousStartPhi(m_x);
        }
        */
    }

    public void postCompute(in double dt)
    {
        int nforces = m_strands.Count;
        for (int i = 0; i < nforces; ++i)
        {
            m_forces[i].postStepScene(dt);
        }
    }

    public Tuple<int, int> getEdge(int edg)
    {
        return m_edges[edg];
    }

    public void applyScript(in double dt)
    {
        int np = getNumParticles();
        for (int i = 0; i < np; ++i)
        {
            if (!isFixed(i)) continue;

            //여기서 물체 전체 움직이는 거에 대한 fixed point 속도 조정
            /*
            int sg_idx = m_script_group[i];
            if (sg_idx < 0 || sg_idx >= m_scripted_translate.Count) continue;

            VectorXs q = m_scripted_rotation[sg_idx];
            VectorXs t = m_scripted_translate[sg_idx];

            VectorXs x0 = new VectorXs(DIM + 1);
            x0.SetSegment(0, DIM + 1, m_base_x.segment(getDof(i), DIM + 1));
            VectorXs xstar = new VectorXs(DIM + 1);
            xstar.SetSegment(0, DIM + 1, m_x.segment(getDof(i), DIM + 1));
            VectorXs p0 = new VectorXs(4);
            p0[0] = x0[0];
            p0[1] = x0[1];
            p0[2] = x0[2];
            p0[3] = 0.0f;
            VectorXs trans_x0 = (q * p0 * q.qinverse()) + t;
            
            m_v.SetSegment(getDof(i), DIM + 1, (trans_x0 - xstar) / dt);
            */

        }

        int nstrand = m_strandEquilibriumParameters.Count;
        for (int i = 0; i < nstrand; ++i)
        {
            if (!m_strandEquilibriumParameters[i].m_valid ||
                !m_strandEquilibriumParameters[i].m_dirty)
            {
                continue;
            }

            updateCurlyHair(m_strandEquilibriumParameters[i].m_dL,
                            ref m_strandEquilibriumParameters[i].m_vertices,
                            m_strandEquilibriumParameters[i].m_curl_radius,
                            m_strandEquilibriumParameters[i].m_curl_density,
                            m_strandEquilibriumParameters[i].m_root_length);
            int nverts = m_strandEquilibriumParameters[i].m_vertices.Count;
            Vectors dof_restshape = new Vectors(nverts * 4 - 1);
            dof_restshape.setZero();

            for (int j = 0; j < nverts; ++j)
            {
                dof_restshape.SetSegment(j * 4, 3,
                    m_strandEquilibriumParameters[i].m_vertices[j]);
            }
            m_strands[i].updateRestShape(dof_restshape, 0);
            m_strandEquilibriumParameters[i].m_dirty = false;
        }
    }

    public void initializeScriptedGroup()
    {
        m_scripted_translate = new List<Vectors>(0);

        applyScript(0);
    }

    public bool isFixed(int particle)
    {
        return m_fixed[particle];
    }

    public void updateCurlyHair(in double dL, ref List<Vectors> vertices,
                     double curl_radius, double curl_density,
                     double root_length)
    {
        int nv = vertices.Count;
        if (nv < 2)
            return;

        // generate an orthonormal frame
        Vectors initnorm = (vertices[1] - vertices[0]).normalized();
        vertices[1] = vertices[0] + initnorm * root_length;

        Vectors p1;
        p1 = CMath.findNormal(initnorm);
        Vectors p2 = p1.cross(initnorm);

        double xa = CMath.M_PI / (curl_density * 4);  // 0 // start curve parameter
        double xb = 0;                          // end curve parameter

        Vectors freepoint = vertices[1];

        for (int j = 2; j < nv; ++j)
        {
            xb = (dL + xa +
                  curl_radius * curl_radius * curl_density * curl_density * xa) /
                 (1 + curl_radius * curl_radius * curl_density *
                          curl_density);  // upate to get length dL along curve
            vertices[j] = freepoint + xb * initnorm +
                          curl_radius * Math.Cos(xb * curl_density) * p2 +
                          curl_radius * Math.Sin(xb * curl_density) * p1;
            xa = xb;  // next...
        }
    }

    public void updateHairConnectivity()
    {
        int nf = m_flows.Count;

        m_particle_to_hairs = new List<int>(getNumParticles());
        for (int i = 0; i < getNumParticles(); i++)
        {
            m_particle_to_hairs.Add(-1);
        }

        m_particle_to_hair_local_indices = new List<int>(getNumParticles());
        for (int i = 0; i < getNumParticles(); i++)
        {
            m_particle_to_hair_local_indices.Add(-1);
        }

        for (int i = 0; i < nf; i++)
        {
            HairFlow flow = m_flows[i];
            var indices = flow.getParticleIndices();
            int nfp = indices.Count;
            for (int j = 0; j < nfp; j++)
            {
                int pidx = indices[j];
                m_particle_to_hairs[pidx] = i;
                m_particle_to_hair_local_indices[pidx] = j;
            }
        }
    }

    public bool isMassSpring()
    {
        return m_massSpringSim;
    }

    public int getNumFlows()
    {
        return m_flows.Count;
    }

    public int getNumDofs()
    {
        return m_x.Size;
    }

    public void preComputeLocal(in Vectors dx, in Vectors dv, double dt)
    {
        int nfl = m_hair_internal_forces.Count;

        if (dx.Size == 0)
        {
            Parallel.For(0, nfl, (i) => {
                foreach (Force f in m_hair_internal_forces[i])
                {
                    if (!f.isPrecomputationParallelized())
                    {
                        f.preCompute(m_x, m_v, m_interpolated_m, dt);
                    }
                }
            });

            for (int i = 0; i < nfl; ++i)
            {
                foreach (Force f in m_hair_internal_forces[i])
                {
                    if (f.isPrecomputationParallelized())
                        f.preCompute(m_x, m_v, m_interpolated_m, dt);
                }
            }

            foreach (Force f in m_external_forces)
            {
                f.preCompute(m_x, m_v, m_interpolated_m, dt);
            }
        }
        else
        {
            Vectors nx = m_x + dx;
            Vectors nv = m_v + dv;
            Parallel.For(0, nfl, (i) => {
                foreach (Force f in m_hair_internal_forces[i])
                {
                    if (!f.isPrecomputationParallelized())
                    {
                        f.preCompute(nx, nv, m_interpolated_m, dt);
                    }
                }
            });

            for (int i = 0; i < nfl; ++i)
            {
                foreach (Force f in m_hair_internal_forces[i])
                {
                    if (f.isPrecomputationParallelized())
                        f.preCompute(nx, nv, m_interpolated_m, dt);
                }
            }

            foreach (Force f in m_external_forces)
            {
                f.preCompute(nx, nv, m_interpolated_m, dt);
            }
        }
    }

    public void updateNumConstraintsLocal(ref int num_constraint_pos_,
        ref int num_constraint_vel_,
        ref int num_J_, ref int num_Jv_,
        ref int num_Jxv_,
        ref int num_tildeK_)
    {
        m_constraint_idx = new VectorXi(6);

        int nfl = getNumFlows();
        for (int i = 0; i < nfl; ++i)
        {
            List<Force> hair_forces = m_hair_internal_forces[i];
            int nhair_forces = hair_forces.Count;

            VectorXi constraint_start = m_constraint_idx.Clone();

            for (int j = 0; j < nhair_forces; ++j)
            {
                int num_pos = hair_forces[j].numConstraintPos();
                int num_vel = hair_forces[j].numConstraintVel();
                int num_J = hair_forces[j].numJ();
                int num_Jv = hair_forces[j].numJv();
                int num_Jxv = hair_forces[j].numJxv();
                int num_TildeK = hair_forces[j].numTildeK();

                hair_forces[j].setInternalIndex(
                    m_constraint_idx[0], m_constraint_idx[1], m_constraint_idx[2],
                    m_constraint_idx[3], m_constraint_idx[4], m_constraint_idx[5]);
                m_constraint_idx[0] += num_pos;
                m_constraint_idx[1] += num_vel;
                m_constraint_idx[2] += num_J;
                m_constraint_idx[3] += num_Jv;
                m_constraint_idx[4] += num_Jxv;
                m_constraint_idx[5] += num_TildeK;
            }

            VectorXi num_constraints = m_constraint_idx - constraint_start;

            m_flows[i].setConstraintParameters(constraint_start, num_constraints);
        }

        num_constraint_pos_ = m_constraint_idx[0];
        num_constraint_vel_ = m_constraint_idx[1];
        num_J_ = m_constraint_idx[2];
        num_Jv_ = m_constraint_idx[3];
        num_Jxv_ = m_constraint_idx[4];
        num_tildeK_ = m_constraint_idx[5];
    }

    [BurstCompile]
    public void localPostPreprocess(ref Vectors lambda, ref Vectors lambda_v,
        ref TripletXs J, ref TripletXs Jv,
        ref TripletXs Jxv, ref TripletXs tildeK,
        ref TripletXs stiffness,
        ref TripletXs damping, ref Vectors Phi,
        ref Vectors Phiv, in Vectors dx,
        in Vectors dv, in double dt)
    {
        int nf = m_hair_internal_forces.Count;

        Vectors local_lambda = lambda;
        Vectors local_lambda_v = lambda_v;
        TripletXs local_J = J;
        TripletXs local_Jv = Jv;
        TripletXs local_Jxv = Jxv;
        TripletXs local_tildeK = tildeK;
        TripletXs local_stiffness = stiffness;
        TripletXs local_damping = damping;
        Vectors local_Phi = Phi;
        Vectors local_Phiv = Phiv;
        Vectors local_dx = dx;
        Vectors local_dv = dv;
        double local_dt = dt;

        if (dx.size() == 0)
        {
            int batchSize = 4;
            int numBatches = (nf + batchSize - 1) / batchSize;

            Parallel.For(0, numBatches, batchIndex => {
                int start = batchIndex * batchSize;
                int end = Math.Min(start + batchSize, nf);

                for (int i = start; i < end; i++)
                {

                    foreach (Force f in m_hair_internal_forces[i])
                    {
                        if (!f.isParallelized())
                        {
                            f.computeIntegrationVars(m_x, m_v, m_interpolated_m, ref local_lambda,
                                                     ref local_lambda_v, ref local_J, ref local_Jv, ref local_Jxv, ref local_tildeK,
                                                     ref local_stiffness, ref local_damping, ref local_Phi, ref local_Phiv, local_dt);
                        }
                    }
                }
            });

            for (int i = 0; i < nf; ++i)
            {
                foreach (Force f in m_hair_internal_forces[i])
                {
                    if (f.isParallelized())
                        f.computeIntegrationVars(m_x, m_v, m_interpolated_m, ref local_lambda,
                                                  ref local_lambda_v, ref local_J, ref local_Jv, ref local_Jxv, ref local_tildeK, ref local_stiffness,
                                                  ref local_damping, ref local_Phi, ref local_Phiv, dt);
                }
            }
        }
        else
        {
            Vectors nx = m_x + dx;
            Vectors nv = m_v + dv;

            int batchSize = 4;
            int numBatches = (nf + batchSize - 1) / batchSize;

            Parallel.For(0, numBatches, batchIndex => {
                int start = batchIndex * batchSize;
                int end = Math.Min(start + batchSize, nf);

                for (int i = start; i < end; i++)
                {

                    foreach (Force f in m_hair_internal_forces[i])
                    {
                        if (!f.isParallelized())
                        {
                            f.computeIntegrationVars(nx, nv, m_interpolated_m, ref local_lambda,
                                                     ref local_lambda_v, ref local_J, ref local_Jv, ref local_Jxv, ref local_tildeK,
                                                     ref local_stiffness, ref local_damping, ref local_Phi, ref local_Phiv, local_dt);
                        }
                    }
                }
            });

            for (int i = 0; i < nf; ++i)
            {
                foreach (Force f in m_hair_internal_forces[i])
                {
                    if (f.isParallelized())
                        f.computeIntegrationVars(nx, nv, m_interpolated_m, ref local_lambda,
                                                  ref local_lambda_v, ref local_J, ref local_Jv, ref local_Jxv, ref local_tildeK, ref local_stiffness,
                                                  ref local_damping, ref local_Phi, ref local_Phiv, dt);
                }
            }
        }

        lambda = local_lambda;
        lambda_v = local_lambda_v;
        J = local_J;
        Jv = local_Jv;
        Jxv = local_Jxv;
        tildeK = local_tildeK;
        stiffness = local_stiffness;
        damping = local_damping;
        Phi = local_Phi;
    }

    public Vectors getInterpolatedM()
    {
        return m_interpolated_m;
    }

    public List<int> getParticleToHairLocalIndices()
    {
        return m_particle_to_hair_local_indices;
    }

    public List<int> getParticleToHairs()
    {
        return m_particle_to_hairs;
    }

    public void accumulateExternalGradU(ref Vectors F, in Vectors dx,
        in Vectors dv)
    {
        if (dx.size() == 0)
        {
            for (int i = 0; i < m_external_forces.Count; i++)
            {
                m_external_forces[i].addGradEToTotal(m_x, m_v, m_interpolated_m, ref F);
            }
        }
        else
        {
            for (int i = 0; i < m_external_forces.Count; i++)
            {
                m_external_forces[i].addGradEToTotal(m_x + dx, m_v + dv, m_interpolated_m, ref F);
            }
        }
    }

    public void storeLambda(in Vectors lambda, in Vectors lambda_v)
    {
        Vectors lambda_ = lambda;
        Vectors lambda_v_ = lambda_v;
        int nf = m_internal_forces.Count;
        Parallel.For(0, nf, i => { 
            m_internal_forces[i].storeLambda(lambda_, lambda_v_);
        });
    }

    public void preComputeInterhair(in Vectors dx, in Vectors dv, in double dt)
    {
        int nf = m_inter_hair_forces.Count;
        if (dx.size() == 0)
        {
            double local_dt = dt;
            Parallel.For(0, nf, i => {
                if (!m_inter_hair_forces[i].isPrecomputationParallelized())
                    m_inter_hair_forces[i].preCompute(m_x, m_v, m_interpolated_m, local_dt);
            });

            for (int i = 0; i < nf; ++i)
            {
                if (m_inter_hair_forces[i].isPrecomputationParallelized())
                    m_inter_hair_forces[i].preCompute(m_x, m_v, m_interpolated_m, dt);
            }
        }
        else
        {
            Vectors nx = m_x + dx;
            Vectors nv = m_v + dv;
            double local_dt = dt;
            Parallel.For(0, nf, i => {
                if (!m_inter_hair_forces[i].isPrecomputationParallelized())
                    m_inter_hair_forces[i].preCompute(nx, nv, m_interpolated_m, local_dt);
            });

            for (int i = 0; i < nf; ++i)
            {
                if (m_inter_hair_forces[i].isPrecomputationParallelized())
                    m_inter_hair_forces[i].preCompute(nx, nv, m_interpolated_m, dt);
            }
        }
    }

    public void updateNumConstraintsInterHair(
        ref int num_constraint_pos_, ref int num_constraint_vel_, ref int num_J_,
    ref int num_Jv_, ref int num_Jxv_, ref int num_tildeK_, ref VectorXi interhair_param,
    ref VectorXi interhair_num)
    {
        interhair_param = m_constraint_idx;

        for (int i = 0; i < m_inter_hair_forces.Count;
             ++i)
        {
            int num_pos = m_inter_hair_forces[i].numConstraintPos();
            int num_vel = m_inter_hair_forces[i].numConstraintVel();
            int num_J = m_inter_hair_forces[i].numJ();
            int num_Jv = m_inter_hair_forces[i].numJv();
            int num_Jxv = m_inter_hair_forces[i].numJxv();
            int num_TildeK = m_inter_hair_forces[i].numTildeK();

            m_inter_hair_forces[i].setInternalIndex(
                m_constraint_idx[0], m_constraint_idx[1], m_constraint_idx[2],
                m_constraint_idx[3], m_constraint_idx[4], m_constraint_idx[5]);
            m_constraint_idx[0] += num_pos;
            m_constraint_idx[1] += num_vel;
            m_constraint_idx[2] += num_J;
            m_constraint_idx[3] += num_Jv;
            m_constraint_idx[4] += num_Jxv;
            m_constraint_idx[5] += num_TildeK;
        }

        num_constraint_pos_ = m_constraint_idx[0];
        num_constraint_vel_ = m_constraint_idx[1];
        num_J_ = m_constraint_idx[2];
        num_Jv_ = m_constraint_idx[3];
        num_Jxv_ = m_constraint_idx[4];
        num_tildeK_ = m_constraint_idx[5];

        interhair_num = m_constraint_idx - interhair_param;
    }

    public void interhairPostPreprocess(
        Vectors lambda, Vectors lambda_v, TripletXs J, TripletXs Jv,
        TripletXs Jxv, TripletXs tildeK, TripletXs stiffness, TripletXs damping,
        Vectors Phi, Vectors Phiv, Vectors dx, Vectors dv,
        double dt)
    {
        int nf = m_inter_hair_forces.Count;
        if (dx.size() == 0)
        {
            // for unparallelized forces
            Parallel.For(0, nf, i => {
                if (!m_inter_hair_forces[i].isParallelized())
                    m_inter_hair_forces[i].computeIntegrationVars(
                        m_x, m_v, m_interpolated_m, ref lambda, ref lambda_v, ref J, ref Jv, ref Jxv, ref tildeK,
                        ref stiffness, ref damping, ref Phi, ref Phiv, dt);
            });

            // for parallelized forces
            for (int i = 0; i < nf; ++i)
            {
                if (m_inter_hair_forces[i].isParallelized())
                    m_inter_hair_forces[i].computeIntegrationVars(
                        m_x, m_v, m_interpolated_m, ref lambda, ref lambda_v, ref J, ref Jv, ref Jxv, ref tildeK,
                        ref stiffness, ref damping, ref Phi, ref Phiv, dt);
            }
        }
        else
        {
            Vectors nx = m_x + dx;
            Vectors nv = m_v + dv;
            // for unparallelized forces
            Parallel.For(0, nf, i => {
                if (!m_inter_hair_forces[i].isParallelized())
                    m_inter_hair_forces[i].computeIntegrationVars(
                        nx, nv, m_interpolated_m, ref lambda, ref lambda_v, ref J, ref Jv, ref Jxv, ref tildeK,
                        ref stiffness, ref damping, ref Phi, ref Phiv, dt);
            });

            // for parallelized forces
            for (int i = 0; i < nf; ++i)
            {
                if (m_inter_hair_forces[i].isParallelized())
                    m_inter_hair_forces[i].computeIntegrationVars(
                        nx, nv, m_interpolated_m, ref lambda, ref lambda_v, ref J, ref Jv, ref Jxv, ref tildeK,
                        ref stiffness, ref damping, ref Phi, ref Phiv, dt);
            }
        }
    }

    public void insertFilmFlow(HairFlow flow)
    {
        m_flows.Add(flow);
    }

    public Vectors getHairRestMass()
    {
        return m_rest_m;
    }

    public void categorizeForces()
    {
        m_external_forces.Clear();
        m_internal_forces.Clear();
        m_inter_hair_forces.Clear();
        m_hair_internal_forces = new List<List<Force>>(m_flows.Count);

        for (int i = 0; i < m_flows.Count; ++i)
        {
            m_hair_internal_forces.Add(new List<Force>(0));
        }

        for (int i = 0; i < m_forces.Count; ++i)
        {
            if (m_forces[i].isExternal())
            {
                m_external_forces.Add(m_forces[i]);
            }
            else
            {
                m_internal_forces.Add(m_forces[i]);

                if (!m_forces[i].isInterHair(new Vectors(0), new Vectors(0)))
                {
                    int hidx = m_forces[i].getAffectedHair(m_particle_to_hairs);
                    if (hidx < 0)
                        continue;
                    m_hair_internal_forces[hidx].Add(m_forces[i]);
                }
                else
                {
                    m_inter_hair_forces.Add(m_forces[i]);
                }
            }
        }
    }

    public void addPositionValue(int startPartNum, int PartNum, Vectors v)
    {
        Parallel.For(startPartNum, PartNum, (i) =>
        {
            m_x.SetSegment(getDof(i), 3, m_x.segment(getDof(i), 3) + v);
        });
    }
}
