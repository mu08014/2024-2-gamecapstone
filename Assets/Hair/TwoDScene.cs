using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

struct Script
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
    VectorXs v;
    VectorXs origin;
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

enum MASS_UPDATE_MODE
{
    MUM_NONE,
    MUM_MASS_ONLY,
    MUM_DIRECT_DIV,
    MUM_MOMENTUM,
    MUM_COUNT
}

class WetHairParameter
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
    public VectorXs gravity;

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
        gravity = new VectorXs(3);
        gravity[0] = 0.0;
        gravity[1] = -981.0;
        gravity[2] = 0.0;
        max_velocity_ratio = 100.0;
    }
};

public class TwoDScene
{
    private VectorXs m_base_x;
    private VectorXs m_x;
    private VectorXs m_v;
    private VectorXs m_m;
    private VectorXs m_rest_m;
    private VectorXs m_radii;
    private VectorXs m_interpolated_m;

    private VectorXs m_fluid_drag_buffer;

    private List<int> m_particle_to_dofs;
    private VectorXi m_dofs_to_component;
    private List<bool> m_is_strand_tip;
    private VectorXi m_dofs_to_particle;
    private List<StrandForce> m_strands;
    private int m_num_strands;
    private List<int> m_edge_to_hair;

    private List<double> m_liquid_on_hair;
    private List<double> m_liquid_free;

    private double m_search_radius;
    private double m_collision_keep_proximity;

    private Vectors m_bb_min;
    private Vectors m_bb_max;

    private bool m_massSpringSim;

    private WetHairParameter m_parameters;

    public int DIM;

    public TwoDScene(in bool isMassSpring)
    {
        m_massSpringSim = isMassSpring;
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

    public VectorXs getRadii()
    {
        return m_radii;
    }

    public double getLiquidDensity()
    {
        return m_parameters.rho;
    }
}
