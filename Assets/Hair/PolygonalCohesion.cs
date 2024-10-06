using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Threading;
using UnityEngine;

public class PolygonalRegion
{
    public int hair_idx_0;
    public int hair_idx_1;
    public List<int> local_indices_0;
    public List<int> local_indices_1;
}

public class EdgeEdgePair
{
    public int base_eidx;
    public int neighbor_hair_idx;
    public double alpha_0;
    public double alpha_1;
    public double neighbor_local_coord_0;
    public double neighbor_local_coord_1;
    public double count_0;
    public double count_1;
}

public class EdgeEdgePairEEC
{
    public int base_eidx;
    public int neighbor_eidx;
    public double alpha_0;
    public double alpha_1;
    public double neighbor_local_coord_0;
    public double neighbor_local_coord_1;

    public double time;
    public double alpha_contack;
    public double neighbor_local_coord_contact;
    public VectorXs avgpos;

    public double count_0;
    public double count_1;

    public bool valid;
    public bool updated;
}

public class ParticleEdgePair
{
    public int pidx;
    public int eidx;
    public double alpha;
    public double dist;
    public double radii_j;
    public double max_dist;
    public double count;
    public bool valid;
    public bool updated;
    public bool should_be_deleted;
    public bool latest;
}

public class PointEdgePair
{
    public int base_eidx;
    public double alpha_point;
    public int neighbor_eidx;
    public double neighbor_alpha_point;
    public double V;
    public double quadrature_weight;
    public double pressure_weight;
    public ulong hash_code;

    public double time;
}

public class ParticleParticlePair
{
    public int[] pidx = new int[2];
    public double d;
    public double r;
}

public class PolygonalCohesion : Force
{
    private const int m_meddling_stencil = 2;
    private const int m_num_quadrature = 1;

    private TwoDScene m_parent;
    private List<Dictionary<int, ParticleEdgePair>> m_adjacency_categorized;
    private List<Dictionary<int, EdgeEdgePairEEC>> m_edge_connections;
    private List<int> m_num_valid_edge_connections;
    private List<int> m_num_edge_connections;
    private List<List<int>> m_particle_to_ppairs;
    private List<int> m_num_adjacency_categorized;
    private List<int> m_counting_valid_adjacency;
    private List<int> m_counting_pp_pair_location;
    private List<ulong> m_counting_pp_pairs;
    private List<ParticleParticlePair> m_particle_particle_pairs;
    private List<PointEdgePair> m_point_edge_pairs;

    private List<List<PointEdgePair>> m_point_edge_pairs_cache;
    private List<int> m_counting_poe_pair_location;

    private VectorXs m_particle_adjacency_hair_size;

    private Dictionary<int, HashSet<int>> m_adjacency_hair_edges_buffer;

    private List<List<int>> m_particle_to_point_edge_pairs;

    private List<HashSet<int>> m_pp_pair_hash;

    private VectorXs m_particle_length;

    private Sorter m_sorter;

    private MatrixXs m_edge_buffer;

    private VectorXs m_gradE;
    private MatrixXs m_hessE;
    private MatrixXs m_hessV;
    private VectorXs m_pair_counts;

    private TripletXs m_hess_buffer;
    private bool m_use_decoupled_force;
    private bool m_compute_particle_poe_mapping;

    private SparseXs m_W_fv_interhair_T;
    private SparseXs m_gradF_global;
    private SparseXs m_gradF_global_T;
    private SparseXs m_gradF_interhair;
    private SparseXs m_gradF_interhair_T;
    private SparseXs m_dir_f_interhair;
    private SparseXs m_dir_f_interhair_T;

    private VectorXs m_G_f_global;
    private VectorXs m_iG_v_global;
    private VectorXs m_G_f_interhair;

    private VectorXs m_u_interhair;
    private VectorXs m_divu_interhair;

    private VectorXs m_area_v_global;
    private VectorXs m_area_e_global;
    private VectorXs m_cur_eta_v_global;
    private VectorXs m_pressure_v_global;

    private MatrixXs m_rhs_offset_v_global;
    private MatrixXs m_cur_hhrr_v_global;

    private VectorXs m_ce_inter_buffer;
    private VectorXs m_ce_inter_short_buffer;
    private VectorXs m_ce_global_buffer;
    private VectorXs m_cv_buffer;

    private List<int> m_pplink_count;

    private Dictionary<ulong, double> m_poep_lambda;

    private Mutex m_pep_mutex;
    private Mutex m_eep_mutex;
    private Mutex m_poep_mutex;

    private CohesionTable m_min_cohesion_table;
    private CohesionTable m_max_cohesion_table;

    Dictionary<ulong, double> m_viscous_start_phi;

    public PolygonalCohesion(ref TwoDScene scene)
    {
        m_parent = scene;
        m_sorter = null;
        m_use_decoupled_force = false;
        m_compute_particle_poe_mapping = true;
    }
}
