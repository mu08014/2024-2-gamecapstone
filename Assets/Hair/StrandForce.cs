using System;
using System.Collections.Generic;

public struct StrandState
{
    public DOFs m_dofs;
    public Edges m_edges;
    public Lengths m_lengths;
    public Tangents m_tangents;
    public ReferenceFrames1 m_referenceFrames1;
    public ReferenceFrames2 m_referenceFrames2;
    public ReferenceTwists m_referenceTwists;
    public Twists m_twists;
    public CurvatureBinormals m_curvatureBinormals;
    public TrigThetas m_trigThetas;
    public MaterialFrames1 m_materialFrames1;
    public MaterialFrames2 m_materialFrames2;
    public Kappas m_kappas;
    public GradKappas m_gradKappas;
    public GradTwists m_gradTwists;
    public GradTwistsSquared m_gradTwistsSquared;
    public HessKappas m_hessKappas;
    public HessTwists m_hessTwists;
    public BendingProducts m_bendingProducts;

    public StrandState(in Vectors initDoFs,
        ref BendingMatrixBase bendingMatrixBase)
    {
        m_dofs = new DOFs(initDoFs);
        m_edges = new Edges(ref m_dofs);
        m_lengths = new Lengths(ref m_edges);
        m_tangents = new Tangents(ref m_edges, ref m_lengths);
        m_referenceFrames1 = new ReferenceFrames1(ref m_tangents);
        m_referenceFrames2 = new ReferenceFrames2(ref m_tangents, ref m_referenceFrames1);
        m_referenceTwists = new ReferenceTwists(ref m_tangents, ref m_referenceFrames1);
        m_twists = new Twists(ref m_referenceTwists, ref m_dofs);
        m_curvatureBinormals = new CurvatureBinormals(ref m_tangents);
        m_trigThetas = new TrigThetas(ref m_dofs);
        m_materialFrames1 = new MaterialFrames1(ref m_trigThetas, ref m_referenceFrames1, ref m_referenceFrames2);
        m_materialFrames2 = new MaterialFrames2(ref m_trigThetas, ref m_referenceFrames1, ref m_referenceFrames2);
        m_kappas = new Kappas(ref m_curvatureBinormals, ref m_materialFrames1, ref m_materialFrames2);
        m_gradKappas = new GradKappas(ref m_lengths, ref m_tangents, ref m_curvatureBinormals,
                   ref m_materialFrames1, ref m_materialFrames2, ref m_kappas);
        m_gradTwists = new GradTwists(ref m_lengths, ref m_curvatureBinormals);
        m_gradTwistsSquared = new GradTwistsSquared(ref m_gradTwists);
        m_hessKappas = new HessKappas(ref m_lengths, ref m_tangents, ref m_curvatureBinormals,
                   ref m_materialFrames1, ref m_materialFrames2, ref m_kappas);
        m_hessTwists = new HessTwists(ref m_tangents, ref m_lengths, ref m_curvatureBinormals);
        m_bendingProducts = new BendingProducts(ref bendingMatrixBase, ref m_gradKappas);
    }
}

public struct StartState
{
    public DOFs m_dofs;
    public Edges m_edges;
    public Lengths m_lengths;
    public Tangents m_tangents;
    public ReferenceFrames1 m_referenceFrames1;
    public ReferenceFrames2 m_referenceFrames2;
    public ReferenceTwists m_referenceTwists;
    public Twists m_twists;
    public CurvatureBinormals m_curvatureBinormals;
    public TrigThetas m_trigThetas;
    public MaterialFrames1 m_materialFrames1;
    public MaterialFrames2 m_materialFrames2;
    public Kappas m_kappas;

    public StartState(in Vectors initDoFs)
    {
        m_dofs = new DOFs(initDoFs);
        m_edges = new Edges(ref m_dofs);
        m_lengths = new Lengths(ref m_edges);
        m_tangents = new Tangents(ref m_edges, ref m_lengths);
        m_referenceFrames1 = new ReferenceFrames1(ref m_tangents);
        m_referenceFrames2 = new ReferenceFrames2(ref m_tangents, ref m_referenceFrames1);
        m_referenceTwists = new ReferenceTwists(ref m_tangents, ref m_referenceFrames1);
        m_twists = new Twists(ref m_referenceTwists, ref m_dofs);
        m_curvatureBinormals = new CurvatureBinormals(ref m_tangents);
        m_trigThetas = new TrigThetas(ref m_dofs);
        m_materialFrames1 = new MaterialFrames1(ref m_trigThetas, ref m_referenceFrames1, ref m_referenceFrames2);
        m_materialFrames2 = new MaterialFrames2(ref m_trigThetas, ref m_referenceFrames1, ref m_referenceFrames2);
        m_twists = new Twists(ref m_referenceTwists, ref m_dofs);
        m_kappas = new Kappas(ref m_curvatureBinormals, ref m_materialFrames1, ref m_materialFrames2);
    }
}

public class StrandForce : Force
{
    public List<int> m_verts;  // in order root to tip
    private int m_globalIndex;         // Global index amongst the hairs
    public StrandParameters m_strandParams;
    private TwoDScene m_scene;
    private bool m_requiresExactForceJacobian;

    // increase memory, reduce re-computation
    private double m_strandEnergyUpdate;
    private Vectors m_strandForceUpdate;
    private TripletXs m_strandHessianUpdate;
    private SparseXs m_hess;

    // Linear Compliant Implicit Euler solve
    private Vectors m_lambda;
    private Vectors m_lambda_v;

    //// Strand State (implicitly the end of timestep state, evolved from rest
    ///config) ////////////////////////
    public StrandState m_strandState;  // future state
    public StartState m_startState;    // current state

    //// Rest shape //////////////////////////////////////////////////////
    public List<double>
        m_restLengths;  // The following four members depend on m_restLengths,
                        // which is why updateEverythingThatDependsOnRestLengths()
                        // must be called
    private double m_totalRestLength;
    public List<double> m_VoronoiLengths;     // rest length around each vertex
    public List<double> m_invVoronoiLengths;  // their inverses
    public List<double> m_vertexMasses;
    public List<Vectors> m_restKappas;
    public List<double> m_restTwists;

    public StrandForce(ref TwoDScene scene, 
        in List<int> consecutiveStrandVertices,
        in int parameterIndex, int globalIndex)
    {
        m_verts = new List<int>(consecutiveStrandVertices.ToArray());
        m_globalIndex = globalIndex;
        m_strandParams = null;
        m_scene = scene;
        m_requiresExactForceJacobian = false;
        m_strandEnergyUpdate = 0;
        m_strandHessianUpdate = new TripletXs();
        m_strandState = new StrandState();
        m_startState = new StartState();
        m_strandParams = m_scene.getStrandParameters(parameterIndex);
        Vectors initDoFs = new Vectors(getNumVertices() * 4 - 1);
        for (int i = 0; i < getNumVertices(); ++i)
        {
             initDoFs.SetSegment(i * 4, 3,
                 m_scene.getX().segment(m_scene.getDof(m_verts[i]), 3));
        }
        var bendingMatirxBase = m_strandParams.getBendingMatrixBase();
        m_strandState =
            new StrandState(initDoFs, ref bendingMatirxBase);
        m_startState = new StartState(initDoFs);
        resizeInternals();
        freezeRestShape(0, getNumEdges(), 0);

        m_lambda = new Vectors(numConstraintNonViscous());
        m_lambda_v = new Vectors(numConstraintViscous());
    }

    public List<Vectors> alterRestKappas()
    {
        return m_restKappas;
    }

    public int getNumVertices() { return m_verts.Count; }

    public int getNumEdges() { return m_verts.Count - 1; }

    public void resizeInternals()
    {
        m_restLengths = new List<double>(getNumEdges());
        m_restKappas = new List<Vectors>(getNumEdges());
        m_restTwists = new List<double>(getNumEdges());
        m_vertexMasses = new List<double>(getNumVertices());
        m_VoronoiLengths = new List<double>(getNumVertices());
        m_invVoronoiLengths = new List<double>(getNumVertices());

        for (int i = 0; i < getNumEdges(); i++)
        {
            m_restLengths.Add(1);
            m_restKappas.Add(new Vectors(2));
            m_restTwists.Add(0);
        }

        for (int i = 0; i < getNumVertices(); i++)
        {
            m_vertexMasses.Add(0);
            m_VoronoiLengths.Add(0);
            m_invVoronoiLengths.Add(0);
        }
    }

    public void updateEverythingThatDependsOnRestLengths()
    {
        m_totalRestLength = 0.0;
        for (ushort vtx = 0; vtx < getNumEdges(); ++vtx)
        {
            m_totalRestLength += m_restLengths[vtx];
        }

        m_VoronoiLengths[0] = 0.5 * m_restLengths[0];
        for (ushort vtx = 1; vtx < getNumEdges(); ++vtx)
        {
            m_VoronoiLengths[vtx] = 0.5 * (m_restLengths[vtx - 1] + m_restLengths[vtx]);
        }
        m_VoronoiLengths[getNumEdges()] = 0.5 * m_restLengths[getNumVertices() - 2];

        for (ushort vtx = 0; vtx < getNumVertices(); ++vtx)
        {
            m_vertexMasses[vtx] = m_strandParams.m_density * m_VoronoiLengths[vtx] *
                                  CMath.M_PI *
                                  m_strandParams.getRadius(vtx, getNumVertices()) *
                                  m_strandParams.getRadius(vtx, getNumVertices());
            m_invVoronoiLengths[vtx] = 1.0 / m_VoronoiLengths[vtx];
        }
    }

    public void updateRestShape(in Vectors dof_restshape, double damping)
    {
        StartState restshape_state = new StartState(dof_restshape);

        int nedges = getNumEdges();
        for (int vtx = 0; vtx < nedges; ++vtx)
        {  // Fix rest lengths
            m_restLengths[vtx] = (1 - damping) * restshape_state.m_lengths[vtx] +
                                 damping * m_restLengths[vtx];
        }
        updateEverythingThatDependsOnRestLengths();

        for (int vtx = 0; vtx < nedges; ++vtx)
        {
            m_restKappas[vtx] = (1 - damping) * restshape_state.m_kappas[vtx] +
                                damping * m_restKappas[vtx];
            m_restTwists[vtx] = (1 - damping) * restshape_state.m_twists[vtx] +
                                damping * m_restTwists[vtx];
        }
    }

    public override int getAffectedHair(in List<int> particle_to_hairs)
    {
        return particle_to_hairs[m_verts[0]];
    }

    public override void getAffectedVars(int pidx, ref HashSet<int> vars)
    {
        throw new System.NotImplementedException();
    }

    public override bool isContained(int pidx)
    {
        throw new System.NotImplementedException();
    }

    public override bool isParallelized()
    {
        return false;
    }

    public override bool isPrecomputationParallelized()
    {
        return false;
    }

    public override string name()
    {
        throw new System.NotImplementedException();
    }

    public override void preCompute(in Vectors x, in Vectors v, in Vectors m, in double dt)
    {
        base.preCompute(x, v, m, dt);
    }

    public void freezeRestShape(int begin, int end, double damping)
    {  // Take the current configuration as rest shape

        for (int vtx = begin; vtx < end; ++vtx)
        {  // Fix rest lengths
            m_restLengths[vtx] = (1 - damping) * m_strandState.m_lengths[vtx] +
                                 damping * m_restLengths[vtx];
        }
        updateEverythingThatDependsOnRestLengths();

        for (int vtx = begin; vtx < end; ++vtx)
        {
            m_restKappas[vtx] = (1 - damping) * m_strandState.m_kappas[vtx] + damping * m_restKappas[vtx];
            m_restTwists[vtx] = (1 - damping) * m_strandState.m_twists[vtx] + damping * m_restTwists[vtx];
        }
    }

    public override void computeIntegrationVars(in Vectors x, in Vectors v, in Vectors m, ref Vectors lambda, ref Vectors lambda_v, ref TripletXs J, ref TripletXs Jv, ref TripletXs Jxv, ref TripletXs tildeK, ref TripletXs stiffness, ref TripletXs damping, ref Vectors Phi, ref Vectors Phiv, in double dt)
    {
        Vectors futureStrandDoFs =
            x.segment(m_scene.getDof(m_verts[0]), getNumVertices() * 4 - 1);
        if (futureStrandDoFs != m_strandState.m_dofs.get())
        {
            m_strandState.m_dofs.set(futureStrandDoFs);
        }

        int ncnv = numConstraintNonViscous();
        lambda.SetSegment(m_internal_index_pos, ncnv, m_lambda);
        if (m_strandParams.m_accumulateWithViscous)
            lambda.SetSegment(m_internal_index_pos + ncnv, numConstraintViscous(),
                m_lambda_v);

        Vectors combinedLambda = m_lambda;
        if (m_strandParams.m_accumulateWithViscous)
        {
            if (!m_strandParams.m_accumulateVisCousOnlyForBendingModes)
            {
                combinedLambda += m_lambda_v;
            }
            else
            {
                combinedLambda.PlusSegment(getNumEdges(), m_lambda_v.size(), m_lambda_v);
            }
        }

        int global_start_dof = m_scene.getDof(m_verts[0]);

        int poffset = 0;
        int jOffset = 0;
        int koffset = 0;

        StretchingForceNonViscous.accumulateIntegrationVars(
            m_internal_index_pos, m_internal_index_J, m_internal_index_tildeK,
            global_start_dof, this, ref combinedLambda, ref J, ref tildeK, ref stiffness, ref Phi,
            poffset);
        poffset += getNumEdges();
        jOffset += getNumEdges() * 6;
        koffset += getNumEdges() * 36;

        TwistingForceNonViscous.accumulateIntegrationVars(
            m_internal_index_pos + poffset, m_internal_index_J + jOffset,
            m_internal_index_tildeK + koffset, global_start_dof, this,
            ref combinedLambda, ref J, ref tildeK, ref stiffness, ref Phi, poffset);
        poffset += (getNumVertices() - 2);
        jOffset += (getNumVertices() - 2) * 11;
        koffset += (getNumVertices() - 2) * 121;

        BendingForceNonViscous.accumulateIntegrationVars(
            m_internal_index_pos + poffset, m_internal_index_J + jOffset,
            m_internal_index_tildeK + koffset, global_start_dof, this,
            ref combinedLambda, ref J, ref tildeK, ref stiffness, ref Phi, poffset);

        poffset += 2 * (getNumVertices() - 2);
        jOffset += 2 * (getNumVertices() - 2) * 11;
        koffset += 2 * (getNumVertices() - 2) * 121;

        if (m_strandParams.m_accumulateWithViscous)
        {
            if (!m_strandParams.m_accumulateVisCousOnlyForBendingModes)
            {
                StretchingForceViscous.accumulateIntegrationVars(
                    m_internal_index_pos + poffset, m_internal_index_J + jOffset,
                    m_internal_index_tildeK + koffset, global_start_dof, this,
                    ref combinedLambda, ref J, ref tildeK, ref stiffness, ref Phi, poffset);
                poffset += getNumEdges();
                jOffset += getNumEdges() * 6;
                koffset += getNumEdges() * 36;
            }
            TwistingForceViscous.accumulateIntegrationVars(
                m_internal_index_pos + poffset, m_internal_index_J + jOffset,
                m_internal_index_tildeK + koffset, global_start_dof, this,
                ref combinedLambda, ref J, ref tildeK, ref stiffness, ref Phi, poffset);
            poffset += (getNumVertices() - 2);
            jOffset += (getNumVertices() - 2) * 11;
            koffset += (getNumVertices() - 2) * 121;

            BendingForceViscous.accumulateIntegrationVars(
                m_internal_index_pos + poffset, m_internal_index_J + jOffset,
                m_internal_index_tildeK + koffset, global_start_dof, this,
                ref combinedLambda, ref J, ref tildeK, ref stiffness, ref Phi, poffset);
        }
    }

    public override void storeLambda(in Vectors lambda, in Vectors lamda_v)
    {
        int ncnv = numConstraintNonViscous();
        m_lambda = lambda.segment(m_internal_index_pos, ncnv);
        if (m_strandParams.m_accumulateWithViscous)
        {
            m_lambda_v =
                lambda.segment(m_internal_index_pos + ncnv, numConstraintViscous());
        }
    }

    public override int numConstraintPos()
    {
        int numConstraintPos = numConstraintNonViscous();
        numConstraintPos += numConstraintViscous();
        return numConstraintPos;
    }

    int numConstraintNonViscous()
    {
        int numConstraintNonViscous = 0;

        numConstraintNonViscous += getNumEdges();  // Stretch

        numConstraintNonViscous += (getNumVertices() - 2);  // Twist

        numConstraintNonViscous += 2 * (getNumVertices() - 2);  // Bend

        return numConstraintNonViscous;
    }

    int numConstraintViscous()
    {
        int numConstraintViscous = 0;
        if (m_strandParams.m_accumulateWithViscous)
        {
            if (!m_strandParams.m_accumulateVisCousOnlyForBendingModes)
            {
                numConstraintViscous += getNumEdges();  // Stretch
            }
            numConstraintViscous += (getNumVertices() - 2);  // Twist
            numConstraintViscous += 2 * (getNumVertices() - 2);  // Bend
        }
        return numConstraintViscous;
    }

    public override int numConstraintVel()
    {
        return 0;
    }

    public override int numJ()
    {
        int numJ = 0;

        numJ += (getNumEdges() * 6);  // Stretch

        numJ += ((getNumVertices() - 2) * 11);  // Twist

        numJ += 2 * ((getNumVertices() - 2) * 11);  // Bend

        if (m_strandParams
                .m_accumulateWithViscous)
        {  // double the active forces for viscous
            if (!m_strandParams.m_accumulateVisCousOnlyForBendingModes)
            {
                numJ += (getNumEdges() * 6);  // Stretch
            }
            numJ += ((getNumVertices() - 2) * 11);  // Twist
            numJ += 2 * ((getNumVertices() - 2) * 11);  // Bend
        }

        return numJ;
    }

    public override int numJv()
    {
        return 0;
    }

    public override int numJxv()
    {
        return numJv();
    }

    public override int numTildeK()
    {
        int numTildeK = 0;
        numTildeK += getNumEdges() * 36;  // Stretch
        numTildeK += (getNumVertices() - 2) * 121;  // Twist
        numTildeK += 2 * (getNumVertices() - 2) * 121;  // Bend

        return numTildeK;
    }

    public void updateStartDoFs(in Vectors x_startOfStep) 
    {
        Vectors currentStrandDoFs = new Vectors(getNumVertices() * 4 - 1);
        currentStrandDoFs.SetSegment(0, getNumVertices() * 4 - 1, x_startOfStep.segment(
        m_scene.getDof(m_verts[0]), getNumVertices() * 4 - 1));
        m_startState.m_dofs.set(currentStrandDoFs);
    }
}

public class Viscous
{
    public static double bendingCoefficient (in StrandForce strand, int vtx) {
        return strand.m_strandParams.viscousBendingCoefficient(
            vtx, strand.getNumVertices());
    }

    public static MatrixXs bendingMatrix(in StrandForce strand, int vtx)
    {
        return strand.m_strandParams.viscousBendingMatrix(vtx,
                                                       strand.getNumVertices());
    }

    public static Vectors kappaBar(in StrandForce strand, int vtx)
    {
        return strand.m_startState.m_kappas[vtx];
    }

    public static double kt(in StrandForce strand, int vtx)
    {
        return strand.m_strandParams.getViscousKt(vtx, strand.getNumVertices());
    }

    public static double thetaBar(in StrandForce strand, int vtx)
    {
        return strand.m_startState.m_twists[vtx];
    }

    public static double ks(in StrandForce strand, int vtx)
    {
        return strand.m_strandParams.getViscousKs(vtx, strand.getNumVertices());
    }

    public static double ellBar(in StrandForce strand, int vtx)
    {
        return strand.m_startState.m_lengths[vtx];
    }
}



public class NonViscous
{
    public static double bendingCoefficient(in StrandForce strand, int vtx) {
        return strand.m_strandParams.bendingCoefficient(vtx,
                                                     strand.getNumVertices());
    }

    public static MatrixXs bendingMatrix(in StrandForce strand, int vtx)
    {
        return strand.m_strandParams.bendingMatrix(vtx, strand.getNumVertices());
    }

    public static Vectors kappaBar(in StrandForce strand, int vtx) {
        return strand.m_restKappas[vtx];
    }

    public static double kt(in StrandForce strand, int vtx)
    {
        return strand.m_strandParams.getKt(vtx, strand.getNumVertices());
    }

    public static double thetaBar(in StrandForce strand, int vtx)
    {
        return strand.m_restTwists[vtx];
    }

    public static double ks(in StrandForce strand, int vtx) {
        return strand.m_strandParams.getKs(vtx, strand.getNumVertices());
    }

    public static double ellBar(in StrandForce strand, int vtx)
    {
        return strand.m_restLengths[vtx];
    }
}

public class BendingForceViscous
{
    public static readonly int s_first = 1;
    public static readonly int s_last = 1;

    public static void accumulateIntegrationVars(in int pos_start, in int j_start,
        in int tildek_start, in int global_start_dof,
        StrandForce strand, ref Vectors lambda, ref TripletXs J, ref TripletXs tildeK,
        ref TripletXs stiffness, ref Vectors Phi, in int lambda_start)
    {
        MatrixXs B = Viscous.bendingMatrix(strand, 1);
        double b = B[0, 0];  // assumes circular rods

        for (int vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx)
        {
            int dfirst = global_start_dof + 4 * (vtx - 1);
            int dsecond = global_start_dof + 4 * vtx;
            int dthird = global_start_dof + 4 * (vtx + 1);

            Vectors kappaBar = Viscous.kappaBar(strand, vtx);
            Vectors kappa = strand.m_strandState.m_kappas[vtx];
            double ilen = strand.m_invVoronoiLengths[vtx];

            int idx_pos = pos_start + 2 * (vtx - 1);
            Phi.SetSegment(idx_pos, 2, kappa - kappaBar);
            stiffness[idx_pos] = new Triplets(idx_pos, idx_pos, b * ilen);
            stiffness[idx_pos + 1] = new Triplets(idx_pos + 1, idx_pos + 1, b * ilen);

            MatrixXs utgk =
                new MatrixXs(strand.m_strandState.m_gradKappas[vtx]).transpose();  //2 x 11
            int idx_j = j_start + 2 * (11 * (vtx - 1));
            for (int r = 0; r < 4; ++r)
            {
                J[idx_j + r] = new Triplets(idx_pos, dfirst + r, utgk[0, r]);
                J[idx_j + 4 + r] = new Triplets(idx_pos, dsecond + r, utgk[0, r + 4]);
                if (r < 3)
                    J[idx_j + 8 + r] = new Triplets(idx_pos, dthird + r, utgk[0, r + 8]);

                J[idx_j + r + 11] = new Triplets(idx_pos + 1, dfirst + r, utgk[1, r]);
                J[idx_j + 4 + r + 11] =
                    new Triplets(idx_pos + 1, dsecond + r, utgk[1, r + 4]);
                if (r < 3)
                    J[idx_j + 8 + r + 11] =
                        new Triplets(idx_pos + 1, dthird + r, utgk[1, r + 8]);
            }
        }
    }
}

public class BendingForceNonViscous
{
    public static readonly int s_first = 1;
    public static readonly int s_last = 1;

    public static void accumulateIntegrationVars(in int pos_start, in int j_start,
        in int tildek_start, in int global_start_dof,
        StrandForce strand, ref Vectors lambda, ref TripletXs J, ref TripletXs tildeK,
        ref TripletXs stiffness, ref Vectors Phi, in int lambda_start)
    {
        MatrixXs B = NonViscous.bendingMatrix(strand, 1);
        double b = B[0, 0];  // assumes circular rods

        for (int vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx)
        {
            int dfirst = global_start_dof + 4 * (vtx - 1);
            int dsecond = global_start_dof + 4 * vtx;
            int dthird = global_start_dof + 4 * (vtx + 1);

            Vectors kappaBar = NonViscous.kappaBar(strand, vtx);
            Vectors kappa = strand.m_strandState.m_kappas[vtx];

            double ilen = strand.m_invVoronoiLengths[vtx];

            int idx_pos = pos_start + 2 * (vtx - 1);
            Phi.SetSegment(idx_pos, 2, kappa - kappaBar);
            stiffness[idx_pos] = new Triplets(idx_pos, idx_pos, b * ilen);
            stiffness[idx_pos + 1] = new Triplets(idx_pos + 1, idx_pos + 1, b * ilen);

            MatrixXs utgk =
                new MatrixXs(strand.m_strandState.m_gradKappas[vtx]).transpose();  //2 x 11
            int idx_j = j_start + 2 * (11 * (vtx - 1));
            for (int r = 0; r < 4; ++r)
            {
                J[idx_j + r] = new Triplets(idx_pos, dfirst + r, utgk[0, r]);
                J[idx_j + 4 + r] = new Triplets(idx_pos, dsecond + r, utgk[0, r + 4]);
                if (r < 3)
                    J[idx_j + 8 + r] = new Triplets(idx_pos, dthird + r, utgk[0, r + 8]);

                J[idx_j + r + 11] = new Triplets(idx_pos + 1, dfirst + r, utgk[1, r]);
                J[idx_j + 4 + r + 11] =
                    new Triplets(idx_pos + 1, dsecond + r, utgk[1, r + 4]);
                if (r < 3)
                    J[idx_j + 8 + r + 11] =
                        new Triplets(idx_pos + 1, dthird + r, utgk[1, r + 8]);
            }

            
            int lidx = lambda_start + (2 * (vtx - 1));
            double weight0 = -lambda[lidx];
            double weight1 = -lambda[lidx + 1];

            Tuple<MatrixXs, MatrixXs> hessKappa =
                strand.m_strandState.m_hessKappas[vtx];
            int idx_tildek = tildek_start + 2 * (121 * (vtx - 1));
            for (int r = 0; r < 11; ++r)
            {
                for (int s = 0; s < 11; ++s)
                {
                    if (r == 3 || r == 7 || s == 3 || s == 7)
                    {
                        tildeK[idx_tildek + r * 11 + s] = new Triplets(
                            dfirst + r, dfirst + s,
                            hessKappa.Item1[r, s] * b * ilen * (kappa[0] - kappaBar[0]));
                        tildeK[idx_tildek + r * 11 + s + 121] = new Triplets(
                            dfirst + r, dfirst + s,
                            hessKappa.Item2[r, s] * b * ilen * (kappa[1] - kappaBar[1]));
                    }
                    else
                    {
                        tildeK[idx_tildek + r * 11 + s] = new Triplets(
                            dfirst + r, dfirst + s, hessKappa.Item1[r, s] * weight0);
                        tildeK[idx_tildek + r * 11 + s + 121] = new Triplets(
                            dfirst + r, dfirst + s, hessKappa.Item2[r, s] * weight1);
                    }
                }
            }
        }
    }
}

public class StretchingForceViscous
{
    public static readonly int s_first = 0;
    public static readonly int s_last = 1;

    public static void accumulateIntegrationVars(
        in int pos_start, in int j_start,
        in int tildek_start, in int global_start_dof,
        StrandForce strand, ref Vectors lambda, ref TripletXs J, ref TripletXs tildeK,
        ref TripletXs stiffness, ref Vectors Phi, in int lambda_start)
    {
        for (int vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx)
        {
            int dfirst = global_start_dof + 4 * vtx;
            int dsecond = global_start_dof + 4 * (vtx + 1);

            double length = strand.m_strandState.m_lengths[vtx];
            int idx_pos = pos_start + vtx;
            Phi[idx_pos] = length - Viscous.ellBar(strand, vtx);
            stiffness[idx_pos] =
                new Triplets(idx_pos, idx_pos,
                         Viscous.ks(strand, vtx) / Viscous.ellBar(strand, vtx));

            Vectors edge = strand.m_strandState.m_tangents[vtx];
            int idx_j = j_start + 6 * vtx;
            for (int r = 0; r < 3; ++r)
            {
                J[idx_j + r] = new Triplets(idx_pos, dfirst + r, -edge[r]);
                J[idx_j + 3 + r] = new Triplets(idx_pos, dsecond + r, edge[r]);
            }
        }
    }
}

public class StretchingForceNonViscous
{
    public static readonly int s_first = 0;
    public static readonly int s_last = 1;

    public static void accumulateIntegrationVars(
        in int pos_start, in int j_start,
        in int tildek_start, in int global_start_dof,
        StrandForce strand, ref Vectors lambda, ref TripletXs J, ref TripletXs tildeK,
        ref TripletXs stiffness, ref Vectors Phi, in int lambda_start)
    {
        for (int vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx)
        {
            int dfirst = global_start_dof + 4 * vtx;
            int dsecond = global_start_dof + 4 * (vtx + 1);

            double length = strand.m_strandState.m_lengths[vtx];
            int idx_pos = pos_start + vtx;
            Phi[idx_pos] = length - NonViscous.ellBar(strand, vtx);
            stiffness[idx_pos] =
                new Triplets(idx_pos, idx_pos,
                         NonViscous.ks(strand, vtx) / NonViscous.ellBar(strand, vtx));

            Vectors edge = strand.m_strandState.m_tangents[vtx];
            int idx_j = j_start + 6 * vtx;
            for (int r = 0; r < 3; ++r)
            {
                J[idx_j + r] = new Triplets(idx_pos, dfirst + r, -edge[r]);
                J[idx_j + 3 + r] = new Triplets(idx_pos, dsecond + r, edge[r]);
            }

            MatrixXs M = (CMath.identity(3) - CMath.matProduct(edge, edge)) / length;
            double weight = -lambda[lambda_start + vtx];

            int idx_tildek = tildek_start + 36 * vtx;
            for (int r = 0; r < 3; ++r)
            {
                for (int s = 0; s < 3; ++s)
                {
                    tildeK[idx_tildek + r * 6 + s] =
                        new Triplets(dfirst + r, dfirst + s, M[r, s] * weight);
                    tildeK[idx_tildek + (r + 3) * 6 + 3 + s] =
                        new Triplets(dsecond + r, dsecond + s, M[r, s] * weight);

                    tildeK[idx_tildek + r * 6 + 3 + s] =
                        new Triplets(dfirst + r, dsecond + s, -M[r, s] * weight);
                    tildeK[idx_tildek + (r + 3) * 6 + s] =
                        new Triplets(dsecond + r, dfirst + s, -M[r, s] * weight);
                }
            }
        }
    }
}


public class TwistingForceViscous
{
    public static readonly int s_first = 1;
    public static readonly int s_last = 1;

    public static void accumulateIntegrationVars(
        in int pos_start, in int j_start,
        in int tildek_start, in int global_start_dof,
        StrandForce strand, ref Vectors lambda, ref TripletXs J, ref TripletXs tildeK,
        ref TripletXs stiffness, ref Vectors Phi, in int lambda_start)
    {
        for (int vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx)
        {
            int dfirst = global_start_dof + 4 * (vtx - 1);
            int dsecond = global_start_dof + 4 * vtx;
            int dthird = global_start_dof + 4 * (vtx + 1);

            double twist = strand.m_strandState.m_twists[vtx];
            double undeformedTwist = Viscous.thetaBar(strand, vtx);
            double kt = Viscous.kt(strand, vtx);
            double ilen = strand.m_invVoronoiLengths[vtx];

            int idx_pos = pos_start + vtx - 1;
            Phi[idx_pos] = twist - undeformedTwist;
            stiffness[idx_pos] = new Triplets(idx_pos, idx_pos, kt * ilen);

            Vectors gt = (strand.m_strandState.m_gradTwists[vtx]);
            int idx_j = j_start + 11 * (vtx - 1);
            for (int r = 0; r < 4; ++r)
            {
                J[idx_j + r] = new Triplets(idx_pos, dfirst + r, gt[r]);
                J[idx_j + 4 + r] = new Triplets(idx_pos, dsecond + r, gt[r + 4]);
                if (r < 3)
                    J[idx_j + 8 + r] = new Triplets(idx_pos, dthird + r, gt[r + 8]);
            }
        }
    }
}

public class TwistingForceNonViscous
{
    public static readonly int s_first = 1;
    public static readonly int s_last = 1;

    public static void accumulateIntegrationVars(
        in int pos_start, in int j_start,
        in int tildek_start, in int global_start_dof,
        StrandForce strand, ref Vectors lambda, ref TripletXs J, ref TripletXs tildeK,
        ref TripletXs stiffness, ref Vectors Phi, in int lambda_start)
    {
        for (int vtx = s_first; vtx < strand.getNumVertices() - s_last; ++vtx)
        {
            int dfirst = global_start_dof + 4 * (vtx - 1);
            int dsecond = global_start_dof + 4 * vtx;
            int dthird = global_start_dof + 4 * (vtx + 1);
            double twist = strand.m_strandState.m_twists[vtx];
            double undeformedTwist = NonViscous.thetaBar(strand, vtx);
            double kt = NonViscous.kt(strand, vtx);
            double ilen = strand.m_invVoronoiLengths[vtx];

            int idx_pos = pos_start + vtx - 1;
            Phi[idx_pos] = twist - undeformedTwist;
            stiffness[idx_pos] = new Triplets(idx_pos, idx_pos, kt * ilen);

            Vectors gt = (strand.m_strandState.m_gradTwists[vtx]);
            int idx_j = j_start + 11 * (vtx - 1);
            for (int r = 0; r < 4; ++r)
            {
                J[idx_j + r] = new Triplets(idx_pos, dfirst + r, gt[r]);
                J[idx_j + 4 + r] = new Triplets(idx_pos, dsecond + r, gt[r + 4]);
                if (r < 3)
                    J[idx_j + 8 + r] = new Triplets(idx_pos, dthird + r, gt[r + 8]);
            }

            double weight = -lambda[lambda_start + (vtx - 1)];
            MatrixXs hessTwist = strand.m_strandState.m_hessTwists[vtx];
            int idx_tildek = tildek_start + 121 * (vtx - 1);
            for (int r = 0; r < 11; ++r)
            {
                for (int s = 0; s < 11; ++s)
                {
                    if (r == 3 || r == 7 || s == 3 || s == 7)
                    {
                        tildeK[idx_tildek + r * 11 + s] = new Triplets(
                        dfirst + r, dfirst + s,
                            hessTwist[r, s] * kt * ilen * (twist - undeformedTwist));
                    }
                    else
                    {
                        tildeK[idx_tildek + r * 11 + s] =
                            new Triplets(dfirst + r, dfirst + s, hessTwist[r, s] * weight);
                    }
                }
            }
        }
    }
}

