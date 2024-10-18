using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;

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
    //GradKappas m_gradKappas;
    //GradTwists m_gradTwists;
    //GradTwistsSquared m_gradTwistsSquared;
    //HessKappas m_hessKappas;
    //HessTwists m_hessTwists;
    //BendingProducts m_bendingProducts;

    public StrandState(in VectorXs initDoFs,
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
        m_twists = new Twists(ref m_referenceTwists, ref m_dofs);
        m_kappas = new Kappas(ref m_curvatureBinormals, ref m_materialFrames1, ref m_materialFrames2);
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

    public StartState(in VectorXs initDoFs)
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
    private VectorXs m_strandForceUpdate;
    private TripletXs m_strandHessianUpdate;
    private SparseXs m_hess;

    // Linear Compliant Implicit Euler solve
    private VectorXs m_lambda;
    private VectorXs m_lambda_v;

    //// Strand State (implicitly the end of timestep state, evolved from rest
    ///config) ////////////////////////
    private StrandState m_strandState;  // future state
    private StartState m_startState;    // current state

    //// Rest shape //////////////////////////////////////////////////////
    public List<double>
        m_restLengths;  // The following four members depend on m_restLengths,
                        // which is why updateEverythingThatDependsOnRestLengths()
                        // must be called
    private double m_totalRestLength;
    private List<double> m_VoronoiLengths;     // rest length around each vertex
    private List<double> m_invVoronoiLengths;  // their inverses
    public List<double> m_vertexMasses;
    private List<Vectors> m_restKappas;
    private List<double> m_restTwists;

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
        VectorXs initDoFs = new VectorXs(getNumVertices() * 4 - 1);
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
        //freezeRestShape(
        //    0, getNumEdges());
    }

    public List<Vectors> alterRestKappas()
    {
        return m_restKappas;
    }

    public int getNumVertices() { return m_verts.Count; }

    public int getNumEdges() { return m_verts.Count - 1; }

    public void resizeInternals()
    {
        m_restLengths = new List<double>(new double[getNumEdges()]);
        m_restKappas = new List<Vectors> (new Vectors[getNumEdges()]);
        m_restTwists = new List<double>(new double[getNumEdges()]);
        m_vertexMasses = new List<double>(new double[getNumVertices()]);
        m_VoronoiLengths = new List<double>(new double[getNumVertices()]);
        m_invVoronoiLengths = new List<double>(new double[getNumVertices()]);
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

    public void updateRestShape(in VectorXs dof_restshape, double damping)
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
            m_restKappas[vtx] = (1 - damping) * restshape_state.m_kappas[vtx].ToVectors() +
                                damping * m_restKappas[vtx];
            m_restTwists[vtx] = (1 - damping) * restshape_state.m_twists[vtx] +
                                damping * m_restTwists[vtx];
        }
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
        throw new System.NotImplementedException();
    }

    public override bool isPrecomputationParallelized()
    {
        throw new System.NotImplementedException();
    }

    public override string name()
    {
        throw new System.NotImplementedException();
    }

    public override int numConstraintPos()
    {
        throw new System.NotImplementedException();
    }

    public override int numConstraintVel()
    {
        throw new System.NotImplementedException();
    }

    public override int numJ()
    {
        throw new System.NotImplementedException();
    }

    public override int numJv()
    {
        throw new System.NotImplementedException();
    }

    public override int numJxv()
    {
        throw new System.NotImplementedException();
    }

    public override int numTildeK()
    {
        throw new System.NotImplementedException();
    }

    public void updateStartDoFs(in VectorXs x_startOfStep) 
    {
        VectorXs currentStrandDoFs = new VectorXs(getNumVertices() * 4 - 1);
        currentStrandDoFs.SetSegment(0, getNumVertices() * 4 - 1, x_startOfStep.segment(
        m_scene.getDof(m_verts[0]), getNumVertices() * 4 - 1));
        m_startState.m_dofs.set(currentStrandDoFs);
    }
}
