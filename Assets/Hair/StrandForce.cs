using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

public struct StrandState
{
    public DOFs m_dofs;
    public Edges m_edges;
    public Lengths m_lengths;
    public Tangents m_tangents;
    //ReferenceFrames1 m_referenceFrames1;
    //ReferenceFrames2 m_referenceFrames2;
    //ReferenceTwists m_referenceTwists;
    //Twists m_twists;
    public CurvatureBinormals m_curvatureBinormals;
    public TrigThetas m_trigThetas;
    //MaterialFrames m_materialFrames1;
    //MaterialFrames m_materialFrames2;
    //Kappas m_kappas;
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

        m_curvatureBinormals = new CurvatureBinormals(ref m_tangents);
        m_trigThetas = new TrigThetas(ref m_dofs);
    }
}

public struct StartState
{
    public DOFs m_dofs;
    public Edges m_edges;
    public Lengths m_lengths;
    public Tangents m_tangents;
    //ReferenceFrames1 m_referenceFrames1;
    //ReferenceFrames2 m_referenceFrames2;
    //ReferenceTwists m_referenceTwists;
    //Twists m_twists;
    public CurvatureBinormals m_curvatureBinormals;
    public TrigThetas m_trigThetas;
    //MaterialFrames<1> m_materialFrames1;
    //MaterialFrames<2> m_materialFrames2;
    //Kappas m_kappas;

    public StartState(in VectorXs initDoFs)
    {
        m_dofs = new DOFs(initDoFs);
        m_edges = new Edges(ref m_dofs);
        m_lengths = new Lengths(ref m_edges);
        m_tangents = new Tangents(ref m_edges, ref m_lengths);

        m_curvatureBinormals = new CurvatureBinormals(ref m_tangents);
        m_trigThetas = new TrigThetas(ref m_dofs);
    }
}

public class StrandForce : Force
{
    private List<int> m_verts;  // in order root to tip
    private int m_globalIndex;         // Global index amongst the hairs
    private StrandParameters m_strandParams;
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
    private List<double>
        m_restLengths;  // The following four members depend on m_restLengths,
                        // which is why updateEverythingThatDependsOnRestLengths()
                        // must be called
    private double m_totalRestLength;
    private List<double> m_VoronoiLengths;     // rest length around each vertex
    private List<double> m_invVoronoiLengths;  // their inverses
    private List<double> m_vertexMasses;
    private List<Vectors> m_restKappas;
    private List<double> m_restTwists;

    public StrandForce(ref TwoDScene scene, 
        in List<int> consecutiveStrandVertices,
        in int parameterIndex, int globalIndex)
    {
        m_verts = consecutiveStrandVertices;
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
            if (m_scene.isTip(m_verts[i]))
                initDoFs.SetSegment(i * 4, 3,
                    m_scene.getX().segment(m_scene.getDof(m_verts[i]), 3));
            else
                initDoFs.SetSegment(i * 4, 4,
                    m_scene.getX().segment(m_scene.getDof(m_verts[i]), 4));
        }
        var bendingMatirxBase = m_strandParams.getBendingMatrixBase();
        m_strandState =
            new StrandState(initDoFs, ref bendingMatirxBase);
        m_startState = new StartState(initDoFs);

        resizeInternals();
        //freezeRestShape(
        //    0, getNumEdges());
    }

    public int getNumVertices() { return m_verts.Count; }

    public int getNumEdges() { return m_verts.Count - 1; }

    public void resizeInternals()
    {
        m_restLengths.Capacity = getNumEdges();
        m_restKappas.Capacity = getNumEdges();
        m_restTwists.Capacity = getNumEdges();
        m_vertexMasses.Capacity = getNumVertices();
        m_VoronoiLengths.Capacity = getNumVertices();
        m_invVoronoiLengths.Capacity = getNumVertices();
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
}
