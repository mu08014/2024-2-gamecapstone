using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;

// Unit : cm for position dofs, no dimension for theta
public class DOFs : DependencyNode<VectorXs>
{
    private ushort m_numEdges;
    public DOFs(in VectorXs dofValues) : base(dofValues)
    {
        m_numEdges = (ushort)((ushort)dofValues.size() / 4);
        setClean();
    }

    public override VectorXs get() { return m_value; }

    public Vectors getVertex(ushort vtx)
    {
        return get().segment(4 * vtx, 3);
    }

    public void setVertex(ushort vtx, Vectors point)
    {
        for (int i = 0; i < 3; i++)
        {
            m_value[4 * vtx + i] = point[i];
        }
        setDependentsDirty();
    }

    public void setDof(ushort i, double val)
    {
        m_value[i] = val;
        setDependentsDirty();
    }

    public double getTheta(ushort vtx)
    {
        return get()[4 * vtx + 3];
    }

    public void setTheta(ushort vtx, double theta)
    {
        m_value[4 * vtx + 3] = theta;
        setDependentsDirty();
    }

    public ushort getNumEdges() { return m_numEdges; }

    public ushort getNumVertices() { return (ushort)(m_numEdges + 1); }

    public string name() { return "DOFs"; }

    protected override void compute()
    {
        throw new System.NotImplementedException();
    }
}

//Unit : cm
public class Edges : DependencyNode<List<Vectors>>
{
    protected DOFs m_dofs;
    public Edges(ref DOFs dofs) : base(0, dofs.getNumEdges())
    {
        m_dofs = dofs;
        m_dofs.addDependent(this);
    }

    public string name() { return "Edges"; }

    protected override void compute()
    {
        m_value.Capacity = m_size;
        VectorXs dofs = m_dofs.get().Clone();
        for (ushort vtx = 0; vtx < m_firstValidIndex; ++vtx)
        {
            m_value[vtx].setZero();
        }
        for (ushort vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            m_value[vtx] = dofs.segment(4 * (vtx + 1), 3) - dofs.segment(4 * vtx, 3);
        }

        setDependentsDirty();
    }
}

//Unit: cm
public class Lengths : DependencyNode<List<double>>
{
    protected Edges m_edges;
    public Lengths(ref Edges edges) : base(0, edges.size())
    {
        m_edges = edges;
        m_edges.addDependent(this);
    }

    public string name() { return "Lengths"; }

    protected override void compute()
    {
        m_value.Capacity = m_size;
        List<Vectors> edges = m_edges.get();
        for (ushort vtx = 0; vtx < m_firstValidIndex; ++vtx)
        {
            m_value[vtx] = 0.0;
        }
        for (ushort vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            m_value[vtx] = edges[vtx].norm();
        }

        setDependentsDirty();
    }
}

// Unit : no dimension
public class Tangents : DependencyNode<List<Vectors>>
{
    protected Edges m_edges;
    protected Lengths m_lengths;

    public Tangents(ref Edges edges, ref Lengths lengths) : base(0, edges.size())
    {
        m_edges = edges;
        m_lengths = lengths;
        m_edges.addDependent(this);
        m_lengths.addDependent(this);
    }

    public string name() { return "Tangents"; }

    protected override void compute()
    {
        m_value.Capacity = m_size;
        List<Vectors> edges = m_edges.get();
        List<double> lengths = m_lengths.get();
        for (ushort vtx = 0; vtx < m_firstValidIndex; ++vtx)
        {
            m_value[vtx].setZero();
        }
        for (ushort vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            m_value[vtx] = edges[vtx] / lengths[vtx];
        }

        setDependentsDirty();
    }
}

public class CurvatureBinormals : DependencyNode<List<Vectors>>
{
    protected Tangents m_tangents;

    public CurvatureBinormals(ref Tangents tangents) : base(1, tangents.size())
    {
        m_tangents = tangents;
        m_tangents.addDependent(this);
    }

    public string name() { return "CurvatureBinormals"; }

    protected override void compute()
    {

    }
}

public class TrigThetas : DependencyNode<Tuple<VectorXs, VectorXs>>
{
    protected DOFs m_dofs;

    public TrigThetas(ref DOFs dofs) : base(new Tuple<VectorXs, VectorXs>(new VectorXs(), new VectorXs()))
    {
        m_dofs = dofs;
        m_dofs.addDependent(this);
    }

    public string name() { return "TrigThetas"; }

    protected override void compute()
    {
        throw new NotImplementedException();
    }
}
