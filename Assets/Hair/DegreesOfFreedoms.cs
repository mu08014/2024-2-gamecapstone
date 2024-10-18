using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;
using static UnityEngine.GraphicsBuffer;

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

    public Vectors this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "Edges"; }

    protected override void compute()
    {
        m_value = new List<Vectors>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new Vectors(3));
        }
        VectorXs dofs = m_dofs.get().Clone();
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

    public double this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "Lengths"; }

    protected override void compute()
    {
        m_value = new List<double>(new double[m_size]);
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

    public Vectors this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "Tangents"; }

    protected override void compute()
    {
        m_value = new List<Vectors>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new Vectors(3));
        }
        List<Vectors> edges = m_edges.get();
        List<double> lengths = m_lengths.get();
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

    public Vectors this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "CurvatureBinormals"; }

    protected override void compute()
    {
        m_value = new List<Vectors>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new Vectors(3));
        }

        List<Vectors> tangents = m_tangents.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            Vectors t1 = tangents[vtx - 1];
            Vectors t2 = tangents[vtx];

            double denominator = 1 + t1.dot(t2);

            if (denominator <= 0)
            {
                if (denominator <= 0)
                {
                    denominator = 1 + t1.normalized().dot(t2.normalized());
                }


                m_value[vtx] = 4 * Math.Tan(5 * Math.Acos(denominator - 1)) *
                                   CMath.findNormal(t1).segment(3, 0);
            }
            else
            {
                m_value[vtx] = 2.0 * t1.cross(t2) / denominator;
            }
        }

        setDependentsDirty();
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

    public void vdSinCos(int n, in VectorXs a, ref Tuple<VectorXs, VectorXs> t)
    {
        for (int i = 0; i < n; i++)
        {
            t.Item1[i] = Math.Sin(a[i]);
            t.Item2[i] = Math.Cos(a[i]);
        }
    }

    protected override void compute()
    {
        VectorXs dofs = m_dofs.get();
        int numThetas = m_dofs.getNumEdges();
        m_value.Item1.resize(numThetas);
        m_value.Item2.resize(numThetas);

        // Extract thetas in their own vector for mkl_vlm
        VectorXs thetasMap = new VectorXs(numThetas);
        for (int i = 0; i < numThetas; i++)
        {
            thetasMap[i] = dofs[i * 4 + 3];
        }

        VectorXs thetaVec = thetasMap.Clone();
        // Compute their sine and cosine
        // assert( typeid(double) == typeid(VecX::scalar) );
        vdSinCos(
            numThetas, thetaVec, ref m_value);  // FIXME this won't compile if scalar != double

        setDependentsDirty();
    }

    public VectorXs getSines() { return get().Item1; }

    public VectorXs getCosines() { return get().Item2; }

}

public class Kappas : DependencyNode<List<VectorXs>>
{
    private CurvatureBinormals m_curvatureBinormals;
    private MaterialFrames1 m_materialFrames1;
    private MaterialFrames2 m_materialFrames2;

    public Kappas(ref CurvatureBinormals curvatureBinormals, 
        ref MaterialFrames1 materialFrames1, ref MaterialFrames2 materialFrames2) : base(new List<VectorXs>(curvatureBinormals.size()), 1, curvatureBinormals.size())
    {
        m_curvatureBinormals = curvatureBinormals;
        m_materialFrames1 = materialFrames1;
        m_materialFrames2 = materialFrames2;
        m_curvatureBinormals.addDependent(this);
        m_materialFrames1.addDependent(this);
        m_materialFrames2.addDependent(this);
    }

    public VectorXs this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "Keppas"; }

    protected override void compute()
    {
        m_value = new List<VectorXs>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new VectorXs(2));
        }
        List<Vectors> curvatureBinormals = m_curvatureBinormals.get();
        List<VectorXs> materialFrames1 = m_materialFrames1.get();
        List<VectorXs> materialFrames2 = m_materialFrames2.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            Vectors kb = curvatureBinormals[vtx];
            VectorXs m1e = materialFrames1[vtx - 1];
            VectorXs m2e = materialFrames2[vtx - 1];
            VectorXs m1f = materialFrames1[vtx];
            VectorXs m2f = materialFrames2[vtx];

            m_value[vtx] = new VectorXs(2);
            m_value[vtx][0] = 0.5 * kb.dot((m2e + m2f).ToVectors());
            m_value[vtx][1] = -0.5 * kb.dot((m1e + m1f).ToVectors());
        }

        setDependentsDirty();
    }
}

public class Twists : DependencyNode<List<double>>
{
    private ReferenceTwists m_refTwists;
    private DOFs m_dofs;

    public Twists(ref ReferenceTwists refTwists, ref DOFs dofs) : base(new List<double>(dofs.getNumEdges()), 1, dofs.getNumEdges())
    {
        m_refTwists = refTwists;
        m_dofs = dofs;
        m_refTwists.addDependent(this);
        m_dofs.addDependent(this);
    }

    public double this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "Twists"; }

    protected override void compute()
    {
        m_value = new List<double>(new double[m_size]);
        List<double> refTwists = m_refTwists.get();
        VectorXs dofs = m_dofs.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            m_value[vtx] = refTwists[vtx] + dofs[4 * vtx + 3] - dofs[4 * vtx - 1];
        }

        setDependentsDirty();
    }
}

public class MaterialFrames1 : DependencyNode<List<VectorXs>>
{
    private TrigThetas m_trigThetas;
    private ReferenceFrames1 m_referencesFrames1;
    private ReferenceFrames2 m_referencesFrames2;

    public MaterialFrames1(ref TrigThetas trigThetas, ref ReferenceFrames1 referencesFrames1,
        ref ReferenceFrames2 referencesFrames2) : base(new List<VectorXs>(referencesFrames1.size()), 0, referencesFrames1.size())
    {
        m_trigThetas = trigThetas;
        m_referencesFrames1 = referencesFrames1;
        m_referencesFrames2 = referencesFrames2;
        m_trigThetas.addDependent(this);
        m_referencesFrames1.addDependent(this);
        m_referencesFrames2.addDependent(this);
    }

    public VectorXs this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "MaterialFrames1"; }

    protected override void compute()
    {
        m_value = new List<VectorXs>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new VectorXs(3));
        }
        List<VectorXs> referenceFrames1 = m_referencesFrames1.get();
        List<VectorXs> referenceFrames2 = m_referencesFrames2.get();
        VectorXs sinThetas = m_trigThetas.getSines();
        VectorXs cosThetas = m_trigThetas.getCosines();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            VectorXs u = referenceFrames1[vtx];
            VectorXs v = referenceFrames2[vtx];
            double s = sinThetas[vtx];
            double c = cosThetas[vtx];

            m_value[vtx] = linearMix(u, v, s, c);
        }

        setDependentsDirty();
    }

    public VectorXs linearMix(in VectorXs u, in VectorXs v, double s, double c)
    {
        return -c * u + s * v;
    }
}

public class MaterialFrames2 : DependencyNode<List<VectorXs>>
{
    private TrigThetas m_trigThetas;
    private ReferenceFrames1 m_referencesFrames1;
    private ReferenceFrames2 m_referencesFrames2;

    public MaterialFrames2(ref TrigThetas trigThetas, ref ReferenceFrames1 referencesFrames1, 
        ref ReferenceFrames2 referencesFrames2) : base(new List<VectorXs>(referencesFrames1.size()), 0, referencesFrames1.size())
    {
        m_trigThetas = trigThetas;
        m_referencesFrames1 = referencesFrames1;
        m_referencesFrames2 = referencesFrames2;
        m_trigThetas.addDependent(this);
        m_referencesFrames1.addDependent(this);
        m_referencesFrames2.addDependent(this);
    }

    public VectorXs this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "MaterialFrames2"; }

    protected override void compute()
    {
        m_value = new List<VectorXs>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new VectorXs(3));
        }
        List<VectorXs> referenceFrames1 = m_referencesFrames1.get();
        List<VectorXs> referenceFrames2 = m_referencesFrames2.get();
        VectorXs sinThetas = m_trigThetas.getSines();
        VectorXs cosThetas = m_trigThetas.getCosines();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            VectorXs u = referenceFrames1[vtx];
            VectorXs v = referenceFrames2[vtx];
            double s = sinThetas[vtx];
            double c = cosThetas[vtx];

            m_value[vtx] = linearMix(u, v, s, c);
        }

        setDependentsDirty();
    }

    public VectorXs linearMix(in VectorXs u, in VectorXs v, double s, double c)
    {
        return -s * u + c * v;
    }
}

public class ReferenceFrames1 : DependencyNode<List<VectorXs>>
{
    private Tangents m_tangents;
    private List<Vectors> m_previousTangents;
    
    public ReferenceFrames1(ref Tangents tangents) : base(new List<VectorXs>(tangents.size()), 0, tangents.size())
    {
        m_tangents = tangents;
        m_tangents.addDependent(this);
        storeInitialFrames(new Vectors(3));
    }

    public VectorXs this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "ReferenceFrame1"; }

    public void storeInitialFrames(in Vectors initRefFrame1) 
    {
        m_value = new List<VectorXs>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new VectorXs(3));
        }
        List<Vectors> tangents = m_tangents.get();
        Vectors tan0 = tangents[0];
        // Do we have an even approximately valid initial reference frame?
        if (initRefFrame1.squaredNorm() > 0.5 &&
            Math.Abs(initRefFrame1.dot(tan0)) < 0.25) {
            // If so, just project it on the plane normal to the tangent vector
            Vectors projectedInitRefFrame1 =
                (initRefFrame1 - initRefFrame1.dot(tan0) * tan0).normalized();
            m_value[0] = projectedInitRefFrame1.ToVectorXs();
        } else  // If a valid initial first reference frame hasn't been provided, use
                // an arbitrary one
        {
            m_value[0].SetSegment(0, 3, CMath.findNormal(tan0).segment(0, 3));
        }
        // Next initial reference frames are obtained by space-parallel transportation
        // along the rod
        for (int vtx = 1; vtx < size(); ++vtx)
        {
            m_value[vtx] = CMath.orthonormalParallelTransport(
                m_value[vtx - 1].ToVectors(), tangents[vtx - 1], tangents[vtx]).ToVectorXs();
            VectorXs v = m_value[vtx];
            CMath.orthoNormalize(ref v, tangents[vtx]);
        }
        // Store tangents backup for time-parallel transport
        m_previousTangents = tangents;

        setClean();
        setDependentsDirty();
    }

    protected override void compute()
    {
        m_value = new List<VectorXs>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new VectorXs(3));
        }
        List<Vectors> tangents = m_tangents.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            Vectors previousTangent = m_previousTangents[vtx];
            Vectors currentTangent = tangents[vtx];

            m_value[vtx] = CMath.orthonormalParallelTransport(m_value[vtx].ToVectors(), previousTangent,
                                                        currentTangent).ToVectorXs();
            VectorXs v = m_value[vtx];
            CMath.orthoNormalize(ref v, currentTangent);

            // Store the current tangent for the next time-parallel transportation
            previousTangent = currentTangent;
        }

        setDependentsDirty();
    }
}

public class ReferenceFrames2 : DependencyNode<List<VectorXs>>
{
    private Tangents m_tangents;
    private ReferenceFrames1 m_referenceFrames1;

    public ReferenceFrames2(ref Tangents tangents, ref ReferenceFrames1 referenceFrames1)
      : base(new List<VectorXs>(tangents.size()), 1, tangents.size())
    {
        m_tangents = tangents;
        m_referenceFrames1 = referenceFrames1;
        m_tangents.addDependent(this);
        m_referenceFrames1.addDependent(this);
    }

    public VectorXs this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "ReferenceFrame2"; }

    protected override void compute()
    {
        m_value = new List<VectorXs>(new VectorXs[m_size]);
        List<Vectors> tangents = m_tangents.get();
        List<VectorXs> referenceFrames1 = m_referenceFrames1.get();
        for (int vtx = 0; vtx < m_firstValidIndex; ++vtx)
        {
            m_value[vtx].setZero();
        }
        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            m_value[vtx] = tangents[vtx].cross(referenceFrames1[vtx].ToVectors()).ToVectorXs();
        }

        setDependentsDirty();
    }
}

public class ReferenceTwists : DependencyNode<List<double>>
{
    private Tangents m_tangents;
    private ReferenceFrames1 m_referenceFrames1;

    public ReferenceTwists(ref Tangents tangents, ref ReferenceFrames1 referenceFrames1)
        : base(new List<double>(tangents.size()), 1, tangents.size())
    {
        m_tangents = tangents;
        m_referenceFrames1 = referenceFrames1;
        m_tangents.addDependent(this);
        m_referenceFrames1.addDependent(this);
    }

    public double this[int index]
    {
        get => m_value[index];
        set => m_value[index] = value;
    }

    public string name() { return "ReferenceTwists"; }

    protected override void compute()
    {
        m_value = new List<double>(new double[m_size]);
        List<Vectors> tangents = m_tangents.get();
        List<VectorXs> referenceFrames1 = m_referenceFrames1.get();
        for (int vtx = 0; vtx < m_firstValidIndex; ++vtx)
        {
            m_value[vtx] = 0.0;
        }
        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            VectorXs u0 = referenceFrames1[vtx - 1];
            VectorXs u1 = referenceFrames1[vtx];
            Vectors tangent = tangents[vtx];

            // transport reference frame to next edge
            Vectors ut = CMath.orthonormalParallelTransport(u0.ToVectors(), tangents[vtx - 1], tangent);

            // rotate by current value of reference twist
            double beforeTwist = m_value[vtx];
            CMath.rotateAxisAngle(ref ut, tangent, beforeTwist);

            // compute increment to reference twist to align reference frames
            m_value[vtx] = beforeTwist + CMath.signedAngle(ut, u1.ToVectors(), tangent);
        }

        setDependentsDirty();
    }
}
