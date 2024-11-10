using JetBrains.Annotations;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
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
        get
        {
            return get()[index];
        }
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
        VectorXs dofs = m_dofs.get();
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
        get
        {
            return get()[index];
        }
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
        get
        {
            return get()[index];
        }
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
        get
        {
            return get()[index];
        }
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

        VectorXs thetaVec = thetasMap;
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
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new VectorXs(2));
        }
        m_curvatureBinormals = curvatureBinormals;
        m_materialFrames1 = materialFrames1;
        m_materialFrames2 = materialFrames2;
        m_curvatureBinormals.addDependent(this);
        m_materialFrames1.addDependent(this);
        m_materialFrames2.addDependent(this);
    }

    public VectorXs this[int index]
    {
        get
        {
            return get()[index];
        }
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

public class GradKappas : DependencyNode<List<MatrixXs>>
{
    protected Lengths m_lengths;
    protected Tangents m_tangents;
    protected CurvatureBinormals m_curvatureBinormals;
    protected MaterialFrames1 m_materialFrames1;
    protected MaterialFrames2 m_materialFrames2;
    protected Kappas m_kappas;

    public GradKappas(ref Lengths lengths, ref Tangents tangents,
        ref CurvatureBinormals curvatureBinormals, ref MaterialFrames1 materialFrames1,
        ref MaterialFrames2 materialFrames2, ref Kappas kappas) : base(new List<MatrixXs>(), 1, curvatureBinormals.size())
    {
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new MatrixXs(11, 2));
        }
        m_lengths = lengths;
        m_tangents = tangents;
        m_curvatureBinormals = curvatureBinormals;
        m_materialFrames1 = materialFrames1;
        m_materialFrames2 = materialFrames2;
        m_kappas = kappas;
        m_lengths.addDependent(this);
        m_tangents.addDependent(this);
        m_curvatureBinormals.addDependent(this);
        m_materialFrames1.addDependent(this);
        m_materialFrames2.addDependent(this);
        m_kappas.addDependent(this);
    }

    public MatrixXs this[int index]
    {
        get
        {
            return get()[index];
        }
        set => m_value[index] = value;
    }

    public string name() { return "GradKeppas"; }

    protected override void compute()
    {
        m_value = new List<MatrixXs>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new MatrixXs(11, 2));
        }
        List<double> lengths = m_lengths.get();
        List<Vectors> tangents = m_tangents.get();
        List<Vectors> curvatureBinormals = m_curvatureBinormals.get();
        List<VectorXs> materialFrames1 = m_materialFrames1.get();
        List<VectorXs> materialFrames2 = m_materialFrames2.get();
        List<VectorXs> kappas = m_kappas.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            MatrixXs gradKappa = m_value[vtx];
            gradKappa.setZero();

            double norm_e = lengths[vtx - 1];
            double norm_f = lengths[vtx];

            Vectors te = tangents[vtx - 1];
            Vectors tf = tangents[vtx];

            VectorXs m1e = materialFrames1[vtx - 1];
            VectorXs m2e = materialFrames2[vtx - 1];
            VectorXs m1f = materialFrames1[vtx];
            VectorXs m2f = materialFrames2[vtx];

            double chi = 1.0 + te.dot(tf);

            //    assert( chi>0 );
            if (chi <= 0)
            {
                chi = 1e-12;
            }

            Vectors tilde_t = (te + tf) / chi;
            Vectors tilde_d1 = ((m1e + m1f) / chi).ToVectors();
            Vectors tilde_d2 = ((m2e + m2f) / chi).ToVectors();

            VectorXs kappa = kappas[vtx];

            Vectors Dkappa0De =
                1.0 / norm_e * (-kappa[0] * tilde_t + tf.cross(tilde_d2));
            Vectors Dkappa0Df =
                1.0 / norm_f * (-kappa[0] * tilde_t - te.cross(tilde_d2));
            Vectors Dkappa1De =
                1.0 / norm_e * (-kappa[1] * tilde_t - tf.cross(tilde_d1));
            Vectors Dkappa1Df =
                1.0 / norm_f * (-kappa[1] * tilde_t + te.cross(tilde_d1));

            gradKappa.Setblock(3, 1, 0, 0, (-1) * Dkappa0De);
            gradKappa.Setblock(3, 1, 4, 0, Dkappa0De - Dkappa0Df);
            gradKappa.Setblock(3, 1, 8, 0, Dkappa0Df);
            gradKappa.Setblock(3, 1, 0, 1, (-1) * Dkappa1De);
            gradKappa.Setblock(3, 1, 4, 1, Dkappa1De - Dkappa1Df);
            gradKappa.Setblock(3, 1, 8, 1, Dkappa1Df);

            Vectors kb = curvatureBinormals[vtx];

            gradKappa[3, 0] = -0.5 * kb.dot(m1e.ToVectors());
            gradKappa[7, 0] = -0.5 * kb.dot(m1f.ToVectors());
            gradKappa[3, 1] = -0.5 * kb.dot(m2e.ToVectors());
            gradKappa[7, 1] = -0.5 * kb.dot(m2f.ToVectors());
            m_value[vtx] = gradKappa;
        }

        setDependentsDirty();
    }
}

public class HessKappas : DependencyNode<List<Tuple<MatrixXs, MatrixXs>>>
{
    protected Lengths m_lengths;
    protected Tangents m_tangents;
    protected CurvatureBinormals m_curvatureBinormals;
    protected MaterialFrames1 m_materialFrames1;
    protected MaterialFrames2 m_materialFrames2;
    protected Kappas m_kappas;

    public HessKappas(ref Lengths lengths, ref Tangents tangents,
             ref CurvatureBinormals curvatureBinormals,
             ref MaterialFrames1 materialFrames1,
             ref MaterialFrames2 materialFrames2, ref Kappas kappas) : base(new List<Tuple<MatrixXs, MatrixXs>>(), 1, curvatureBinormals.size())
    {
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new Tuple<MatrixXs, MatrixXs>(new MatrixXs(11, 11), new MatrixXs(11, 11)));
        }
        m_lengths = lengths;
        m_tangents = tangents;
        m_curvatureBinormals = curvatureBinormals;
        m_materialFrames1 = materialFrames1;
        m_materialFrames2 = materialFrames2;
        m_kappas = kappas;
        m_lengths.addDependent(this);
        m_tangents.addDependent(this);
        m_curvatureBinormals.addDependent(this);
        m_materialFrames1.addDependent(this);
        m_materialFrames2.addDependent(this);
        m_kappas.addDependent(this);
    }

    public Tuple<MatrixXs, MatrixXs> this[int index]
    {
        get
        {
            return get()[index];
        }
        set => m_value[index] = value;
    }

    public string name() { return "HessKeppas"; }

    protected override void compute()
    {
        m_value = new List<Tuple<MatrixXs, MatrixXs>>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new Tuple<MatrixXs, MatrixXs>(new MatrixXs(11, 11), new MatrixXs(11, 11)));
        }
        List<double> lengths = m_lengths.get();
        List<Vectors> tangents = m_tangents.get();
        List<Vectors> curvatureBinormals = m_curvatureBinormals.get();
        List<VectorXs> materialFrames1 = m_materialFrames1.get();
        List<VectorXs> materialFrames2 = m_materialFrames2.get();
        List<VectorXs> kappas = m_kappas.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            Tuple<MatrixXs, MatrixXs> HessKappa = m_value[vtx];

            MatrixXs DDkappa1 = HessKappa.Item1;
            MatrixXs DDkappa2 = HessKappa.Item2;

            DDkappa1.setZero();
            DDkappa2.setZero();

            double norm_e = lengths[vtx - 1];
            double norm_f = lengths[vtx];
            double norm2_e = norm_e * norm_e;  // That's bloody stupid, taking the square of a square root.
            double norm2_f = norm_f * norm_f;

            Vectors te = tangents[vtx - 1];
            Vectors tf = tangents[vtx];

            VectorXs m1e = materialFrames1[vtx - 1];
            VectorXs m2e = materialFrames2[vtx - 1];
            VectorXs m1f = materialFrames1[vtx];
            VectorXs m2f = materialFrames2[vtx];

            double chi = 1.0 + te.dot(tf);

            //    assert( chi>0 );
            if (chi <= 0)
            {
                chi = 1e-12;
            }

            Vectors tilde_t = (te + tf) / chi;
            Vectors tilde_d1 = ((m1e + m1f) / chi).ToVectors();
            Vectors tilde_d2 = ((m2e + m2f) / chi).ToVectors();

            VectorXs kappa = kappas[vtx];

            Vectors kb = curvatureBinormals[vtx];

            MatrixXs tt_o_tt = CMath.outerProd(tilde_t, tilde_t);
            MatrixXs tf_c_d2t_o_tt = CMath.outerProd(tf.cross(tilde_d2), tilde_t);
            MatrixXs tt_o_tf_c_d2t = tf_c_d2t_o_tt.transpose();
            MatrixXs kb_o_d2e = CMath.outerProd(kb, m2e.ToVectors());
            MatrixXs d2e_o_kb = kb_o_d2e.transpose();

            MatrixXs Id = CMath.identity(3);

            {
                MatrixXs D2kappa1De2 =
                    1.0 / norm2_e *
                        (2 * kappa[0] * tt_o_tt - (tf_c_d2t_o_tt + tt_o_tf_c_d2t)) -
                    kappa[0] / (chi * norm2_e) * (Id - CMath.outerProd(te, te)) +
                    1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

                MatrixXs te_c_d2t_o_tt = CMath.outerProd(te.cross(tilde_d2), tilde_t);
                MatrixXs tt_o_te_c_d2t = te_c_d2t_o_tt.transpose();
                MatrixXs kb_o_d2f = CMath.outerProd(kb, m2f.ToVectors());
                MatrixXs d2f_o_kb = kb_o_d2f.transpose();

                MatrixXs D2kappa1Df2 =
                    1.0 / norm2_f *
                        (2 * kappa[0] * tt_o_tt + (te_c_d2t_o_tt + tt_o_te_c_d2t)) -
                    kappa[0] / (chi * norm2_f) * (Id - CMath.outerProd(tf, tf)) +
                    1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

                MatrixXs D2kappa1DeDf =
                    -kappa[0] / (chi * norm_e * norm_f) * (Id + CMath.outerProd(te, tf)) +
                    1.0 / (norm_e * norm_f) *
                        (2 * kappa[0] * tt_o_tt - tf_c_d2t_o_tt + tt_o_te_c_d2t -
                         CMath.crossMat(tilde_d2));
                MatrixXs D2kappa1DfDe = D2kappa1DeDf.transpose();

                double D2kappa1Dthetae2 = -0.5 * kb.dot(m2e.ToVectors());
                double D2kappa1Dthetaf2 = -0.5 * kb.dot(m2f.ToVectors());
                Vectors D2kappa1DeDthetae =
                    1.0 / norm_e *
                    (0.5 * kb.dot(m1e.ToVectors()) * tilde_t - 1.0 / chi * tf.cross(m1e.ToVectors()));
                Vectors D2kappa1DeDthetaf =
                    1.0 / norm_e *
                    (0.5 * kb.dot(m1f.ToVectors()) * tilde_t - 1.0 / chi * tf.cross(m1f.ToVectors()));
                Vectors D2kappa1DfDthetae =
                    1.0 / norm_f *
                    (0.5 * kb.dot(m1e.ToVectors()) * tilde_t + 1.0 / chi * te.cross(m1e.ToVectors()));
                Vectors D2kappa1DfDthetaf =
                    1.0 / norm_f *
                    (0.5 * kb.dot(m1f.ToVectors()) * tilde_t + 1.0 / chi * te.cross(m1f.ToVectors()));

                DDkappa1.Setblock(3, 3, 0, 0, D2kappa1De2);
                DDkappa1.Setblock(3, 3, 0, 4, (-1) * D2kappa1De2 + D2kappa1DeDf);
                DDkappa1.Setblock(3, 3, 4, 0, (-1) * D2kappa1De2 + D2kappa1DfDe);
                DDkappa1.Setblock(3, 3, 4, 4,
                    D2kappa1De2 - (D2kappa1DeDf + D2kappa1DfDe) + D2kappa1Df2);
                DDkappa1.Setblock(3, 3, 0, 8, -1 * D2kappa1DeDf);
                DDkappa1.Setblock(3, 3, 8, 0, -1 * D2kappa1DfDe);
                DDkappa1.Setblock(3, 3, 4, 8, D2kappa1DeDf - D2kappa1Df2);
                DDkappa1.Setblock(3, 3, 8, 4, D2kappa1DfDe - D2kappa1Df2);
                DDkappa1.Setblock(3, 3, 8, 8, D2kappa1Df2);
                DDkappa1[3, 3] = D2kappa1Dthetae2;
                DDkappa1[7, 7] = D2kappa1Dthetaf2;
                DDkappa1[3, 7] = DDkappa1[7, 3] = 0;
                DDkappa1.Setblock(3, 1, 0, 3, -1 * D2kappa1DeDthetae);
                DDkappa1.Setblock(1, 3, 3, 0, DDkappa1.block(3, 1, 0, 3).transpose());
                DDkappa1.Setblock(3, 1, 4, 3, D2kappa1DeDthetae - D2kappa1DfDthetae);
                DDkappa1.Setblock(1, 3, 3, 4, DDkappa1.block(3, 1, 4, 3).transpose());
                DDkappa1.Setblock(3, 1, 8, 3, D2kappa1DfDthetae);
                DDkappa1.Setblock(1, 3, 3, 8, DDkappa1.block(3, 1, 8, 3).transpose());
                DDkappa1.Setblock(3, 1, 0, 7, -1 * D2kappa1DeDthetaf);
                DDkappa1.Setblock(1, 3, 7, 0, DDkappa1.block(3, 1, 0, 7).transpose());
                DDkappa1.Setblock(3, 1, 4, 7, D2kappa1DeDthetaf - D2kappa1DfDthetaf);
                DDkappa1.Setblock(1, 3, 7, 4, DDkappa1.block(3, 1, 4, 7).transpose());
                DDkappa1.Setblock(3, 1, 8, 7, D2kappa1DfDthetaf);
                DDkappa1.Setblock(1, 3, 7, 8, DDkappa1.block(3, 1, 8, 7).transpose());
            }

            {
                MatrixXs tf_c_d1t_o_tt = CMath.outerProd(tf.cross(tilde_d1), tilde_t);
                MatrixXs tt_o_tf_c_d1t = tf_c_d1t_o_tt.transpose();
                MatrixXs kb_o_d1e = CMath.outerProd(kb, m1e.ToVectors());
                MatrixXs d1e_o_kb = kb_o_d1e.transpose();

                MatrixXs D2kappa2De2 =
                    1.0 / norm2_e *
                        (2 * kappa[1] * tt_o_tt + (tf_c_d1t_o_tt + tt_o_tf_c_d1t)) -
                    kappa[1] / (chi * norm2_e) * (Id - CMath.outerProd(te, te)) -
                    1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);

                MatrixXs te_c_d1t_o_tt = CMath.outerProd(te.cross(tilde_d1), tilde_t);
                MatrixXs tt_o_te_c_d1t = te_c_d1t_o_tt.transpose();
                MatrixXs kb_o_d1f = CMath.outerProd(kb, m1f.ToVectors());
                MatrixXs d1f_o_kb = kb_o_d1f.transpose();

                MatrixXs D2kappa2Df2 =
                    1.0 / norm2_f *
                        (2 * kappa[1] * tt_o_tt - (te_c_d1t_o_tt + tt_o_te_c_d1t)) -
                    kappa[1] / (chi * norm2_f) * (Id - CMath.outerProd(tf, tf)) -
                    1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb);

                MatrixXs D2kappa2DeDf =
                    -kappa[1] / (chi * norm_e * norm_f) * (Id + CMath.outerProd(te, tf)) +
                    1.0 / (norm_e * norm_f) *
                        (2 * kappa[1] * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t +
                         CMath.crossMat(tilde_d1));
                MatrixXs D2kappa2DfDe = D2kappa2DeDf.transpose();

                double D2kappa2Dthetae2 = 0.5 * kb.dot(m1e.ToVectors());
                double D2kappa2Dthetaf2 = 0.5 * kb.dot(m1f.ToVectors());
                Vectors D2kappa2DeDthetae =
                    1.0 / norm_e *
                    (0.5 * kb.dot(m2e.ToVectors()) * tilde_t - 1.0 / chi * tf.cross(m2e.ToVectors()));
                Vectors D2kappa2DeDthetaf =
                    1.0 / norm_e *
                    (0.5 * kb.dot(m2f.ToVectors()) * tilde_t - 1.0 / chi * tf.cross(m2f.ToVectors()));
                Vectors D2kappa2DfDthetae =
                    1.0 / norm_f *
                    (0.5 * kb.dot(m2e.ToVectors()) * tilde_t + 1.0 / chi * te.cross(m2e.ToVectors()));
                Vectors D2kappa2DfDthetaf =
                    1.0 / norm_f *
                    (0.5 * kb.dot(m2f.ToVectors()) * tilde_t + 1.0 / chi * te.cross(m2f.ToVectors()));

                DDkappa2.Setblock(3, 3, 0, 0, D2kappa2De2);
                DDkappa2.Setblock(3, 3, 0, 4, (-1) * D2kappa2De2 + D2kappa2DeDf);
                DDkappa2.Setblock(3, 3, 4, 0, (-1) * D2kappa2De2 + D2kappa2DfDe);
                DDkappa2.Setblock(3, 3, 4, 4,
                    D2kappa2De2 - (D2kappa2DeDf + D2kappa2DfDe) + D2kappa2Df2);
                DDkappa2.Setblock(3, 3, 0, 8, -1 * D2kappa2DeDf);
                DDkappa2.Setblock(3, 3, 8, 0, -1 * D2kappa2DfDe);
                DDkappa2.Setblock(3, 3, 4, 8, D2kappa2DeDf - D2kappa2Df2);
                DDkappa2.Setblock(3, 3, 8, 4, D2kappa2DfDe - D2kappa2Df2);
                DDkappa2.Setblock(3, 3, 8, 8, D2kappa2Df2);
                DDkappa2[3, 3] = D2kappa2Dthetae2;
                DDkappa2[7, 7] = D2kappa2Dthetaf2;
                DDkappa2[3, 7] = DDkappa2[7, 3] = 0;
                DDkappa2.Setblock(3, 1, 0, 3, -1 * D2kappa2DeDthetae);
                DDkappa2.Setblock(1, 3, 3, 0, DDkappa2.block(3, 1, 0, 3).transpose());
                DDkappa2.Setblock(3, 1, 4, 3, D2kappa2DeDthetae - D2kappa2DfDthetae);
                DDkappa2.Setblock(1, 3, 3, 4, DDkappa2.block(3, 1, 4, 3).transpose());
                DDkappa2.Setblock(3, 1, 8, 3, D2kappa2DfDthetae);
                DDkappa2.Setblock(1, 3, 3, 8, DDkappa2.block(3, 1, 8, 3).transpose());
                DDkappa2.Setblock(3, 1, 0, 7, -1 * D2kappa2DeDthetaf);
                DDkappa2.Setblock(1, 3, 7, 0, DDkappa2.block(3, 1, 0, 7).transpose());
                DDkappa2.Setblock(3, 1, 4, 7, D2kappa2DeDthetaf - D2kappa2DfDthetaf);
                DDkappa2.Setblock(1, 3, 7, 4, DDkappa2.block(3, 1, 4, 7).transpose());
                DDkappa2.Setblock(3, 1, 8, 7, D2kappa2DfDthetaf);
                DDkappa2.Setblock(1, 3, 7, 8, DDkappa2.block(3, 1, 8, 7).transpose());
            }

            HessKappa = new Tuple<MatrixXs, MatrixXs>(DDkappa1, DDkappa2);
            m_value[vtx] = HessKappa;
        }

        setDependentsDirty();
    }
}

public class Twists : DependencyNode<List<double>>
{
    private ReferenceTwists m_refTwists;
    private DOFs m_dofs;

    public Twists(ref ReferenceTwists refTwists, ref DOFs dofs) : base(new List<double>(new double[dofs.getNumEdges()]), 1, dofs.getNumEdges())
    {
        m_refTwists = refTwists;
        m_dofs = dofs;
        m_refTwists.addDependent(this);
        m_dofs.addDependent(this);
    }

    public double this[int index]
    {
        get
        {
            return get()[index];
        }
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

public class GradTwists : DependencyNode<List<Vectors>>
{
    protected Lengths m_lengths;
    protected CurvatureBinormals m_curvatureBinormals;

    public GradTwists(ref Lengths lengths, ref CurvatureBinormals curvatureBinormals) : base(new List<Vectors>(), 1, lengths.size())
    {
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new Vectors(11));
        }
        m_lengths = lengths;
        m_curvatureBinormals = curvatureBinormals;
        m_lengths.addDependent(this);
        m_curvatureBinormals.addDependent(this);
    }

    public Vectors this[int index]
    {
        get
        {
            return get()[index];
        }
        set => m_value[index] = value;
    }

    public string name() { return "GradTwists"; }

    protected override void compute()
    {
        m_value = new List<Vectors>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new Vectors(11));
        }

        List<Vectors> curvatureBinormals = m_curvatureBinormals.get();
        List<double> lengths = m_lengths.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            VectorXs Dtwist = m_value[vtx].ToVectorXs();
            Dtwist.setZero();

            Vectors kb = curvatureBinormals[vtx];

            Dtwist.SetSegment(0, 3, -0.5 / lengths[vtx - 1] * kb);
            Dtwist.SetSegment(8, 3, 0.5 / lengths[vtx] * kb);
            Dtwist.SetSegment(4, 3, (-1) * (Dtwist.segment(0, 3) + Dtwist.segment(8, 3)));
            Dtwist[3] = -1;
            Dtwist[7] = 1;
            m_value[vtx] = Dtwist.ToVectors();
        }

        setDependentsDirty();
    }
}

public class HessTwists : DependencyNode<List<MatrixXs>>
{
    protected Tangents m_tangents;
    protected Lengths m_lengths;
    protected CurvatureBinormals m_curvatureBinormals;

    public HessTwists(ref Tangents tangents, ref Lengths lengths,
             ref CurvatureBinormals curvatureBinormals) : base(new List<MatrixXs>(), 1, lengths.size())
    {
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new MatrixXs(11, 11));
        }
        m_tangents = tangents;
        m_lengths = lengths;
        m_curvatureBinormals = curvatureBinormals;
        m_tangents.addDependent(this);
        m_lengths.addDependent(this);
        m_curvatureBinormals.addDependent(this);
    }

    public MatrixXs this[int index]
    {
        get
        {
            return get()[index];
        }
        set => m_value[index] = value;
    }

    public string name() { return "HessTwists"; }

    protected override void compute()
    {
        m_value = new List<MatrixXs>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new MatrixXs(11, 11));
        }
        List<Vectors> tangents = m_tangents.get();
        List<double> lengths = m_lengths.get();
        List<Vectors> curvatureBinormals = m_curvatureBinormals.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            MatrixXs DDtwist = m_value[vtx];

            DDtwist.setZero();

            Vectors te = tangents[vtx - 1];
            Vectors tf = tangents[vtx];
            double norm_e = lengths[vtx - 1];
            double norm_f = lengths[vtx];
            Vectors kb = curvatureBinormals[vtx];

            double chi = 1 + te.dot(tf);

            //    assert( chi>0 );
            if (chi <= 0)
            {
                chi = 1e-12;
            }

            Vectors tilde_t = 1.0 / chi * (te + tf);

            MatrixXs D2mDe2 =
                -0.25 / (norm_e * norm_e) *
                (CMath.outerProd(kb, te + tilde_t) + CMath.outerProd(te + tilde_t, kb));
            MatrixXs D2mDf2 =
                -0.25 / (norm_f * norm_f) *
                (CMath.outerProd(kb, tf + tilde_t) + CMath.outerProd(tf + tilde_t, kb));
            MatrixXs D2mDeDf =
                0.5 / (norm_e * norm_f) *
                (2.0 / chi * CMath.crossMat(te) - CMath.outerProd(kb, tilde_t));
            MatrixXs D2mDfDe = D2mDeDf.transpose();

            DDtwist.Setblock(3, 3, 0, 0, D2mDe2);
            DDtwist.Setblock(3, 3, 0, 4, -1 * D2mDe2 + D2mDeDf);
            DDtwist.Setblock(3, 3, 4, 0, -1 * D2mDe2 + D2mDfDe);
            DDtwist.Setblock(3, 3, 4, 4, D2mDe2 - (D2mDeDf + D2mDfDe) + D2mDf2);
            DDtwist.Setblock(3, 3, 0, 8, -1 * D2mDeDf);
            DDtwist.Setblock(3, 3, 8, 0, -1 * D2mDfDe);
            DDtwist.Setblock(3, 3, 8, 4, D2mDfDe - D2mDf2);
            DDtwist.Setblock(3, 3, 4, 8, D2mDeDf - D2mDf2);
            DDtwist.Setblock(3, 3, 8, 8, D2mDf2);

            m_value[vtx] = DDtwist;
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
        get
        {
            return get()[index];
        }
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
        get
        {
            return get()[index];
        }
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
        get
        {
            return get()[index];
        }
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
            m_value[vtx] = v;

            // Store the current tangent for the next time-parallel transportation
            m_previousTangents[vtx] = currentTangent;
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
        get
        {
            return get()[index];
        }
        set => m_value[index] = value;
    }

    public string name() { return "ReferenceFrame2"; }

    protected override void compute()
    {
        m_value = new List<VectorXs>(m_size);
        for (int vtx = 0; vtx < m_size; ++vtx)
        {
            m_value.Add(new VectorXs(3));
        }
        List<Vectors> tangents = m_tangents.get();
        List<VectorXs> referenceFrames1 = m_referenceFrames1.get();

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
        get
        {
            return get()[index];
        }
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

public class BendingProducts : DependencyNode<List<MatrixXs>>
{
    private BendingMatrixBase m_bendingMatrixBase;
    private GradKappas m_gradKappas;


    public BendingProducts(ref BendingMatrixBase bendingMatrixBase, ref GradKappas gradKappas) : base(1, gradKappas.size())
    {
        m_bendingMatrixBase = bendingMatrixBase;
        m_gradKappas = gradKappas;
        m_bendingMatrixBase.addDependent(this);
        m_gradKappas.addDependent(this);

        compute();
    }

    public MatrixXs this[int index]
    {
        get
        {
            return get()[index];
        }
        set => m_value[index] = value;
    }

    public string name() { return "BendingProducts"; }

    protected override void compute()
    {
        m_value = new List<MatrixXs>(m_size);
        for (int i = 0; i < m_size; i++)
        {
            m_value.Add(new MatrixXs(11, 11));
        }

        MatrixXs bendingMatrix = m_bendingMatrixBase.get();
        List<MatrixXs> gradKappas = m_gradKappas.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            CMath.symBProduct(11, m_value[vtx], bendingMatrix, gradKappas[vtx]);
        }

        setDependentsDirty();
    }
}

public class GradTwistsSquared : DependencyNode<List<MatrixXs>>
{
    private GradTwists m_gradTwists;
    public GradTwistsSquared(ref GradTwists gradTwists) : base(1, gradTwists.size())
    {
        m_gradTwists = gradTwists;
        m_gradTwists.addDependent(this);

        compute();
    }

    public string name() { return "GradTwistsSquared"; }

    protected override void compute()
    {
        m_value = new List<MatrixXs>(m_size);
        for (int i = 0; i < m_size; ++i)
        {
            m_value.Add(new MatrixXs(11, 11));
        }

        List<Vectors> gradTwists = m_gradTwists.get();

        for (int vtx = m_firstValidIndex; vtx < size(); ++vtx)
        {
            Vectors gradTwist = gradTwists[vtx];
            m_value[vtx] = CMath.matProduct(gradTwist, gradTwist);
        }

        setDependentsDirty();
    }
}
