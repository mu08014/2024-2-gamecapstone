using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

// Unit : cm
public class PhysicalRadius : DependencyNode<double>
{
    public PhysicalRadius(double radius) : base(radius)
    {
        setClean();
    }

    public virtual string name() { return "PhyscialRadius"; }

    protected override void compute() { }
}

// Unit: No dimension
public class BaseRotation : DependencyNode<double>
{
    public BaseRotation(double radius) : base(radius)
    {
        setClean();
    }
    public virtual string name() { return "BaseRotation"; }

    protected override void compute() { }
}

//Unit: cm^4
public class BendingMatrixBase : DependencyNode<MatrixXs>
{
    protected PhysicalRadius m_physicalRadius;
    protected BaseRotation m_baseRotation;

    public BendingMatrixBase(PhysicalRadius rad, BaseRotation baseRotation) : base(new MatrixXs(2, 2))
    {
        m_physicalRadius = rad;
        m_baseRotation = baseRotation;
        m_value.setZero();
        m_physicalRadius.addDependent(this);
        m_baseRotation.addDependent(this);
    }

    public virtual string name() { return "BendingMatrixBase"; }

    protected override void compute()
    {
        double radius = m_physicalRadius.get();
        double baseRotation = m_baseRotation.get();

        MatrixXs B = m_value.Clone();
        B[0, 0] = CMath.PI_4 * radius * CMath.cube(radius);
        B[1, 1] = CMath.PI_4 * radius * CMath.cube(radius);

        MatrixXs rot = new Rotation2D(baseRotation).toRotationMatrix();
        B = (rot * B * rot.transpose()).Clone();
        B[0, 1] = B[1, 0] = 
            0.5 * (B[0, 1] + B[1, 0]);

        setDependentsDirty();
    }
}

// Unit: dPa = g cm^-1 s^-2
public class YoungsModulus : DependencyNode<double>
{
    public YoungsModulus(double youngsModulus) : base(youngsModulus)
    {
        setClean();
    }

    public virtual string name() { return "YoungsModulus"; }

    protected override void compute() { }
}

// Unit: dPa = g cm^-1 s^-2
public class ShearModulus : DependencyNode<double>
{
    public ShearModulus(double shearModulus) : base(shearModulus)
    {
        setClean();
    }

    public virtual string name() { return "ShearModulus"; }

    protected override void compute() { }
}

// Unit : 10^-5 N = g cm s^-2
public class ElasticKs : DependencyNode<double>
{
    protected PhysicalRadius m_physicalRadius;
    protected YoungsModulus m_youngsModulus;

    public ElasticKs(PhysicalRadius rad, YoungsModulus ym) : base(double.NaN)
    {
        m_physicalRadius = rad;
        m_youngsModulus = ym;
        m_physicalRadius.addDependent(this);
        m_youngsModulus.addDependent(this);
    }

    public virtual string name() { return "ElasticKs"; }

    protected override void compute()
    {
        double radius = m_physicalRadius.get();
        double youngsModulus = m_youngsModulus.get();

        m_value = CMath.M_PI * radius * radius * youngsModulus;

        setDependentsDirty();
    }
}

// Unit: 10^-5 cm^2 N = g cm^3 s^-2
public class ElasticKt : DependencyNode<double>
{
    protected PhysicalRadius m_physicalRadius;
    protected ShearModulus m_shearModulus;

    public ElasticKt(PhysicalRadius rad, ShearModulus sm) : base(double.NaN)
    {
        m_physicalRadius = rad;
        m_shearModulus = sm;
        m_physicalRadius.addDependent(this);
        m_shearModulus.addDependent(this);
    }

    public virtual string name() { return "ElasticKt"; }

    protected override void compute()
    {
        double radius = m_physicalRadius.get();
        double shearModulus = m_shearModulus.get();

        m_value = CMath.PI_4 * radius * radius * (radius * radius + radius * radius) * shearModulus;
        setDependentsDirty();
    }
}