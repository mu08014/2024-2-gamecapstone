using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class StrandParameters
{
    public double m_density;
    public double m_viscosity;

    public double m_viscousBendingCoefficientBase;
    public double m_viscousKit;
    public double m_viscousKs;

    public PhysicalRadius m_physicalRadius;
    public BaseRotation m_baseRotation;
    public BendingMatrixBase m_bendingMatrixBase;
    public YoungsModulus m_youngsModulus;
    public ShearModulus m_shearModulus;
    public double m_stretchingMultiplier;
    public ElasticKs m_ks;
    public ElasticKt m_kt;

    StrandParameters(double radius, double YoungsModulus, double shearModulus,
        double StretchingMultiplier, double density,
        double viscosity, double baseRotation, double dt, Vector3 color,
        bool accumViscous = true, bool accumViscousBend = true,
        bool variableRadiusHair = false, double straightHairs = 1
        )
    {

    }
}

public class StrandEquilibriumParameters
{

}
