using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class StrandEquilibriumParameters
{
    public List<VectorXs> m_vertices;
    public double m_curl_radius;
    public double m_curl_density;
    public double m_dL;
    public double m_root_length;
    public bool m_valid;
    public bool m_dirty;

    StrandEquilibriumParameters(List<VectorXs> vertices, double curl_radius, double curl_density, double dL, double root_length, bool valid)
    {
        m_vertices = vertices;
        m_curl_radius = curl_radius;
        m_curl_density = curl_density;
        m_dL = dL;
        m_root_length = root_length;
        m_valid = valid;
        m_dirty = false;
    }
}

public class StrandParameters
{
    public double m_density;
    public double m_viscosity;

    public double m_viscousBendingCoefficientBase;
    public double m_viscousKt;
    public double m_viscousKs;

    public PhysicalRadius m_physicalRadius;
    public BaseRotation m_baseRotation;
    public BendingMatrixBase m_bendingMatrixBase;
    public YoungsModulus m_youngsModulus;
    public ShearModulus m_shearModulus;
    public double m_stretchingMultiplier;
    public ElasticKs m_ks;
    public ElasticKt m_kt;

    public Vector3 m_color;

    public bool m_accumulateWithViscous;
    public bool m_accumulateVisCousOnlyForBendingModes;
    public bool m_variableRadiusHair;
    public double m_straightHairs;

    StrandParameters(double radius, double YoungsModulus, double shearModulus,
        double stretchingMultiplier, double density,
        double viscosity, double baseRotation, double dt, Vector3 color,
        bool accumViscous = true, bool accumViscousBend = true,
        bool variableRadiusHair = false, double straightHairs = 1
        )
    {
        m_density = density;
        m_viscosity = viscosity;
        m_physicalRadius = new PhysicalRadius(radius);
        m_baseRotation = new BaseRotation(baseRotation);
        m_bendingMatrixBase = new BendingMatrixBase(m_physicalRadius, m_baseRotation);
        m_youngsModulus = new YoungsModulus(YoungsModulus);
        m_shearModulus = new ShearModulus(shearModulus);
        m_stretchingMultiplier = stretchingMultiplier;
        m_ks = new ElasticKs(m_physicalRadius, m_youngsModulus);
        m_kt = new ElasticKt(m_physicalRadius, m_shearModulus);
        m_accumulateWithViscous = accumViscous;
        m_accumulateVisCousOnlyForBendingModes = accumViscousBend;
        m_variableRadiusHair = variableRadiusHair;
        m_straightHairs = straightHairs;
        m_color = color;
    }

    public double interpolatedRadiusMultiplier(int vtx, int numVertices)
    {
        if (m_variableRadiusHair)
        {
            double s = vtx / (numVertices - 1);
            return (Math.Exp(-3.4612 * s)) * m_straightHairs +
                (1 - m_straightHairs);
        }
        return 1;
    }

    public Vector3 getColor() { return m_color; }

    public double getKs(int vtx, int numVertices)
    {
        double interpol = interpolatedRadiusMultiplier(vtx, numVertices);
        return interpol * interpol * m_ks.get() * m_stretchingMultiplier;
    }

    public double getKt(int vtx, int numVertices)
    {
        double interpol = interpolatedRadiusMultiplier(vtx, numVertices);
        double interpol2 = interpol * interpol;
        return interpol2 * interpol2 * m_kt.get();
    }

    public double getRadius(int vtx, int numVertices)
    {
        return interpolatedRadiusMultiplier(vtx, numVertices) *
            m_physicalRadius.get();
    }

    public double bendingCoefficient(int vtx, int numVertices)
    {
        double interpol = interpolatedRadiusMultiplier(vtx, numVertices);
        double interpol2 = interpol * interpol;

        return interpol2 * interpol2 * m_youngsModulus.get();
    }

    public BendingMatrixBase getBendingMatrixBase() { return m_bendingMatrixBase; }

    public MatrixXs bendingMatrix(int vtx, int numVertices)
    {
        return bendingCoefficient(vtx, numVertices) * m_bendingMatrixBase.get();
    }

    public void computeViscousForceCoefficients(double dt)
    {
        double m_radius = m_physicalRadius.get();

        m_viscousKs = CMath.M_PI * m_radius * m_radius * 3 * m_viscosity / dt;
        m_viscousKt = CMath.M_PI_4 * m_radius * m_radius *
            (m_radius * m_radius + m_radius * m_radius) * m_viscosity /
            dt;
        m_viscousBendingCoefficientBase = 3 * m_viscosity / dt;
    }

    public double viscousBendingCoefficient(int vtx, int numVertices)
    {
        double interpol = interpolatedRadiusMultiplier(vtx, numVertices);
        double interpol2 = interpol * interpol;

        return interpol2 * interpol2 * m_viscousBendingCoefficientBase;
    }

    public MatrixXs viscousBendingMatrix(int vtx, int numVertices)
    {
        return viscousBendingCoefficient(vtx, numVertices) *
            m_bendingMatrixBase.get();
    }

    public double getViscousKs(int vtx, int numVertices)
    {
        double interpol = interpolatedRadiusMultiplier(vtx, numVertices);

        return interpol * interpol * m_viscousKs;
    }

    public double getViscousKt(int vtx, int numVertices)
    {
        double interpol = interpolatedRadiusMultiplier(vtx, numVertices);
        double interpol2 = interpol * interpol;

        return interpol2 * interpol2 * m_viscousKt;
    }
}
