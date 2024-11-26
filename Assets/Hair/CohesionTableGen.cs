using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

public class CohesionTable
{
    private double m_sigma;
    private double m_theta;
    private double m_radii;
    private double m_max_alpha;
    private double m_max_d0;

    private double m_min_d0;
    private double m_min_d0_planar;

    private const double m_ang_epsilon = 0.008;

    private int m_discretization;

    private MatrixXs m_A_table;
    private MatrixXs m_alpha_table;
    private MatrixXs m_dEdd_table;

    private MatrixXs m_A_planar_table;
    private MatrixXs m_alpha_planar_table;
    private MatrixXs m_dEdd_planar_table;

    private Vectors m_max_As;
    private Vectors m_min_As;

    private Vectors m_max_A_planars;
    private Vectors m_min_A_planars;

    private double m_radius_multiplier;
    private double m_radius_multiplier_planar;

    private double m_collision_stiffness;
    private double m_collision_stiffness_planar;

    public CohesionTable(double radius_multiplier,
        double collision_stiffness,
        double radius_multiplier_planar,
        double collision_stiffness_planar)
    {
        m_sigma = 72.0;
        m_theta = CMath.M_PI / 4;
        m_radii = 0.004;
        m_max_alpha = CMath.M_PI - m_theta;
        m_max_d0 = 0.2;
        m_min_d0 = 2.0 * m_radii;
        m_min_d0_planar = m_radii;
        m_radius_multiplier = radius_multiplier;
        m_collision_stiffness = collision_stiffness;
        m_collision_stiffness_planar = collision_stiffness_planar;
        m_radius_multiplier_planar = radius_multiplier_planar;
        m_discretization = 256;
    }

    public double interpolate_table(in double A, in double d0,
        in MatrixXs mat,
        in double dmin)
    {
        double d_inc = m_max_d0 / (double)m_discretization;
        double p = Math.Max(0, d0 - dmin) / d_inc;
        double fp = Math.Min(p - Math.Floor(p), 1.0);
        int ip0 = Math.Max(0, Math.Min(m_discretization - 2, (int)p));
        int ip1 = ip0 + 1;

        double a_inc0 = (m_max_As[ip0] - m_min_As[ip0]) / m_discretization;
        double a_inc1 = (m_max_As[ip1] - m_min_As[ip1]) / m_discretization;

        if (a_inc0 == 0.0 || a_inc1 == 0.0)
            return 0.0;

        double q0 = (Math.Max(m_min_As[ip0], Math.Min(m_max_As[ip0], A)) - m_min_As[ip0]) / a_inc0;
        double q1 = (Math.Max(m_min_As[ip1], Math.Min(m_max_As[ip1], A)) - m_min_As[ip1]) / a_inc1;

        double fq0 = Math.Min(q0 - Math.Floor(q0), 1.0);
        double fq1 = Math.Min(q1 - Math.Floor(q1), 1.0);

        int iq00 = Math.Max(0, Math.Min(m_discretization - 2, (int)q0));
        int iq10 = Math.Max(0, Math.Min(m_discretization - 2, (int)q1));
        int iq01 = iq00 + 1;
        int iq11 = iq10 + 1;

        double dEdd0;
        double dEdd1;

        double v00 = mat[iq00, ip0];
        double v01 = mat[iq01, ip0];
        double v10 = mat[iq10, ip1];
        double v11 = mat[iq11, ip1];

        dEdd0 = CMath.lerp(v00, v01, fq0);
        dEdd1 = CMath.lerp(v10, v11, fq1);

        return CMath.lerp(dEdd0, dEdd1, fp);
    }

    public double interpolate_table_grad(in double A, in double d0,
        in MatrixXs mat, in double dmin)
    {
        double d_inc = m_max_d0 / (double)m_discretization;
        double p = Math.Max(0, d0 - dmin) / d_inc;
        int ip0 = Math.Max(0, Math.Min(m_discretization - 2, (int)p));
        int ip1 = ip0 + 1;

        double a_inc0 = (m_max_As[ip0] - m_min_As[ip0]) / m_discretization;
        double a_inc1 = (m_max_As[ip1] - m_min_As[ip1]) / m_discretization;

        if (a_inc0 == 0.0 || a_inc1 == 0.0)
            return 0.0;

        double q0 = (Math.Max(m_min_As[ip0], Math.Min(m_max_As[ip0], A)) - m_min_As[ip0]) / a_inc0;
        double q1 = (Math.Max(m_min_As[ip1], Math.Min(m_max_As[ip1], A)) - m_min_As[ip1]) / a_inc1;

        double fq0 = Math.Min(q0 - Math.Floor(q0), 1.0);
        double fq1 = Math.Min(q1 - Math.Floor(q1), 1.0);

        int iq00 = Math.Max(0, Math.Min(m_discretization - 2, (int)q0));
        int iq10 = Math.Max(0, Math.Min(m_discretization - 2, (int)q1));
        int iq01 = iq00 + 1;
        int iq11 = iq10 + 1;

        double dEdd0;
        double dEdd1;

        double v00 = mat[iq00, ip0];
        double v01 = mat[iq01, ip0];
        double v10 = mat[iq10, ip1];
        double v11 = mat[iq11, ip1];

        dEdd0 = CMath.lerp(v00, v01, fq0);
        dEdd1 = CMath.lerp(v10, v11, fq1);

        return CMath.lerp(dEdd0, dEdd1, d_inc);
    }

    public double interpolate_dEdd(in double A, in double d0)
    {
        return interpolate_table(A, d0, m_dEdd_table, 0.0);
    }

    public double interpolate_d2Edd2(in double A, in double d0)
    {
        return interpolate_table_grad(A, d0, m_dEdd_table, 0.0);
    }

    public double interpolate_alpha(in double A, in double d0)
    {
        return interpolate_table(A, d0, m_alpha_table, m_min_d0);
    }

    public double interpolate_dEdd_planar(in double A, in double d0)
    {
        return interpolate_table(A, d0, m_dEdd_planar_table, 0.0);
    }

    public double interpolate_d2Edd2_planar(in double A, in double d0)
    {
        return interpolate_table_grad(A, d0, m_dEdd_planar_table, 0.0);
    }

    public double interpolate_alpha_planar(in double A, in double d0)
    {
        return interpolate_table(A, d0, m_alpha_planar_table, m_min_d0_planar);
    }

    public void setParameter(in double sigma, in double theta, in double radii, in double max_d0, in int disc)
    {
        m_sigma = sigma;
        m_theta = theta;
        m_radii = radii;
        m_max_d0 = max_d0;
        m_discretization = disc;
        m_min_d0 = radii * 2.0;
        m_max_alpha = CMath.M_PI - m_theta;
    }

    public double getRadii()
    {
        return m_radii;
    }

    public double computeH(in double R, in double alpha)
    {
        return m_radii * Math.Sin(alpha) - R * (1.0 - Math.Sin(m_theta + alpha));
    }

    public double computeR(in double alpha, in double d0)
    {
        return (d0 - 2.0 * m_radii * Math.Cos(alpha)) / (2.0 * Math.Cos(m_theta + alpha));
    }

    public double computeA(in double R, in double alpha)
    {
        return 2.0 * R * R *
            (alpha + m_theta - CMath.M_PI / 2 + 0.5 * Math.Sin(2.0 * (alpha + m_theta))) +
            2.0 * m_radii * R * (Math.Sin(2.0 * alpha + m_theta) - Math.Sin(m_theta)) -
            m_radii * m_radii * (2.0 * alpha - Math.Sin(2.0 * alpha));
    }

    public double computeApproxA(in double alpha, in double d0)
    {
        double gamma = alpha + m_theta;
        double t2 = m_radii * m_radii;
        double t3 = d0 * d0;
        double t4 = Math.Cos(m_theta);
        double t5 = t4 * t4;
        double t6 = Math.Sin(m_theta);
        return -t2 * Math.Sin(m_theta * 2.0) + t2 * CMath.M_PI * (1.0 / 3.0) -
         t3 * CMath.M_PI * (1.0 / 6.0) - gamma * t2 * (8.0 / 3.0) +
         gamma * t3 * (1.0 / 3.0) + t2 * m_theta * 2.0 -
         t2 * t5 * CMath.M_PI * (4.0 / 3.0) + d0 * m_radii * t4 * 2.0 +
         gamma * t2 * t5 * (8.0 / 3.0) +
         d0 * gamma * m_radii * t6 * (2.0 / 3.0) -
         d0 * m_radii * t6 * CMath.M_PI * (1.0 / 3.0);
    }

    public double computeApproxdEdd(in double alpha, in double d0)
    {
        double gamma = alpha + m_theta;

        double t2 = Math.Sin(m_theta);
        double t3 = m_radii * m_radii;
        double t4 = Math.Cos(m_theta);
        double t5 = d0 * d0;
        double t6 = d0 * m_radii * t2 * 2.0;
        double t7 = m_theta * 2.0;
        double t8 = Math.Sin(t7);
        return (m_sigma *
                (t3 * -8.0 + t5 + t6 + t3 * (t4 * t4) * 8.0 + t3 * t8 * CMath.M_PI * 2.0 -
                 gamma * t3 * t8 * 4.0 - d0 * gamma * m_radii * t4 * 2.0 +
                 d0 * m_radii * t4 * CMath.M_PI) *
                2.0) /
               (t5 + t6 - (t2 * t2) * t3 * 8.0);
    }

    public double computedEddPlanar(in double R, in double alpha)
    {
        if (R == 0.0)
        {
            return 0.0;
        }
        else
        {
            double t2 = m_theta * 2.0;
            double t3 = Math.Sin(alpha);
            double t4 = alpha + m_theta;
            double t5 = Math.Sin(t4);
            double t6 = alpha + t2;
            double t7 = m_radii * m_radii;
            double t8 = Math.Sin(t6);
            double t9 = R * R;
            double t10 = alpha * 2.0;
            double t11 = Math.Sin(t2);
            double t12 = t2 + t10;
            double t13 = Math.Sin(t12);
            double t14 = m_theta * 3.0;
            double t15 = alpha + t14;
            double t16 = Math.Cos(t10);
            double t17 = Math.Cos(t12);
            double t18 = Math.Cos(m_theta);
            double t19 = t10 + m_theta;
            double t20 = Math.Cos(t19);
            return (m_sigma *
                    (t7 * CMath.M_PI * -2.0 - t9 * CMath.M_PI * 2.0 + alpha * t7 * 2.0 +
                     alpha * t9 * 2.0 + t3 * t7 * 4.0 + t3 * t9 * 4.0 + t7 * t8 * 2.0 +
                     t8 * t9 * 4.0 - t7 * t11 * 2.0 + t7 * t13 * 2.0 - t9 * t11 +
                     t9 * t13 * 3.0 + t7 * m_theta * 4.0 + t9 * m_theta * 4.0 +
                     t7 * Math.Sin(t10) * 2.0 + t7 * Math.Sin(alpha - t2) * 2.0 +
                     t9 * Math.Sin(t10 + m_theta * 4.0) - R * m_radii * Math.Sin(t14) +
                     R * m_radii * Math.Sin(t15) * 2.0 + R * m_radii * Math.Sin(t19) * 5.0 -
                     R * m_radii * Math.Sin(m_theta) * 3.0 +
                     R * m_radii * Math.Sin(alpha - m_theta) * 6.0 + t7 * t16 * CMath.M_PI * 2.0 +
                     t9 * t17 * CMath.M_PI * 2.0 + R * m_radii * t5 * 8.0 -
                     alpha * t7 * t16 * 2.0 - alpha * t9 * t17 * 2.0 -
                     t7 * t16 * m_theta * 4.0 - t9 * t17 * m_theta * 4.0 +
                     R * m_radii * Math.Sin(t10 + t14) * 3.0 +
                     R * m_radii * t18 * m_theta * 8.0 -
                     R * m_radii * t20 * m_theta * 8.0 -
                     R * m_radii * t18 * CMath.M_PI * 4.0 + R * m_radii * t20 * CMath.M_PI * 4.0 +
                     R * alpha * m_radii * t18 * 4.0 -
                     R * alpha * m_radii * t20 * 4.0)) /
                   (R * (m_radii * 2.0 + R * t18 * 4.0 + R * Math.Cos(t4) * 3.0 +
                         R * Math.Cos(t15) + m_radii * Math.Cos(alpha) * 2.0 +
                         m_radii * Math.Cos(t2) * 2.0 + m_radii * Math.Cos(t6) * 2.0 -
                         R * t5 * CMath.M_PI * 2.0 + R * alpha * t5 * 2.0 -
                         m_radii * t3 * CMath.M_PI * 2.0 + R * t5 * m_theta * 4.0 +
                         alpha * m_radii * t3 * 2.0 + m_radii * t3 * m_theta * 4.0));
        }
    }

    public double computedEdd(in double R, in double alpha)
    {
        if (R == 0.0)
        {
            return 0.0;
        }
        else
        {
            double t2 = Math.Sin(alpha);
            double t3 = alpha + m_theta;
            double t4 = Math.Sin(t3);
            double t5 = R * R;
            double t6 = m_radii * m_radii;
            double t7 = m_theta * 2.0;
            double t8 = alpha * 2.0;
            double t9 = t7 + t8;
            double t10 = Math.Sin(t9);
            double t11 = Math.Cos(t8);
            double t12 = Math.Cos(t9);
            double t13 = Math.Cos(m_theta);
            double t14 = t8 + m_theta;
            double t15 = Math.Cos(t14);
            return (m_sigma *
                    (-t5 * CMath.M_PI - t6 * CMath.M_PI + alpha * t5 * 2.0 + alpha * t6 * 2.0 +
                     t5 * t10 * 2.0 + t6 * t10 + t5 * m_theta * 2.0 +
                     t6 * m_theta * 2.0 - t6 * Math.Sin(t7) + t6 * Math.Sin(t8) +
                     R * m_radii * Math.Sin(t14) * 3.0 - R * m_radii * Math.Sin(m_theta) * 2.0 +
                     R * m_radii * Math.Sin(t8 + m_theta * 3.0) + t5 * t12 * CMath.M_PI +
                     t6 * t11 * CMath.M_PI - alpha * t5 * t12 * 2.0 - alpha * t6 * t11 * 2.0 -
                     t5 * t12 * m_theta * 2.0 - t6 * t11 * m_theta * 2.0 +
                     R * m_radii * t13 * m_theta * 4.0 -
                     R * m_radii * t15 * m_theta * 4.0 -
                     R * m_radii * t13 * CMath.M_PI * 2.0 + R * m_radii * t15 * CMath.M_PI * 2.0 +
                     R * alpha * m_radii * t13 * 4.0 -
                     R * alpha * m_radii * t15 * 4.0)) /
                   (R * (m_radii * Math.Cos(alpha + t7) + R * Math.Cos(t3) * 2.0 +
                         m_radii * Math.Cos(alpha) - R * t4 * CMath.M_PI + R * alpha * t4 * 2.0 -
                         m_radii * t2 * CMath.M_PI + R * t4 * m_theta * 2.0 +
                         alpha * m_radii * t2 * 2.0 + m_radii * t2 * m_theta * 2.0));
        }
    }

    public double computeRPlanar(in double alpha, in double d0)
    {
        return (d0 - m_radii * Math.Cos(alpha)) / (Math.Cos(m_theta + alpha) + Math.Cos(m_theta));
    }

    public double computeAPlanar(in double R, in double alpha)
    {
        return 2.0 * (0.5 * m_radii * m_radii * Math.Sin(alpha) * Math.Cos(alpha) +
                m_radii * Math.Sin(alpha) * R * Math.Cos(m_theta + alpha) +
                0.5 * R * R * Math.Sin(m_theta + alpha) * Math.Cos(m_theta + alpha)) +
         2.0 *
             (R * Math.Sin(m_theta + alpha) - R * Math.Sin(m_theta) +
              m_radii * Math.Sin(alpha)) *
             R * Math.Cos(m_theta) +
         R * R * Math.Sin(m_theta) * Math.Cos(m_theta) -
         (alpha * m_radii * m_radii + R * R * (CMath.M_PI - 2.0 * m_theta - alpha));
    }

    public double computeHPlanar(in double R, in double alpha)
    {
        return computeH(R, alpha);
    }

    public double computeApproxAPlanar(in double alpha, in double d0)
    {
        if (m_theta == 0.0)
            return 0.0;
        else
        {
            double gamma = alpha + m_theta * 2.0;
            double t2 = m_theta * 3.0;
            double t3 = Math.Cos(t2);
            double t4 = m_radii * m_radii;
            double t5 = d0 * d0;
            double t6 = Math.Cos(m_theta);
            double t7 = CMath.M_PI * CMath.M_PI;
            double t8 = m_theta * 5.0;
            double t9 = Math.Cos(t8);
            double t10 = gamma * gamma;
            double t11 = Math.Sin(m_theta);
            double t12 = Math.Sin(t2);
            double t13 = Math.Sin(t8);
            return 1.0 / (t11 * t11 * t11) *
                   (t3 * t4 * 6.0 - t3 * t5 * 1.2E1 + t5 * t6 * 1.2E1 - t4 * t9 * 6.0 -
                    t4 * t11 * CMath.M_PI * 4.0 - t4 * t12 * CMath.M_PI * 1.6E1 -
                    t5 * t11 * CMath.M_PI * 3.2E1 + t4 * t13 * CMath.M_PI * 4.0 -
                    d0 * m_radii * t3 * 2.4E1 + d0 * m_radii * t6 * 2.4E1 -
                    gamma * t4 * t11 * 3.2E1 + gamma * t4 * t12 * 2.8E1 +
                    gamma * t5 * t11 * 3.2E1 - gamma * t4 * t13 * 4.0 -
                    t3 * t4 * t7 * 3.0 - t3 * t4 * t10 * 3.0 + t4 * t6 * t7 * 2.2E1 +
                    t5 * t6 * t7 * 2.0E1 + t4 * t6 * t10 * 2.2E1 + t4 * t7 * t9 +
                    t5 * t6 * t10 * 2.0E1 + t4 * t9 * t10 + t4 * t11 * m_theta * 7.2E1 -
                    t4 * t12 * m_theta * 2.4E1 + d0 * gamma * m_radii * t11 * 4.0E1 +
                    d0 * gamma * m_radii * t12 * 8.0 + d0 * m_radii * t6 * t7 * 4.0E1 +
                    d0 * m_radii * t6 * t10 * 4.0E1 -
                    d0 * m_radii * t11 * CMath.M_PI * 4.0E1 -
                    d0 * m_radii * t12 * CMath.M_PI * 8.0 + gamma * t3 * t4 * CMath.M_PI * 6.0 -
                    gamma * t4 * t6 * CMath.M_PI * 4.4E1 - gamma * t5 * t6 * CMath.M_PI * 4.0E1 -
                    gamma * t4 * t9 * CMath.M_PI * 2.0 -
                    d0 * gamma * m_radii * t6 * CMath.M_PI * 8.0E1) *
                   (1.0 / 4.8E1);
        }
    }

    public double computeApproxdEddPlanar(in double alpha, in double d0)
    {
        double gamma = alpha + m_theta * 2.0;

        if (m_theta == 0.0)
        {
            double t2 = d0 * d0;
            double t3 = m_radii * m_radii;
            return m_sigma * 1.0 / Math.Pow(d0 + m_radii, 2.0) *
                   (t2 * CMath.M_PI * 4.0 + t3 * CMath.M_PI * 4.0 - gamma * t2 * 4.0 -
                    gamma * t3 * 4.0 + d0 * m_radii * CMath.M_PI * 8.0 -
                    d0 * gamma * m_radii * 8.0) *
                   (1.0 / 2.0);
        }
        else
        {
            double t2 = m_theta * 2.0;
            double t3 = Math.Cos(t2);
            double t4 = d0 * d0;
            double t5 = m_radii * m_radii;
            double t6 = t3 * t3;
            double t7 = CMath.M_PI * CMath.M_PI;
            double t8 = gamma * gamma;
            double t9 = Math.Sin(t2);
            double t10 = m_theta * 4.0;
            double t11 = Math.Sin(t10);
            return (m_sigma * 1.0 / Math.Pow(d0 + m_radii * t3, 2.0) *
                    (t4 * -2.0 + t3 * t4 * 2.0 + t4 * t7 - t5 * t6 * 2.0 + t4 * t8 -
                     t5 * t7 * 2.0 - t5 * t8 * 2.0 - gamma * t4 * CMath.M_PI * 2.0 +
                     gamma * t5 * CMath.M_PI * 4.0 - t4 * t9 * CMath.M_PI * 2.0 - t5 * t11 * CMath.M_PI -
                     d0 * m_radii * t3 * 4.0 + d0 * m_radii * t6 * 4.0 -
                     d0 * m_radii * t7 - d0 * m_radii * t8 + gamma * t4 * t9 * 2.0 +
                     gamma * t5 * t11 - t3 * t4 * t7 + t3 * t5 * t6 * 2.0 -
                     t3 * t4 * t8 + t3 * t5 * t7 + t3 * t5 * t8 + t5 * t6 * t7 +
                     t5 * t6 * t8 + d0 * gamma * m_radii * t9 * 2.0 +
                     d0 * gamma * m_radii * t11 + d0 * m_radii * t6 * t7 +
                     d0 * m_radii * t6 * t8 + d0 * gamma * m_radii * CMath.M_PI * 2.0 -
                     d0 * m_radii * t9 * CMath.M_PI * 2.0 - d0 * m_radii * t11 * CMath.M_PI +
                     gamma * t3 * t4 * CMath.M_PI * 2.0 - gamma * t3 * t5 * CMath.M_PI * 2.0 -
                     gamma * t5 * t6 * CMath.M_PI * 2.0 -
                     d0 * gamma * m_radii * t6 * CMath.M_PI * 2.0) *
                    (-1.0 / 2.0)) /
                   Math.Sin(m_theta);
        }
    }

    public double computedEddAreaDist(in double A_target, in double d0)
    {
        if (d0 < getDStar())
            return computedEddAreaDist(A_target, getDStar());

        if (d0 < 2.0 * Math.Sqrt(A_target / CMath.M_PI + 2.0 * m_radii * m_radii) - m_radii * 2.0)
            return 0.0;

        double alpha = interpolate_alpha(A_target, d0);
        double gamma = alpha + m_theta;

        double dEdd;

        if (gamma < Math.PI / 2 + m_ang_epsilon && gamma > Math.PI / 2 - m_ang_epsilon)
            dEdd = computeApproxdEdd(alpha, d0);
        else
        {
            double R_target = computeR(alpha, d0);
            dEdd = computedEdd(R_target, alpha);
        }

        return Math.Max(0.0, dEdd);
    }

    public double computedEddAreaDistPlanar(in double A_target, in double d0)
    {
        if (d0 < getDStarPlanar())
            return computedEddAreaDistPlanar(A_target, getDStarPlanar());

        if (m_theta > 0.0 && d0 < Math.Sqrt((1.0 - Math.Cos(m_theta)) * (1.0 - Math.Cos(m_theta)) *
            (A_target + CMath.M_PI * m_radii * m_radii) /
            (m_theta - 0.5 * Math.Sin(2.0 * m_theta))) -
            m_radii)
            return 0.0;
        
        double alpha = interpolate_alpha(A_target, d0);

        double gamma = alpha + m_theta * 2.0;

        double dEdd;

        if (gamma < CMath.M_PI + m_ang_epsilon && gamma > CMath.M_PI - m_ang_epsilon)
        {
            dEdd = computeApproxdEddPlanar(alpha, d0);
        }
        else
        {
            double R_target = computeRPlanar(alpha, d0);
            dEdd = computedEddPlanar(R_target, alpha);
        }

        return Math.Max(0.0, dEdd);
    }

    public void construct_alpha_table()
    {
        double alpha_inc = m_max_alpha / (double)m_discretization;
        double d_inc = m_max_d0 / (double)m_discretization;

        m_alpha_table.resize(m_discretization, m_discretization);
        m_A_table.resize(m_discretization, m_discretization);
        m_dEdd_table.resize(m_discretization, m_discretization);
        m_max_As.resize(m_discretization);
        m_min_As.resize(m_discretization);

        double dmin = getDMin();

        for (int i = 0; i < m_discretization; ++i)
        {
            double d0 = d_inc * (double)i + dmin;

            double maxA = -1e+20;
            double minA = 1e+20;
            for (int j = 0; j < m_discretization; ++j)
            {
                double alpha = alpha_inc * (double)j;

                double gamma = alpha + m_theta;

                double A;

                if (gamma < Math.PI / 2 + m_ang_epsilon && gamma > Math.PI / 2 - m_ang_epsilon)
                {
                    A = computeApproxA(alpha, d0);
                }
                else
                {
                    double R = computeR(alpha, d0);
                    A = computeA(R, alpha);
                }

                m_A_table[j, i] = A;
                maxA = Math.Max(maxA, A);
                minA = Math.Min(minA, A);
            }

            m_max_As[i] = maxA;
            m_min_As[i] = minA;
        }

        for (int i = 0; i < m_discretization; ++i)
        {
            Vectors v = m_A_table.col(i);

            double A_inc = (m_max_As[i] - m_min_As[i]) / m_discretization;

            for (int j = 0; j < m_discretization; ++j)
            {
                double A_target = A_inc * (double)j + m_min_As[i];

                int retidx = Math.Min(CMath.bipart_closest(v, A_target),
                                      m_discretization - 1);
                int otheridx = Math.Min(m_discretization - 1, retidx + 1);

                double A0 = m_A_table[retidx, i];
                double A1 = m_A_table[otheridx, i];

                double alpha_target;

                if (A1 == A0)
                    alpha_target = retidx * alpha_inc;
                else
                    alpha_target =
                        ((A_target - A0) / (A1 - A0) * (double)(otheridx - retidx) +
                         retidx) *
                        alpha_inc;

                if (alpha_target > m_max_alpha)
                {
                    alpha_target = m_max_alpha;
                }
                else if (alpha_target < 0)
                {
                    alpha_target = 0.0;
                }

                m_alpha_table[j, i] = alpha_target;
            }
        }

        for (int j = 0; j < m_discretization; ++j)
        {
            for (int i = 0; i < m_discretization; ++i)
            {
                double A_inc = (m_max_As[i] - m_min_As[i]) / m_discretization;
                double A_target = Math.Max(m_min_As[i], A_inc * (double)j + m_min_As[i]);

                double d0 = d_inc * (double)i;
                m_dEdd_table[j, i] = computedEddAreaDist(A_target, d0);
            }
        }
    }

    public void construct_planar_alpha_table()
    {
        double alpha_inc = m_max_alpha / (double)m_discretization;
        double d_inc = m_max_d0 / (double)m_discretization;

        m_alpha_planar_table.resize(m_discretization, m_discretization);
        m_A_planar_table.resize(m_discretization, m_discretization);
        m_dEdd_planar_table.resize(m_discretization, m_discretization);
        m_max_A_planars.resize(m_discretization);
        m_min_A_planars.resize(m_discretization);

        double dmin = getDMinPlanar();

        for (int i = 0; i < m_discretization; ++i)
        {
            double d0 = d_inc * (double)i + dmin;

            double maxA = -1e+20;
            double minA = 1e+20;
            for (int j = 0; j < m_discretization; ++j)
            {
                double alpha = alpha_inc * (double)j;

                double gamma = alpha + m_theta * 2.0;

                double A;

                if (gamma < CMath.M_PI + m_ang_epsilon && gamma > CMath.M_PI - m_ang_epsilon)
                {
                    A = computeApproxAPlanar(alpha, d0);
                }
                else
                {
                    double R = computeRPlanar(alpha, d0);
                    A = computeAPlanar(R, alpha);
                }

                m_A_planar_table[j, i] = A;
                maxA = Math.Max(maxA, A);
                minA = Math.Min(minA, A);
            }

            m_max_A_planars[i] = maxA;
            m_min_A_planars[i] = minA;
        }

        for (int i = 0; i < m_discretization; ++i)
        {
            Vectors v = m_A_planar_table.col(i);

            double A_inc = (m_max_A_planars[i] - m_min_A_planars[i]) / m_discretization;

            for (int j = 0; j < m_discretization; ++j)
            {
                double A_target = A_inc * (double)j + m_min_A_planars[i];

                int retidx = Math.Min(CMath.bipart_closest(v, A_target),
                                      m_discretization - 1);
                int otheridx = Math.Min(m_discretization - 1, retidx + 1);

                double A0 = m_A_planar_table[retidx, i];
                double A1 = m_A_planar_table[otheridx, i];

                double alpha_target;

                if (A1 == A0)
                    alpha_target = retidx * alpha_inc;
                else
                    alpha_target =
                        ((A_target - A0) / (A1 - A0) * (double)(otheridx - retidx) +
                         retidx) *
                        alpha_inc;

                if (alpha_target > m_max_alpha)
                {
                    alpha_target = m_max_alpha;
                }
                else if (alpha_target < 0)
                {
                    alpha_target = 0.0;
                }

                m_alpha_planar_table[j, i] = alpha_target;
            }
        }

        for (int j = 0; j < m_discretization; ++j)
        {
            for (int i = 0; i < m_discretization; ++i)
            {
                double A_inc =
                    (m_max_A_planars[i] - m_min_A_planars[i]) / m_discretization;
                double A_target =
                    Math.Max(m_min_A_planars[i], A_inc * (double)j + m_min_A_planars[i]);

                double d0 = d_inc * (double)i;
                m_dEdd_planar_table[j, i] = computedEddAreaDistPlanar(A_target, d0);
            }
        }
    }

    public double getStiffness(in double d0, in double A_target,
        in double pressure_weight)
    {
        double d_star = getDStar();
        double d_hat = getDHat();
        double dmin = getDMin();
        double dry_criterion = 1e-12;

        if (d0 < d_hat)
        {
            if (d0 < dmin)
            {
                return m_collision_stiffness;
            }
            else
            {
                double dEdd_d_hat =
                    A_target < dry_criterion
                        ? 0.0
                        : (interpolate_dEdd(A_target, d_hat) * pressure_weight);

                double stiffness_d_hat = dEdd_d_hat / (d_hat - d_star);

                if (d0 < d_star)
                {
                    // linear interpolate
                    double stiffness =
                        ((stiffness_d_hat - m_collision_stiffness) * d0 +
                         m_collision_stiffness * d_star - stiffness_d_hat * dmin) /
                        (d_star - dmin);
                    return stiffness;
                }
                else
                {
                    // extrapolate the stiffness at d_hat
                    return stiffness_d_hat;
                }
            }
        }
        else
        {
            double dEdd_d0 = A_target < dry_criterion
                                 ? 0.0
                                 : (interpolate_dEdd(A_target, d0) * pressure_weight);
            double stiffness_d0 = dEdd_d0 / (d0 - d_star);
            return stiffness_d0;
        }
    }

    public double getStiffnessPlanar(in double d0, in double A_target, in double pressure_weight)
    {
        double d_star = getDStarPlanar();
        double d_hat = getDHatPlanar();
        double dmin = getDMinPlanar();
        double dry_criterion = 1e-12;

        if (d0 < d_hat)
        {
            if (d0 < dmin)
            {
                return m_collision_stiffness_planar;
            }
            else
            {
                double dEdd_d_hat =
                    A_target < dry_criterion
                        ? 0.0
                        : (interpolate_dEdd_planar(A_target, d_hat) * pressure_weight);

                double stiffness_d_hat = dEdd_d_hat / (d_hat - d_star);

                if (d0 < d_star)
                {
                    // linear interpolate
                    double stiffness =
                        ((stiffness_d_hat - m_collision_stiffness_planar) * d0 +
                         m_collision_stiffness_planar * d_star - stiffness_d_hat * dmin) /
                        (d_star - dmin);
                    return stiffness;
                }
                else
                {
                    // extrapolate the stiffness at d_hat
                    return stiffness_d_hat;
                }
            }
        }
        else
        {
            double dEdd_d0 =
                A_target < dry_criterion
                    ? 0.0
                    : (interpolate_dEdd_planar(A_target, d0) * pressure_weight);
            double stiffness_d0 = dEdd_d0 / (d0 - d_star);
            return stiffness_d0;
        }
    }

    public double getDMinPlanar() 
    {
        return m_min_d0_planar;
    }

    public double getDHatPlanar()
    {
        return getDMinPlanar() * (2.0 * m_radius_multiplier_planar - 1.0);
    }

    public double getDStarPlanar() 
    {
        return getDMinPlanar() * m_radius_multiplier_planar;
    }

    public double getDMin()
    {
        return m_min_d0;
    }

    public double getDHat()
    {
        return getDMin() * (2.0 * m_radius_multiplier - 1.0);
    }

    public double getDStar()
    {
        return getDMin() * m_radius_multiplier;
    }

    public double getRadiusMultiplierPlanar()
    {
        return m_radius_multiplier_planar;
    }

    public double getCollisionStiffnessPlanar()
    {
        return m_collision_stiffness_planar;
    }

    public double getRadiusMultiplier()
    {
        return m_radius_multiplier;
    }

    public double getCollisionStiffness()
    {
        return m_collision_stiffness;
    }
}
