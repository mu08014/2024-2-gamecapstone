using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SceneStepper
{
    protected List<double> m_timing_statistics;

    protected Vectors m_old_v;
    protected Vectors m_a;
    protected Vectors m_next_x;

    protected int DIM;

    public SceneStepper(int dim = 3)
    {
        m_timing_statistics = new List<double>();
        m_old_v = new Vectors();
        m_a = new Vectors();
        m_next_x = new Vectors();
        DIM = dim;
    }

    public Vectors getAcceleration()
    {
        return m_a;
    }

    public virtual void accept(ref TwoDScene scene, double dt)
    {
        Vectors x = scene.getX();
        Vectors v = scene.getV();

        for (int i = 0; i < v.Size; i++)
        {
            v[i] = (m_next_x[i] - x[i]) / dt;
        }
        for (int i = 0; i < x.Size; i++)
        {
            x[i] = m_next_x[i];
        }
        //v = (m_next_x - x) / dt;
        //x = m_next_x;
        m_a = (v - m_old_v) / dt;
        m_old_v = v.Clone();
    }

    public void setNextX(in Vectors nextx)
    {
        m_next_x = nextx.Clone();
    }

    public Vectors getNextX()
    {
        return m_next_x;
    }

    public virtual void PostStepScene(ref TwoDScene scene, double dt)
    {
        scene.postCompute(dt);
    }

    public void init(ref TwoDScene scene)
    {
        m_old_v = scene.getV().Clone();
        m_a.resize(m_old_v.size());
        m_a.setZero();
    }

    public List<double> getTimingStatics()
    {
        return m_timing_statistics;
    }
}
