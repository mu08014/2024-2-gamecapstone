using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SceneStepper
{
    protected List<double> m_timing_statistics;

    protected VectorXs m_old_v;
    protected VectorXs m_a;
    protected VectorXs m_next_x;

    protected int DIM;

    public SceneStepper(int dim = 3)
    {
        m_timing_statistics = new List<double>();
        m_old_v = new VectorXs();
        m_a = new VectorXs();
        m_next_x = new VectorXs();
        DIM = dim;
    }

    public VectorXs getAcceleration()
    {
        return m_a;
    }

    public void accept(ref TwoDScene scene, double dt)
    {
        VectorXs x = scene.getX();
        VectorXs v = scene.getV();

        v = (m_next_x - x) / dt;
        x = m_next_x;
        m_a = (v - m_old_v) / dt;
        m_old_v = v.Clone();
    }

    public void setNextX(in VectorXs nextx)
    {
        m_next_x = nextx.Clone();
    }

    public VectorXs getNextX()
    {
        return m_next_x;
    }

    public void PostStepScene(ref TwoDScene scene, double dt)
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
