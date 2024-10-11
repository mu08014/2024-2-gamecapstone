using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;


public abstract class DependencyBase
{

    private bool m_dirty;

    protected List<DependencyBase> m_dependents;
   
    protected void setDirtyWithoutPropagating()
    {
        m_dirty = true;
    }

    protected abstract void compute();

    public DependencyBase()
    {
        m_dirty = true;
    }

    public bool isDirty()
    {
        return m_dirty;
    }

    public void setDirty()
    {
        if (!m_dirty)
        {
            m_dirty = true;
            setDependentsDirty();
        }
    }

    public void setDependentsDirty()
    {
        foreach (var dep in m_dependents)
        {
            dep.setDirty();
        }
    }

    public void setClean()
    {
        m_dirty = false;
    }

    public void addDependent(DependencyBase depentent)
    {
        m_dependents.Add(depentent);
    }
}

public abstract class DependencyNode<ValueT> : DependencyBase
{
    protected ValueT m_value;
    public DependencyNode(ValueT value)
    {
        m_value = value;
    }

    public virtual ValueT get()
    {
        if (isDirty())
        {
            compute();
            setClean();
        }
        return m_value;
    }

    public virtual void set(ValueT value)
    {
        setDependentsDirty();
        m_value = value;
    }

    public void free()
    {
        setDirtyWithoutPropagating();
    }
}
