using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

[Serializable]
public class HairParticle 
{

    public Vector3 position;
    public Vector3 velocity;
    public Vector3 buf0;
    public MatrixXs c;
    List<int> bridges;


    double radii;
    double fresh;
    double pressure;
    double edge_alpha;

    public bool IsFixed;
    public int strandID;

    public int edge_idx;
    public bool deceased;

    public HairParticle(Vector3 position, Vector3 velocity, bool IsFixed)
    {
        this.position = position;
        this.velocity = velocity;
        this.IsFixed = IsFixed;
        this.strandID = -1;
    }
    public HairParticle(Vector3 position)
    {
        this.position = position;
        this.velocity = Vector3.zero;
        this.IsFixed = false;
        this.strandID = -1;
    }

    public HairParticle(Vector3 position, int strandID)
    {
        this.position = position;
        this.velocity = Vector3.zero;
        this.IsFixed = false;
        this.strandID = strandID;
    }
    public HairParticle(double x, double y, double z, double vx, double vy, double vz)
    {
        position = new Vector3((float)x, (float)y, (float)z);
        velocity = new Vector3((float)x, (float)y, (float)z);
        IsFixed= false;
    }

}

public class HairParticleTestMode : MonoBehaviour
{
    private Vector3 m_Position;
    private Vector3 m_Velocity;
    private bool m_Fix;

    public Vector3 Position
    {
        get
        {
            return m_Position;
        }
        set
        {
            m_Position = value;
        }
    }

    public Vector3 Velocity
    {
        get
        {
            return m_Velocity;
        }
        set
        {
            m_Velocity = value;
        }
    }

    public bool Fix
    {
        get
        {
            return m_Fix;
        }
        set
        {
            m_Fix = value;
        }
    }

    public HairParticleTestMode(double x, double y, double z, double vx, double vy, double vz)
    {
        m_Position = new Vector3((float)x, (float)y, (float)z);
        m_Velocity = new Vector3((float)x, (float)y, (float)z);
        m_Fix = false;
    }

    void Awake()
    {
        transform.position = m_Position;
    }

    // Update is called once per frame
    void Update()
    {
        transform.position = m_Position;
    }
}
