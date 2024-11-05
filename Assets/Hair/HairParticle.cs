using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class HairParticle : MonoBehaviour
{
    private Vector3 m_Position;
    private Vector3 m_Velocity;

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

    public HairParticle(double x, double y, double z, double vx, double vy, double vz)
    {
        m_Position = new Vector3((float)x, (float)y, (float)z);
        m_Velocity = new Vector3((float)x, (float)y, (float)z);
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
