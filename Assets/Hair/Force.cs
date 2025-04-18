using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using System.Security.Cryptography;
using UnityEngine;
using UnityEngine.Pool;

public class Force
{
    protected int m_internal_index_pos;
    protected int m_internal_index_vel;

    protected int m_internal_index_J;
    protected int m_internal_index_Jv;
    protected int m_internal_index_Jxv;
    protected int m_internal_index_tildeK;

    public Force() { }

    public Force(Force other)
    {
        this.m_internal_index_pos = other.m_internal_index_pos;
        this.m_internal_index_vel = other.m_internal_index_vel;
        this.m_internal_index_J = other.m_internal_index_J;
        this.m_internal_index_Jv = other.m_internal_index_Jv;
        this.m_internal_index_Jxv = other.m_internal_index_Jxv;
        this.m_internal_index_tildeK = other.m_internal_index_tildeK;
    }

    public Force Clone()
    {
        return new Force(this);
    }

    virtual public void addEnergyToTotal(in Vectors x, in Vectors v, in Vectors m, double E) { }

    virtual public void addGradEToTotal(in Vectors x, in Vectors v, in Vectors m, ref Vectors gradE) { }

    virtual public void addHessXToTotal(in Vectors x, in Vectors v, in Vectors m, ref TripletXs hessE) { }

    virtual public void addHessVToTotal(in Vectors x, in Vectors v, in Vectors m, ref TripletXs hessE) { }

    virtual public void addGradEToTotal(in Vectors x, in Vectors v, in Vectors m, ref Vectors gradE, int pidx) { }

    virtual public void addHessXToTotal(in Vectors x, in Vectors v, in Vectors m, ref Vectors hessE, int pidx) { }

    virtual public void addHessVToTotal(in Vectors x, in Vectors v, in Vectors m, ref Vectors hessE, int pidx) { }

    virtual public void computeIntegrationVars(in Vectors x, in Vectors v, in Vectors m, ref Vectors lambda,
        ref Vectors lambda_v, ref TripletXs J, ref TripletXs Jv, ref TripletXs Jxv, ref TripletXs tildeK, ref TripletXs stiffness,
        ref TripletXs damping, ref Vectors Phi, ref Vectors Phiv, in double dt) { }

    virtual public int numConstraintPos() { return 0; }

    virtual public int numConstraintVel() { return 0; }

    virtual public int numJ() { return 0; }

    virtual public int numJv() { return 0; }

    virtual public int numJxv() { return 0; }

    virtual public int numTildeK() { return 0; }

    virtual public bool isParallelized() { return false; }

    virtual public bool isPrecomputationParallelized() { return false; }

    virtual public string name() { return ""; }

    virtual public void getAffectedVars(int pidx, ref HashSet<int> vars) { }

    virtual public bool isContained(int pidx) { return false; }


    virtual public void preCompute(in Vectors x, in Vectors v, in Vectors m, in double dt) { }

    virtual public bool isInterHair(in Vectors lambda, in Vectors lambda_v)
    {
        return false;
    }

    virtual public void storeLambda(in Vectors lambda, in Vectors lamda_v) { }

    virtual public void postStepScene(in double dt) { }

    virtual public int getAffectedHair(in List<int> particle_to_hairs){ return -1;}

    virtual public bool isExternal() { return false; }

    virtual public void setInternalIndex(int index_pos, int index_vel, int index_J, int index_Jv, int index_Jxv, int index_tildeK)
    {
        m_internal_index_pos = index_pos;
        m_internal_index_vel = index_vel;
        m_internal_index_J = index_J;
        m_internal_index_Jv = index_Jv;
        m_internal_index_Jxv = index_Jxv;
        m_internal_index_tildeK = index_tildeK;
    }

}