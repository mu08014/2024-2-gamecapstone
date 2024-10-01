using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using System.Security.Cryptography;
using UnityEngine;
using UnityEngine.Pool;

abstract public class Force
{
    protected int m_internal_index_pos;
    protected int m_internal_index_vel;

    protected int m_internal_index_J;
    protected int m_internal_index_Jv;
    protected int m_internal_index_Jxv;
    protected int m_intermal_index_tildeK;

    virtual public void addEnergyToTotal(in VectorXs x, in VectorXs v, in VectorXs m, double E) { }

    virtual public void addGradEToTotal(in VectorXs x, in VectorXs v, in VectorXs m, ref VectorXs gradE) { }

    virtual public void addHessXToTotal(in VectorXs x, in VectorXs v, in VectorXs m, ref TripletXs hessE) { }

    virtual public void addHessVToTotal(in VectorXs x, in VectorXs v, in VectorXs m, ref TripletXs hessE) { }

    virtual public void addGradEToTotal(in VectorXs x, in VectorXs v, in VectorXs m, ref VectorXs gradE, int pidx) { }

    virtual public void addHessXToTotal(in VectorXs x, in VectorXs v, in VectorXs m, ref VectorXs hessE, int pidx) { }

    virtual public void addHessVToTotal(in VectorXs x, in VectorXs v, in VectorXs m, ref VectorXs hessE, int pidx) { }

    virtual public void computeIntegrationVars(in VectorXs x, in VectorXs v, in VectorXs m, ref VectorXs lambda,
        ref VectorXs lambda_v, ref TripletXs J, ref TripletXs Jv, ref TripletXs Jxv, ref TripletXs tildeK, ref TripletXs stiffness,
        ref TripletXs damping, ref VectorXs Phi, ref VectorXs Phiv, in double dt) { }

    abstract public int numConstraintPos();

    abstract public int numConstraintVel();

    abstract public int numJ();

    abstract public int numJv();

    abstract public int numJxv();

    abstract public int numTildeK();

    abstract public bool isParallelized();

    abstract public bool isPrecomputationParallelized();

    abstract public string name();

    abstract public void getAffectedVars(int pidx, ref HashSet<int> vars);

    abstract public bool isContained(int pidx);


    virtual public void preCompute(in VectorXs x, in VectorXs v, in VectorXs m, in double dt) { }

    virtual public bool isInterHair(in VectorXs lambda, in VectorXs lambda_v)
    {
        return false;
    }

    virtual public void storeLambda(in VectorXs lambda, in VectorXs lamda_v) { }

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
        m_intermal_index_tildeK = index_tildeK;
    }

}