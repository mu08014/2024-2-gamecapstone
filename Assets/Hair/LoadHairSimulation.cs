using System;
using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using UnityEngine;

public class LoadHairSimulation
{
    public TwoDScene scene;
    public SceneStepper scene_stepper;

    public void loadDERSimulation()
    {
        scene = new TwoDScene(false);

        //List<Script> scripts
        //loadScripts(scripts)

        //loadLiquidParameter()
        loadStrandParameters(ref scene);
        loadStrandParticleEdges(scene);

        //loadSpringForces()
        //loadDragDampingForces()
        //loadLinearSpringForces()
        //loadLinearBenginForces()

        //loadFluidSim()

        //PolygonalCohesion
        //scene->insertForce()
        //scene->setPolygonalCohesion


        scene_stepper = new SceneStepper();
        loadIntegrator(ref scene, ref scene_stepper); //여기서 scene_stepper 생성 , StrandCompliantManager
        scene_stepper.init(ref scene);

        scene.updateHairConnectivity();
        scene.categorizeForces();
        //scene.initializesScriptedGroup(scripts);
    }

    public void loadStrandParameters(ref TwoDScene scene)
    {
        double radius = 0.0025;
        double YoungsModulus = 1e10;
        double shearModulus = 3.4e9;
        double density = 1.3;
        double viscosity = 1e3;
        double stretchingMultiplier = 1.0;
        double baseRotation = 0;
        bool accumulateWithViscous = true;
        bool accumulateViscousOnlyForBendingModes = true;
        bool variableRadiusHair = false;
        double straightHairs = 1;

        StrandParameters newParam = new StrandParameters(
            radius, YoungsModulus, shearModulus, stretchingMultiplier, density,
            viscosity, baseRotation, scene.getDt() / scene.getHairSteps(),
            accumulateWithViscous, accumulateViscousOnlyForBendingModes,
            variableRadiusHair, straightHairs);
        scene.insertStrandParameters(ref newParam);
    }

    public void loadStrandParticleEdges(TwoDScene scene)
    {
        int numstrands = 0;
        int numparticles = 0;
        int numedges = 0;

        List<int> haircount = GameObject.Find("HairInterface").GetComponent<HairInterface>().getHairsCount();

        foreach (int i in haircount)
        {
            numparticles += i;
            numedges += i - 1;
            numstrands++;
        }
        scene.resizeSystem(numparticles, numedges, numstrands);
        List<int> vert_to_dof = new List<int>(numparticles); 
        vert_to_dof.AddRange(new int[numparticles]);
        VectorXi dofVars = new VectorXi(numparticles * 4 - numstrands + 1);
        VectorXi dofVerts = new VectorXi(numparticles * 4 - numstrands + 1);
        VectorXi dofs = new VectorXi(4);
        List<bool> tipVerts = new List<bool>(numparticles);
        tipVerts.AddRange(new bool[numparticles]);
        for (int i = 0; i < 4; i++)
        {
            dofs[i] = i;
        }
        int dof = 0;
        numparticles = 0;
        foreach (int i in haircount)
        {
            for (int j = 0; j < i; j++)
            {
                dofVars.SetSegment(dof, 4, dofs);
                dofVerts.SetSegment(dof, 4, dofVerts.segment(dof, 4).setConstant(numparticles));
                vert_to_dof[numparticles++] = dof;
                dof += 4;
            }
            --dof;
            tipVerts[numparticles - 1] = true;
        }
        dofVars.conservativeResize(numparticles * 4 - numstrands);
        dofVerts.conservativeResize(numparticles * 4 - numstrands);
        scene.setVertToDoFMap(vert_to_dof, dofVars, tipVerts, dofVerts);
        List<VectorXs> particle_pos = new List<VectorXs>();
        VectorXs pos = new VectorXs(3);
        VectorXs vel = new VectorXs(3);
        int vtx = 0;
        int edx = 0;
        int paramsIndex = 0;
        int globalStrandID = 0;

        List<List<int>> particle_indices_vec = new List<List<int>>();
        List<List<double>> particle_eta_vec = new List<List<double>>();
        List<List<bool>> particle_state_vec = new List<List<bool>>();

        foreach(var hair in GameObject.Find("HairInterface").GetComponent<HairInterface>().Hairs)
        {
            List<double> particle_eta = new List<double>();
            List<bool> particle_state = new List<bool>();
            List<int> particle_indices = new List<int>();
            int globalvtx = vtx;

            foreach (var hairparticle in hair)
            {
                particle_indices.Add(vtx);
                double eta = 0;
                particle_eta.Add(eta);
                bool isSource = false;
                particle_state.Add(isSource);
                pos[0] = hairparticle.GetComponent<HairParticle>().Position.x;
                pos[1] = hairparticle.GetComponent<HairParticle>().Position.y;
                pos[2] = hairparticle.GetComponent<HairParticle>().Position.z;

                scene.setPosition(vtx, pos.ToVectors());
                particle_pos.Add(pos);

                vel[0] = hairparticle.GetComponent<HairParticle>().Velocity.x;
                vel[1] = hairparticle.GetComponent<HairParticle>().Velocity.y;
                vel[2] = hairparticle.GetComponent<HairParticle>().Velocity.z;
                scene.setVelocity(vtx, vel.ToVectors());

                bool fixe = false;
                scene.setFixed(vtx, fixe);

                if (vtx !=  globalvtx)
                {
                    Tuple<int, int> newedge = new Tuple<int, int>(vtx - 1, vtx);
                    scene.setEdge(edx, newedge, 0.055);
                    scene.setEdgeRestLength(edx, (scene.getPosition(newedge.Item1) -
                                      scene.getPosition(newedge.Item2))
                                         .norm());
                    edx++;
                }

                ++vtx;
            }
            particle_indices_vec.Add(particle_indices);
            particle_eta_vec.Add(particle_eta);
            particle_state_vec.Add(particle_state);
            Force newSF = new StrandForce(ref scene, particle_indices, paramsIndex, globalStrandID++);
            StrandEquilibriumParameters newSEP = new StrandEquilibriumParameters(particle_pos, 0, 0, 0, 0, false);
            scene.insertForce(ref newSF);
            scene.insertStrandEquilibriumParameters(ref newSEP);
        }

        scene.computeMassesAndRadiiFromStrands();

        List<HairFlow> hairflow = scene.getFilmFlows();
        for (int i = 0; i < particle_indices_vec.Count; ++i)
        {
            hairflow.Add(null);
        }

        Parallel.For(0, particle_indices_vec.Count, p => {
            scene.getFilmFlows()[p] = new CylindricalShallowFlow(
            ref scene, particle_indices_vec[p],
            new VectorXs(particle_eta_vec[p], 0,
                                 particle_eta_vec[p].Count),
            particle_state_vec[p]);
        });
    }

    public void loadIntegrator(ref TwoDScene twodscene, ref SceneStepper scenestepper)
    {
        int max_iters = 50;
        int max_newton = 0;
        double criterion = 1e-8;
        bool compute_interhair = true;
        bool use_preconditioner = true;
        scenestepper = new StrandCompliantManager(
            ref twodscene, max_newton, max_iters, criterion, compute_interhair,
            use_preconditioner);
        //twodscene.notifyGlobalIntegrator();
        //twodscene.notifyFullIntegrator(); 둘 다 PolygonalCohesion 관련 함수
    }

    // WetHairCore의 stepSystem()
    public void hairUpdate(double dt)
    {
        WetHairParameter parameter = scene.getParameter();
        double hairsubstep = dt;

        scene.applyScript(dt);

        scene.updateStrandParamsTimestep(hairsubstep);

        for (int i = 0; i < parameter.hairsteps; ++i)
        {
            //scene.updatePolygonalStructure(dt); Polygonal Cohesion 구현 후 작업
            Debug.Log("vx : " + scene.getVelocity(0)[0] + " vy : " + scene.getVelocity(0)[1] + " vz : " + scene.getVelocity(0)[2]);
            scene.updateStrandStartStates();
            Debug.Log("vx : " + scene.getVelocity(0)[0] + " vy : " + scene.getVelocity(0)[1] + " vz : " + scene.getVelocity(0)[2]);
            ((StrandCompliantManager)scene_stepper).stepScene(ref scene, hairsubstep, true);
            Debug.Log("vx : " + scene.getVelocity(0)[0] + " vy : " + scene.getVelocity(0)[1] + " vz : " + scene.getVelocity(0)[2]);
            ((StrandCompliantManager)scene_stepper).accept(ref scene, hairsubstep);
            Debug.Log("vx : " + scene.getVelocity(0)[0] + " vy : " + scene.getVelocity(0)[1] + " vz : " + scene.getVelocity(0)[2]);
            ((StrandCompliantManager)scene_stepper).PostStepScene(ref scene, hairsubstep);
            Debug.Log("vx : " + scene.getVelocity(0)[0] + " vy : " + scene.getVelocity(0)[1] + " vz : " + scene.getVelocity(0)[2]);
        }
        Debug.Log("vx : " + scene.getVelocity(0)[0] + " vy : " + scene.getVelocity(0)[1] + " vz : " + scene.getVelocity(0)[2]);
    }
}
