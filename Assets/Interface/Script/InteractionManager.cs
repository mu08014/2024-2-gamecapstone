using System;
using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using Unity.VisualScripting;
using UnityEngine;

public class InteractionManager : MonoBehaviour
{
    private List<LiquidComponent> LiquidComponents;

    private List<HairComponent> HairComponents;

    private StrandParameters[] allStrandParameters;

    public TwoDScene libwethairScene;

    public SceneStepper scene_stepper;


    private bool IsFluidSimul = false;
    private bool IsHairSimul = false;
    public float executeSpeed = 1.0f;


    [SerializeField]
    private WetHairParameter SimulationParameter;

    [SerializeField]
    private bool isMassSpring = false;

    [SerializeField]
    private int numOfAllHairParticles;



    private void Awake()
    {
        libwethairScene = new TwoDScene(isMassSpring);

        LiquidComponents = new List<LiquidComponent>();
        HairComponents = new List<HairComponent>();
        allStrandParameters = new StrandParameters[0];
        var Liquids = FindObjectsOfType(typeof(LiquidComponent)); // need Assertion Check
        
        foreach (var obj in Liquids)
        {
            LiquidComponents.Add(obj.GetComponent<LiquidComponent>());
        }
        if (Liquids.Length > 0) IsFluidSimul = true;

        var Hairs = FindObjectsOfType(typeof(HairComponent));
        HairComponent tempComp;
        List<StrandParameters> tempParam = new List<StrandParameters>() ;
        foreach (var obj in Hairs)
        {
            tempComp = obj.GetComponent<HairComponent>();
            HairComponents.Add(tempComp);
            foreach (var param in tempComp.GetStrandHairParameters()) {
                tempParam.Add(param);
            }
        }
        allStrandParameters = tempParam.ToArray();
        if (Hairs.Length > 0) IsHairSimul=true;
        


        InitSimulScene();

        scene_stepper = new SceneStepper();

        loadIntegrator(ref libwethairScene, ref scene_stepper); //여기서 scene_stepper 생성 , StrandCompliantManager
        scene_stepper.init(ref libwethairScene);


        libwethairScene.updateHairConnectivity();
        libwethairScene.categorizeForces();
        libwethairScene.initializeScriptedGroup();
    }
    void Start()
    {
        Application.targetFrameRate = 30;
        StartCoroutine(HairParticleUpdate());
    }


    public void Clear()
    {

    }
    // Update is called once per frame
    private void FixedUpdate()
    {
        
    }

    private void InitSimulScene() //Awake()에서 호출, haircomponents에서 m, v, fixed 추출
    {

        for (int i = 0; i < allStrandParameters.Length;i++)
        {
            libwethairScene.insertStrandParameters(ref allStrandParameters[i]);
        }
        if (IsFluidSimul) { }

        if (IsHairSimul)
        {

            int numparticles = 0;
            int numstrands = 0;
            int numedges = 0;  

            for (int i = 0; i < HairComponents.Count; i++)
            {

                numparticles += HairComponents[i].hairParticles.Count;
                numstrands += HairComponents[i].numofStrands;
                numedges += HairComponents[i].hairParticles.Count - HairComponents[i].numofStrands;

            }
            numOfAllHairParticles = numparticles;
            int[] strandEnd = new int[numstrands];
            int strandCount = 0;
            for (int i  = 0; i < HairComponents.Count; i++)
            {
                int start = 0;
                int end = 0;

                for (int j = 0; j < HairComponents[i].hairParticles.Count; j++)
                {
                    if (HairComponents[i].m_hair_fixed[j])
                    {
                        end = j;

                        if (j == 0) continue;
                        if (strandCount == 0)
                        {
                            strandEnd[strandCount] = end - start- 1;
                        }
                        else
                        {
                            strandEnd[strandCount] = strandEnd[strandCount - 1] + end - start;
                        }
                        start = j;
                        strandCount++;
                    }
                }
                strandEnd[strandCount] = strandEnd[strandCount-1] +  HairComponents[i].hairParticles.Count - start;
            }

            libwethairScene.resizeSystem
                (numparticles,
                numedges,
                numstrands);
            Debug.Log("Numparticles = "+numparticles);
            Debug.Log("Numstrands = " + numstrands);
            Debug.Log("Numedges = " + numedges);


            List<int> vert_to_dof = new List<int>();
            VectorXi dofVars = new VectorXi(numparticles * 4 - numstrands + 1);
            VectorXi dofVerts = new VectorXi(numparticles * 4 - numstrands + 1);
            
            VectorXi dofs = new VectorXi(4);
            List<bool> tipVerts = new List<bool>(numparticles);
            tipVerts.AddRange(new bool[numparticles]);


            for (int k = 0; k < 4; k++)
            {
                dofs[k] = k;
            }

            int dof = 0;
            strandCount = 0;
            for (int i = 0; i < numparticles; i++)
            {
                dofVars.SetSegment(dof,4, dofs);
                dofVerts.SetSegment(dof, 4, dofVerts.segment(dof, 4).setConstant(i));
                vert_to_dof.Add(dof);
                dof += 4;
                if (i  == strandEnd[strandCount]) 
                {
                    --dof; 
                    strandCount++;
                    tipVerts[i] = true;
                }
            }
            dofVars.conservativeResize(numparticles * 4 - numstrands);
            dofVerts.conservativeResize(numparticles * 4 - numstrands);
            libwethairScene.setVertToDoFMap(vert_to_dof, dofVars, tipVerts, dofVerts);
            List<Vectors> particle_pos = new List<Vectors>();
            int vtx = 0;
            int edx = 0;
            int paramsIndex = 0;
            int globalStrandID = 0;

            List<List<int>> particle_indices_vec = new List<List<int>>();
            List<List<double>> particle_eta_vec = new List<List<double>>();
            List<List<bool>> particle_state_vec = new List<List<bool>>();

            strandCount = 0;
            List<double> particle_eta = new List<double>();
            List<bool> particle_state = new List<bool>();
            List<int> particle_indices = new List<int>();
            int globalvtx = vtx;

            for (int i = 0; i < HairComponents.Count; i++)
            {
                for (int j = 0; j < HairComponents[i].m_hair_x.Length; j++)
                {
                    particle_indices.Add(vtx);
                    double eta = 0;
                    particle_eta.Add(eta);
                    bool isSource = false;
                    particle_state.Add(isSource);
                    libwethairScene.setPosition(vtx, HairComponents[i].m_hair_x[j]);
                    particle_pos.Add(HairComponents[i].m_hair_x[j]);

                    libwethairScene.setVelocity(vtx, HairComponents[i].m_hair_v[j]);
                    libwethairScene.setFixed(vtx, HairComponents[i].m_hair_fixed[j]);


                    if(vtx != globalvtx)
                    {
                        Tuple<int, int> newedge = new Tuple<int, int>(vtx - 1, vtx);
                        libwethairScene.setEdge(edx, newedge, 0.055);
                        libwethairScene.setEdgeRestLength(edx, (libwethairScene.getPosition(newedge.Item1) -
                                          libwethairScene.getPosition(newedge.Item2))
                                             .norm());
                        edx++;
                        
                    }


                    //Debug.Log("numofedges" + edx);
                    if (vtx == strandEnd[strandCount])
                    {
                        particle_indices_vec.Add(particle_indices);
                        particle_eta_vec.Add(particle_eta);
                        particle_state_vec.Add(particle_state);

                        Force newSF = new StrandForce(ref libwethairScene, particle_indices, paramsIndex, globalStrandID++);
                        StrandEquilibriumParameters newSEP = new StrandEquilibriumParameters(particle_pos, 0, 0, 0, 0, false);
                        libwethairScene.insertForce(ref newSF);
                        libwethairScene.insertStrandEquilibriumParameters(ref newSEP);

                        particle_eta = new List<double>();
                        particle_state = new List<bool>();
                        particle_indices = new List<int>();
                        globalvtx = ++vtx;

                        strandCount++;
                        
                    }
                    else vtx++;

                }

                Debug.Log("A : " + particle_indices_vec.Count);
                Debug.Log("B : " + particle_eta_vec.Count);
                Debug.Log("C : " + particle_state_vec.Count);
            }
            libwethairScene.computeMassesAndRadiiFromStrands();

            List<HairFlow> hairflow = libwethairScene.getFilmFlows();
            for (int i = 0; i < particle_indices_vec.Count; ++i)
            {
                hairflow.Add(null);
            }

            Parallel.For(0, particle_indices_vec.Count, p => {
                libwethairScene.getFilmFlows()[p] = new CylindricalShallowFlow(
                ref libwethairScene, particle_indices_vec[p],
                new Vectors(particle_eta_vec[p], 0,
                                     particle_eta_vec[p].Count),
                particle_state_vec[p]);
            });

        }
        //for (int i = 0; i < (libwethairScene.getX().size() / 3); i++)
        //    Debug.Log("X position of " + i + "th particle is " + libwethairScene.getPosition(i)[0]);
        //Debug.Log("Completed!");
    }
    IEnumerator HairParticleUpdate()
    {
        while (true)
        {
            
            WetHairParameter parameter = libwethairScene.getParameter();
            double hairsubstep = Time.fixedDeltaTime * executeSpeed;

            libwethairScene.applyScript(hairsubstep);

            libwethairScene.updateStrandParamsTimestep(hairsubstep);
            HairComponent hairComp = null;

            int particleCount = 0;
            for (int i = 0; i < HairComponents.Count; i++)
            {
                hairComp = HairComponents[i];
                Vectors basePos = (Vectors)hairComp.gameObject.transform.position - (Vectors)hairComp.anchorPos;
                libwethairScene.addPositionValue(particleCount, hairComp.hairParticles.Count, basePos); // 큐브 위치랑 헤어 위치 동기화
                hairComp.anchorPos += (Vector3)basePos;
                particleCount += hairComp.hairParticles.Count;
            }


            for (int i = 0; i < parameter.hairsteps; i++)
            {
                //scene.updatePolygonalStructure(dt); Polygonal Cohesion 구현 후 작업
                libwethairScene.updateStrandStartStates();

                ((StrandCompliantManager)scene_stepper).stepScene(ref libwethairScene, hairsubstep, true);

                ((StrandCompliantManager)scene_stepper).accept(ref libwethairScene, hairsubstep);

                ((StrandCompliantManager)scene_stepper).PostStepScene(ref libwethairScene, hairsubstep);

            }

            particleCount = 0;
            for (int i = 0; i < HairComponents.Count; ++i)
            {
                foreach (HairParticle p in HairComponents[i].hairParticles)
                {
                    Vectors x_vec = libwethairScene.getPosition(particleCount).Clone();
                    Vectors v_vec = libwethairScene.getVelocity(particleCount).Clone();
                    p.position = new Vector3((float)x_vec[0], (float)x_vec[1], (float)x_vec[2]);
                    p.velocity = new Vector3((float)v_vec[0], (float)v_vec[1], (float)v_vec[2]);

                    particleCount++;
                } 
            }
            yield return new WaitForSeconds(Time.fixedDeltaTime);
        }
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
}



