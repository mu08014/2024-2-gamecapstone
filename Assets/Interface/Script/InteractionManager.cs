using System;
using System.Collections;
using System.Collections.Generic;
using TMPro.EditorUtilities;
using Unity.VisualScripting;
using UnityEngine;

public class InteractionManager : MonoBehaviour
{
    private List<LiquidComponent> LiquidComponents;

    private List<HairComponent> HairComponents;

    private List<StrandParameters> allStrandParameters;

    private TwoDScene libwethairScene;
    [SerializeField]
    private WetHairParameter SimulationParameter;

    [SerializeField]
    private bool isMassSpring = false;

    private bool IsFluidSimul = false;
    private bool IsHairSimul = false;
    // Start is called before the first frame update

    private void Awake()
    {

        LiquidComponents = new List<LiquidComponent>();
        HairComponents = new List<HairComponent>();
        allStrandParameters = new List<StrandParameters>();
        var Liquids = FindObjectsOfType(typeof(LiquidComponent)); // need Assertion Check
        
        foreach (var obj in Liquids)
        {
            Debug.Log(obj);
            LiquidComponents.Add(obj.GetComponent<LiquidComponent>());
        }
        if (Liquids.Length > 0) IsFluidSimul = true;

        var Hairs = FindObjectsOfType(typeof(HairComponent));
        foreach (var obj in Hairs)
        {
            Debug.Log(obj.name);
            HairComponent tempComp = obj.GetComponent<HairComponent>();
            HairComponents.Add(tempComp);
            foreach (var param in tempComp.GetStrandHairParameters()) {
                allStrandParameters.Add(param);
            }
            
        }
        if (Hairs.Length > 0) IsHairSimul=true;

        InitSimulScene();
    }
    void Start()
    {
    }

    // Update is called once per frame
    private void FixedUpdate()
    {
        
    }

    private void InitSimulScene() //Awake()에서 호출, haircomponents에서 m, v, fixed 추출
    {
        libwethairScene = new TwoDScene(isMassSpring);

        if (IsFluidSimul) { }

        Debug.Log(IsHairSimul);
        if (IsHairSimul)
        {

            int numparticles = 0;
            int numstrands = 0;

            for (int i = 0; i < HairComponents.Count; i++)
            {

                numparticles += HairComponents[i].hairParticles.Count;
                numstrands += HairComponents[i].numofStrands;

            }
            List<bool> tipVerts = new List<bool>(numparticles);
            int[] strandEnd = new int[numstrands];
            int strandCount = 0;
            for (int i = 0; i < HairComponents.Count; i++)
            {
                for (int j = 0; j < HairComponents[i].StartofEachParam.Length; j++)
                {
                    if (strandCount == 0)
                    {
                        strandEnd[strandCount] =  HairComponents[i].EndofEachParam[j] + HairComponents[i].StartofEachParam[j];
                    }
                    else
                    {
                        strandEnd[strandCount] = strandEnd[strandCount - 1] + 1 + HairComponents[i].EndofEachParam[j] + HairComponents[i].StartofEachParam[j];
                    }
                    strandCount++;
                }
            }

            libwethairScene.resizeSystem
                (numparticles,
                numparticles - 1,
                numstrands);

            List<int> vert_to_dof = new List<int>();
            VectorXi dofVars = new VectorXi(numparticles * 4 - numstrands + 1);
            VectorXi dofVerts = new VectorXi(numparticles * 4 - numstrands + 1);
            VectorXi dofs = new VectorXi(4);

            
            for (int k = 0; k < 4; k++)
            {
                dofs[k] = k;
            }

            int dof = 0;
            strandCount = 0;
            for (int i = 0; i < numparticles; i++)
            {
                dofVars.SetSegment(dof,4, dofs);
                dofVerts.SetSegment(dof, 4, dofVerts.segment(dof, 4).setConstant(numparticles));
                vert_to_dof.Add(dof);
                dof += 4;
                if (i == strandEnd[strandCount]) 
                {
                    dof-- ; 
                    strandCount++;
                    tipVerts[i] = true;
                }
            }
            dofVars.conservativeResize(numparticles * 4 - numstrands);
            dofVerts.conservativeResize(numparticles * 4 - numstrands);
            libwethairScene.setVertToDoFMap(vert_to_dof, dofVars, tipVerts, dofVerts);

            int partNum = 0;
            for (int i = 0; i < HairComponents.Count; i++)
            {
                for (int j = 0; j < HairComponents[i].m_hair_x.Length; j++)
                {
                    libwethairScene.setPosition(partNum, HairComponents[i].m_hair_x[j]);
                    libwethairScene.setVelocity(partNum, HairComponents[i].m_hair_v[j]);
                    libwethairScene.setFixed(partNum, HairComponents[i].m_hair_fixed[j]);
                    partNum++;
                }
                
            }
                
        }
        for (int i = 0; i < (libwethairScene.getX().size() / 3); i++)
            Debug.Log("X position of " + i + "th particle is " + libwethairScene.getPosition(i)[0]);
        Debug.Log("Completed!");
    }
   
}
