using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class InteractionManager : MonoBehaviour
{
    private List<LiquidComponent> LiquidComponents;

    private List<HairComponent> HairComponents;
    private bool IsAnySimulComponent = false;
    // Start is called before the first frame update

    private void Awake()
    {
        GameObject[] Liquids = GameObject.FindGameObjectsWithTag("Liquid"); // need Assertion Check
        foreach ( GameObject obj in Liquids)
        {
            LiquidComponents.Add(GetComponent<LiquidComponent>());
        }
        GameObject[] Hairs = GameObject.FindGameObjectsWithTag("Hair");
        foreach (GameObject obj in Hairs)
        {
            HairComponents.Add(GetComponent<HairComponent>());
        }
        if (Liquids.Length > 0 && Hairs.Length > 0) IsAnySimulComponent = true;

    }
    void Start()
    {
    }

    // Update is called once per frame
    private void FixedUpdate()
    {
        
    }
}
