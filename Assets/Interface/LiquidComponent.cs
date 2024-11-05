using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class LiquidComponent : MonoBehaviour
{
    /*
     * public class WetHairParameter -> [Serializable] public class WetHairParameter
     * Default Constructor in TwoDScene.cs needed to be changed 
     *dt -> Time.fixedDeltaTime : change deltatime which is specified for each .xml sample to Unity deltaTime
     *
     *math.cs line71 need to be fixed  : elements[i] -> elements.Add(i)
     *https://stackoverflow.com/questions/4236629/argumentoutofrangeexception-on-initialized-list
     */
    [SerializeField]
    private WetHairParameter LiquidParameter; 
    
    
    [SerializeField]
    private Vector3Int LiquidVolume;

    private int NumOfLiquidParticle
    {
        get
        {
            return LiquidVolume.x * LiquidVolume.y * LiquidVolume.z;
        }
    }

    private void Awake()
    {
        LiquidVolume = new Vector3Int(10, 10, 10);

}
// Start is called before the first frame update
void Start()
    {

    }

    // Update is called once per frame
    void Update()
    {

    }
}
