using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class HairComponent : MonoBehaviour {

    /*
     * change needed to StrandParameters Class at StrandParameters
     *
     * StrandParameters class  -> [Serializable] public StrandParameters : needed for Serialization
     * 
     * StrandParameters default Constructor needed:
     * StrandParameters(double radius = 0.085, double YoungsModulus = 2.4e10, double shearModulus=1.5e10,
     *  double stretchingMultiplier, double density=1.13,
     *  double viscosity=1e8, double baseRotation=0.0, double dt=Time.fixedDeltaTime, Vector3 color,
     *  bool accumViscous = true, bool accumViscousBend = true,
     *  bool variableRadiusHair = false, double straightHairs = 1
     *   )
     * 
     *Vector3 Color -> Color.black : Vector3 class to Unity Color class for only managing Col
     */
    [SerializeField]
    private List<StrandParameters> StrandHairs;   //Parameter list for A strand hair set

    private bool IsAnyHairExist = false;

    public List<StrandParameters> GetStrandHairParameters()
    {
        return StrandHairs;
    }


    private void Awake()
    {
        //Call AddFur 
        //Mesh Hair = # get mesh after every strand hair attached to parents' mesh
        //if (mesh.GetVertices() > 0 )
        //  {
        //      IsAnyHairExist = true;
        //   }
    }
    // Start is called before the first frame update
    void Start()
    {     

    }

    // Update is called once per frame
    void FixedUpdate()
    {
        if (IsAnyHairExist)
        {
            //Call WetHairCore::StepSystem ( Needed to be modified? )
        }
    }

    void Update()
    {
        //call Render Function : ex) LineMeshTest
    }
}
