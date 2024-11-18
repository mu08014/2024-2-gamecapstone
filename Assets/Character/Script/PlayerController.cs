using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PlayerController : MonoBehaviour
{

    private Movement3D mMovement;
    
    // Start is called before the first frame update
    void Start()
    {
        mMovement = GetComponent<Movement3D>();
    }

    // Update is called once per frame
    void Update()
    {


    }
}
