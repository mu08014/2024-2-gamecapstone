using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FloorBugFix : MonoBehaviour
{
    // Update is called once per frame
    private Vector3 pos;

    private void Start()
    {
        pos = new Vector3();
    }

    void Update()
    {
        pos = gameObject.transform.position;
        //Debug.Log(pos);
        if (pos.y < -0.3f)
        {
            gameObject.transform.position = new Vector3(pos.x, -0.3f, pos.z);
        }
    }
}
