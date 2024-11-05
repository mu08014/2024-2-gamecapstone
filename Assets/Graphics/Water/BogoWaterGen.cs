using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using UnityEngine;

/// <summary>
/// 물 그래픽 테스트용
/// </summary>
public class BogoWaterGen : MonoBehaviour
{
    public List<Transform> trs = new();
    public WaterGraphics water;

    private void Start()
    {
        int cc = transform.childCount;
        for (int i = 0; i < cc; i++)
        {
            transform.GetChild(i).GetComponent<MeshRenderer>().enabled = false;
            trs.Add(transform.GetChild(i));
        }
        water.InitBuffer(trs.Count);
    }

    void Update()
    {
        List<Vector3> vec = new();
        foreach (var tr in trs)
        {
            vec.Add(tr.position);
        }
        water.UpdateWaterPoses(vec.ToArray());
    }
}
