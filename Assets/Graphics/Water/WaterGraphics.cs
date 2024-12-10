using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// 물 그래픽 표현에 사용되는 쉐이더를 관리
/// </summary>
public class WaterGraphics : MonoBehaviour
{
    public LBM lbm = null;

    public MeshRenderer targetMesh;
    Material _material;
    Material material
    {
        get => _material != null ? _material : _material = targetMesh.material;
    }

    int bufferID;
    int countID;
    ComputeBuffer buffer = null;

    private void Start()
    {
        bufferID = Shader.PropertyToID("_ParticlePoses");
        countID = Shader.PropertyToID("_ParticleCount");

        if (lbm == null)
        {
            Debug.Log("Please connect the lbm script");
        } else
        {
            StartCoroutine(LateStart());    
        }
    }

    private IEnumerator LateStart()
    {
        yield return null;
        InitBuffer(lbm.spheres.Length);
    }

    public void InitBuffer(int size)
    {
        buffer = new ComputeBuffer(size, sizeof(float) * 3);
    }

    public void UpdateWaterPoses(Vector3[] poses)
    {
        if (buffer == null) return;
        buffer.SetData(poses);
        material.SetBuffer(bufferID, buffer);
        material.SetInt("_ParticleCount", poses.Length);
    }

    private void Update()
    {
        if (lbm != null)
        {
            List<Vector3> poses = new();
            var spheres = lbm.spheres;
            foreach (var i in spheres)
            {
                if (i != null)
                    poses.Add(i.transform.position);
            }
            if (poses.Count != 0)
                UpdateWaterPoses(poses.ToArray());
        }
    }
}
