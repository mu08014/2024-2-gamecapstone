using UnityEngine;

/// <summary>
/// 물 그래픽 표현에 사용되는 쉐이더를 관리
/// </summary>
public class WaterGraphics : MonoBehaviour
{
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
    }

    public void InitBuffer(int size)
    {
        buffer = new ComputeBuffer(size, sizeof(float) * 3);
    }

    public void UpdateWaterPoses(Vector3[] poses)
    {
        buffer.SetData(poses);
        material.SetBuffer(bufferID, buffer);
        material.SetInt("_ParticleCount", poses.Length);
    }
}
