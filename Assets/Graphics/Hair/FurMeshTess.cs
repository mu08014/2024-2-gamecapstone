using System.Collections.Generic;
using System.Linq;
using UnityEditor;
using UnityEngine;

/// <summary>
/// Mesh of hairs aka line mesh with spline method
/// </summary>
public class FurMeshTess : FurMesh
{
    private List<Vector2> idx = new();
    private ComputeBuffer posBuffer = null;
    private ComputeBuffer tangentBuffer = null;
    private Material mat = null;

    public void UpdateMeshTess(int level, int rod, int triCount)
    {

        particles = new List<HairParticle>();
        foreach (var hair in _hairs.Select((value, index) => (value, index)))
        {
            foreach (var dot in hair.value._positions.Select((value2, index2) => (value2, index2)))
            {
                HairParticle particle = new HairParticle(dot.value2);
                particles.Add(particle);
                if (dot.index2 == 0)
                {
                    particle.IsFixed = true;
                }
            }
        }

        mat = mat != null ? mat : GetComponent<MeshRenderer>().sharedMaterial;

        positions.Clear();
        normals.Clear();
        uvs.Clear();
        idx.Clear();
        for (int k = 0; k < _hairs.Count; k++)
        {
            Hair hair = _hairs[k];

            for (int j = 0; j < hair._positions.Length; j++)
            {
                positions.Add(hair._positions[j]);
                normals.Add(hair._normal);
                uvs.Add(hair._uv);
                idx.Add(new Vector2(k, j));
            }
        }

        indices.Clear();
        for (int i = 0; i < triCount; i++)
        {
            int offset = i * (rod + 1) * (level + 2) * (level + 1) / 2;
            for (int r = 0; r <= rod; r++)
            {
                int start = 0;
                for (int j = level; j > 0; j--)
                {
                    for (int k = 0; k < j; k++)
                    {
                        indices.Add((start + k) * (rod + 1) + r + offset);
                        indices.Add((start + k + 1) * (rod + 1) + r + offset);
                        indices.Add((start + k + j + 1) * (rod + 1) + r + offset);
                    }
                    for (int k = 0; k < j - 1; k++)
                    {
                        indices.Add((start + k + 1) * (rod + 1) + r + offset);

                        indices.Add((start + k + j + 2) * (rod + 1) + r + offset);
                        indices.Add((start + k + j + 1) * (rod + 1) + r + offset);
                    }
                    start += j + 1;
                }
            }
        }

        hairMesh = meshFilter.sharedMesh;
        if (hairMesh == null)
        {
            hairMesh = new Mesh();
            meshFilter.mesh = hairMesh;
        }

        hairMesh.Clear();
        hairMesh.SetVertices(positions.ToArray());
        hairMesh.SetUVs(2, normals.ToArray());
        hairMesh.SetUVs(3, uvs.ToArray());
        hairMesh.SetUVs(4, idx.ToArray());
        hairMesh.SetIndices(indices.ToArray(), MeshTopology.Triangles, 0);
        hairMesh.RecalculateBounds();

        //posBuffer.SetData(positions.ToArray());
        //tangentBuffer.SetData(tangents.ToArray());

        GetComponentInParent<HairComponent>().SetHairInfo(this.Parent, particles);

        float ratio = (float)Screen.width / Screen.height;
        meshRenderer.sharedMaterial.SetFloat("_AspectRatio", ratio);

        //UpdateMesh();
    }

    public override void UpdateHairPos(Vector3[] pos)
    {
        posBuffer ??= new ComputeBuffer(pos.Length, sizeof(float) * 3);
        tangentBuffer ??= new ComputeBuffer(pos.Length, sizeof(float) * 3);
        posBuffer.SetData(pos);
        var rend = GetComponentInParent<MeshRenderer>();
        hairMesh.SetVertices(pos);
        rend.sharedMaterial.SetBuffer("_VertexPosition", posBuffer);

        float ratio = (float)Screen.width / Screen.height;
        meshRenderer.sharedMaterial.SetFloat("_AspectRatio", ratio);
    }

    private void OnDestroy()
    {
        posBuffer?.Release();
    }
}
