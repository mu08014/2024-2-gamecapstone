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
    private List<Vector3> emptyPos = new();
    private ComputeBuffer posBuffer = null;
    private ComputeBuffer tangentBuffer = null;
    private Material mat = null;

    public void UpdateMeshTess(int level, int rod)
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

        UpdateMesh();
    }

    public override void UpdateMesh()
    {
        // init compute buffer
        posBuffer ??= new ComputeBuffer(_hairs.Count * 15, sizeof(float) * 3);
        tangentBuffer ??= new ComputeBuffer(_hairs.Count * 15, sizeof(float) * 3);
        mat = mat != null ? mat : GetComponent<MeshRenderer>().material;

        hairMesh = meshFilter.sharedMesh;
        if (hairMesh == null)
        {
            hairMesh = new Mesh();
            meshFilter.mesh = hairMesh;
        }
        
        tangents.Clear();
        uvs.Clear();
        normals.Clear();
        positions.Clear();
        indices.Clear();
        idx.Clear();

        int iin = 0;
        for (int k = 0; k < _hairs.Count; k++)
        {
            Hair hair = _hairs[k];

            for (int j = 0; j < hair._positions.Length; j++)
            {
                emptyPos.Add(Vector3.zero); // just add dummy position
                normals.Add(hair._normal);
                uvs.Add(hair._uv);
                indices.Add(iin++);
                idx.Add(new Vector2(k, j));
            }
        }

        particles = new List<HairParticle>();
        foreach (var hair in _hairs.Select((value, index) => (value, index)))
        {
            foreach (var dot in hair.value._positions.Select((value2, index2) => (value2,index2)))
            {
                HairParticle particle = new HairParticle(dot.value2);
                particles.Add(particle);
                if (dot.index2 == 0)
                {
                    particle.IsFixed = true;
                }
            }
        }

        hairMesh.Clear();
        hairMesh.SetVertices(emptyPos.ToArray());
        //hairMesh.SetUVs(1, tangents.ToArray());
        hairMesh.SetUVs(2, normals.ToArray());
        hairMesh.SetUVs(3, uvs.ToArray());
        hairMesh.SetUVs(4, idx.ToArray());
        hairMesh.SetIndices(indices.ToArray(), MeshTopology.Points, 0);
        hairMesh.RecalculateBounds();

        posBuffer.SetData(positions.ToArray());
        //tangentBuffer.SetData(tangents.ToArray());

        GetComponentInParent<HairComponent>().SetHairInfo(this.Parent, particles);
        
        float ratio = (float)Screen.width / Screen.height;
        meshRenderer.material.SetFloat("_AspectRatio", ratio);
        meshRenderer.material.SetBuffer("_VertexPosition", posBuffer);
    }

    public override void UpdateHairPos(Vector3[] pos)
    {
        posBuffer ??= new ComputeBuffer(pos.Length, sizeof(float) * 3);
        tangentBuffer ??= new ComputeBuffer(pos.Length, sizeof(float) * 3);
        posBuffer.SetData(pos);
        var rend = GetComponentInParent<MeshRenderer>();
        rend.material.SetBuffer("_VertexPosition", posBuffer);
        //rend.material.SetBuffer("_Tengents", tangentBuffer);

        //idx.Clear();
        //for (int k = 0; k < pos.Length / 3; k++)
        //{
        //    indices.Add(k);
        //    idx.Add(new Vector2(k, 0));
        //}
        //hairMesh.SetIndices(indices.ToArray(), MeshTopology.Points, 0);
        //hairMesh.SetUVs(4, idx.ToArray());
    }

    private void OnDestroy()
    {
        posBuffer?.Release();
    }
}
