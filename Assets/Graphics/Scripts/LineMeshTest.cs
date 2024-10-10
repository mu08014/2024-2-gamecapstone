using System.Collections;
using System.Collections.Generic;
using Unity.Collections;
using UnityEngine;
using UnityEngine.Rendering;

public class LineMeshTest : MonoBehaviour
{
    [System.Runtime.InteropServices.StructLayout(System.Runtime.InteropServices.LayoutKind.Sequential)]
    struct VertexData
    {
        public Vector3 position;
        public Vector3 tangent;
    }

    public int lineCount = 10;
    public float length = 1.0f;

    private MeshFilter _meshFilter;

    //https://docs.unity3d.com/ScriptReference/Mesh.SetVertexBufferData.html
    public void Start()
    {
        var mesh = new Mesh();

        var layout = new[]
        {
            new VertexAttributeDescriptor(VertexAttribute.Position, VertexAttributeFormat.Float32, 3),
            new VertexAttributeDescriptor(VertexAttribute.Position, VertexAttributeFormat.Float32, 3)
        };
        var vertexCount = lineCount * 2;
        mesh.SetVertexBufferParams(vertexCount, layout);

        //var verts = new NativeArray<VertexData>(vertexCount, Allocator.Temp);
        //var vindc = 0;

        var vertices = new List<Vector3>();
        var indices = new List<int>();
        var tangents = new List<Vector3>();
        int indc = 0;

        for (int i = 0; i < lineCount; ++i)
        {
            tangents.Add(new Vector3(0, 0.0f, 0.5f));
            tangents.Add(new Vector3(0, 0.0f, 0.5f));
            tangents.Add(new Vector3(0, 0.5f, 0.5f));
            tangents.Add(new Vector3(0, 0.5f, 0.0f));

            float z_offset = -Mathf.Sin(1.0f / lineCount * i * Mathf.PI) / 10;

            vertices.Add(new Vector3(1.0f / lineCount * i, 1.0f, 1.0f + z_offset));
            vertices.Add(new Vector3(1.0f / lineCount * i, 1.0f, 0.5f + z_offset));
            vertices.Add(new Vector3(1.0f / lineCount * i, 0.5f, 0.0f + z_offset));
            vertices.Add(new Vector3(1.0f / lineCount * i, 0.0f, 0.0f + z_offset));
            
            indices.Add(indc);
            indices.Add(indc + 1);

            indices.Add(indc + 1);
            indices.Add(indc + 2);
            
            indices.Add(indc + 2);
            indices.Add(indc + 3);

            indc += 4;
        }

        mesh.SetVertices(vertices.ToArray());
        mesh.SetIndices(indices.ToArray(), MeshTopology.Lines, 0);
        mesh.SetUVs(1, tangents.ToArray());

        //mesh.SetVertexBufferData(verts, 0, 0, vertexCount);

        // recaculate mesh's bounds to prevent object culling by unity
        mesh.RecalculateBounds();
        //mesh.RecalculateNormals();
        //mesh.RecalculateTangents();

        _meshFilter = GetComponent<MeshFilter>();
        _meshFilter.mesh = mesh;
    }
}
