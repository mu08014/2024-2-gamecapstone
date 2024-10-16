using System.Collections.Generic;
using UnityEditor;
using UnityEngine;
using UnityEngine.Assertions;

/// <summary>
/// Mesh of hairs aka line mesh
/// </summary>
public class FurMesh : MonoBehaviour
{
    private MeshFilter _meshFilter = null;
    public MeshFilter meshFilter
    {
        get {
            if (_meshFilter == null)
                _meshFilter = GetComponent<MeshFilter>();
            return _meshFilter;
        }
    }

    /// <summary>
    /// Represetaion of one hair.
    /// </summary>
    public class Hair
    {
        public Vector3[] _positions;
        public Vector3 _normal;
        public Vector2 _uv;

        public Hair(List<Vector3> positions, Vector3 normal, Vector2 uv)
        {
            // There's no dot hair.
            Assert.IsTrue(positions.Count >= 2);
            _positions = positions.ToArray();
            _normal = normal;
            _uv = uv;
        }
    }

    private List<Hair> _hairs = new();

    public void AddHair(Hair hair)
    {
        _hairs.Add(hair);
    }

    public void Clear()
    {
        _hairs.Clear();
    }

    private List<Vector3> tangets = new();
    private List<Vector2> uvs = new();
    private List<Vector3> normals = new();
    private List<Vector3> positions = new();
    private List<int> indices = new();

    public void UpdateMesh()
    {
        Mesh mesh = meshFilter.mesh;
        mesh.Clear();
        tangets.Clear();
        uvs.Clear();
        normals.Clear();
        positions.Clear();
        indices.Clear();

        foreach (Hair hair in _hairs)
        {
            int head = positions.Count;

            tangets.Add(hair._positions[0] - hair._positions[1]);
            for (int i = 1; i < hair._positions.Length; i++)
                tangets.Add(hair._positions[i - 1] - hair._positions[i]);

            for (int i = 0; i < hair._positions.Length; i++)
                positions.Add(hair._positions[i]);

            for (int i = 0; i < hair._positions.Length; i++)
                normals.Add(hair._normal);

            for (int i = 0; i < hair._positions.Length; i++)
                uvs.Add(hair._uv);

            for (int i = 0; i < hair._positions.Length - 1; i++)
            {
                indices.Add(head + i);
                indices.Add(head + i + 1);
            }
        }

        mesh.SetVertices(positions.ToArray());
        mesh.SetUVs(3, uvs.ToArray());
        mesh.SetUVs(1, tangets.ToArray());
        mesh.SetUVs(2, normals.ToArray());
        mesh.SetIndices(indices.ToArray(), MeshTopology.Lines, 0);
        mesh.RecalculateBounds();
        //meshFilter.mesh.Clear();
        //meshFilter.mesh = mesh;
    }
}
