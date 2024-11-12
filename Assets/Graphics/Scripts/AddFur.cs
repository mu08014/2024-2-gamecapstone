using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Assertions;
using UnityEngine.UI;
using UnityEngine.UIElements;

public class AddFur : MonoBehaviour
{
    [Tooltip("Iteration of fur generation")]
    [Range(1, 3)]
    public int Iter = 1;

    [Tooltip("Length of generated hairs")]
    public float Length = 1.0f;

    private FurMesh _furmesh;
    public FurMesh furmesh
    {
        get
        {
            if (_furmesh == null)
                _furmesh = GetComponent<FurMesh>();
            return _furmesh;
        }
        set
        {
            _furmesh = value;
        }
    }
    private Mesh _mesh;
    private MeshFilter _meshFilter;

    private class Vertex {
        public Vector3 position;
        public Vector3 normal;
        public Vector2 uv;

        public Vertex(Vector3 pos, Vector3 norm, Vector2 u)
        {
            position = pos; normal = norm.normalized; uv = u;
        }

        static public Vertex operator+(Vertex a, Vertex b)
        {
            return new Vertex(a.position + b.position, a.normal + b.normal, a.uv + b.uv);
        }

        static public Vertex operator/(Vertex a, float b)
        {
            return new Vertex(a.position / b, a.normal / b, a.uv / b);
        }
    }

    void MakeFur(int level, Vertex v1, Vertex v2, Vertex v3)
    {
        var av = (v1 + v2 + v3) / 3;
        _furmesh.AddHair(new FurMesh.Hair(
            new List<Vector3> { av.position, av.position + av.normal * Length }, av.normal, av.uv));

        if (level < Iter)
        {
            MakeFur(level + 1, v1, (v1 + v2) / 2, (v1 + v3) / 2);
            MakeFur(level + 1, (v1 + v2) / 2, v2, (v2 + v3) / 2);
            MakeFur(level + 1, (v1 + v3) / 2, (v2 + v3) / 2, v3);
        }
    }

    public void initMesh(GameObject _furmeshPrefab)
    {

        _meshFilter = GetComponent<MeshFilter>();
        if (_meshFilter == null)
        {
            _mesh = GetComponent<SkinnedMeshRenderer>().sharedMesh;
        }
        else
        {
            _mesh = _meshFilter.sharedMesh;
        }

        Debug.Assert( _mesh != null );
        Mesh baseMesh = _mesh;//컴포넌트  존재여부 확인 필요
        _furmesh = gameObject.GetComponentInChildren<FurMesh>();
        if (_furmesh == null)
        {
            _furmesh = Instantiate(_furmeshPrefab, transform).GetComponent<FurMesh>();
            Debug.Log("A new furmesh object created!");
        }
        _furmesh.Clear();
        _furmesh.Parent = this;

        var vertices = baseMesh.vertices;
        var normals = baseMesh.normals;
        var uvs = baseMesh.uv;
        var tris = baseMesh.GetTriangles(0);

        for (int i = 0; i < tris.Length; i += 3)
        {
            var v1 = new Vertex(vertices[tris[i]], normals[tris[i]], uvs[tris[i]]);
            var v2 = new Vertex(vertices[tris[i+1]], normals[tris[i+1]], uvs[tris[i+1]]);
            var v3 = new Vertex(vertices[tris[i+2]], normals[tris[i+2]], uvs[tris[i+2]]);
            MakeFur(1, v1, v2, v3);
        }

        _furmesh.UpdateMesh();
    }

}
