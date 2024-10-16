using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;

public class AddFur : MonoBehaviour
{
    [Tooltip("Iteration of fur generation")]
    [Range(1, 3)]
    public int Iter = 1;

    [Tooltip("Length of generated hairs")]
    public float Length = 1.0f;

    [SerializeField] GameObject _furmeshprefab;
    private FurMesh _furmesh;
    private MeshFilter _mesh;

    void Start()
    {
        _mesh = GetComponent<MeshFilter>();
        _furmesh = Instantiate(_furmeshprefab, transform).GetComponent<FurMesh>();
        initMesh();
    }

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

    void initMesh()
    {
        _furmesh.Clear();

        var vertices = _mesh.mesh.vertices;
        var normals = _mesh.mesh.normals;
        var uvs = _mesh.mesh.uv;
        var tris = _mesh.mesh.GetTriangles(0);

        for (int i = 0; i < tris.Length; i += 3)
        {
            var v1 = new Vertex(vertices[tris[i]], normals[tris[i]], uvs[tris[i]]);
            var v2 = new Vertex(vertices[tris[i+1]], normals[tris[i+1]], uvs[tris[i+1]]);
            var v3 = new Vertex(vertices[tris[i+2]], normals[tris[i+2]], uvs[tris[i+2]]);
            MakeFur(1, v1, v2, v3);
        }

        _furmesh.UpdateMesh();
    }

    void Update()
    {
        initMesh();
    }
}
