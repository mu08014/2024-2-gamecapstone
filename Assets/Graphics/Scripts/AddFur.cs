using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
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
    public float NoisePower = 0.1f;

    [Range(1, 3)]
    public int Rod = 1;

    [SerializeField] GameObject _furmeshprefab;
    [SerializeField] Material _furmaterial = null;
    [SerializeField] Texture2D _lengthTexture;
    private FurMesh _furmesh;
    private MeshFilter _mesh;

    void Start()
    {
        _mesh = GetComponent<MeshFilter>();
        _furmesh = Instantiate(_furmeshprefab, transform).GetComponent<FurMesh>();
        //if (_furmesh != null)
        //    _furmesh.GetComponent<MeshRenderer>().material = _furmaterial;
        initMesh();

        //StartCoroutine(HairUpdateCo());
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
            return new Vertex(a.position / b, a.normal, a.uv / b);
        }
    }

    void MakeFur(int level, Vertex v1, Vertex v2, Vertex v3)
    {
        var av = (v1 + v2 + v3) / 3;

        var noiseLength = 0.0f;
        if (_lengthTexture != null)
            noiseLength = (_lengthTexture.GetPixelBilinear(av.uv.x, av.uv.y).r - 0.5f) * NoisePower;

        List<Vector3> poses = new();
        for (int i = 0; i <= Rod; i++)
            poses.Add(av.position + (Length + noiseLength) * i / Rod * av.normal);

        _furmesh.AddHair(new FurMesh.Hair(
                poses, av.normal, av.uv));

        if (level < Iter)
        {
            MakeFur(level + 1, v1, (v1 + v2) / 2, (v1 + v3) / 2);
            MakeFur(level + 1, (v1 + v2) / 2, v2, (v2 + v3) / 2);
            MakeFur(level + 1, (v1 + v3) / 2, (v2 + v3) / 2, v3);
        }
    }

    public void initMesh()
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
        //initMesh();
    }

    IEnumerator HairUpdateCo()
    {
        var wait = new WaitForSeconds(1.0f);
        while (true)
        {
            initMesh();
            yield return wait;
        }
    }
}
