using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Assertions;
using UnityEngine.UI;
using UnityEngine.UIElements;

public class AddFurTess : AddFur
{
    [SerializeField] GameObject _furmeshprefab;
    [SerializeField] Material _furmaterial = null;
    [SerializeField] Texture2D _lengthTexture;


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

        static public Vertex operator /(Vertex a, float b)
        {
            return new Vertex(a.position / b, a.normal, a.uv / b);
        }

        static public Vertex operator *(Vertex a, float b)
        {
            return new Vertex(a.position * b, a.normal, a.uv * b);
        }
    }

    void GenHair(Vertex av)
    {
        var noiseLength = 0.0f;
        if (_lengthTexture != null)
            noiseLength = (_lengthTexture.GetPixelBilinear(av.uv.x, av.uv.y).r - 0.5f) * NoisePower;

        List<Vector3> poses = new();
        for (int i = 0; i <= Rod; i++)
            poses.Add(av.position + (Length + noiseLength) * i / Rod * av.normal);

        furmesh.AddHair(new FurMesh.Hair(
                poses, av.normal, av.uv));
    }

    void MakeFur(int level, Vertex v1, Vertex v2, Vertex v3)
    {
        for (int i = 0; i < level; i++)
        {
            Vertex left = (v3 * i + v1 * (level - i)) / level;
            Vertex right = (v3 * i + v2 * (level - i)) / level;
            for (int j = 0; j <= level - i; j++)
            {
                Vertex av = (right * j + left * (level - i - j)) / (level - i);
                GenHair(av);
            }
        }
        GenHair(v3);
    }


    public override void initMesh(GameObject _furmeshPrefab)

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
        furmesh = gameObject.GetComponentInChildren<FurMesh>();
        if (furmesh == null)
        {
            furmesh = Instantiate(_furmeshPrefab, transform).GetComponent<FurMesh>();
            Debug.Log("A new furmesh object created!");
        }
        furmesh.Clear();
        furmesh.Parent = this;

        var vertices = baseMesh.vertices;
        var normals = baseMesh.normals;
        var uvs = baseMesh.uv;
        var tris = baseMesh.GetTriangles(0);

        for (int i = 0; i < tris.Length; i += 3)
        {
            var v1 = new Vertex(vertices[tris[i]], normals[tris[i]], uvs[tris[i]]);
            var v2 = new Vertex(vertices[tris[i+1]], normals[tris[i+1]], uvs[tris[i+1]]);
            var v3 = new Vertex(vertices[tris[i+2]], normals[tris[i+2]], uvs[tris[i+2]]);
            MakeFur(Iter, v1, v2, v3);
        }

        //furmesh.UpdateMesh();
        if (furmesh is FurMeshTess mesh)
            mesh.UpdateMeshTess(Iter, Rod);
        else
            furmesh.UpdateMesh();
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
            initMesh(_furmeshprefab);
            yield return wait;
        }
    }

}
