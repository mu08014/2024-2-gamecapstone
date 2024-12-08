using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;
using Unity.Collections;
using Unity.VisualScripting;
using UnityEditor;
using UnityEngine;
using UnityEngine.Assertions;

/// <summary>
/// Mesh of hairs aka line mesh
/// </summary>
public class FurMesh : MonoBehaviour
{
    private MeshFilter _meshFilter;
    public MeshFilter meshFilter
    {
        get {
            if (_meshFilter == null)
                _meshFilter = GetComponent<MeshFilter>();
            return _meshFilter;
        }
    }
    [SerializeField, HideInInspector]
    public Mesh hairMesh;


    private AddFur _parent; // 유니티 계층구조 parent와 별개. furmesh 인스턴스를 생성시킨 addfur 저장하기 위해 사용

    public AddFur Parent { get { return _parent; } set { _parent = value; } }

    private MeshRenderer _meshRenderer;
    public MeshRenderer meshRenderer
    {
        get
        {
            if (_meshRenderer == null)
                _meshRenderer = GetComponent<MeshRenderer>();
            return _meshRenderer;
        }
    }
    
    public List<HairParticle> particles;

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

    public List<Hair> _hairs = new();

    public void AddHair(Hair hair)
    {
        _hairs.Add(hair);
    }

    public void Clear()
    {
        _hairs.Clear();
    }
    
    protected List<Vector3> tangents = new();
    protected List<Vector2> uvs = new();
    protected List<Vector3> normals = new();
    protected List<Vector3> positions = new();
    protected List<int> indices = new();

    public Vector3[] getHairParticlePosArray()
    {
        List<Vector3> answer = new();
        foreach (var hair in _hairs)
        {
            answer.AddRange(hair._positions);
        }
        return answer.ToArray();
    }

    public virtual void UpdateMesh()
    {

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


        for (int k = 0; k < _hairs.Count; k++)
        {

            Hair hair = _hairs[k];

            int head = positions.Count;

            tangents.Add(hair._positions[0] - hair._positions[1]);
            for (int i = 1; i < hair._positions.Length; i++)
                tangents.Add(hair._positions[i - 1] - hair._positions[i]);

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

        particles = new List<HairParticle>();
        foreach (var hair in _hairs.Select((value, index) => (value, index)))
        {
            foreach (var dot in hair.value._positions.Select((value2, index2) => (value2,index2)))
            {
                Debug.Log(dot.value2 + ", " + transform.TransformPoint(dot.value2));
                HairParticle particle = new HairParticle(transform.TransformPoint(dot.value2));
                particles.Add(particle);
                if (dot.index2 == 0)
                {
                    particle.IsFixed = true;
                    //Debug.Log(particles.Count +"th particle is Fixed");
                }
            }
            //Debug.Log("Number of particles in " + hair.index + "th hair is " + particles.Count / (hair.index + 1));
        }

        //Debug.Log("Number of Hairs is " + _hairs.Count);

        hairMesh.Clear();

        hairMesh.SetVertices(positions.ToArray());
        hairMesh.SetUVs(3, uvs.ToArray());
        hairMesh.SetUVs(1, tangents.ToArray());
        hairMesh.SetUVs(2, normals.ToArray());
        hairMesh.SetIndices(indices.ToArray(), MeshTopology.Lines, 0);
        hairMesh.RecalculateBounds();
        this._parent.GetComponent<HairComponent>().SetHairInfo(this._parent, particles);


        
        float ratio = (float)Screen.width / Screen.height;
        //meshRenderer.material.SetFloat("_AspectRatio", ratio);

        /*
        tangents.Dispose();
        uvs.Dispose();
        normals.Dispose();
        positions.Dispose();
        indices.Dispose();
        */
    }   

    public virtual void UpdateHairPos(Vector3[] pos)
    {
        hairMesh.SetVertices(pos);
    }

}
