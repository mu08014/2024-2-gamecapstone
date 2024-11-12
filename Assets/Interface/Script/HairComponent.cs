using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using Unity.VisualScripting.Antlr3.Runtime.Tree;
using UnityEditor.SceneManagement;
using UnityEngine;
using UnityEngine.UIElements;

[ExecuteAlways]
public class HairComponent : MonoBehaviour {

    [SerializeField]
    private StrandParameters[] strandHairs;   //Parameter list for A strand hair set

    [SerializeField] 
    GameObject _furmeshprefab;

    [HideInInspector]
    public List<HairParticle> hairParticles;

    private AddFur[] _addFur;

    [HideInInspector]
    public Vectors[] m_hair_x;

    [ HideInInspector]
    public Vectors[] m_hair_v;

    [HideInInspector]
    public List<bool> m_hair_fixed;

    [HideInInspector]
    public int numofStrands = 0;

    [HideInInspector]
    public int[] StartofEachParam; //strandparameter별로 생성한 particle을 m_hair_x, m_hair_v, m_hair fixed에 몰아넣으므로, 각 strandparameter의 particle이 어디부터 시작인지 알기 위한 index 값을 저장

    [HideInInspector]
    public int[] EndofEachParam;

    public ref StrandParameters[] GetStrandHairParameters()
    {
        return ref strandHairs;
    }

    public int GetNumOfHairParam()
    {
        return strandHairs.Length;
    }
    public bool IsAnyHairExist()
    {
        return (strandHairs.Length > 0);
    }

    public void Clear()
    {
        hairParticles = null;
        m_hair_x = null;
        m_hair_v = null;   
        m_hair_fixed = null;
        numofStrands = 0;   
    }

    public void WhenGenerationButtonClicked() //커스텀 에디터 버튼 전용 함수
    {
        if (!IsAnyHairExist())
        {
            Debug.LogError("There's no StrandHairParameter Instance!");
            return;
        }
        this.Clear();
        
        _addFur = GetComponents<AddFur>(); // Addfur component를 haircomponent에 생성된 strandparameter 수만큼 초기화;
        if (_addFur.Length !=  GetNumOfHairParam())
        {
            var temp = new AddFur[GetNumOfHairParam()];
            for (int i = 0; i < GetNumOfHairParam(); i++)
            {
                if (i < _addFur.Length)
                {
                    temp[i] = _addFur[i];
                }
                else
                {
                    temp[i] = gameObject.AddComponent<AddFur>();
                }

            }
            _addFur = temp;
        }

        StartofEachParam = new int[] { GetNumOfHairParam() };
        EndofEachParam = new int[GetNumOfHairParam()];  

        foreach (var fur in _addFur) // hair의 패러미터 종류 수만큼 fur 생성
        {
            fur.initMesh(_furmeshprefab);
            Debug.Log("Hair Mesh" + /* index*/ "Created!");
        }
    }
    
    public void SetHairInfo(AddFur fur,  List<HairParticle> particles) //Furmesh에서 생성된 mesh의 정점 정보를 파티클로 변환해 저장
    {
        if (hairParticles == null)
        {
            hairParticles = new List<HairParticle>();
        }
        for (int i = 0; i < _addFur.Length; i++)
        {
            if (fur == _addFur[i])
            {
                Debug.Log("It's from" + i + "th _addfur component!");
                Mesh tempMesh = _addFur[i].furmesh.meshFilter.sharedMesh;
                int OriginLength = hairParticles.Count;
                hairParticles.AddRange(particles);
                if (m_hair_x == null) 
                {
                    m_hair_x = new Vectors[hairParticles.Count];
                    m_hair_v = new Vectors[hairParticles.Count];
                    m_hair_fixed = new List<bool>();
                }


                Vectors[] new_x = new Vectors[hairParticles.Count];
                if (m_hair_x.Length == new_x.Length) { StartofEachParam[i] = 0; EndofEachParam[i] = new_x.Length - 1; }
                else { StartofEachParam[i] = m_hair_x.Length; EndofEachParam[i] = new_x.Length - 1; }
                Debug.Log("Start particle index of " + i + "th strandparameter is " + StartofEachParam[i]);

                for (int j = 0; j < m_hair_x.Length; j++)
                {
                    new_x[j] = m_hair_x[j];
                }
                foreach (var vert in tempMesh.vertices.Select((value, index) => (value, index)))
                {
                    new_x[StartofEachParam[i] + vert.index] = new Vectors(3);
                    new_x[StartofEachParam[i] + vert.index][0] = vert.value.x;
                    new_x[StartofEachParam[i] + vert.index][1] = vert.value.y;
                    new_x[StartofEachParam[i] + vert.index][2] = vert.value.z;
                }
                m_hair_x = new_x;

                m_hair_v = new Vectors[hairParticles.Count];
                for (int j = 0; j < hairParticles.Count; j++)
                {
                    m_hair_v[j] = new Vectors(3);
                    m_hair_v[j].setZero();
                }
                int strandId = -1;
                int fixedCount = 0;
                for (int j = 0; j < hairParticles.Count; j++)
                {
                    if (strandId != hairParticles[j].strandID || j == OriginLength)
                    {
                        m_hair_fixed.Add(true);
                        strandId = hairParticles[j].strandID;
                        fixedCount++;
                    }
                    else m_hair_fixed.Add(false);
                }
                Debug.Log("a is " + m_hair_fixed.Count);

                numofStrands += fixedCount;
                Debug.Log("Number of Fixed Particle is " + numofStrands);
            }

        }
        
    }
    
    // Start is called before the first frame update
    void Start()
    {
        Debug.Log("Number of m_hair_x are " + m_hair_x.Length);

        Debug.Log("Number of hairparticles are " + hairParticles.Count);
    }



    
}


   
