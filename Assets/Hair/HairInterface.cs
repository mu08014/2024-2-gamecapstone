using System.Collections;
using System.Collections.Generic;
using System.IO;
using Unity.VisualScripting;
using UnityEditor.Overlays;
using UnityEngine;

public class HairInterface : MonoBehaviour
{
    private List<List<GameObject>> hairs = new List<List<GameObject>>();

    public List<List<GameObject>> Hairs
    {
        get
        {
            return hairs;
        }
    }

    private void Start()
    {
        MakeHair();

        LoadHairSimulation loadHairSimulation = new LoadHairSimulation();
        loadHairSimulation.loadDERSimulation();
        Debug.Log("0");
        VectorXs v = loadHairSimulation.scene.getX().Clone();
        Debug.Log("0");
        for (int i = 0; i < v.Size; i += 3)
        {
            Debug.Log("x: " + v[i].ToString() + ", y: " + v[i + 1].ToString() + ", z: " + v[i + 2].ToString());
        }
    }

    private void Update()
    {

    }

    public List<int> getHairsCount()
    {
        List<int> result = new List<int>();
        foreach (var hair in hairs)
        {
            result.Add(hair.Count);
        }
        return result;
    }

    void MakeHair()
    {
        var filePath = Path.Combine(Application.dataPath, "Hair/HairData.txt");
        try
        {
            if (File.Exists(filePath))
            {
                List<GameObject> list = new List<GameObject>();
                foreach (string line in File.ReadLines(filePath))
                {
                    if (line.Contains("END"))
                    {
                        return;
                    }
                    else if (string.IsNullOrEmpty(line))
                    {
                        hairs.Add(list);
                        list = new List<GameObject>();
                        continue;
                    }
                    string[] data = line.Split(' ');
                    float[] pos = new float[3];
                    float[] vel = new float[3];
                    foreach (string value in data)
                    {
                        if (value == "")
                        {
                            continue;
                        }
                        else if (value[..2].Equals("x="))
                        {
                            string[] ps = value[2..].Split(',');
                            int i = 0;
                            foreach (string position in ps)
                            {
                                pos[i++] = float.Parse(position);
                            }

                        }
                        else if (value[..2].Equals("v="))
                        {
                            string[] vs = value[2..].Split(',');
                            int i = 0;
                            foreach (string velocity in vs)
                            {
                                vel[i++] = float.Parse(velocity);
                            }
                        }
                    }
                    GameObject newHairObject = new GameObject("hairparticle");
                    newHairObject.AddComponent<HairParticle>();
                    newHairObject.GetComponent<HairParticle>().Position = new Vector3(pos[0], pos[1], pos[2]);
                    newHairObject.GetComponent<HairParticle>().Velocity = new Vector3(vel[0], vel[1], vel[2]);
                    newHairObject.AddComponent<MeshFilter>();
                    newHairObject.GetComponent<MeshFilter>().mesh = Resources.GetBuiltinResource<Mesh>("Sphere.fbx");
                    newHairObject.AddComponent<MeshRenderer>();
                    newHairObject.GetComponent<Renderer>().material = new Material(Shader.Find("Standard"));
                    newHairObject.transform.localScale = new Vector3(0.1f, 0.1f, 0.1f);

                    list.Add(newHairObject);
                    Debug.Log("create hair");
                }
            }
            else
            {
                Debug.Log("파일을 찾을 수 없습니다: " + filePath);
            }
        }
        catch (System.Exception e)
        {
            Debug.Log("파일을 읽는 중 오류 발생: " + e.Message);
        }
    }

}


