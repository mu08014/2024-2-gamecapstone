using System.Collections;
using System.Collections.Generic;
using System.IO;
using Unity.VisualScripting;
using UnityEditor.Overlays;
using UnityEngine;

public class HairInterface : MonoBehaviour
{
    private List<List<GameObject>> hairs = new List<List<GameObject>>();
    private LoadHairSimulation loadHairSimulation = new LoadHairSimulation();
    private int hairParticleCount;

    public float hairUpdateTime = 0.1f;
    public float executeSpeed = 1f;

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
        hairParticleCount = getHairParticleCount();

        loadHairSimulation.loadDERSimulation();
        for (int i = 0; i < hairParticleCount; i++)
        {
            Debug.Log("Particle " + i);
            Vectors v = loadHairSimulation.scene.getPosition(i).Clone();
            Debug.Log("x: " + v[0].ToString() + " y: " + v[1].ToString() + " z: " + v[2].ToString());

            Vectors vel = loadHairSimulation.scene.getVelocity(i);
            Debug.Log("vx : " + vel[0] + " vy : " + vel[1] + " vz : " + vel[2]);
        }

        StartCoroutine(HairUpdateCoroutine());
    }

    private void Update()
    {
        //HairUpdate();
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

    public int getHairParticleCount()
    {
        int result = 0;
        foreach (var hair in hairs)
        {
            foreach (var particle in hair)
                result++;
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
                    bool isFix = false;
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
                        else if (value[..2].Equals("f="))
                        {
                            string fs = value[2..];
                            if (fs[0] == '1')
                            {
                                isFix = true;
                            }
                        }
                    }
                    GameObject newHairObject = new GameObject("hairparticle");
                    newHairObject.AddComponent<HairParticle>();
                    newHairObject.GetComponent<HairParticle>().Position = new Vector3(pos[0], pos[1], pos[2]);
                    newHairObject.GetComponent<HairParticle>().Velocity = new Vector3(vel[0], vel[1], vel[2]);
                    newHairObject.GetComponent<HairParticle>().Fix = isFix;
                    newHairObject.AddComponent<MeshFilter>();
                    newHairObject.GetComponent<MeshFilter>().mesh = Resources.GetBuiltinResource<Mesh>("Sphere.fbx");
                    newHairObject.AddComponent<MeshRenderer>();
                    newHairObject.GetComponent<Renderer>().material = new Material(Shader.Find("Standard"));
                    newHairObject.transform.localScale = new Vector3(0.3f, 0.3f, 0.3f);

                    list.Add(newHairObject);
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

    public void HairUpdate()
    {
        loadHairSimulation.hairUpdate(Time.deltaTime);

        int row = 0;
        int col = 0;
        for (int i = 0; i < hairParticleCount; i++)
        {
            Vectors x_vec = loadHairSimulation.scene.getPosition(i).Clone();
            Vectors v_vec = loadHairSimulation.scene.getVelocity(i).Clone();
            hairs[row][col].GetComponent<HairParticle>().Position = new Vector3((float)x_vec[0], (float)x_vec[1], (float)x_vec[2]);
            hairs[row][col++].GetComponent<HairParticle>().Velocity = new Vector3((float)v_vec[0], (float)v_vec[1], (float)v_vec[2]);
            if (col >= hairs[row].Count)
            {
                row++;
                col = 0;
            }
        }
    }

    IEnumerator HairUpdateCoroutine()
    {
        while (true)
        {
            loadHairSimulation.hairUpdate(hairUpdateTime * executeSpeed);

            int row = 0;
            int col = 0;
            for (int i = 0; i < hairParticleCount; i++)
            {
                Vectors x_vec = loadHairSimulation.scene.getPosition(i).Clone();
                Vectors v_vec = loadHairSimulation.scene.getVelocity(i).Clone();
                hairs[row][col].GetComponent<HairParticle>().Position = new Vector3((float)x_vec[0], (float)x_vec[1], (float)x_vec[2]);
                hairs[row][col++].GetComponent<HairParticle>().Velocity = new Vector3((float)v_vec[0], (float)v_vec[1], (float)v_vec[2]);
                if (col >= hairs[row].Count)
                {
                    row++;
                    col = 0;
                }
            }
            yield return new WaitForSeconds(hairUpdateTime);
        }
    }

}


