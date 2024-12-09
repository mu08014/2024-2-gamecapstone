using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class LBM : MonoBehaviour
{
    public GameObject spherePrefab;

    [Header("Grid Settings")]
    public Vector3Int gridSize;
    public Vector3 spawnBoxCenter;
    public Vector3Int spawnBox;
    public float simulationSpeed = 1.0f;
    public float boundValue = -0.35f;
    public float particleRadius;
    public float bounceMIN = 2.0f;
    public float repulsionStrength = 150.0f;

    [Header("Physics Settings")]
    public Vector3 gravity = new Vector3(0, -9.81f, 0);
    public float tau = 1.0f;
    [Range(0f, 1f)]
    public float velocityBlendFactor = 0.1f;

    [Header("Collision Settings")]
    public float collisionForce = 300;

    float[,,,] f; //grid
    float[,,,] feq; //평형분포
    float[,,,] fNew; //스트리밍 단계 다음 분포
    internal GameObject[,,] spheres; //물 입자의 좌표
    internal Vector3[,,] velocities; //물 입자의 속도

    const int NL = 19;

    public static readonly float[] Vx = { 0f, 1f, -1f, 0f, 0f, 0f, 0f, 1f, -1f, 1f, -1f, 1f, -1f, 1f, -1f, 0f, 0f, 0f, 0f };
    public static readonly float[] Vy = { 0f, 0f, 0f, 1f, -1f, 0f, 0f, 1f, -1f, -1f, 1f, 0f, 0f, 0f, 0f, 1f, -1f, 1f, -1f };
    public static readonly float[] Vz = { 0f, 0f, 0f, 0f, 0f, 1f, -1f, 0f, 0f, 0f, 0f, 1f, -1f, -1f, 1f, 1f, -1f, -1f, 1f };

    public static readonly float[] Weights = {
        1.0f / 3.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f, 1.0f / 18.0f,
        1.0f / 18.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f,
        1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f, 1.0f / 36.0f
    };

    private void OnValidate()
    {
        gridSize.x = Mathf.Max(gridSize.x, 0);
        gridSize.y = Mathf.Max(gridSize.y, 0);
        gridSize.z = Mathf.Max(gridSize.z, 0);

        spawnBox.x = Mathf.Clamp(spawnBox.x, 0, gridSize.x);
        spawnBox.y = Mathf.Clamp(spawnBox.y, 0, gridSize.y);
        spawnBox.z = Mathf.Clamp(spawnBox.z, 0, gridSize.z);

        Vector3 halfSpawnBox = new Vector3(spawnBox.x / 2f, spawnBox.y / 2f, spawnBox.z / 2f);

        float minX = -gridSize.x / 2f + halfSpawnBox.x;
        float maxX = gridSize.x / 2f - halfSpawnBox.x;
        spawnBoxCenter.x = Mathf.Clamp(spawnBoxCenter.x, minX, maxX);

        float minY = -gridSize.y / 2f + halfSpawnBox.y;
        float maxY = gridSize.y / 2f - halfSpawnBox.y;
        spawnBoxCenter.y = Mathf.Clamp(spawnBoxCenter.y, minY, maxY);

        float minZ = -gridSize.z / 2f + halfSpawnBox.z;
        float maxZ = gridSize.z / 2f - halfSpawnBox.z;
        spawnBoxCenter.z = Mathf.Clamp(spawnBoxCenter.z, minZ, maxZ);
    }

    void OnDrawGizmos()
    {
        Gizmos.color = Color.blue;
        Gizmos.DrawWireCube(Vector3.zero, gridSize);

        Gizmos.color = Color.cyan;
        Gizmos.DrawWireCube(spawnBoxCenter, spawnBox);

        //ㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡㅡ

        //if (spheres == null || velocities == null)
        //    return;

        //Gizmos.color = Color.red;
        //for (int x = 0; x < gridSize.x; x++)
        //{
        //    for (int y = 0; y < gridSize.y; y++)
        //    {
        //        for (int z = 0; z < gridSize.z; z++)
        //        {
        //            GameObject sphere = spheres[x, y, z];
        //            if (sphere != null)
        //            {
        //                // 속도 벡터를 시각화
        //                Gizmos.DrawLine(sphere.transform.position, sphere.transform.position + velocities[x, y, z]);
        //            }
        //        }
        //    }
        //}

    }

    void Init()
    {
        f = new float[gridSize.x, gridSize.y, gridSize.z, NL];
        feq = new float[gridSize.x, gridSize.y, gridSize.z, NL];
        fNew = new float[gridSize.x, gridSize.y, gridSize.z, NL];
        velocities = new Vector3[gridSize.x, gridSize.y, gridSize.z];
        spheres = new GameObject[gridSize.x, gridSize.y, gridSize.z];

        float initialDensity = 1.0f;
        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    for (int i = 0; i < NL; i++)
                    {
                        f[x, y, z, i] = initialDensity * Weights[i];
                    }

                    // 랜덤 속도
                    velocities[x, y, z] = new Vector3(
                        Random.Range(-1.0f, 1.0f),
                        Random.Range(-1.0f, 1.0f),
                        Random.Range(-1.0f, 1.0f)
                    );
                }
            }
        }
    }

    void Start()
    {
        Init();
        CreateSpheres();
    }

    float CalDensity(int x, int y, int z)
    {
        float density = 0;
        for (int i = 0; i < NL; i++)
        {
            density += f[x, y, z, i];
        }
        return density;
    }

    Vector3 CalVelocity(int x, int y, int z, float density)
    {
        Vector3 velocity = Vector3.zero;
        for (int i = 0; i < NL; i++)
        {
            velocity.x += f[x, y, z, i] * Vx[i];
            velocity.y += f[x, y, z, i] * Vy[i];
            velocity.z += f[x, y, z, i] * Vz[i];
        }
        return velocity / density;
    }

    void CalEquilibrium(int x, int y, int z, float density, Vector3 velocity)
    {
        for (int i = 0; i < NL; i++)
        {
            float cu = Vx[i] * velocity.x + Vy[i] * velocity.y + Vz[i] * velocity.z;
            float usq = velocity.sqrMagnitude;
            feq[x, y, z, i] = Weights[i] * density * (1.0f + 3.0f * cu + 4.5f * cu * cu - 1.5f * usq);
        }
    }

    //평형상태
    void Equilibrium()
    {
        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    float density = CalDensity(x, y, z);
                    Vector3 velocity = CalVelocity(x, y, z, density);
                    CalEquilibrium(x, y, z, density, velocity);

                    for (int i = 0; i < NL; i++)
                    {
                        f[x, y, z, i] += (feq[x, y, z, i] - f[x, y, z, i]) / tau;
                    }
                }
            }
        }
    }

    //평형 단계의 값을 물 입자에게 전달
    void Streaming()
    {
        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    for (int i = 0; i < NL; i++)
                    {
                        int targetX = (x + (int)Vx[i] + gridSize.x) % gridSize.x;
                        int targetY = (y + (int)Vy[i] + gridSize.y) % gridSize.y;
                        int targetZ = (z + (int)Vz[i] + gridSize.z) % gridSize.z;

                        fNew[targetX, targetY, targetZ, i] = f[x, y, z, i];
                    }
                }
            }
        }

        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    for (int i = 0; i < NL; i++)
                    {
                        f[x, y, z, i] = fNew[x, y, z, i];
                    }
                }
            }
        }
    }

    void CreateSpheres()
    {
        Vector3 centerspawn = spawnBoxCenter - new Vector3(spawnBox.x / 2f, spawnBox.y / 2f, spawnBox.z / 2f);

        for (int x = 0; x < spawnBox.x; x++)
        {
            for (int y = 0; y < spawnBox.y; y++)
            {
                for (int z = 0; z < spawnBox.z; z++)
                {
                    GameObject sphere = Instantiate(spherePrefab, new Vector3(x, y, z) + centerspawn, Quaternion.identity);

                    sphere.transform.localScale = Vector3.one * particleRadius * 2f;

                    Collision collisionScript = sphere.GetComponent<Collision>();
                    if (collisionScript != null)
                    {
                        collisionScript.lbmScript = this;
                    }
                    spheres[x, y, z] = sphere;
                }
            }
        }
    }

    void AddGravity()
    {
        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    if (spheres[x, y, z] != null)
                    {
                        velocities[x, y, z] += gravity * Time.fixedDeltaTime;
                    }
                }
            }
        }
    }

    void Boundary()
    {
        Vector3 topRight = new Vector3(gridSize.x / 2f, gridSize.y / 2f, gridSize.z / 2f);
        Vector3 bottomLeft = -topRight;

        for (int x = 0; x < spawnBox.x; x++)
        {
            for (int y = 0; y < spawnBox.y; y++)
            {
                for (int z = 0; z < spawnBox.z; z++)
                {
                    GameObject sphere = spheres[x, y, z];
                    if (sphere == null) continue;

                    Vector3 position = sphere.transform.position;
                    Vector3 velocity = velocities[x, y, z];

                    if (position.x - particleRadius < bottomLeft.x)
                    {
                        velocity.x *= boundValue;
                        position.x = bottomLeft.x + particleRadius;
                    }

                    if (position.x + particleRadius > topRight.x)
                    {
                        velocity.x *= boundValue;
                        position.x = topRight.x - particleRadius;
                    }

                    if (position.y - particleRadius < bottomLeft.y)
                    {
                        if (Mathf.Abs(velocity.y) < bounceMIN)
                        {
                            velocity.y = 0f;
                            position.y = bottomLeft.y + particleRadius;
                        }
                        else
                        {
                            velocity.y *= boundValue;
                            position.y = bottomLeft.y + particleRadius;
                        }
                    }

                    if (position.y + particleRadius > topRight.y)
                    {

                        velocity.y *= boundValue;
                        position.y = topRight.y - particleRadius;
                    }

                    if (position.z - particleRadius < bottomLeft.z)
                    {
                        velocity.z *= boundValue;
                        position.z = bottomLeft.z + particleRadius;
                    }

                    if (position.z + particleRadius > topRight.z)
                    {
                        velocity.z *= boundValue;
                        position.z = topRight.z - particleRadius;
                    }

                    sphere.transform.position = position;
                    velocities[x, y, z] = velocity;
                }
            }
        }
    }

    void Overlap()
    {
        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    GameObject particle1 = spheres[x, y, z];
                    if (particle1 == null) continue;

                    Vector3 pos1 = particle1.transform.position;

                    //탐색
                    for (int i = Mathf.Max(x - 5, 0); i <= Mathf.Min(x + 5, gridSize.x - 1); i++)
                    {
                        for (int j = Mathf.Max(y - 5, 0); j <= Mathf.Min(y + 5, gridSize.y - 1); j++)
                        {
                            for (int k = Mathf.Max(z - 5, 0); k <= Mathf.Min(z + 5, gridSize.z - 1); k++)
                            {
                                if (i == x && j == y && k == z) continue;

                                GameObject particle2 = spheres[i, j, k];
                                if (particle2 == null) continue;

                                Vector3 pos2 = particle2.transform.position;

                                // 두 입자 간 거리 계산
                                Vector3 direction = pos2 - pos1;
                                float distance = direction.magnitude;

                                if (distance < 2 * particleRadius)
                                {
                                    // 반발력 계산
                                    Vector3 repulsionForce = direction.normalized * (2 * particleRadius - distance) * repulsionStrength;

                                    velocities[x, y, z] -= repulsionForce * Time.fixedDeltaTime * 0.5f;
                                    velocities[i, j, k] += repulsionForce * Time.fixedDeltaTime * 0.5f;

                                    Vector3 correction = direction.normalized * (2 * particleRadius - distance) * 0.1f;
                                    spheres[x, y, z].transform.position -= correction * 0.5f;
                                    spheres[i, j, k].transform.position += correction * 0.5f;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void UpdateLBMVelocities()
    {
        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    if (spheres[x, y, z] != null)
                    {
                        float density = CalDensity(x, y, z);
                        Vector3 velocityFromLBM = CalVelocity(x, y, z, density);

                        // 기존 속도와 LBM 속도 조합
                        velocities[x, y, z] = Vector3.Lerp(velocities[x, y, z], velocityFromLBM, velocityBlendFactor);
                    }
                }
            }
        }
    }

    void Update()
    {
        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    if (spheres[x, y, z] != null)
                    {
                        spheres[x, y, z].transform.position += velocities[x, y, z] * Time.fixedDeltaTime * simulationSpeed;
                    }
                }
            }
        }
    }

    void LogVelocities()
    {
        for (int x = 0; x < gridSize.x; x++)
        {
            for (int y = 0; y < gridSize.y; y++)
            {
                for (int z = 0; z < gridSize.z; z++)
                {
                    if (spheres[x, y, z] != null)
                    {
                        Vector3 velocity = velocities[x, y, z];
                        Debug.Log(spheres[x, y, z].transform.position);
                        Debug.Log(velocity);
                    }
                }
            }
        }
    }

    void FixedUpdate()
    {
        Equilibrium();
        Streaming();
        UpdateLBMVelocities();
        AddGravity();
        Overlap();
        Update();
        Boundary();
        LogVelocities();
    }
}