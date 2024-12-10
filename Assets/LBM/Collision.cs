using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Collision : MonoBehaviour
{
    public LBM lbmScript; // LBM 스크립트를 참조

    private void OnTriggerStay(Collider other)
    {
        for (int x = 0; x < lbmScript.gridSize.x; x++)
        {
            for (int y = 0; y < lbmScript.gridSize.y; y++)
            {
                for (int z = 0; z < lbmScript.gridSize.z; z++)
                {
                    GameObject sphere = lbmScript.spheres[x, y, z];
                    if (sphere == null) continue;

                    Vector3 spherePos = sphere.transform.position;
                    Vector3 closestPoint = other.ClosestPoint(spherePos);
                    float dist = Vector3.Distance(spherePos, closestPoint);

                    float boundarySize = Mathf.Min(other.bounds.extents.x, other.bounds.extents.y, other.bounds.extents.z);

                    if (dist < lbmScript.particleRadius) // 겹침 발생 기준을 입자 반지름으로 설정
                    {
                        Vector3 direction = (spherePos - closestPoint).normalized;
                        float overlapDistance = lbmScript.particleRadius - dist;

                        // 위치 보정
                        sphere.transform.position += direction * overlapDistance;

                        // 반발력 계산
                        float scale = Mathf.Clamp01((boundarySize - dist) / boundarySize);
                        float forceStrength = lbmScript.collisionForce * scale;
                        Vector3 repulsionForce = direction * forceStrength * Time.deltaTime;

                        // 속도 보정
                        lbmScript.velocities[x, y, z] = Vector3.Lerp(lbmScript.velocities[x, y, z], repulsionForce, 0.5f);

                        // 디버그용 로그 출력
                        //Debug.Log($"Collision corrected at ({x}, {y}, {z}). Overlap: {overlapDistance}, Force: {repulsionForce}");
                    }
                }
            }
        }
    }
}
