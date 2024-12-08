using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Collision : MonoBehaviour
{
    public LBM lbmScript;

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

                    if (dist < lbmScript.particleRadius)
                    {
                        Vector3 direction = (spherePos - closestPoint).normalized;
                        float overlapDistance = lbmScript.particleRadius - dist;

                        sphere.transform.position += direction * overlapDistance;

                        float scale = Mathf.Clamp01((boundarySize - dist) / boundarySize);
                        float forceStrength = lbmScript.collisionForce * scale;
                        Vector3 repulsionForce = direction * forceStrength * Time.deltaTime;

                        lbmScript.velocities[x, y, z] = Vector3.Lerp(lbmScript.velocities[x, y, z], repulsionForce, 0.5f);
                    }
                }
            }
        }
    }
}
