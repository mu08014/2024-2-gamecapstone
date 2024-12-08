using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;

public class RigidbodyMovement : MonoBehaviour
{

    private Rigidbody m_RigidBody;
    private MeshCollider m_MeshCollider;
    public float speed = 1.0f;
    public float jumpPower = 5.0f;
    public float rayLen = 1.0f;
    public bool isGrounded = true;

    private void Awake()
    {
        m_RigidBody = GetComponent<Rigidbody>();
    }

    // Start is called before the first frame update
    // Update is called once per frame
    void FixedUpdate()
    {   
        Vector3 center = transform.position;
        float moveHorizontal = Input.GetAxis("Horizontal");
        float moveVertical = Input.GetAxis("Vertical");


        Vector3 movement = new Vector3(moveHorizontal, 0, moveVertical);
        m_RigidBody.AddForce(movement * speed);


        if(Physics.Raycast(center, Vector3.down, rayLen))
        {
            Debug.DrawRay(center, Vector3.down * rayLen, Color.red);
            isGrounded = true;
        }
        else
        {
            Debug.DrawRay(center, Vector3.down * rayLen, Color.green);
            isGrounded= false;
        }
        if (isGrounded && Input.GetButtonDown("Jump"))
        {
            m_RigidBody.AddForce(Vector3.up * jumpPower, ForceMode.Impulse);
        }

    }
}
