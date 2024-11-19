using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UIElements;

public class Movement3D : MonoBehaviour
{

    [SerializeField]
    private float movespeed = 30.0f;
    [SerializeField]
    private float rotationSpeed = 30.0f;
    private Vector3 moveDirection;
    private float decelration = 0.8f;
   

    private CharacterController mCharacterController;
    // Start is called before the first frame update

    private void Awake()
    {
        mCharacterController = GetComponent<CharacterController>();
        moveDirection = Vector3.zero;
    }
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        if (mCharacterController.isGrounded)
        {
            float x = Input.GetAxisRaw("Horizontal");
            float z = Input.GetAxisRaw("Vertical");
            if (z != 0) { decelration = 0.8f; }
            else { decelration = 1.0f; }
            moveDirection = new Vector3(x,0,z);
            if (moveDirection != Vector3.zero) 
            {
                transform.forward = Vector3.Lerp(transform.forward, moveDirection, Time.deltaTime * rotationSpeed);
            }
            if (Input.GetAxisRaw("Jump") > 0)
            {
                Debug.Log("jump");
                moveDirection.y = 4.0f;
            }
        }
        mCharacterController.Move(moveDirection * movespeed * Time.deltaTime*decelration);

    }

    private void FixedUpdate()
    {
        moveDirection.y += Physics.gravity.y * Time.deltaTime * 0.8f;
    }
}
