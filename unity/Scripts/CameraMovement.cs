using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CameraMovement : MonoBehaviour
{
    public float turnSpeed;
    public float moveSpeed;
    public float minTurnAngle;
    public float maxTurnAngle;
    public float minFov;
    public float maxFov;
    
    private float rotX;

    private void Start()
    {
        Cursor.visible = false;

        turnSpeed = 0.01f;
        moveSpeed = 1000f;
        minTurnAngle = -90.0f;
        maxTurnAngle = 90.0f;
        minFov = 1.0f;
        maxFov = 120.0f;
    }

    void Update()
    {
        MouseAiming();
        KeyboardMovement();
    }
    void MouseAiming()
    {
        float fov = gameObject.GetComponent<Camera>().fieldOfView;
        Debug.Log("FOV " + fov);
        // get the mouse inputs
        float y = Input.GetAxis("Mouse X") * turnSpeed * fov;
        rotX += Input.GetAxis("Mouse Y") * turnSpeed * fov;
        // clamp the vertical rotation
        rotX = Mathf.Clamp(rotX, minTurnAngle, maxTurnAngle);
        // rotate the camera
        transform.eulerAngles = new Vector3(-rotX, transform.eulerAngles.y + y, 0);

        // zoom the camera
        float newFOV = fov * Mathf.Pow(1.05f, -Input.mouseScrollDelta.y);
        gameObject.GetComponent<Camera>().fieldOfView = Mathf.Clamp(newFOV, minFov, maxFov);
        
    }
    void KeyboardMovement()
    {
        Vector3 dir = new Vector3(0, 0, 0);
        dir.x = Input.GetAxis("Horizontal");
        dir.y = (Input.GetKey(KeyCode.Space) ? 1f : 0f) - (Input.GetKey(KeyCode.LeftControl) ? 1f : 0f);
        dir.z = Input.GetAxis("Vertical");
        transform.Translate(dir * moveSpeed * Time.deltaTime);
    }
}
