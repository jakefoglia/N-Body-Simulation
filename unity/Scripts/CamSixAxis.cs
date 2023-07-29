
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CamSixAxis : MonoBehaviour
{
	CharacterController cc;

	public float speed;
	public float lookspeed;
	public float rollspeed;
	public bool useRoll;
	public bool useJoystick;

	bool isAxisEnabled(string axisName)
	{
		try
		{
			Input.GetAxis(axisName);
			return true;
		}
		catch (UnityException e)
		{
			Debug.Log("CameraRig.cs requires Roll axis to be set up in Input Manager to incorporate Roll feature.");
			Debug.Log("CameraRig.cs requires Joystick X and Joystick Y axes to be set up in Input Manager to incorporate controller-based look.");
			return false;
		}
	}

	// Use this for initialization
	void Start()
	{
		speed = 20.0f;
		lookspeed = 80.0f;
		rollspeed = 40.0f;
		useRoll = true;
		useJoystick = true;

		cc = GetComponent<CharacterController>();
		useRoll = isAxisEnabled("Roll");
		useJoystick = isAxisEnabled("Joystick X") && isAxisEnabled("Joystick Y");
	}

	// Update is called once per frame
	void Update()
	{
		float forward = Input.GetAxis("Vertical");
		float strafe = Input.GetAxis("Horizontal");
		float jump = Input.GetAxis("Jump");
		float lookside = Input.GetAxis("Mouse X");
		float lookup = -Input.GetAxis("Mouse Y");

		cc.Move(
			speed * Time.deltaTime * (
				forward * transform.forward
				+
				strafe * transform.right
				+
				jump * transform.up
			)
		);

		if (Input.GetMouseButton(1))
		{
			Vector3 rotAxis = (lookside * transform.up + lookup * transform.right);
			float rotValue = rotAxis.magnitude;
			transform.Rotate(rotAxis.normalized, rotValue * lookspeed * Time.deltaTime, Space.World);
		}
		else if (useJoystick)
		{
			lookside = Input.GetAxis("Joystick X");
			lookup = Input.GetAxis("Joystick Y");
			Vector3 rotAxis = (lookside * transform.up + lookup * transform.right);
			float rotValue = rotAxis.magnitude;
			transform.Rotate(rotAxis.normalized, rotValue * lookspeed * Time.deltaTime, Space.World);
		}
		if (useRoll)
		{
			float roll = -Input.GetAxis("Roll");
			transform.Rotate(transform.forward, rollspeed * roll * Time.deltaTime, Space.World);
		}
	}
}