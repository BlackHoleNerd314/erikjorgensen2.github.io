using UnityEngine;

public class ObjectMover : MonoBehaviour
{
    [Header("Velocity in World Frame")]
    public Vector3 velocityWorld = new Vector3(0.0f, 0f, 0f);

    [Header("Initial Conditions")]
    public Vector3 initialPosition = Vector3.zero;

    //public float c = 1f;//30.44f;
    float c = CGHscale.c;

    public Rigidbody rb;

    public Vector3 GetVelocity()
    {
        return velocityWorld * c;
    }

    public float GetC()
    {
        return c;
    }

    void Start()
    {
        rb.position = initialPosition;
        float v = velocityWorld.magnitude;
        float gamma = 1f / Mathf.Sqrt(1f - (v * v) / (c * c));
        rb.linearVelocity = velocityWorld * gamma;
        Time.timeScale = 1f;
    }


}

