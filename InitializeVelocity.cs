using UnityEngine;

public class InitializeVelocity : MonoBehaviour
{
    public Rigidbody rb;
    public Vector3 velocity;
    // Start is called once before the first execution of Update after the MonoBehaviour is created
    void Start()
    {
        rb.linearVelocity = velocity;
    }

}
