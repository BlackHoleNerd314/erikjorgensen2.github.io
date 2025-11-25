using UnityEngine;

public class PlayerMover : MonoBehaviour
{
    [Header("Velocity in World Frame")]
    public Vector3 velocityWorld = new Vector3(0f, 0f, 0f);

    [Header("Initial Conditions")]
    public Vector3 initialPosition = Vector3.zero;

    [Header("Speed of Light")]
    float c = CGHscale.c;
    //float c = 1f;//30.44f;
    float G = CGHscale.G;

    float acceleration = 0.0216f; // Proper acceleration (in player's rest frame)
    public Rigidbody player;
    public Transform X;
    public Transform Y;
    public Transform Z;
    private Vector3 deltaV = new Vector3(0, 0, 0);

    public Vector3 GetVelocity()
    {
        return velocityWorld * c;
    }

    public float GetC()
    {
        return c;
    }

    private void Start()
    {
        player.position = initialPosition;
    }

    Vector3 acceleration0;

    void Update()
    {

        //
        
        //Attractor script = GetComponent<Attractor>();

        acceleration0 = Vector3.zero;
        foreach (Rigidbody source in Object.FindObjectsByType<Rigidbody>(FindObjectsSortMode.None))
        {
            Vector3 direction = (player.position - source.position);
            float r = direction.magnitude;
            if (r > 0.01f)
            {
                float a = G * source.mass / (r * r);
                acceleration0 += (direction.normalized * a);
            }
            
        }


        //

        float dt = Time.deltaTime;

        Vector3 position = player.position;
        Vector3 X0 = (X.position - position).normalized;
        Vector3 Y0 = (Y.position - position).normalized;
        Vector3 Z0 = (Z.position - position).normalized;

        deltaV = Vector3.zero;
        // Transform input to world frame
        if (Input.GetKey(KeyCode.W))
        {
            deltaV = acceleration * Z0;
        }
        if (Input.GetKey(KeyCode.S))
        {
            deltaV = -acceleration * Z0;
        }
        if (Input.GetKey(KeyCode.D))
        {
            deltaV = acceleration * X0;
        }
        if (Input.GetKey(KeyCode.A))
        {
            deltaV = -acceleration * X0;
        }
        if (Input.GetKey(KeyCode.R))
        {
            deltaV = acceleration * Y0;
        }
        if (Input.GetKey(KeyCode.F))
        {
            deltaV = -acceleration * Y0;
        }
        // Apply relativistic velocity addition
        Vector3 deltaV0 = deltaV - acceleration0 * dt;
        velocityWorld = RelativisticVelocityAddition(velocityWorld, deltaV0, c);
        float v = velocityWorld.magnitude;
        float gamma = 1f / Mathf.Sqrt(1f - (v * v) / (c * c));
        player.linearVelocity = velocityWorld * gamma;
        //Time.timeScale = 1f;


    }

    Vector3 RelativisticVelocityAddition(Vector3 u, Vector3 v, float c)
    {
        float c2 = c * c;
        float u2 = u.sqrMagnitude;
        float gamma_u = 1.0f / Mathf.Sqrt(1.0f - u2 / c2);
        float dot = Vector3.Dot(u, v);

        float denom = 1.0f + dot / c2;
        Vector3 term1 = u;
        Vector3 term2 = v / gamma_u;
        Vector3 term3 = (gamma_u / (gamma_u + 1.0f)) * (dot / c2) * u;

        return (term1 + term2 + term3) / denom;
    }



}

