using UnityEngine;

public class QuantumParticle : MonoBehaviour
{

    float h = CGHscale.h;
    public Rigidbody rb;
    Vector3 deltaPosition;

    void Update()
    {
        
        float lengthDeltaScale = Mathf.Sqrt(h * Time.deltaTime / rb.mass);
        float dx = Random.Range(-1f, 1f);
        float dy = Random.Range(-1f, 1f);
        float dz = Random.Range(-1f, 1f);
        deltaPosition = new Vector3(dx, dy, dz) * lengthDeltaScale;
        rb.position += deltaPosition;
        //attractor.Translate(deltaPosition, Space.World);
    }

//    public static float RandomStandardNormal()
//    {
//        float u1 = 1.0f - Random.value; // avoid 0
//        float u2 = 1.0f - Random.value;
//        float randStdNormal = Mathf.Sqrt(-2.0f * Mathf.Log(u1)) * Mathf.Sin(2.0f * Mathf.PI * u2);
//        return randStdNormal;
//    }


}
