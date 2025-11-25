using UnityEngine;

public class ColorChanger : MonoBehaviour
{
    public Material material;
    public Rigidbody Player;
    public Rigidbody Object;
    public Material surface;
    //float c = //1f;//30.44f;
    float c = CGHscale.c;
    void Update()
    {// Velocity and Lorentz inverse boost
        Vector3 vPlayer0 = -Player.linearVelocity;//.GetComponent<PlayerMover>().velocityWorld;
        Vector3 vObject0 = -Object.linearVelocity;//.GetComponent<ObjectMover>().velocityWorld;
        float vPlayer2 = vPlayer0.x * vPlayer0.x + vPlayer0.y * vPlayer0.y + vPlayer0.z * vPlayer0.z;
        float EPlayer1 = Mathf.Sqrt(1 + vPlayer2);
        float vObject2 = vObject0.x * vObject0.x + vObject0.y * vObject0.y + vObject0.z * vObject0.z;
        float EObject1 = Mathf.Sqrt(1 + vObject2);
        Vector3 vPlayer = vPlayer0 / EPlayer1;
        Vector3 vObject = vObject0 / EObject1;
        Matrix4x4 Lplayer = LorentzMatrix.GetInverseBoost(vPlayer, c);
        Matrix4x4 Lplayer_inv = LorentzMatrix.GetInverseBoost(-vPlayer, c);
        Matrix4x4 Lobject = LorentzMatrix.GetInverseBoost(vObject, c);
        Matrix4x4 Lobject_inv = LorentzMatrix.GetInverseBoost(-vObject, c);
        // Compute Relative Velocity
        Vector4 objectVel = new Vector4(1f, 0f, 0f, 0f);
        Vector4 worldVel = Lobject * objectVel;
        Vector4 playerVel = Lplayer_inv * worldVel;
        // Find Relative Position
        Vector3 RelPos = Object.position - Player.position;
        Vector4 objectRelPos = new Vector4(- RelPos.magnitude,RelPos.x,RelPos.y,RelPos.z);
        // Calculate Doppler Shift(in player frame)
        Vector4 worldRelPos = Lobject * objectRelPos;
        Vector4 playerRelPos = Lplayer_inv * worldRelPos;
        float TimeRetarded = - playerRelPos.x;
        float Distance = RelPos.magnitude;
        float ColorTemp = 5778f * (TimeRetarded/Distance);
        // Use Color
        Color color = BlackbodyToRGB(ColorTemp);
        //
        Color surfaceColor = color; // or "_Color" in Standard
        Color Tint = surface.color;
        Color Emission = surfaceColor * Tint;
        GetComponent<Renderer>().material.color = Emission;
        material.SetColor("_EmissionColor", Emission);
        material.EnableKeyword("_EMISSION");



    }

    // Approximate blackbody to RGB (normalized)
    Color BlackbodyToRGB(float kelvin)
    {
        kelvin = Mathf.Clamp(kelvin, 1000f, 40000f) / 100f;

        float r, g, b;

        // Red
        r = kelvin <= 66f ? 1.0f : Mathf.Clamp01(1.292936186062745f * Mathf.Pow(kelvin - 60f, -0.1332047592f));

        // Green
        if (kelvin <= 66f)
        {
            g = Mathf.Clamp01(0.3900815787690196f * Mathf.Log(kelvin) - 0.6318414437886275f);
        }
        else
        {
            g = Mathf.Clamp01(1.129890860895294f * Mathf.Pow(kelvin - 60f, -0.0755148492f));
        }

        // Blue
        if (kelvin >= 66f)
        {
            b = 1.0f;
        }
        else if (kelvin <= 19f)
        {
            b = 0.0f;
        }
        else
        {
            b = Mathf.Clamp01(0.5432067891101961f * Mathf.Log(kelvin - 10f) - 1.19625408914f);
        }

        return new Color(r, g, b);
    }
}
