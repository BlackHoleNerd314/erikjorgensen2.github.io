using UnityEngine;

public static class LorentzMatrix
{
    /// <summary>
    /// Returns the inverse Lorentz boost matrix for a given 3-velocity.
    /// Boosts from the world frame to the rest frame of the moving observer.
    /// </summary>
    public static Matrix4x4 GetInverseBoost(Vector3 v, float c)//30.44f)
    {
        float v2 = v.sqrMagnitude;
        float beta2 = v2 / (c * c);

        if (beta2 >= 1f)
        {
            Debug.LogWarning("Speed must be less than c. Clamping.");
            beta2 = 0.999999f;
        }

        float gamma = 1f / Mathf.Sqrt(1f - beta2);
        Vector3 n = (v2 > 0f) ? v.normalized : Vector3.zero;

        float vx = v.x / c;
        float vy = v.y / c;
        float vz = v.z / c;

        Matrix4x4 L = Matrix4x4.identity;

        // Time-time component
        L[0, 0] = gamma;

        // Time-space components
        L[0, 1] = -gamma * vx;
        L[0, 2] = -gamma * vy;
        L[0, 3] = -gamma * vz;

        // Space-time components
        L[1, 0] = -gamma * vx;
        L[2, 0] = -gamma * vy;
        L[3, 0] = -gamma * vz;


        if (beta2 < 1e-8f)
        {
            return Matrix4x4.identity;
        }

        // Spatial-spatial components
        L[1, 1] = 1f + (gamma - 1f) * (vx * vx) / beta2;
        L[1, 2] = (gamma - 1f) * (vx * vy) / beta2;
        L[1, 3] = (gamma - 1f) * (vx * vz) / beta2;

        L[2, 1] = (gamma - 1f) * (vy * vx) / beta2;
        L[2, 2] = 1f + (gamma - 1f) * (vy * vy) / beta2;
        L[2, 3] = (gamma - 1f) * (vy * vz) / beta2;

        L[3, 1] = (gamma - 1f) * (vz * vx) / beta2;
        L[3, 2] = (gamma - 1f) * (vz * vy) / beta2;
        L[3, 3] = 1f + (gamma - 1f) * (vz * vz) / beta2;

        return L;
    }
}
