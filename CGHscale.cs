using UnityEngine;

public class CGHscale : MonoBehaviour
{
    static float c0 = 299792458f;
    static float G0 = 6.67408e-11f;
    static float h0 = 6.62607015e-34f;
    static float kB0 = 1.380649e-23f;
    public static float ToyScale = 0.175f;
    public static float c = Mathf.Pow(c0, ToyScale);
    public static float G = Mathf.Pow(G0, ToyScale);
    public static float h = Mathf.Pow(h0, ToyScale);
    public static float kB = Mathf.Pow(kB0, ToyScale);
}
