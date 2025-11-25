using UnityEngine;

public class PlayerScale : MonoBehaviour
{
    public Transform playerTR;
    public Rigidbody playerRB;
    // Update is called once per frame
    void Update()
    {
        if (Input.GetKeyDown(KeyCode.T))
        {
            Time.timeScale *= 2f;
            Time.fixedDeltaTime = 0.02f * Time.timeScale;
        }
        if (Input.GetKeyDown(KeyCode.B))
        {
            Time.timeScale *= 0.5f;
            Time.fixedDeltaTime = 0.02f * Time.timeScale;
        }
        if (Input.GetKeyDown(KeyCode.Y))
        {
            playerTR.localScale *= 2f;
        }
        if (Input.GetKeyDown(KeyCode.N))
        {
            playerTR.localScale *= 0.5f;
        }
        if (Input.GetKeyDown(KeyCode.U))
        {
            playerRB.mass *= 2f;
        }
        if (Input.GetKeyDown(KeyCode.M))
        {
            playerRB.mass *= 0.5f;
        }

        if (Input.GetKeyDown(KeyCode.G))
        {
            Time.timeScale = 1f;
            Time.fixedDeltaTime = 0.02f;
        }
        if (Input.GetKeyDown(KeyCode.H))
        {
            playerTR.localScale = new Vector3(1f, 1f, 1f);
        }
        if (Input.GetKeyDown(KeyCode.J))
        {
            playerRB.mass = 1f;
        }

    }
}
