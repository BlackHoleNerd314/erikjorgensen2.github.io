using UnityEngine;

public class CameraRotate : MonoBehaviour
{

    public Transform Player;
    public Transform X;
    public Transform Y; 
    public Transform Z; 
    private Vector3 torq;
    private Vector3 AngMoment;
    void Start()
    {
        torq = new Vector3(0, 0, 0);
        AngMoment = new Vector3(0, 0, 0);
    }
    // Update is called once per frame
    void Update()
    {
        Vector3 X0 = new Vector3(1, 0, 0);//(X.position - transform.position).normalized;
        Vector3 Y0 = new Vector3(0, 1, 0);//(Y.position - transform.position).normalized;
        Vector3 Z0 = new Vector3(0, 0, 1);//(Z.position - transform.position).normalized;
        if (Input.GetKey(KeyCode.Q))
        {
            torq += Z0;
        }
        if (Input.GetKey(KeyCode.E))
        {
            torq += -Z0;
        }

        if (Input.GetKey(KeyCode.Z))
        {
            torq += Y0;
        }
        if (Input.GetKey(KeyCode.X))
        {
            torq += -Y0;
        }

        if (Input.GetKey(KeyCode.C))
        {
            torq += X0;
        }
        if (Input.GetKey(KeyCode.V))
        {
            torq += -X0;
        }

        Vector3 AngMoment = torq/1000f;
        float angleInDegrees = AngMoment.magnitude;
        //Player.rotation *= Quaternion.Euler(AngMoment.x,AngMoment.y,AngMoment.z);
        Player.rotation *= Quaternion.AngleAxis(angleInDegrees, AngMoment.normalized);
    }
}
