using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class MeshManager : MonoBehaviour
{
    public Slider simplificationSlider; // ���Ƽ򻯳̶�
    public Button remeshButton;         // ������������
    public Text infoText;

    private SkinnedMeshRenderer skinnedMeshRenderer;
    private Mesh bakedMesh;             // ��̬�決������
    private Mesh originalMesh;          // ԭʼ����

    void Start()
    {
        skinnedMeshRenderer = GetComponent<SkinnedMeshRenderer>();
        if (skinnedMeshRenderer == null )
        {
            Debug.LogError("SkinnedMeshRenderer not found!");
            return;
        }

        // ����ԭʼ����
        originalMesh = skinnedMeshRenderer.sharedMesh;
        // ��ʼ��Slider
        simplificationSlider.onValueChanged.AddListener(OnSimplificationSliderValueChanged);
        // ��ʼ���������񻯰�ť
        remeshButton.onClick.AddListener(RemeshMesh);
        // ��ʾ��ʼ��Ϣ
        UpdateInfoText();
    }

    void Update()
    {
        
    }

    void OnSimplificationSliderValueChanged(float value)
    {
        // ���ݻ�������ֵ��������򻯳̶�
        SimplifyMesh(value);

        // ������ʾ��Ϣ
        UpdateInfoText();
    }

    void SimplifyMesh(float simplificationFactor)
    {
        if (skinnedMeshRenderer == null || originalMesh == null)
            return;

        // �決��ǰ����
        bakedMesh = new Mesh();
        skinnedMeshRenderer.BakeMesh(bakedMesh);   // �ѵ�ǰrenderer�����mesh����bakedMesh��

        // ��ȡbakedMesh�Ķ����������
        Vector3[] verticies = bakedMesh.vertices;
        int[] triangles = bakedMesh.triangles;     // �Զ����������ɵ�int���飬����OpenGL����������

        // Ŀ�궥������
        int targetVertexCount = Mathf.RoundToInt(verticies.Length * (1 - simplificationFactor));
        targetVertexCount = Mathf.Max(targetVertexCount, 3);    // ��֤�򻯺�������һ����������

        // �򻯶���
        Vector3[] simplifiedVerticies = new Vector3[targetVertexCount];
        for (int i = 0; i < targetVertexCount; i++)
        {
            simplifiedVerticies[i] = verticies[i];
        }

        // ��������������
        int targetTriangleCount = (targetVertexCount / 3) * 3;  // ��ȡ�򻯺�����������ε���Ч������
        int[] simplifiedTriangles = new int[targetTriangleCount];
        for (int i = 0; i < targetTriangleCount; i++)
        {
            simplifiedTriangles[i] = triangles[i];  // �򵥵ش�ǰ����ȡ����������ȷ��һ�����Թ���������
        }

        // ��������
        bakedMesh.vertices = simplifiedVerticies;
        bakedMesh.triangles = simplifiedTriangles;
        bakedMesh.RecalculateNormals();     // �������κͶ������¼�������ķ���

        // Ӧ��bakedMesh��skinnedmeshrenderer
        skinnedMeshRenderer.sharedMesh = bakedMesh;
    }

    void RemeshMesh()
    {
        // TODO:�������»�
    }

    void UpdateInfoText()
    {
        if (infoText != null && bakedMesh != null)
        {
            int vertexCount = bakedMesh.vertexCount;
            int triangleCount = bakedMesh.triangles.Length / 3;
            infoText.text = $"Verticies: {vertexCount}\nTriangles: {triangleCount}";
        }
    }
}
