using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class MeshManager : MonoBehaviour
{
    public Slider simplificationSlider; // 控制简化程度
    public Button remeshButton;         // 控制重新网格化
    public Text infoText;

    private SkinnedMeshRenderer skinnedMeshRenderer;
    private Mesh bakedMesh;             // 动态烘焙的网格
    private Mesh originalMesh;          // 原始网格

    void Start()
    {
        skinnedMeshRenderer = GetComponent<SkinnedMeshRenderer>();
        if (skinnedMeshRenderer == null )
        {
            Debug.LogError("SkinnedMeshRenderer not found!");
            return;
        }

        // 保存原始网格
        originalMesh = skinnedMeshRenderer.sharedMesh;
        // 初始化Slider
        simplificationSlider.onValueChanged.AddListener(OnSimplificationSliderValueChanged);
        // 初始化重新网格化按钮
        remeshButton.onClick.AddListener(RemeshMesh);
        // 显示初始信息
        UpdateInfoText();
    }

    void Update()
    {
        
    }

    void OnSimplificationSliderValueChanged(float value)
    {
        // 根据滑动条的值更新网格简化程度
        SimplifyMesh(value);

        // 更新显示信息
        UpdateInfoText();
    }

    void SimplifyMesh(float simplificationFactor)
    {
        if (skinnedMeshRenderer == null || originalMesh == null)
            return;

        // 烘焙当前网格
        bakedMesh = new Mesh();
        skinnedMeshRenderer.BakeMesh(bakedMesh);   // 把当前renderer里面的mesh存入bakedMesh中

        // 获取bakedMesh的顶点和三角形
        Vector3[] verticies = bakedMesh.vertices;
        int[] triangles = bakedMesh.triangles;     // 对顶点的索引组成的int数组，类似OpenGL的索引数组

        // 目标顶点数量
        int targetVertexCount = Mathf.RoundToInt(verticies.Length * (1 - simplificationFactor));
        targetVertexCount = Mathf.Max(targetVertexCount, 3);    // 保证简化后至少有一个三角形面

        // 简化顶点
        Vector3[] simplifiedVerticies = new Vector3[targetVertexCount];
        for (int i = 0; i < targetVertexCount; i++)
        {
            simplifiedVerticies[i] = verticies[i];
        }

        // 更新三角形索引
        int targetTriangleCount = (targetVertexCount / 3) * 3;  // 获取简化后能组成三角形的有效顶点数
        int[] simplifiedTriangles = new int[targetTriangleCount];
        for (int i = 0; i < targetTriangleCount; i++)
        {
            simplifiedTriangles[i] = triangles[i];  // 简单地从前往后取顶点索引，确保一定可以构成三角形
        }

        // 更新网格
        bakedMesh.vertices = simplifiedVerticies;
        bakedMesh.triangles = simplifiedTriangles;
        bakedMesh.RecalculateNormals();     // 从三角形和顶点重新计算网格的法线

        // 应用bakedMesh到skinnedmeshrenderer
        skinnedMeshRenderer.sharedMesh = bakedMesh;
    }

    void RemeshMesh()
    {
        // TODO:网格重新化
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
