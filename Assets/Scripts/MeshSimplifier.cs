using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

/// <summary>
/// 顶点数据
/// </summary>
class VData
{
    public Vector3 pos;              // 顶点坐标
    public List<int> tris = new();   // 该顶点参与的三角形的索引
    public Matrix4x4 Q = Matrix4x4.zero;
}

/// <summary>
/// 每个三角形的数据：指向3个顶点索引——(指向 VData 数组)
/// </summary>
struct TData
{
    public int i0, i1, i2;
    public TData(int a, int b, int c) { i0 = a; i1 = b; i2 = c; }
}

/// <summary>
/// 边数据：两个顶点索引，合并的代价，最优合并点pos
/// </summary>
struct EData
{
    public int vA, vB;          // 顶点索引（保证 vA < vB）
    public float cost;          // 误差代价
    public Vector3 optimalPos;  // 合并后新顶点的位置

    public EData(int a, int b, float c, Vector3 p)
    { vA = a; vB = b; cost = c; optimalPos = p; }
}

[RequireComponent(typeof(MeshFilter))]
public class MeshSimplifier : MonoBehaviour
{
    public Text statsText; 

    private MeshFilter meshFilter;
    private Mesh workingMesh;

    //private Vector3[] vertices;
    //private int[] triangles;

    private List<VData> vDatas;   // 所有顶点
    private List<TData> tDatas;   // 所有三角形
    private List<EData> eDatas;   // 所有边
    private PriorityQueue<int, float> edgePQ;   // 存储边索引和cost

    [Header("曲面简化控制参数")]
    public int targetTriangleCount = 1000;  // 简化后的面数
    public int collapsesPerFrame = 10;      // 每帧坍缩边数量
    public KeyCode startKey = KeyCode.S;    // 开启简化的按键

    private bool hasStart = false;


    private void Awake()
    {
        meshFilter = GetComponent<MeshFilter>();

        // 拿取mesh副本
        workingMesh = Instantiate(meshFilter.sharedMesh);
        meshFilter.sharedMesh = workingMesh;    // 副本放回filter保证不修改原mesh资源

        //vertices = workingMesh.vertices;
        //triangles = workingMesh.triangles;

        BuildVT();
        BuildEdges();
        BuildEdgePriority();

        UpdateStatsUI();
    }

    private void Update()
    {
        if (!hasStart && Input.GetKeyDown(startKey))
        {
            hasStart = true;
            StartCoroutine(SimplifyRoutine());
        }
    }

    private IEnumerator SimplifyRoutine()
    {
        while (tDatas.Count > targetTriangleCount && edgePQ.Count > 0)
        {
            int step = Mathf.Min(collapsesPerFrame, edgePQ.Count);
            for (int i = 0; i < step; i++)
            {
                if (!CollapseCheapestEdge()) 
                    break;
            }

            // 刷新 UI, 让一帧结束
            UpdateStatsUI();
            yield return null;
        }
        Debug.Log($"[MeshSimplifier] 简化完成：顶点 {vDatas.Count}  三角形 {tDatas.Count}");
    }

    /// <summary>
    /// 边坍缩过程
    /// </summary>
    /// <returns></returns>
    private bool CollapseCheapestEdge()
    {
        while (edgePQ.Count > 0)
        {
            // 1. 取出最小代价边索引
            int eIndex = edgePQ.Dequeue();
            if (eIndex >= eDatas.Count)          // 边已失效
                continue;

            var edge = eDatas[eIndex];

            // 如果该边对应的顶点的三角形索引列表为0，说明该顶点是之前坍缩中被删除的
            if (vDatas[edge.vA].tris.Count == 0 || vDatas[edge.vB].tris.Count == 0)
                continue;

            if (edge.cost < 0) edge.cost = 0f;                 // 负零修正

            // 2. 顶点索引
            int ia = edge.vA;
            int ib = edge.vB;
            VData va = vDatas[ia];
            VData vb = vDatas[ib];

            // 3. 生成新顶点
            VData vNew = new VData()
            {
                pos = edge.optimalPos,
                Q = AddMatrix(va.Q, vb.Q)
            };
            int newIndex = vDatas.Count;
            vDatas.Add(vNew);

            //------------------------------------------------------------------
            // 4. 更新所有三角形，把 ia / ib 替换成 newIndex
            //------------------------------------------------------------------
            List<int> affectedTris = new List<int>();           // 收集需要检查的三角形
            affectedTris.AddRange(va.tris);
            affectedTris.AddRange(vb.tris);

            foreach (int tIdx in affectedTris)
            {
                if (tIdx >= tDatas.Count)   // 索引已无效，跳过
                    continue;

                TData td = tDatas[tIdx];

                // 用 newIndex 替换
                if (td.i0 == ia || td.i0 == ib) td.i0 = newIndex;
                if (td.i1 == ia || td.i1 == ib) td.i1 = newIndex;
                if (td.i2 == ia || td.i2 == ib) td.i2 = newIndex;

                // 若出现重复顶点 -> 退化三角形，标记删除
                if (td.i0 == td.i1 || td.i1 == td.i2 || td.i2 == td.i0)
                {
                    //tDatas[tIdx] = tDatas[^1];   // 与最后一个交换
                    //tDatas.RemoveAt(tDatas.Count - 1);

                    RemoveTriangle(tIdx);
                }
                else
                {
                    tDatas[tIdx] = td;           // 保存修改
                }
            }

            //------------------------------------------------------------------
            // 5. 失效旧顶点（不从列表删除，置空其索引的三角形）
            //------------------------------------------------------------------
            vDatas[ia].tris.Clear();
            vDatas[ib].tris.Clear();

            //------------------------------------------------------------------
            // 6. 重新为 newIndex 累加 Q   &  构建它的 tris 列表
            //------------------------------------------------------------------
            vNew.tris = new List<int>();
            for (int t = 0; t < tDatas.Count; t++)
            {
                var td = tDatas[t];
                if (td.i0 == newIndex || td.i1 == newIndex || td.i2 == newIndex)
                    vNew.tris.Add(t);
            }
            vDatas[newIndex] = vNew;   // 写回

            //------------------------------------------------------------------
            // 7. 为受影响的边重新计算 cost 并入堆
            //------------------------------------------------------------------
            RebuildEdgesAround(newIndex);

            return true;
        }
        // 优先队列空或无合法边
        return false;
    }

    /// <summary>
    /// 删除 tDatas[tIdx]，并保持所有顶点的 tris 列表一致
    /// </summary>
    /// <param name="tIdx"></param>
    private void RemoveTriangle(int tIdx)
    {
        int lastIdx = tDatas.Count - 1;
        var dead = tDatas[tIdx];   // 要删除的三角形
        var lastTri = tDatas[lastIdx];

        // 1) 先从 dead 的三个顶点列表里移除 tIdx
        vDatas[dead.i0].tris.Remove(tIdx);
        vDatas[dead.i1].tris.Remove(tIdx);
        vDatas[dead.i2].tris.Remove(tIdx);

        // 2) 如果 dead 不是最后一个，拿 lastTri 填补空位
        if (tIdx != lastIdx)
        {
            tDatas[tIdx] = lastTri;   // 覆盖

            // 把 lastTri 在三个顶点的 tris 列表里把 lastIdx 改成 tIdx
            ReplaceTriIndex(vDatas[lastTri.i0].tris, lastIdx, tIdx);
            ReplaceTriIndex(vDatas[lastTri.i1].tris, lastIdx, tIdx);
            ReplaceTriIndex(vDatas[lastTri.i2].tris, lastIdx, tIdx);
        }

        // 3) 真正删掉最后一个元素
        tDatas.RemoveAt(lastIdx);
    }

    /// <summary>
    /// 帮助函数：把 list 中的 oldIdx 替换为 newIdx
    /// </summary>
    /// <param name="list"></param>
    /// <param name="oldIdx"></param>
    /// <param name="newIdx"></param>
    private static void ReplaceTriIndex(List<int> list, int oldIdx, int newIdx)
    {
        int p = list.IndexOf(oldIdx);
        if (p >= 0) list[p] = newIdx;
    }

    /// <summary>
    /// 扫描新顶点所关联的所有三角形，把相关边重新计算cost并加入优先队列
    /// </summary>
    /// <param name="newIndex"></param>
    private void RebuildEdgesAround(int newIndex)
    {
        VData vn = vDatas[newIndex];

        // 遍历新点所在所有三角形
        foreach (var tIdx in vn.tris)
        {
            TData td = tDatas[tIdx];
            int[] vertexIdx = { td.i0, td.i1, td.i2 };
            for (int k = 0; k < 3; k++)
            {
                int a = vertexIdx[k];
                int b = vertexIdx[(k + 1) % 3];
                if (a == b)
                    continue;
                if (a > b)
                    (a, b) = (b, a);

                // 计算该边的最优合并点及其cost
                var (pos, c) = ComputeEdgeCost(vDatas[a], vDatas[b]);

                int edgeIdx = eDatas.Count;
                eDatas.Add(new EData(a, b, c, pos));
                edgePQ.Enqueue(edgeIdx, c);
            }
        }
    }

    /// <summary>
    /// 把workingMesh的数据转换为VData和TData数据
    /// 计算每个顶点的误差矩阵Q
    /// </summary>
    private void BuildVT()
    {
        Vector3[] verts = workingMesh.vertices;
        int[] tris = workingMesh.triangles;

        int vCount = verts.Length;
        int tCount = tris.Length / 3;

        // 初始化VData
        vDatas = new List<VData>(vCount);
        for (int i = 0; i < vCount; i++)
        {
            vDatas.Add(new VData { pos = verts[i] });
        }

        // 初始化TData，把三角形索引给到顶点
        tDatas = new List<TData>(tCount);
        for (int t = 0; t < tCount; t++)
        {
            int i0 = tris[t * 3 + 0];
            int i1 = tris[t * 3 + 1];
            int i2 = tris[t * 3 + 2];

            tDatas.Add(new TData(i0, i1, i2));
            vDatas[i0].tris.Add(t);
            vDatas[i1].tris.Add(t);
            vDatas[i2].tris.Add(t);
        }

        // 为每个三角形求平面方程，累加到三个顶点的Q矩阵
        for (int t = 0; t < tCount; t++)
        {
            var td = tDatas[t];
            Vector3 p0 = vDatas[td.i0].pos;
            Vector3 p1 = vDatas[td.i1].pos;
            Vector3 p2 = vDatas[td.i2].pos;

            // 计算平面方程(a,b,c,d)
            Vector3 n = Vector3.Cross(p1 - p0, p2 - p0).normalized;
            float d = -Vector3.Dot(n, p0);
            Vector4 plane = new(n.x, n.y, n.z, d);

            // 构造4x4误差矩阵 K = plane * plane^T
            Matrix4x4 K = OuterProduct(plane);

            // 累加K到三个顶点中
            vDatas[td.i0].Q = AddMatrix(vDatas[td.i0].Q, K);
            vDatas[td.i1].Q = AddMatrix(vDatas[td.i1].Q, K);
            vDatas[td.i2].Q = AddMatrix(vDatas[td.i2].Q, K);
        }

        Debug.Log($"[MeshSimplifier] BuildVT finish：vertex {vCount}，triangle {tCount}");
    }
    
    /// <summary>
    /// 扫描tData，生成无重复的边并计算cost
    /// </summary>
    private void BuildEdges()
    {
        // 确保没有重复边
        var edgeSet = new HashSet<(int, int)>();
        eDatas = new List<EData>();

        for (int t = 0; t < tDatas.Count; t++)
        {
            var td = tDatas[t];
            AddEdgeUnique(td.i0, td.i1, edgeSet);
            AddEdgeUnique(td.i1, td.i2, edgeSet);
            AddEdgeUnique(td.i2, td.i0, edgeSet);
        }

        Debug.Log($"[MeshSimplifier] edge build：{eDatas.Count}");
    }

    /// <summary>
    /// 添加边
    /// </summary>
    /// <param name="i2"></param>
    /// <param name="i0"></param>
    /// <param name="edgeSet"></param>
    private void AddEdgeUnique(int a, int b, HashSet<(int, int)> set)
    {
        if (a == b)
            return;
        if (a > b)
            (a, b) = (b, a);    // 保证a < b
        if (!set.Add((a, b)))   // 已存在就return
            return;

        // 计算cost
        var v1 = vDatas[a];
        var v2 = vDatas[b];
        (Vector3 pos, float cost) = ComputeEdgeCost(v1, v2);

        // 添加到边列表
        eDatas.Add(new EData(a, b, cost, pos));
    }

    /// <summary>
    /// 初始化边-cost优先队列
    /// </summary>
    private void BuildEdgePriority()
    {
        edgePQ = new PriorityQueue<int, float>();
        for (int i = 0; i < eDatas.Count; i++)
        {
            edgePQ.Enqueue(i, eDatas[i].cost);
        }
    }

    /// <summary>
    /// 根据 Q1+Q2 计算最优位置和误差代价
    /// </summary>
    /// <param name="v1"></param>
    /// <param name="v2"></param>
    /// <returns></returns>
    private (Vector3 pos, float cost) ComputeEdgeCost(VData v1, VData v2)
    {
        Matrix4x4 Q = AddMatrix(v1.Q, v2.Q);

        // 拆3x3的A 和 3x1的b
        Matrix4x4 A = new Matrix4x4();
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                A[r, c] = Q[r, c];
        Vector3 b = new Vector3(Q[0, 3], Q[1, 3], Q[2, 3]);

        // 解线性方程Ax = -b，若不可逆则直接选取中点作为合并点
        bool invertible = Invert3x3(A, out Matrix4x4 Ainv);
        Vector3 pos;
        if (invertible)
            pos = Ainv.MultiplyVector(-b);   // x = -A^-1 · b
        else
            pos = 0.5f * (v1.pos + v2.pos);  // 不可逆 -- 用中点

        // 代价 = v_h^T Q v_h  合并点对于误差矩阵计算出来的二次度量误差
        Vector4 vh = new Vector4(pos.x, pos.y, pos.z, 1f);
        float cost = Vector4.Dot(vh, Q * vh);

        return (pos, cost);
    }

    /// <summary>
    /// 求3x3矩阵逆，不可逆返回false
    /// </summary>
    /// <param name="m"></param>
    /// <param name="inv"></param>
    /// <returns></returns>
    private static bool Invert3x3(Matrix4x4 m, out Matrix4x4 inv)
    {
        inv = Matrix4x4.zero;
        float det =
            m.m00 * (m.m11 * m.m22 - m.m12 * m.m21) -
            m.m01 * (m.m10 * m.m22 - m.m12 * m.m20) +
            m.m02 * (m.m10 * m.m21 - m.m11 * m.m20);

        if (Mathf.Abs(det) < 1e-12f) return false;

        float id = 1f / det;

        inv.m00 = (m.m11 * m.m22 - m.m12 * m.m21) * id;
        inv.m01 = -(m.m01 * m.m22 - m.m02 * m.m21) * id;
        inv.m02 = (m.m01 * m.m12 - m.m02 * m.m11) * id;

        inv.m10 = -(m.m10 * m.m22 - m.m12 * m.m20) * id;
        inv.m11 = (m.m00 * m.m22 - m.m02 * m.m20) * id;
        inv.m12 = -(m.m00 * m.m12 - m.m02 * m.m10) * id;

        inv.m20 = (m.m10 * m.m21 - m.m11 * m.m20) * id;
        inv.m21 = -(m.m00 * m.m21 - m.m01 * m.m20) * id;
        inv.m22 = (m.m00 * m.m11 - m.m01 * m.m10) * id;

        return true;
    }

    /// <summary>
    /// 返回p * p^T 的矩阵
    /// </summary>
    /// <param name="p"></param>
    /// <returns></returns>
    private static Matrix4x4 OuterProduct(Vector4 p)
    {
        Matrix4x4 m = new Matrix4x4();
        for (int r = 0; r < 4; r++)
        {
            for (int c = 0; c < 4; c++)
            {
                m[r, c] = p[r] * p[c];
            }
        }
        return m;
    }

    public static Matrix4x4 AddMatrix(Matrix4x4 a, Matrix4x4 b)
    {
        Matrix4x4 res = new Matrix4x4();
        for (int r = 0;r < 4;r++)
        {
            for (int c = 0;c < 4;c++)
            {
                res[r, c] = a[r, c] + b[r, c];
            }
        }
        return res;
    }

    /// <summary>
    /// 更新UI显示
    /// </summary>
    private void UpdateStatsUI()
    {
        if (statsText != null)
        {
            statsText.text = $"Vertex: {vDatas.Count}\nTriangles: {tDatas.Count}";
        }

        // 把简化后的几何写回 Mesh（只位置/索引，法线可后算）
        Vector3[] newVerts = new Vector3[vDatas.Count];
        for (int i = 0; i < vDatas.Count; i++) 
            newVerts[i] = vDatas[i].pos;
        int[] newTris = new int[tDatas.Count * 3];
        for (int t = 0; t < tDatas.Count; t++)
        {
            var td = tDatas[t];
            newTris[t * 3 + 0] = td.i0;
            newTris[t * 3 + 1] = td.i1;
            newTris[t * 3 + 2] = td.i2;
        }

        workingMesh.Clear();
        workingMesh.vertices = newVerts;
        workingMesh.triangles = newTris;
        workingMesh.RecalculateNormals();
    }
}