using System;
using System.Collections;
using System.Collections.Generic;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UI;

// 抄一个优先队列来用
public class PriorityQueue<TElement, TPriority> : IEnumerable<(TElement element, TPriority priority)>
    where TPriority : IComparable<TPriority>
{
    private List<(TElement element, TPriority priority)> heap = new List<(TElement, TPriority)>();

    public int Count => heap.Count;

    // 清空元素
    public void Clear()
    {
        heap.Clear(); // 直接清空内部列表
    }

    // 插入元素
    public void Enqueue(TElement element, TPriority priority)
    {
        heap.Add((element, priority));
        HeapifyUp(heap.Count - 1);
    }

    // 移除并返回优先级最高的元素
    public TElement Dequeue()
    {
        if (heap.Count == 0) throw new InvalidOperationException("The priority queue is empty.");
        TElement result = heap[0].element;

        // 用最后一个元素替代根元素
        heap[0] = heap[heap.Count - 1];
        heap.RemoveAt(heap.Count - 1);

        if (heap.Count > 0)
            HeapifyDown(0);

        return result;
    }

    // 获取优先级最高的元素但不移除
    public TElement Peek()
    {
        if (heap.Count == 0) throw new InvalidOperationException("The priority queue is empty.");
        return heap[0].element;
    }

    // 移除特定元素
    public bool Remove(TElement element)
    {
        int index = heap.FindIndex(x => EqualityComparer<TElement>.Default.Equals(x.element, element));
        if (index == -1) return false;

        heap[index] = heap[heap.Count - 1];
        heap.RemoveAt(heap.Count - 1);

        if (index < heap.Count)
        {
            HeapifyUp(index);
            HeapifyDown(index);
        }

        return true;
    }

    // 实现 IEnumerable 接口，支持 foreach 遍历
    public IEnumerator<(TElement element, TPriority priority)> GetEnumerator()
    {
        return heap.GetEnumerator();
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    // 堆化：从下往上
    private void HeapifyUp(int index)
    {
        while (index > 0)
        {
            int parentIndex = (index - 1) / 2;
            if (heap[index].priority.CompareTo(heap[parentIndex].priority) >= 0)
                break;

            Swap(index, parentIndex);
            index = parentIndex;
        }
    }

    // 堆化：从上往下
    private void HeapifyDown(int index)
    {
        while (true)
        {
            int leftChild = 2 * index + 1;
            int rightChild = 2 * index + 2;
            int smallest = index;

            if (leftChild < heap.Count && heap[leftChild].priority.CompareTo(heap[smallest].priority) < 0)
                smallest = leftChild;

            if (rightChild < heap.Count && heap[rightChild].priority.CompareTo(heap[smallest].priority) < 0)
                smallest = rightChild;

            if (smallest == index)
                break;

            Swap(index, smallest);
            index = smallest;
        }
    }

    // 交换堆中的两个元素
    private void Swap(int i, int j)
    {
        var temp = heap[i];
        heap[i] = heap[j];
        heap[j] = temp;
    }
}


// 边的结构体
public struct Edge : IComparable<Edge>
{
    public int v1;  // 边的第一个顶点索引
    public int v2;  // 边的第二个顶点索引

    public Edge(int v1, int v2)
    {
        this.v1 = v1;
        this.v2 = v2;
    }

    // 实现比较器，用于优先队列按代价排序
    public int CompareTo(Edge other)
    {
        // 根据顶点索引比较
        if (v1 != other.v1)
            return v1.CompareTo(other.v1);
        else
            return v2.CompareTo(other.v2);
    }

    public override bool Equals(object obj)
    {
        if (obj is Edge other)  // 转换obj为Edge类型，如果可以
        {
            // 检查两个顶点是否相同（无序）
            return (v1 == other.v1 && v2 == other.v2) || (v1 == other.v2 && v2 == other.v1);
        }
        return false;
    }

    public override int GetHashCode()
    {
        // 生成无序哈希码，确保 {v1, v2} 和 {v2, v1} 被视为同一条边
        return v1.GetHashCode() ^ v2.GetHashCode();
    }
}

// 误差Q结构体
public struct Quadric
{
    public Matrix4x4 Q;     // 误差矩阵

    public Quadric(Matrix4x4 q)
    {
        this.Q = q;
    }

    // 计算某个顶点的位置误差
    public float ComputeError(Vector3 vertex)
    {
        Vector4 v = new Vector4(vertex.x, vertex.y, vertex.z, 1); // 转为齐次坐标
        return Vector4.Dot(v, Q * v);   // 计算误差————距离平方转为矩阵形式得到
    }

    // 定义Matrix4x4加法
    public static Matrix4x4 AddMatrix4x4(Matrix4x4 m1, Matrix4x4 m2)
    {
        Matrix4x4 res = new Matrix4x4();
        for (int row = 0; row < 4; row++)
        {
            for (int col = 0; col < 4; col++)
            {
                res[row, col] = m1[row, col] + m2[row, col];
            }
        }
        return res;
    }

    // 合并两个Quadric
    public static Quadric operator +(Quadric q1, Quadric q2)
    {
        return new Quadric(AddMatrix4x4(q1.Q, q2.Q));
    }
}

public class QEMSimplifier : MonoBehaviour
{
    public Slider simplificationSlider;     // 滑动条，用于控制简化程度
    public Text infoText;                   // 显示顶点数和三角形数的文本
    private Mesh originalMesh;              // 原始网格
    private Mesh workingMesh;               // 当前正在简化的网格
    private Dictionary<int, Quadric> quadrics;  // 存储顶点的误差矩阵
    private PriorityQueue<Edge, float> edgeQueue; // 优先队列，存储边和其代价

    // Start is called before the first frame update
    void Start()
    {
        // 获取模型原始网格
        originalMesh = GetComponent<MeshFilter>().mesh;
        workingMesh = Instantiate(originalMesh); // 创建副本用于操作
        quadrics = new Dictionary<int, Quadric>();
        edgeQueue = new PriorityQueue<Edge, float>();

        // 初始化顶点误差矩阵和边队列
        InitializeQuadrics(workingMesh);
        InitializeEdgeQueue(workingMesh);

        // 初始化滑动条事件监听
        simplificationSlider.onValueChanged.AddListener(OnSliderValueChanged);

        // 更新显示信息
        UpdateInfoText();
    }

    void OnSliderValueChanged(float value)
    {
        // 根据滑动条值计算目标顶点数
        int targetVertexCount = Mathf.RoundToInt(originalMesh.vertexCount * (1 - value));

        if (targetVertexCount >= originalMesh.vertexCount)
        {
            // 恢复原始网格
            GetComponent<MeshFilter>().mesh = originalMesh;
        }
        else
        {
            // 执行网格简化
            SimplifyMesh(targetVertexCount);
        }

        // 更新显示信息
        UpdateInfoText();
    }

    // 边坍缩网格简化的总函数
    void SimplifyMesh(int targetVertexCount)
    {
        // 创建新的网格
        workingMesh = Instantiate(originalMesh);

        // 清空并且再次初始化优先队列,初始化误差矩阵
        edgeQueue.Clear();
        InitializeEdgeQueue(workingMesh);
        InitializeQuadrics(workingMesh);

        // 获取顶点和三角形
        Vector3[] vertices = workingMesh.vertices;
        int[] triangles = workingMesh.triangles;

        // 持续网格简化直到顶点数量不足目标简化数量————不断坍缩边以更新顶点位置和三角形
        while (vertices.Length > targetVertexCount && edgeQueue.Count > 0)
        {
            // 获取代价最小的边
            Edge edge = edgeQueue.Dequeue();
            // 计算边坍缩后新点的位置
            Vector3 optimalPoint;
            float cost = ComputeEdgeCost(edge, vertices, out optimalPoint);

            // 更新坍缩边的顶点位置
            vertices[edge.v1] = optimalPoint;

            // 更新三角形的顶点索引
            UpdateTriangles(triangles, edge.v1, edge.v2);

            // 更新误差矩阵——这里edge.v1存储的误差本来还是坍缩前v1的误差
            quadrics[edge.v1] = quadrics[edge.v1] + quadrics[edge.v2];

            // 更新边的优先队列
            UpdateEdgeQueue(edge.v1, vertices);
        }

        // 重新构建Mesh的顶点和索引
        RebuildMesh(vertices, triangles);

        // 应用模型
        GetComponent<MeshFilter>().mesh = workingMesh;
    }

    // 重新构建Mesh的顶点和索引
    void RebuildMesh(Vector3[] vertices, int[] triangles)
    {
        // 创建旧顶点到新顶点索引的映射
        Dictionary<int, int> indexMap = new Dictionary<int, int>();
        List<Vector3> newVertices = new List<Vector3>();
        List<int> newTriangles = new List<int>();

        // 遍历三角形，收集被使用的顶点
        for (int i = 0; i < triangles.Length; i += 3)
        {
            int v0 = triangles[i];
            int v1 = triangles[i + 1];
            int v2 = triangles[i + 2];

            // 跳过退化三角形（所有顶点索引相同或为无效索引）
            if (v0 == v1 || v1 == v2 || v0 == v2 || v0 < 0 || v1 < 0 || v2 < 0)
                continue;

            // 确保每个顶点都映射到新顶点数组
            if (!indexMap.ContainsKey(v0))
            {
                indexMap[v0] = newVertices.Count;
                newVertices.Add(vertices[v0]);
            }
            if (!indexMap.ContainsKey(v1))
            {
                indexMap[v1] = newVertices.Count;
                newVertices.Add(vertices[v1]);
            }
            if (!indexMap.ContainsKey(v2))
            {
                indexMap[v2] = newVertices.Count;
                newVertices.Add(vertices[v2]);
            }

            // 添加新三角形的索引
            newTriangles.Add(indexMap[v0]);
            newTriangles.Add(indexMap[v1]);
            newTriangles.Add(indexMap[v2]);
        }

        // 检查是否有未更新的索引（冗余数据）
        foreach (var triangle in newTriangles)
        {
            if (triangle < 0 || triangle >= newVertices.Count)
            {
                Debug.LogError($"Invalid triangle index detected: {triangle}");
                return;
            }
        }

        // 更新网格数据
        workingMesh.Clear(); // 清空旧的网格数据，防止残留数据干扰
        workingMesh.vertices = newVertices.ToArray();
        workingMesh.triangles = newTriangles.ToArray();
        workingMesh.RecalculateNormals(); // 重新计算法线
    }

    // 边坍缩后更新边的优先队列
    void UpdateEdgeQueue(int updatedVertex, Vector3[] vertices)
    {
        List<Edge> affectedEdges = new List<Edge>();

        // 查找受影响的边（与 updatedVertex 相连的所有边）
        foreach (var edge in edgeQueue)
        {
            if (edge.element.v1 == updatedVertex || edge.element.v2 == updatedVertex)
            {
                affectedEdges.Add(edge.element);
            }
        }

        // 从优先队列中移除受影响的边
        foreach (var edge in affectedEdges)
        {
            edgeQueue.Remove(edge);
        }

        // 重新计算受影响边的代价并加入队列
        foreach (var edge in affectedEdges)
        {
            Vector3 optimalPoint;
            float cost = ComputeEdgeCost(edge, vertices, out optimalPoint);
            edgeQueue.Enqueue(edge, cost);
        }
    }

    // 初始化每个顶点的误差矩阵Q
    void InitializeQuadrics(Mesh mesh)
    {
        Vector3[] vertices = mesh.vertices;
        int[] triangles = mesh.triangles;

        for (int i = 0; i < triangles.Length; i += 3)
        {
            // 获取三角形中存储的顶点的索引
            int v0 = triangles[i];
            int v1 = triangles[i + 1];
            int v2 = triangles[i + 2];

            // 根据索引在顶点数组取出顶点（位置坐标）
            Vector3 p0 = vertices[v0];
            Vector3 p1 = vertices[v1];
            Vector3 p2 = vertices[v2];

            // 计算面法线
            Vector3 normal = Vector3.Cross(p1 - p0, p2 - p0).normalized;

            // d是平面的偏移量，计算公式为d = −n dot p，是通过点法式的平面方程化简系数得到的
            float d = -Vector3.Dot(normal, p0);

            // 计算平面误差矩阵
            Matrix4x4 planeQ = ComputePlaneQuadric(normal, d);

            // 为每个顶点累加该平面的误差矩阵
            if (!quadrics.ContainsKey(v0)) quadrics[v0] = new Quadric(planeQ);
            else quadrics[v0] = quadrics[v0] + new Quadric(planeQ);

            if (!quadrics.ContainsKey(v1)) quadrics[v1] = new Quadric(planeQ);
            else quadrics[v1] = quadrics[v1] + new Quadric(planeQ);

            if (!quadrics.ContainsKey(v2)) quadrics[v2] = new Quadric(planeQ);
            else quadrics[v2] = quadrics[v2] + new Quadric(planeQ);
        }
    }

    // 初始化边代价为权的优先队列
    void InitializeEdgeQueue(Mesh mesh)
    {
        HashSet<Edge> processEdges = new HashSet<Edge>();
        Vector3[] vertices = mesh.vertices;
        int[] triangles = mesh.triangles;

        for (int i = 0; i < triangles.Length; i += 3)
        {
            // 按顺序获取每个三角形中存储的顶点的索引
            int v0 = triangles[i];
            int v1 = triangles[i + 1];
            int v2 = triangles[i + 2];

            // 按顺序构建边
            Edge[] edges = { new Edge(v0,v1), new Edge(v1,v2), new Edge(v2,v0) };

            foreach (var edge in edges)
            {
                // 如果没有处理过该边
                if (!processEdges.Contains(edge))
                {
                    processEdges.Add(edge);
                    float cost = ComputeEdgeCost(edge, vertices, out _);
                    edgeQueue.Enqueue(edge, cost);
                }
            }
        }
    }

    // 计算坍缩一条指定边的代价，out出坍缩后顶点的位置
    float ComputeEdgeCost(Edge edge, Vector3[] vertices, out Vector3 optimalPoint)
    {
        Quadric q1 = quadrics[edge.v1];
        Quadric q2 = quadrics[edge.v2];

        // 合并误差矩阵
        Matrix4x4 Q = Quadric.AddMatrix4x4(q1.Q, q2.Q);

        // 提取算法论文中的Q矩阵系数A和b
        Matrix4x4 A = ExtractTopLeft3x3(Q);
        Vector3 b = ExtractColumn(Q, 3);

        // 如果A可逆，则可以数值计算出边坍缩后最优点的位置
        if (IsInvertible(A))
        {
            optimalPoint = -(A.inverse * b);
            // 计算坍缩边的cost代价————近似为新顶点与原顶点的误差矩阵的误差计算加和得到
            return q1.ComputeError(optimalPoint) + q2.ComputeError(optimalPoint);
        }
        else
        {
            // 如果A不可逆，那么近似就取该边的中点位置作为合并后的新顶点
            optimalPoint = (vertices[edge.v1] + vertices[edge.v2]) / 2;
            return q1.ComputeError(vertices[edge.v1]) + q2.ComputeError(vertices[edge.v2]);
        }
    }

    // 计算一个三角形面的误差的系数矩阵————用于给其上三个顶点计算误差使用
    Matrix4x4 ComputePlaneQuadric(Vector3 normal, float d)
    {
        return new Matrix4x4(
            new Vector4(normal.x * normal.x, normal.x * normal.y, normal.x * normal.z, normal.x * d),
            new Vector4(normal.x * normal.y, normal.y * normal.y, normal.y * normal.z, normal.y * d),
            new Vector4(normal.x * normal.z, normal.y * normal.z, normal.z * normal.z, normal.z * d),
            new Vector4(normal.x * d, normal.y * d, normal.z * d, d * d)
        );
    }

    // 提取误差矩阵Q中的左上角A 系数矩阵
    Matrix4x4 ExtractTopLeft3x3(Matrix4x4 Q)
    {
        return new Matrix4x4(
            Q.GetRow(0),
            Q.GetRow(1),
            Q.GetRow(2),
            new Vector4(0, 0, 0, 1)
        );
    }

    // 提取误差矩阵Q中的第四列b 系数向量，这里column应该都是取3获得第四列
    Vector3 ExtractColumn(Matrix4x4 Q, int column)
    {
        return new Vector3(Q[0, column], Q[1, column], Q[2, column]);
    }

    // 4x4矩阵是否可逆
    bool IsInvertible(Matrix4x4 A)
    {
        // 矩阵可逆 = 满秩 = 行列式不为0
        return Mathf.Abs(A.determinant) > Mathf.Epsilon;
    }

    // 边坍缩操作后更新Mesh的三角形顶点索引
    void UpdateTriangles(int[] triangles, int v1, int v2)
    {
        for (int i = 0; i < triangles.Length; i++)
        {
            // 如果当前顶点索引等于v2（被移除的顶点索引），将其替换为v1（保留的顶点索引）
            if (triangles[i] == v2)
                triangles[i] = v1;
        }
    }

    // 更新UI的text信息显示
    void UpdateInfoText()
    {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        infoText.text = $"Vertices: {mesh.vertexCount}\nTriangles: {mesh.triangles.Length / 3}";
    }
}
