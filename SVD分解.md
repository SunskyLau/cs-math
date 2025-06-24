# Tensor

## ✅ 第一部分：线性代数基础

这部分是整个课程的基础，涉及向量空间、线性变换、矩阵运算、子空间等。虽然你可能已经遗忘了一些概念，但我们将从头开始，用简单易懂的方式重新建立这些知识体系。

---

### 1. 向量空间（Vector Space）

#### 📌 定义：
- 一个集合，里面的元素称为**向量**。
- 可以进行两种操作：
  - **加法**：$ v + w \in V $
  - **标量乘法**：$ c \cdot v \in V $, 其中 $ c \in \mathbb{R} $

#### 📌 常见例子：
- $ \mathbb{R}^n $：所有长度为 n 的实数向量
- 比如：$ \begin{bmatrix} 1 \\ 2 \end{bmatrix} \in \mathbb{R}^2 $

#### 📌 向量空间性质：
- 封闭性（加法、数乘后仍在空间内）
- 存在零向量
- 每个向量都有负向量
- 分配律、结合律成立

---

### 2. 线性变换（Linear Map）

#### 📌 定义：
- 函数 $ L: V \rightarrow W $，满足：
  - 加法保持不变：$ L(v + w) = L(v) + L(w) $
  - 标量乘法保持不变：$ L(cv) = cL(v) $

#### 📌 举例：
- 设 $ V = \mathbb{R}^2 $, $ W = \mathbb{R}^2 $, 线性变换 $ L(x) = Ax $，其中：
  $$
  A = \begin{bmatrix}
  1 & 2 \\
  3 & 4
  \end{bmatrix}, \quad x = \begin{bmatrix} a \\ b \end{bmatrix} \Rightarrow L(x) = Ax = \begin{bmatrix} a+2b \\ 3a+4b \end{bmatrix}
  $$

#### 📌 线性变换与矩阵的关系：
- 每个线性变换都可以表示为一个矩阵乘法
- 矩阵记录了基向量的变化

---

### 3. 矩阵的四个基本子空间

#### 📌 给定一个矩阵 $ A \in \mathbb{R}^{m \times n} $，它有以下四个重要子空间：

| 名称 | 数学定义 | 物理意义 |
|------|----------|----------|
| 列空间 $ C(A) $ | $ \{Ax \mid x \in \mathbb{R}^n\} $ | 所有列向量的线性组合 |
| 零空间 $ N(A) $ | $ \{x \mid Ax = 0\} $ | 使矩阵作用为零的所有向量 |
| 行空间 $ R(A) $ | $ C(A^T) $ | 所有行向量的线性组合 |
| 左零空间 $ N(A^T) $ | $ \{x \mid A^T x = 0\} $ | 使转置矩阵作用为零的向量 |

#### 📌 正交关系：
- 列空间 ⊥ 左零空间
- 行空间 ⊥ 零空间

---

### 4. 投影与最小二乘问题（Least Squares）

#### 📌 目标：
- 解不一致方程组 $ Ax = b $，其中没有精确解
- 找到最接近的解，使得误差最小

#### 📌 方法：
- 最小化目标函数：
  $$
  \min_x \|Ax - b\|^2
  $$
- 解为：
  $$
  x = (A^TA)^{-1} A^T b
  $$

#### 📌 示例：
设：
$$
A = \begin{bmatrix}
1 & 2 \\
3 & 4 \\
5 & 6
\end{bmatrix}, \quad
b = \begin{bmatrix}
1 \\
2 \\
3
\end{bmatrix}
$$

计算：
- $ A^T A = \begin{bmatrix} 35 & 44 \\ 44 & 56 \end{bmatrix} $
- $ A^T b = \begin{bmatrix} 38 \\ 49 \end{bmatrix} $
- 解得：
  $$
  x = (A^TA)^{-1} A^T b = \begin{bmatrix} -1.0 \\ 1.0 \end{bmatrix}
  $$

---

## ✅ 第二部分：矩阵分解方法详解（SVD、QR、NMF、交替最小化）

这一部分是重点中的重点，尤其 SVD 和 QR 是张量分解的基础。

---

### 1. 奇异值分解（Singular Value Decomposition, SVD）

#### 📌 形式：
$$
A = U \Sigma V^T
$$
- $ U \in \mathbb{R}^{m \times m} $：正交矩阵（左奇异向量）
- $ \Sigma \in \mathbb{R}^{m \times n} $：对角矩阵，对角线上是奇异值
- $ V \in \mathbb{R}^{n \times n} $：正交矩阵（右奇异向量）

#### 📌 算法步骤：
1. 计算 $ A^T A $，得到右奇异向量矩阵 $ V $
2. 计算 $ AA^T $，得到左奇异向量矩阵 $ U $
3. 奇异值为 $ \sigma_i = \sqrt{\lambda_i} $，即 $ A^T A $ 的特征值开根号
4. 构造 $ \Sigma $，将奇异值按降序排列

#### 📌 示例：
给定：
$$
A = \begin{bmatrix}
1 & 2 \\
3 & 4
\end{bmatrix}
$$

- $ A^T A = \begin{bmatrix} 10 & 14 \\ 14 & 20 \end{bmatrix} $
- 特征值为 $ \lambda_1 = 29.87, \lambda_2 = 0.13 $，奇异值为 $ \sigma_1 = \sqrt{29.87}, \sigma_2 = \sqrt{0.13} $
- 得到：
  $$
  U = \begin{bmatrix}
  -0.4046 & -0.9145 \\
  -0.9145 & 0.4046
  \end{bmatrix}, \quad
  \Sigma = \begin{bmatrix}
  5.465 & 0 \\
  0 & 0.366
  \end{bmatrix}, \quad
  V = \begin{bmatrix}
  -0.5760 & 0.8174 \\
  -0.8174 & -0.5760
  \end{bmatrix}
  $$

---

### 2. QR 分解

#### 📌 形式：
$$
A = QR
$$
- $ Q \in \mathbb{R}^{m \times n} $：正交矩阵（列正交）
- $ R \in \mathbb{R}^{n \times n} $：上三角矩阵

#### 📌 算法步骤（Gram-Schmidt 正交化）：
1. 对矩阵 $ A $ 的每一列做 Gram-Schmidt 正交化
2. 得到一组正交向量 $ q_1, q_2, ..., q_n $
3. 构造矩阵 $ Q = [q_1, q_2, ..., q_n] $
4. $ R $ 的元素为 $ r_{ij} = q_i^T a_j $

#### 📌 示例：
$$
A = \begin{bmatrix}
1 & 1 \\
1 & 0 \\
0 & 1
\end{bmatrix}
\Rightarrow
Q = \begin{bmatrix}
\frac{1}{\sqrt{2}} & \frac{1}{\sqrt{6}} \\
\frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{6}} \\
0 & \frac{2}{\sqrt{6}}
\end{bmatrix}, \quad
R = \begin{bmatrix}
\sqrt{2} & \frac{1}{\sqrt{2}} \\
0 & \frac{\sqrt{6}}{2}
\end{bmatrix}
$$

---

### 3. 非负矩阵分解（Nonnegative Matrix Factorization, NMF）

#### 📌 形式：
$$
M \approx AW, \quad M \in \mathbb{R}_{\geq 0}^{m \times n}, A \in \mathbb{R}_{\geq 0}^{m \times r}, W \in \mathbb{R}_{\geq 0}^{r \times n}
$$

#### 📌 目标函数：
$$
\min_{A,W} \|M - AW\|_F^2
$$

#### 📌 算法步骤（交替最小化 Alternating Minimization）：
1. 初始化 $ A $ 或 $ W $（非负）
2. 固定 $ A $，求最优 $ W $（使用非负最小二乘）
3. 固定 $ W $，求最优 $ A $
4. 重复步骤 2~3，直到收敛或达到最大迭代次数

#### 📌 更新公式（乘法更新规则）：
- $ W_{ik} \leftarrow W_{ik} \cdot \frac{(A^T M)_{ik}}{(A^T A W)_{ik}} $
- $ A_{ik} \leftarrow A_{ik} \cdot \frac{(M W^T)_{ik}}{(A W W^T)_{ik}} $

---

### 4. 交替最小化（Alternating Minimization）

#### 📌 思想：
- 将复杂优化问题拆成两个子问题交替求解
- 每次固定一个变量，优化另一个

#### 📌 应用场景：
- NMF
- 张量分解（PARAFAC）
- 聚类（K-means）

#### 📌 优点：
- 实现简单
- 收敛速度快

#### 📌 缺点：
- 只能保证局部最优
- 对初始化敏感

---
