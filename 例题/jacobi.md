我们来求一个 **3×3 矩阵的 Jacobi 迭代矩阵** 和 **Gauss-Seidel (GS) 迭代矩阵**。

---

## 一、背景知识

给定线性方程组：

$$
A\mathbf{x} = \mathbf{b}
$$

其中 $ A \in \mathbb{R}^{n \times n} $，我们可以将矩阵 $ A $ 分解为：

$$
A = D - L - U
$$

- $ D $：对角矩阵（主对角线上元素）
- $ -L $：严格下三角部分（不含对角线）
- $ -U $：严格上三角部分（不含对角线）

### Jacobi 迭代法迭代矩阵：
$$
J = D^{-1}(L + U)
$$

### Gauss-Seidel 迭代法迭代矩阵：
$$
G = (D - L)^{-1} U
$$

---

## 二、具体计算示例

我们选取一个具体的 3×3 矩阵 $ A $ 来演示：

设：

$$
A = 
\begin{bmatrix}
4 & -1 & 0 \\
-1 & 4 & -1 \\
0 & -1 & 4
\end{bmatrix}
$$

这个矩阵是对称正定的，适合用 Jacobi 和 GS 方法求解。

---

### 步骤 1：分解矩阵

#### 对角部分 $ D $

$$
D = 
\begin{bmatrix}
4 & 0 & 0 \\
0 & 4 & 0 \\
0 & 0 & 4
\end{bmatrix}
$$

#### 严格下三角部分 $ -L $

$$
-L = 
\begin{bmatrix}
0 & 0 & 0 \\
-1 & 0 & 0 \\
0 & -1 & 0
\end{bmatrix}, \quad
L = 
\begin{bmatrix}
0 & 0 & 0 \\
1 & 0 & 0 \\
0 & 1 & 0
\end{bmatrix}
$$

#### 严格上三角部分 $ -U $

$$
-U = 
\begin{bmatrix}
0 & -1 & 0 \\
0 & 0 & -1 \\
0 & 0 & 0
\end{bmatrix}, \quad
U = 
\begin{bmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
0 & 0 & 0
\end{bmatrix}
$$

---

### 步骤 2：Jacobi 迭代矩阵

$$
J = D^{-1}(L + U)
$$

先求 $ D^{-1} $：

$$
D^{-1} = 
\begin{bmatrix}
\frac{1}{4} & 0 & 0 \\
0 & \frac{1}{4} & 0 \\
0 & 0 & \frac{1}{4}
\end{bmatrix}
$$

再算 $ L + U $：

$$
L + U = 
\begin{bmatrix}
0 & 1 & 0 \\
1 & 0 & 1 \\
0 & 1 & 0
\end{bmatrix}
$$

所以 Jacobi 迭代矩阵：

$$
J = 
\begin{bmatrix}
0 & \frac{1}{4} & 0 \\
\frac{1}{4} & 0 & \frac{1}{4} \\
0 & \frac{1}{4} & 0
\end{bmatrix}
$$

---

### 步骤 3：Gauss-Seidel 迭代矩阵

$$
G = (D - L)^{-1} U
$$

先计算 $ D - L $：

$$
D - L = 
\begin{bmatrix}
4 & 0 & 0 \\
-1 & 4 & 0 \\
0 & -1 & 4
\end{bmatrix}
$$

接下来需要计算它的逆矩阵：

$$
(D - L)^{-1} = ?
$$

为了简化，我们可以使用数值方法或公式来求。这里给出结果（你可以验证）：

$$
(D - L)^{-1} = 
\begin{bmatrix}
\frac{1}{4} & 0 & 0 \\
\frac{1}{16} & \frac{1}{4} & 0 \\
\frac{1}{64} & \frac{1}{16} & \frac{1}{4}
\end{bmatrix}
$$

然后乘以 $ U $：

$$
U = 
\begin{bmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
0 & 0 & 0
\end{bmatrix}
$$

计算 $ G = (D - L)^{-1} U $：

$$
G = 
\begin{bmatrix}
0 & \frac{1}{4} & \frac{1}{16} \\
0 & \frac{1}{16} & \frac{1}{4} \\
0 & \frac{1}{64} & \frac{1}{16}
\end{bmatrix}
$$

---

## 三、最终结果总结

对于矩阵：

$$
A = 
\begin{bmatrix}
4 & -1 & 0 \\
-1 & 4 & -1 \\
0 & -1 & 4
\end{bmatrix}
$$

对应的：

### Jacobi 迭代矩阵：

$$
J = 
\begin{bmatrix}
0 & \frac{1}{4} & 0 \\
\frac{1}{4} & 0 & \frac{1}{4} \\
0 & \frac{1}{4} & 0
\end{bmatrix}
$$

### Gauss-Seidel 迭代矩阵：

$$
G = 
\begin{bmatrix}
0 & \frac{1}{4} & \frac{1}{16} \\
0 & \frac{1}{16} & \frac{1}{4} \\
0 & \frac{1}{64} & \frac{1}{16}
\end{bmatrix}
$$

---

如果你有特定的矩阵 $ A $，我可以帮你代入并重新计算相应的 Jacobi 和 GS 迭代矩阵。欢迎继续提问！