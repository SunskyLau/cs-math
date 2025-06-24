# QR分解
## **1. 最小二乘问题（Least Square Approximation）**
### **1.1 问题背景**
- **超定方程组**：当方程数量 $ m $ 大于未知数数量 $ n $（即 $ A \in \mathbb{R}^{m \times n}, m > n $）时，通常没有精确解。
- **目标**：找到向量 $ x \in \mathbb{R}^n $，使得残差 $ \|b - Ax\|_2 $ 最小。

**例子**：  
假设我们有以下超定方程组：
$$
\begin{cases}
x + y = 1 \\
2x + y = 2 \\
x - y = 0 \\
\end{cases}
\Rightarrow A = \begin{bmatrix} 1 & 1 \\ 2 & 1 \\ 1 & -1 \end{bmatrix}, \quad b = \begin{bmatrix} 1 \\ 2 \\ 0 \end{bmatrix}
$$
由于三个方程无法同时满足，我们需要找到 $ x, y $ 使得 $ \|b - Ax\|_2 $ 最小。

### **1.2 正规方程（Normal Equations）**
- **公式**：$ A^T A x = A^T b $
- **推导**：对目标函数 $ f(x) = \|b - Ax\|_2^2 $ 求导并令导数为零，得到上述方程。
$$  \|v\|_2 = \sqrt{v_1^2 + v_2^2 + \cdots + v_m^2} $$

**例子**：  
计算 $ A^T A $ 和 $ A^T b $：
$$
A^T A = \begin{bmatrix} 1 & 2 & 1 \\ 1 & 1 & -1 \end{bmatrix} \begin{bmatrix} 1 & 1 \\ 2 & 1 \\ 1 & -1 \end{bmatrix} = \begin{bmatrix} 6 & 4 \\ 4 & 3 \end{bmatrix}, \quad A^T b = \begin{bmatrix} 1 & 2 & 1 \\ 1 & 1 & -1 \end{bmatrix} \begin{bmatrix} 1 \\ 2 \\ 0 \end{bmatrix} = \begin{bmatrix} 5 \\ 3 \end{bmatrix}
$$
解方程组：
$$
\begin{cases}
6x + 4y = 5 \\
4x + 3y = 3 \\
\end{cases}
\Rightarrow x = 1.5, \quad y = -1.0
$$

### **1.3 正规方程的缺陷**
- **数值不稳定性**：$ A^T A $ 的条件数是 $ A $ 的平方，导致计算误差放大。
- **解决方法**：使用 QR 分解代替正规方程。

---
## **2. QR 分解**
### **2.1 基本概念**
- **定义**：将矩阵 $ A \in \mathbb{R}^{m \times n} $ 分解为 $ A = QR $，其中：
   - $ Q \in \mathbb{R}^{m \times m} $ 是正交矩阵（$ Q^T Q = I $）。
   - $ R \in \mathbb{R}^{m \times n} $ 是上三角矩阵。
- **用途**：用于最小二乘问题，避免 $ A^T A $ 的数值问题。

### **2.2 Gram-Schmidt 正交化**
- **步骤**：
   1. 从矩阵 $ A $ 的列向量 $ a_1, a_2, \dots, a_n $ 出发。
   2. **逐步正交化并归一化**，得到一组标准正交基 $ q_1, q_2, \dots, q_n $。

> ✅ 所谓"正交化"，是指通过减去已有正交向量在当前向量上的投影，使得新的向量与之前所有向量垂直；  
> "归一化"则是将其长度变为 1，形成单位向量。

**例子**：  
对矩阵  
$$
A = \begin{bmatrix} 1 & 1 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}
$$  
进行 Gram-Schmidt 分解，得到分解结果：$ A = QR $

其中：
- $ Q $ 是一个列向量为标准正交基的矩阵（即正交矩阵）
- $ R $ 是一个上三角矩阵

#### **第一步：正交化第一个向量**

我们从第一个列向量 $ a_1 $ 开始：

$$
a_1 = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}
$$

由于是第一步，没有其他向量可以影响它，所以直接设：

$$
v_1 = a_1 = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}
$$

然后进行归一化（使其长度为 1）：

$$
\|v_1\| = \sqrt{1^2 + 0^2 + 0^2} = 1,\quad
q_1 = \frac{v_1}{\|v_1\|} = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}
$$

#### **第二步：正交化第二个向量**

现在处理第二个列向量：

$$
a_2 = \begin{bmatrix} 1 \\ 1 \\ 0 \end{bmatrix}
$$

我们要做的是：**从 $ a_2 $ 中减去它在已有的正交向量 $ q_1 $ 上的投影**，从而得到一个新的向量 $ v_2 $，它与 $ q_1 $ 正交。

**投影公式回顾**：

$$
\text{proj}_{q_1}(a_2) = (q_1^T a_2)\, q_1
$$

计算内积：

$$
q_1^T a_2 = 1 \cdot 1 + 0 \cdot 1 + 0 \cdot 0 = 1
$$

所以：

$$
\text{proj}_{q_1}(a_2) = 1 \cdot \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}
$$

于是：

$$
v_2 = a_2 - \text{proj}_{q_1}(a_2) = \begin{bmatrix} 1 \\ 1 \\ 0 \end{bmatrix} - \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}
$$

再归一化：

$$
\|v_2\| = \sqrt{0^2 + 1^2 + 0^2} = 1,\quad
q_2 = \frac{v_2}{\|v_2\|} = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}
$$

#### **构造矩阵 $ Q $ 和 $ R $**

我们现在得到了两个标准正交向量 $ q_1 $ 和 $ q_2 $，它们组成矩阵 $ Q $：

$$
Q = \begin{bmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}
$$

接下来构造上三角矩阵 $ R $，其元素由以下关系决定：

$$
A = QR
$$

我们可以反推 $ R $：

$$
R = Q^T A
$$

先计算 $ Q^T $：

$$
Q^T = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{bmatrix}
$$

再计算：

$$
R = Q^T A = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{bmatrix}
\begin{bmatrix} 1 & 1 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}
= \begin{bmatrix} 1 & 1 \\ 0 & 1 \end{bmatrix}
$$

**最终结果**：

$$
A = QR =
\begin{bmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}
\cdot
\begin{bmatrix} 1 & 1 \\ 0 & 1 \end{bmatrix}
=
\begin{bmatrix} 1 & 1 \\ 0 & 1 \\ 0 & 0 \end{bmatrix}
$$
完美还原原始矩阵！


### **2.3 Householder 变换**
- **思想**：通过反射变换将矩阵 $ A $ 逐步变为上三角矩阵。
- **反射矩阵构造**：对向量 $ x \in \mathbb{R}^m $，构造反射向量 $ u = \|x\| e_1 - x $，反射矩阵为：
  $$
  F = I - 2 \frac{u u^T}{u^T u}
  $$

**例子**：  
对向量 $ x = \begin{bmatrix} 12 \\ 6 \\ -4 \end{bmatrix} $ 构造反射矩阵。

1. 计算 $ \|x\| = \sqrt{12^2 + 6^2 + (-4)^2} = 14 $。
2. 构造 $ u = \|x\| e_1 - x = \begin{bmatrix} 14 \\ 0 \\ 0 \end{bmatrix} - \begin{bmatrix} 12 \\ 6 \\ -4 \end{bmatrix} = \begin{bmatrix} 2 \\ -6 \\ 4 \end{bmatrix} $。
3. 计算 $ u^T u = 2^2 + (-6)^2 + 4^2 = 56 $。
4. 构造反射矩阵：
   $$
   F = I - 2 \cdot \frac{uu^T}{56} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix} - \frac{2}{56} \begin{bmatrix} 4 & -12 & 8 \\ -12 & 36 & -24 \\ 8 & -24 & 16 \end{bmatrix} = \begin{bmatrix} \frac{6}{7} & \frac{3}{7} & -\frac{2}{7} \\ \frac{3}{7} & -\frac{2}{7} & \frac{6}{7} \\ -\frac{2}{7} & \frac{6}{7} & \frac{3}{7} \end{bmatrix}
   $$

### **2.4 Givens 旋转**
- **思想**：通过平面旋转变换零化特定元素。
- **旋转矩阵**：对第 $ i $ 行和第 $ k $ 行，构造：
  $$
  G(i,k,\theta) = \begin{bmatrix} 1 & \cdots & 0 & \cdots & 0 & \cdots & 0 \\ \vdots & \ddots & \vdots & & \vdots & & \vdots \\ 0 & \cdots & c & \cdots & -s & \cdots & 0 \\ \vdots & & \vdots & \ddots & \vdots & & \vdots \\ 0 & \cdots & s & \cdots & c & \cdots & 0 \\ \vdots & & \vdots & & \vdots & \ddots & \vdots \\ 0 & \cdots & 0 & \cdots & 0 & \cdots & 1 \end{bmatrix}
  $$
  其中 $ c = \cos\theta, s = \sin\theta $。

**例子**：  
将矩阵 $ A = \begin{bmatrix} 6 & 5 \\ 5 & 1 \\ 0 & 4 \end{bmatrix} $ 的第2行第1列元素变为0。

1. 取子向量 $ [6, 5] $，计算旋转参数：
   $$
   r = \sqrt{6^2 + 5^2} = 7.8102, \quad c = \frac{6}{r} = 0.7682, \quad s = -\frac{5}{r} = -0.6402
   $$
2. 构造 Givens 矩阵：
   $$
   G = \begin{bmatrix} c & -s & 0 \\ s & c & 0 \\ 0 & 0 & 1 \end{bmatrix} = \begin{bmatrix} 0.7682 & 0.6402 & 0 \\ -0.6402 & 0.7682 & 0 \\ 0 & 0 & 1 \end{bmatrix}
   $$
3. 应用变换：
   $$
   GA = \begin{bmatrix} 7.8102 & 4.4813 \\ 0 & -2.4327 \\ 0 & 4 \end{bmatrix}
   $$

---
## **3. 最小二乘问题的 QR 解法**
### **3.1 步骤**
1. 对矩阵 $ A $ 进行 QR 分解：$ A = QR $。
2. 将残差 $ \|b - Ax\|_2 $ 转化为：
   $$
   \|b - Ax\|_2 = \|Q^T b - Rx\|_2
   $$
3. 分解 $ Q^T b = \begin{bmatrix} c_1 \\ c_2 \end{bmatrix} $，其中 $ c_1 \in \mathbb{R}^n $。
4. 解上三角方程 $ R_1 x = c_1 $，其中 $ R_1 $ 是 $ R $ 的前 $ n $ 行。

**例子**：  
用 QR 分解求解之前的超定方程组：
$$
A = \begin{bmatrix} 1 & 1 \\ 2 & 1 \\ 1 & -1 \end{bmatrix}, \quad b = \begin{bmatrix} 1 \\ 2 \\ 0 \end{bmatrix}
$$

1. 对 $ A $ 进行 QR 分解：
   $$
   Q = \begin{bmatrix} \frac{1}{\sqrt{6}} & \frac{1}{\sqrt{2}} \\ \frac{2}{\sqrt{6}} & 0 \\ \frac{1}{\sqrt{6}} & -\frac{1}{\sqrt{2}} \end{bmatrix}, \quad R = \begin{bmatrix} \sqrt{6} & \frac{2}{\sqrt{6}} \\ 0 & \sqrt{2} \end{bmatrix}
   $$

2. 计算 $ Q^T b = \begin{bmatrix} \frac{1 \cdot 1 + 2 \cdot 2 + 1 \cdot 0}{\sqrt{6}} \\ \frac{1 \cdot 1 + 0 \cdot 2 + (-1) \cdot 0}{\sqrt{2}} \end{bmatrix} = \begin{bmatrix} \frac{5}{\sqrt{6}} \\ \frac{1}{\sqrt{2}} \end{bmatrix} $

3. 解方程 $ R x = Q^T b $：
   $$
   \begin{cases}
   \sqrt{6} x + \frac{2}{\sqrt{6}} y = \frac{5}{\sqrt{6}} \\
   \sqrt{2} y = \frac{1}{\sqrt{2}} \\
   \end{cases}
   \Rightarrow y = 0.5, \quad x = \frac{5 - 2 \cdot 0.5}{6} = 0.6667
   $$

---

---
## **4. 三种 QR 分解方法对比**
| 方法          | 优点                          | 缺点                  |
|---------------|-------------------------------|-----------------------|
| Gram-Schmidt  | 易于理解，适合教学            | 数值稳定性较差        |
| Householder   | 数值稳定性好，适合稠密矩阵    | 实现复杂度较高        |
| Givens 旋转   | 适合稀疏矩阵，可选择性零化元素| 计算量较大，需多次迭代|

---
## **5. 期末考试重点总结**
1. **最小二乘问题**：理解正规方程的推导和缺陷。
2. **QR 分解**：掌握 Gram-Schmidt 和 Householder 的步骤，能手动计算简单矩阵的分解。
3. **QR 解最小二乘**：能通过 QR 分解求解 $ x $，并理解其几何意义（投影到列空间）。
4. **方法对比**：了解三种 QR 分解方法的适用场景和优缺点。

**考试题型预测**：
- 给定矩阵 $ A $ 和向量 $ b $，要求用 QR 分解求最小二乘解。
- 判断某种 QR 分解方法是否适用于特定矩阵（如稀疏矩阵）。
- 推导正规方程或 QR 分解的某一步骤。


---
## **FINAL: 最小二乘问题的QR解法完整总结**

### **1. 问题定义**
给定超定线性方程组 $A \mathbf{x} = \mathbf{b}$（$A \in \mathbb{R}^{m \times n}, m > n$），最小二乘问题求解：
$$
\min_{\mathbf{x}} \| A \mathbf{x} - \mathbf{b} \|_2^2
$$

---

### **2. QR分解的两种形式**
#### **(1) 完全QR分解（Full QR）**
- 分解形式：
   $$
   A = Q R = \begin{bmatrix} Q_1 & Q_2 \end{bmatrix} \begin{bmatrix} R_1 \\ 0 \end{bmatrix}
   $$
   - $Q \in \mathbb{R}^{m \times m}$ 是正交矩阵（$Q^T Q = I_m$），
   - $Q_1 \in \mathbb{R}^{m \times n}$ 是 $A$ 的列空间的正交基，
   - $Q_2 \in \mathbb{R}^{m \times (m-n)}$ 是 $A$ 的左零空间的正交基，
   - $R_1 \in \mathbb{R}^{n \times n}$ 是上三角矩阵。

#### **(2) 瘦QR分解（Skinny QR）**
- 分解形式：
   $$
   A = Q R
   $$
   - $Q \in \mathbb{R}^{m \times n}$ 列正交（$Q^T Q = I_n$），
   - $R \in \mathbb{R}^{n \times n}$ 是上三角矩阵。

---

### **3. 残差转化的数学推导**
#### **(1) 完全QR分解下的残差转化**
$$
\| A \mathbf{x} - \mathbf{b} \|_2^2 = \left\| \begin{bmatrix} Q_1 & Q_2 \end{bmatrix} \begin{bmatrix} R_1 \\ 0 \end{bmatrix} \mathbf{x} - \mathbf{b} \right\|_2^2
$$
利用正交矩阵性质 $Q^T Q = I$：
$$
= \| Q_1 R_1 \mathbf{x} - \mathbf{b} \|_2^2 = \| R_1 \mathbf{x} - Q_1^T \mathbf{b} \|_2^2 + \| Q_2^T \mathbf{b} \|_2^2
$$
最小化残差等价于解：
$$
R_1 \mathbf{x} = Q_1^T \mathbf{b}
$$

#### **(2) 瘦QR分解下的残差转化**
$$
\| A \mathbf{x} - \mathbf{b} \|_2^2 = \| Q R \mathbf{x} - \mathbf{b} \|_2^2
$$
利用 $Q^T Q = I_n$ 和投影矩阵性质：
$$
= \| Q (R \mathbf{x} - Q^T \mathbf{b}) + (I - Q Q^T) \mathbf{b} \|_2^2
$$
由于 $Q (R \mathbf{x} - Q^T \mathbf{b})$ 与 $(I - Q Q^T) \mathbf{b}$ 正交：
$$
= \| R \mathbf{x} - Q^T \mathbf{b} \|_2^2 + \| (I - Q Q^T) \mathbf{b} \|_2^2
$$
最小化残差等价于解：
$$
R \mathbf{x} = Q^T \mathbf{b}
$$

---

### **4. 算法步骤**
#### **(1) 完全QR分解解法**
1. **计算完全QR分解**：
    $$
    A = Q R = \begin{bmatrix} Q_1 & Q_2 \end{bmatrix} \begin{bmatrix} R_1 \\ 0 \end{bmatrix}
    $$
    （使用Householder变换或Givens旋转）。
2. **计算 $Q_1^T \mathbf{b}$**：
    $$
    \mathbf{c}_1 = Q_1^T \mathbf{b}
    $$
3. **解上三角系统**：
    $$
    R_1 \mathbf{x} = \mathbf{c}_1
    $$
4. **残差计算**（可选）：
    $$
    \text{残差} = \| Q_2^T \mathbf{b} \|_2
    $$

#### **(2) 瘦QR分解解法**
1. **计算瘦QR分解**：
    $$
    A = Q R
    $$
    （使用Gram-Schmidt或Householder变换的紧凑形式）。
2. **计算 $Q^T \mathbf{b}$**：
    $$
    \mathbf{c} = Q^T \mathbf{b}
    $$
3. **解上三角系统**：
    $$
    R \mathbf{x} = \mathbf{c}
    $$
4. **残差计算**（可选）：
    $$
    \text{残差} = \| \mathbf{b} - Q R \mathbf{x} \|_2
    $$

---

### **5. 数值稳定性与计算效率**
| **方法**               | **稳定性**          | **适用场景**                     |
|------------------------|---------------------|----------------------------------|
| **Householder QR**     | 最优（推荐）         | 通用问题，尤其是稠密矩阵         |
| **Gram-Schmidt (MGS)** | 较好（修正后稳定）   | 需要显式 $Q$ 时              |
| **Givens QR**         | 高（适合稀疏矩阵）   | 稀疏矩阵或结构化矩阵             |

---

### **6. 关键点总结**
1. **残差转化的普适性**：
    - 无论完全QR还是瘦QR分解，最小二乘问题均可转化为上三角系统的求解。
    - 瘦QR分解中，$Q$ 的 $m \times n$ 形式不影响残差最小化的正确性。
2. **正交性保证**：
    - 完全QR分解的 $Q$ 是方阵，瘦QR分解的 $Q$ 列正交。
    - 瘦QR的 $Q Q^T$ 是投影矩阵，$I - Q Q^T$ 计算残差的垂直分量。
3. **实际应用选择**：
    - 显式需要 $Q$ 时（如迭代算法），用瘦QR（MGS）。
    - 仅求最小二乘解时，用Householder QR（无需显式 $Q$）。

---

### **7. 示例验证**
以 $A = \begin{bmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \end{bmatrix}, \mathbf{b} = \begin{bmatrix} 2 \\ 3 \\ 4 \end{bmatrix}$ 为例：
1. **瘦QR分解**（Gram-Schmidt）：
    - $Q = \begin{bmatrix} \frac{1}{\sqrt{3}} & -\frac{1}{\sqrt{2}} \\ \frac{1}{\sqrt{3}} & 0 \\ \frac{1}{\sqrt{3}} & \frac{1}{\sqrt{2}} \end{bmatrix}, R = \begin{bmatrix} \sqrt{3} & 2\sqrt{3} \\ 0 & \sqrt{2} \end{bmatrix}$
    - 解 $R \mathbf{x} = Q^T \mathbf{b}$ 得 $\mathbf{x} = [1, 1]^T$（精确解）。
2. **完全QR分解**：
    - 计算 $Q_1^T \mathbf{b}$ 和 $R_1$，结果一致。

---

### **最终结论**
QR解法通过正交分解将最小二乘问题转化为上三角系统的求解，**无论是完全QR还是瘦QR分解，均能保证残差最小化的正确性**。实际应用中：
- 优先选择数值稳定的Householder QR；
- 瘦QR分解节省存储且适用于列满秩矩阵；
- 残差计算时需注意投影分量的处理。