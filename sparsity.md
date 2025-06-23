# Lecture 5: Sparsity（稀疏性）详细讲解

---

## 一、什么是稀疏性（Sparsity）？

**定义：**
一个向量是稀疏的，如果它只有“少数”非零元素。  
例如：
- 向量 [1, 0, 0, 0, 2] 是稀疏的，因为只有两个非零值。
- 向量 [1, 2, 3, 4, 5] 不是稀疏的，因为所有元素都非零。

**为什么稀疏性重要？**
在信号处理、图像压缩、基因分析等领域中，很多数据本身具有稀疏特性，或者可以在某些变换域（如傅里叶变换、小波变换等）下变得稀疏。利用这种性质可以简化问题、减少计算复杂度，并提高模型性能。

---

## 二、L0、L1 和 L2 范数（Norms）

### 1. **L0 范数**
- 定义：向量中非零元素的数量。
- 例子：
  - x = [1, 0, 0, 2] → ||x||₀ = 2
- 特点：直接衡量稀疏性，但它是非凸且 NP-hard 的优化目标。

### 2. **L1 范数**
- 定义：向量各元素绝对值之和。
- 公式：$$ \|x\|_1 = \sum_{i=1}^n |x_i| $$
- 例子：
  - x = [1, -2, 3] → ||x||₁ = |1| + |-2| + |3| = 6
- 特点：用于近似 L0 范数，是凸函数，可以用线性规划求解。

### 3. **L2 范数**
- 定义：向量各元素平方和的开方。
- 公式：$$ \|x\|_2 = \sqrt{\sum_{i=1}^n x_i^2} $$
- 例子：
  - x = [3, 4] → ||x||₂ = √(3² + 4²) = 5
- 特点：常用于最小二乘法，能防止过拟合。

---

## 三、线性逆问题（Linear Inverse Problems）

### 1. 满秩系统（Full-determined system）
- 方程个数 = 未知数个数
- 解唯一（若矩阵 A 可逆）：
  $$ x^* = A^{-1}b $$

**例题：**
设 A = [[1, 2], [3, 4]], b = [5, 6]

$$
A^{-1} = \frac{1}{(1)(4)-(2)(3)} \begin{bmatrix} 4 & -2 \\ -3 & 1 \end{bmatrix} = \frac{1}{-2} \begin{bmatrix} 4 & -2 \\ -3 & 1 \end{bmatrix}
= \begin{bmatrix} -2 & 1 \\ 1.5 & -0.5 \end{bmatrix}
$$

$$
x = A^{-1}b = \begin{bmatrix} -2 & 1 \\ 1.5 & -0.5 \end{bmatrix} \begin{bmatrix} 5 \\ 6 \end{bmatrix}
= \begin{bmatrix} (-2)(5) + (1)(6) \\ (1.5)(5) + (-0.5)(6) \end{bmatrix}
= \begin{bmatrix} -4 \\ 4.5 \end{bmatrix}
$$

---

### 2. 超定系统（Over-determined system）
- 方程个数 > 未知数个数
- 最小二乘解（Least Squares Solution）：
  $$ x^* = (A^TA)^{-1}A^Tb $$

**例题：**
设 A = [[1, 2], [3, 4], [5, 6]], b = [7, 8, 9]

先算 Gram 矩阵：
$$
A^T = \begin{bmatrix} 1 & 3 & 5 \\ 2 & 4 & 6 \end{bmatrix}, \quad
A^TA = \begin{bmatrix} 1^2+3^2+5^2 & 1×2+3×4+5×6 \\
2×1+4×3+6×5 & 2^2+4^2+6^2 \end{bmatrix}
= \begin{bmatrix} 35 & 44 \\ 44 & 56 \end{bmatrix}
$$

再求逆：
$$
(A^TA)^{-1} = \frac{1}{(35)(56)-(44)^2} \begin{bmatrix} 56 & -44 \\ -44 & 35 \end{bmatrix}
= \frac{1}{(1960 - 1936)} \begin{bmatrix} 56 & -44 \\ -44 & 35 \end{bmatrix}
= \frac{1}{24} \begin{bmatrix} 56 & -44 \\ -44 & 35 \end{bmatrix}
$$

最后：
$$
A^Tb = \begin{bmatrix} 1×7 + 3×8 + 5×9 \\ 2×7 + 4×8 + 6×9 \end{bmatrix}
= \begin{bmatrix} 7 + 24 + 45 \\ 14 + 32 + 54 \end{bmatrix}
= \begin{bmatrix} 76 \\ 100 \end{bmatrix}
$$

$$
x^* = (A^TA)^{-1}A^Tb = \frac{1}{24} \begin{bmatrix} 56×76 - 44×100 \\ -44×76 + 35×100 \end{bmatrix}
= \frac{1}{24} \begin{bmatrix} 4256 - 4400 \\ -3344 + 3500 \end{bmatrix}
= \frac{1}{24} \begin{bmatrix} -144 \\ 156 \end{bmatrix}
= \begin{bmatrix} -6 \\ 6.5 \end{bmatrix}
$$

---

### 3. 欠定系统（Under-determined system）
- 方程个数 < 未知数个数
- 有无穷多解，需加正则化项选择最优解

#### （1）最小 L2 范数解（Minimum L2 Norm Solution）：
$$
x^* = A^T(AA^T)^{-1}b
$$

**例题：**
设 A = [[1, 2, 3], [4, 5, 6]], b = [7, 8]

$$
AA^T = \begin{bmatrix} 1^2+2^2+3^2 & 1×4+2×5+3×6 \\
4×1+5×2+6×3 & 4^2+5^2+6^2 \end{bmatrix}
= \begin{bmatrix} 14 & 32 \\ 32 & 77 \end{bmatrix}
$$

求逆：
$$
(AA^T)^{-1} = \frac{1}{14×77 - 32^2} \begin{bmatrix} 77 & -32 \\ -32 & 14 \end{bmatrix}
= \frac{1}{1078 - 1024} \begin{bmatrix} 77 & -32 \\ -32 & 14 \end{bmatrix}
= \frac{1}{54} \begin{bmatrix} 77 & -32 \\ -32 & 14 \end{bmatrix}
$$

$$
A^T = \begin{bmatrix} 1 & 4 \\ 2 & 5 \\ 3 & 6 \end{bmatrix}, \quad
A^Tb = \begin{bmatrix} 1×7 + 4×8 \\ 2×7 + 5×8 \\ 3×7 + 6×8 \end{bmatrix}
= \begin{bmatrix} 7 + 32 \\ 14 + 40 \\ 21 + 48 \end{bmatrix}
= \begin{bmatrix} 39 \\ 54 \\ 69 \end{bmatrix}
$$

$$
x^* = A^T(AA^T)^{-1}b = \frac{1}{54} \begin{bmatrix} 1 & 4 \\ 2 & 5 \\ 3 & 6 \end{bmatrix}
\begin{bmatrix} 77×7 - 32×8 \\ -32×7 + 14×8 \end{bmatrix}
= \cdots
$$

---

## 四、稀疏优化方法（Sparse Optimization Approaches）

### 1. 匹配追踪（Matching Pursuit）

匹配追踪算法是一种用于求解稀疏表示问题的迭代贪心算法。其核心思想是在每一步迭代中，从给定的字典中选择最能匹配当前残差的原子（基向量），从而逐步构建信号的稀疏表示。

#### 算法目标
给定信号 $\mathbf{b} \in \mathbb{R}^m$ 和字典矩阵 $\mathbf{A} = [\mathbf{a}_1, \mathbf{a}_2, ..., \mathbf{a}_n] \in \mathbb{R}^{m \times n}$（其中每列 $\mathbf{a}_i$ 为标准化的基向量），找到一个稀疏向量 $\mathbf{x}$，使得 $\mathbf{Ax} \approx \mathbf{b}$。

#### 详细算法步骤
1. **初始化**：
   - 设置残差 $\mathbf{r}_0 = \mathbf{b}$
   - 设置迭代计数器 $k = 0$
   - 初始化系数向量 $\mathbf{x} = \mathbf{0}$

2. **迭代过程**：
   - 计算所有字典原子与当前残差的内积：
     $$ \text{correlation}_j = |\langle \mathbf{r}_k, \mathbf{a}_j \rangle| = |\mathbf{a}_j^\top \mathbf{r}_k|, \quad j = 1,\ldots,n $$
   - 选择相关性最大的原子：
     $$ j^* = \arg\max_j |\langle \mathbf{r}_k, \mathbf{a}_j \rangle| $$
   - 计算投影系数：
     $$ \alpha_k = \langle \mathbf{r}_k, \mathbf{a}_{j^*} \rangle $$
   - 更新对应位置的系数：
     $$ x_{j^*} = x_{j^*} + \alpha_k $$
   - 更新残差：
     $$ \mathbf{r}_{k+1} = \mathbf{r}_k - \alpha_k \mathbf{a}_{j^*} $$
   - $k = k + 1$

3. **终止条件**（满足任一即可）：
   - 残差范数小于预设阈值：$\|\mathbf{r}_k\|_2 < \varepsilon$
   - 达到最大迭代次数：$k > K_{\text{max}}$
   - 最大相关性小于阈值：$\max_{j} |\langle \mathbf{r}_k, \mathbf{a}_j \rangle| < \delta$

#### 算法示例
考虑一个简单的信号分解问题：

**输入**：
- 信号向量 $\mathbf{b} = [3, 2, 1]^T$
- 字典矩阵 $\mathbf{A} = \begin{bmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 0 \end{bmatrix}$
- 停止阈值 $\varepsilon = 0.1$

**迭代过程**：

1. **初始化**：
   - $\mathbf{r}_0 = \mathbf{b} = [3, 2, 1]^T$
   - $\mathbf{x} = [0, 0, 0]^T$
   - $k = 0$

2. **第一次迭代**：
   - 计算相关性公式：$|\langle \mathbf{r}_k, \mathbf{a}_j \rangle| = |\mathbf{a}_j^\top \mathbf{r}_k|$
   - 代入计算：
     $|\langle \mathbf{r}_0, \mathbf{a}_1 \rangle| = |3\times1 + 2\times0 + 1\times1| = |3+0+1| = 4$
     $|\langle \mathbf{r}_0, \mathbf{a}_2 \rangle| = |3\times0 + 2\times1 + 1\times1| = |0+2+1| = 3$
     $|\langle \mathbf{r}_0, \mathbf{a}_3 \rangle| = |3\times1 + 2\times1 + 1\times0| = |3+2+0| = 5$
   - 选择 $j^* = 3$（最大相关性）
   - 投影系数公式：$\alpha_k = \langle \mathbf{r}_k, \mathbf{a}_{j^*} \rangle$
   - 代入计算：$\alpha_0 = \mathbf{a}_3^\top \mathbf{r}_0 = [1, 1, 0] \cdot [3, 2, 1]^T = 1\times3 + 1\times2 + 0\times1 = 3+2+0 = 5$
   - 更新系数公式：$x_{j^*} = x_{j^*} + \alpha_k$
   - 代入计算：$x_3 = 0 + 5 = 5$
   - 残差更新公式：$\mathbf{r}_{k+1} = \mathbf{r}_k - \alpha_k \mathbf{a}_{j^*}$
   - 代入计算：$\mathbf{r}_1 = \mathbf{r}_0 - 5\mathbf{a}_3 = [3, 2, 1]^T - 5[1, 1, 0]^T = [3-5, 2-5, 1-0]^T = [-2, -3, 1]^T$

3. **第二次迭代**：
   - 使用相关性公式计算：
     $|\langle \mathbf{r}_1, \mathbf{a}_1 \rangle| = |-2\times1 + (-3)\times0 + 1\times1| = |-2+0+1| = |-1| = 1$
     $|\langle \mathbf{r}_1, \mathbf{a}_2 \rangle| = |-2\times0 + (-3)\times1 + 1\times1| = |0-3+1| = |-2| = 2$
     $|\langle \mathbf{r}_1, \mathbf{a}_3 \rangle| = |-2\times1 + (-3)\times1 + 1\times0| = |-2-3+0| = |-5| = 5$
   - 选择 $j^* = 3$ (最大相关性)
   - 使用投影系数公式计算：$\alpha_1 = \langle \mathbf{r}_1, \mathbf{a}_3 \rangle = [-2, -3, 1] \cdot [1, 1, 0]^T = -2\times1 + (-3)\times1 + 1\times0 = -2-3+0 = -5$
   - 使用系数更新公式：$x_3 = x_3 + \alpha_1 = 5 + (-5) = 0$
   - 使用残差更新公式计算：$\mathbf{r}_2 = \mathbf{r}_1 - \alpha_1 \mathbf{a}_3 = [-2, -3, 1]^T - (-5)[1, 1, 0]^T = [-2, -3, 1]^T + 5[1, 1, 0]^T = [-2+5, -3+5, 1+0]^T = [3, 2, 1]^T$

4. **终止**：
   残差范数计算公式：$\|\mathbf{r}_k\|_2 = \sqrt{\sum_{i=1}^n r_i^2}$
   代入计算得 $\|\mathbf{r}_2\|_2 = \|[3, 2, 1]^T\|_2 = \sqrt{3^2 + 2^2 + 1^2} = \sqrt{9+4+1} = \sqrt{14} \approx 3.74$
   因为 $\|\mathbf{r}_2\|_2 \approx 3.74 \not< \varepsilon = 0.1$，所以算法不应终止，需要继续迭代。

**输出**：
最终得到的稀疏表示为 $\mathbf{x} = [0, 0, 0]^T$ (此结果为当前迭代步骤的输出，算法尚未收敛)

---

### 2. 平滑重构技巧（Smooth Reformulation Tricks）

#### （1）正负拆分法（Positive-Negative Split Trick）
将变量 x 分成正部分 p 和负部分 n：
$$
x = p - n,\quad p,n \geq 0
$$

这样可以把 L1 范数转换为线性规划问题。

**例题：**
最小化：
$$
\min_x \|x\|_1 \quad s.t.\ Ax = b
$$

转化为：
$$
\min_{p,n} \mathbf{1}^T(p+n),\quad A(p-n)=b,\quad p,n \geq 0
$$

---

### 3. 字典学习（Dictionary Learning）

目标是找到一个字典 A 和稀疏表示 X，使得 AX ≈ B，且 X 尽可能稀疏。

**算法步骤：**
1. 固定 A，求 X（LASSO 问题）
2. 固定 X，求 A（最小二乘）
3. 交替迭代直至收敛

**应用：图像去噪、超分辨率重建等**

---

## 五、稀疏性的实际应用举例

### 1. 图像去噪（Image Denoising）
- 模型：$y \approx Ax + e$，其中 $x$ 是稀疏的
- 目标：
$$ \min_x \|Ax - y\|_2^2 + \lambda \|x\|_1 $$

### 2. 图像修复（Image Inpainting）
- 缺失像素用掩码 $R$ 表示
- 模型：
$$ \min_x \|RAx - y\|_2^2 + \lambda \|x\|_1 $$

### 3. 鲁棒人脸识别（Robust Face Recognition）
- 模型：$b = Ax + e$，其中 $x$ 和 $e$ 都是稀疏的
- 利用稀疏性识别不同人脸

### 4. 压缩感知（Compressive Sensing）
- 通过少量测量恢复原始信号
- 条件：信号在某个基下是稀疏的

---

## 六、补充知识点：矩阵微积分（Matrix Calculus）

### 1. 跟踪迹（Trace）运算
- tr(A) = sum of diagonal elements
- tr(AB) = tr(BA)

### 2. 最小二乘推导
- 最小化：
$$
\|Ax - b\|_2^2 = (Ax - b)^T(Ax - b)
= x^T A^T A x - 2x^T A^T b + b^T b
$$

对 x 求导：
$$
\nabla_x \|Ax - b\|^2 = 2A^T A x - 2A^T b
\Rightarrow x^* = (A^T A)^{-1} A^T b
$$

---

## 七、总结复习建议

| 主题 | 核心公式 | 应用 |
|------|----------|------|
| L0 范数 | 非零元素数量 | 衡量稀疏性 |
| L1 范数 | 绝对值和 | LASSO、稀疏建模 |
| L2 范数 | 平方和开根号 | 正则化、最小二乘 |
| 最小二乘 | x* = (AᵀA)⁻¹Aᵀb | 回归、插值 |
| 最小 L2 解 | x* = Aᵀ(AAᵀ)⁻¹b | 欠定系统 |
| 字典学习 | min‖AX - B‖F² + λ‖X‖₁ | 图像处理、语音识别 |
| 匹配追踪 | 贪心选基 | 信号逼近 |