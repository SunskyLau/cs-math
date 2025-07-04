# **基于LU分解的线性方程组求解** 
## **1. 线性方程组与LU分解的核心思想**  
### **1.1 问题描述**  
给定线性方程组 $ Ax = b $，其中 $ A \in \mathbb{R}^{n \times n} $ 非奇异 **(矩阵可逆，行列式不为0，det(A)≠0)**，目标是高效求解 $ x $。  
- **直接法（如LU分解）**：通过矩阵分解将 $ A $ 拆解为更易求解的形式（如三角矩阵），适用于中小规模稠密矩阵。  
- **迭代法**：适用于大规模稀疏矩阵，但此处重点讨论直接法。  

### **1.2 LU分解的定义**  
矩阵 $ A $ 可分解为 $ A = LU $，其中：  
- $ L \in \mathbb{R}^{n \times n} $ 是**单位下三角矩阵**（对角线为1，上方全为0）。  
- $ U \in \mathbb{R}^{n \times n} $ 是**上三角矩阵**（下方全为0）。  

**分解目的**：将原问题 $ Ax = b $ 转化为两个三角形方程组：  
1. 解 $ Lz = b $（前向替换）。  
2. 解 $ Ux = z $（后向替换）。  

---

## **2. 前向替换与后向替换**  
### **2.1 前向替换（解 $ Lz = b $）**  
**原理**：从第一个方程开始，逐步代入求解。  
**公式**：  
$$
z_i = \frac{b_i - \sum_{j=1}^{i-1} L_{ij} z_j}{L_{ii}} \quad (i=1,2,\dots,n)
$$
**示例**：  
设  
$$
L = \begin{bmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 3 & 4 & 1 \end{bmatrix}, \quad b = \begin{bmatrix} 1 \\ 4 \\ 11 \end{bmatrix}
$$
求解步骤：  
1. $ z_1 = \frac{1}{1} = 1 $  
2. $ z_2 = \frac{4 - 2 \cdot z_1}{1} = 2 $  
3. $ z_3 = \frac{11 - 3 \cdot z_1 - 4 \cdot z_2}{1} = 0 $  
结果：$ z = [1, 2, 0]^T $  

### **2.2 后向替换（解 $ Ux = z $）**  
**原理**：从最后一个方程开始，反向代入求解。  
**公式**：  
$$
x_i = \frac{z_i - \sum_{j=i+1}^n U_{ij} x_j}{U_{ii}} \quad (i=n,n-1,\dots,1)
$$
**示例**：  
设  
$$
U = \begin{bmatrix} 1 & 2 & 3 \\ 0 & 4 & 5 \\ 0 & 0 & 6 \end{bmatrix}, \quad z = \begin{bmatrix} 1 \\ 2 \\ 0 \end{bmatrix}
$$
求解步骤：  
1. $ x_3 = \frac{0}{6} = 0 $  
2. $ x_2 = \frac{2 - 5 \cdot 0}{4} = 0.5 $  
3. $ x_1 = \frac{1 - 2 \cdot 0.5 - 3 \cdot 0}{1} = 0 $  
结果：$ x = [0, 0.5, 0]^T $  

---

## **3. 高斯消元法与LU分解**  
### **3.1 核心思想**  
通过**初等变换**将 $ A $ 转化为上三角矩阵 $ U $，并记录消元因子构造 $ L $。  
**步骤**：  
1. 第 $ k $ 步消元：用第 $ k $ 行消去第 $ k+1 $ 到 $ n $ 行的第 $ k $ 列元素。  
2. 记录消元因子 $ l_{ik} = \frac{a_{ik}^{(k-1)}}{a_{kk}^{(k-1)}} $ 作为 $ L(i, k) $ 的元素。  

### **3.2 示例：3阶矩阵的LU分解**  
给定矩阵  
$$
A = \begin{bmatrix} 2 & 4 & 6 \\ 1 & 3 & 5 \\ 3 & 7 & 10 \end{bmatrix}
$$
**步骤1**：消去第2、3行的第1列：  
- 消元因子 $ l_{21} = \frac{1}{2} $，$ l_{31} = \frac{3}{2} $。  
- 更新矩阵：  
$$
A^{(1)} = \begin{bmatrix} 2 & 4 & 6 \\ 0 & 1 & 2 \\ 0 & 1 & 1 \end{bmatrix}
$$
**步骤2**：消去第3行的第2列：  
- 消元因子 $ l_{32} = \frac{1}{1} = 1 $。  
- 更新矩阵：  
$$
U = \begin{bmatrix} 2 & 4 & 6 \\ 0 & 1 & 2 \\ 0 & 0 & -1 \end{bmatrix}, \quad L = \begin{bmatrix} 1 & 0 & 0 \\ 0.5 & 1 & 0 \\ 1.5 & 1 & 1 \end{bmatrix}
$$

---

## **4. 特殊矩阵的分解方法**  
### **4.1 Cholesky分解（对称正定矩阵）**  
若 $ A $ 对称正定，则可分解为 $ A = LL^T $，其中 $ L $ 为下三角矩阵，对角线元素为正。  
**公式推导**：  
$$
A_{ij} = \sum_{k=1}^j L_{ik} L_{jk} \quad (i \geq j)
$$
**示例**：  
设  
$$
A = \begin{bmatrix} 4 & 2 & 2 \\ 2 & 5 & 2 \\ 2 & 2 & 6 \end{bmatrix}
$$
分解步骤：  
我们从矩阵乘法的角度出发。设：

$$
A = LL^T
$$

其中 $ L $ 是一个下三角矩阵，形式如下：

$$
L = \begin{bmatrix}
l_{11} & 0 & 0 \\
l_{21} & l_{22} & 0 \\
l_{31} & l_{32} & l_{33}
\end{bmatrix}, \quad
L^T = \begin{bmatrix}
l_{11} & l_{21} & l_{31} \\
0 & l_{22} & l_{32} \\
0 & 0 & l_{33}
\end{bmatrix}
$$

那么：

$$
LL^T = \begin{bmatrix}
l_{11}^2 & l_{11}l_{21} & l_{11}l_{31} \\
l_{21}l_{11} & l_{21}^2 + l_{22}^2 & l_{21}l_{31} + l_{22}l_{32} \\
l_{31}l_{11} & l_{31}l_{21} + l_{32}l_{22} & l_{31}^2 + l_{32}^2 + l_{33}^2
\end{bmatrix}
$$

与原矩阵 $ A $ 比较，就可以逐步求出每个 $ l_{ij} $ 的值。
1. $ L_{11} = \sqrt{4} = 2 $。  
2. $ L_{21} = \frac{2}{2} = 1 $，$ L_{22} = \sqrt{5 - 1^2} = 2 $。  
3. $ L_{31} = 1 $，$ L_{32} = \frac{2 - 1 \cdot 1}{2} = 0.5 $，$ L_{33} = \sqrt{6 - 1^2 - 0.5^2} = \sqrt{4.75} $。  

结果：  
$$
L = \begin{bmatrix} 2 & 0 & 0 \\ 1 & 2 & 0 \\ 1 & 0.5 & \sqrt{4.75} \end{bmatrix}
$$

---

## **5. 条件数与稳定性分析**  
### **5.1 条件数定义**  
矩阵 $ A $ 的条件数 $ \kappa(A) = \|A\| \cdot \|A^{-1}\| $，衡量方程组对输入扰动的敏感度：  
- $ \kappa(A) \approx 1 $：良态（稳定）。  
- $ \kappa(A) \gg 1 $：病态（不稳定）。  

### **5.2 条件数与误差界**  
若 $ A + \delta A $ 和 $ b + \delta b $ 有扰动，则解的相对误差满足：  
$$
\frac{\|x - \hat{x}\|}{\|x\|} \leq \kappa(A) \left( \frac{\|\delta A\|}{\|A\|} + \frac{\|\delta b\|}{\|b\|} \right)
$$
**示例**：  
设 $ A = \begin{bmatrix} 1 & 1 \\ 1 & 1.0001 \end{bmatrix} $，其条件数 $ \kappa(A) \approx 4 \times 10^4 $，即使扰动很小，解也可能剧烈变化。  

---

## **6. 实际应用与注意事项**  
### **6.1 部分选主元（GEPP）**  
为避免小主元导致数值不稳定，每步消元选择当前列中绝对值最大的元素作为主元，并交换行。  
**示例**：若 $ A = \begin{bmatrix} 0.001 & 1 \\ 1 & 0 \end{bmatrix} $，直接消元会导致除以0.001，误差放大1000倍；选主元后交换行，避免问题。  

### **6.2 病态矩阵的判断**  
- **行列式陷阱**：即使 $ \det(A) \approx 0 $，矩阵可能并非病态（如 $ D_n = \text{diag}(0.1, \dots, 0.1) $）。  
- **条件数优先**：用 $ \kappa(A) $ 直接判断。  

---

## **7. 总结与复习要点**  
1. **LU分解**：掌握分解步骤、存储方式（L和U复用A的存储空间）。  
2. **三角求解**：前向/后向替换的公式与实现（时间复杂度 $ O(n^2) $）。  
3. **稳定性**：理解条件数的意义、GEPP的作用。  
4. **特殊分解**：Cholesky分解的条件与应用。  

**附录：常见矩阵范数与条件数关系**  
- $ \|A\|_2 = \sigma_{\max}(A) $，$ \|A\|_F = \sqrt{\sum_{i,j} a_{ij}^2} $。  
- $ \kappa_2(A) = \frac{\sigma_{\max}(A)}{\sigma_{\min}(A)} $。  