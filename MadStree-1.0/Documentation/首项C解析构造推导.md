# 首项 C 解析构造推导（indicial 幂次与商公式）

本文件给出 `MSBoundaryData` 奇异边界首项向量 C 的显式构造推导，
并把每一步对应到代码位置。对象是 $k_0\to\infty$ 边界曲线拉回后的 dlog DE 在 $t=0$
处的 Frobenius 首项数据；手册正文见 `用户手册.pdf` 的 generic-sector 边界一节。

## 1. 边界曲线与留数矩阵

`msGenericBoundaryCurveData` 按所选 strict chart 构造分级 Laurent 曲线：

$$E_v(t)=-i\sigma_v B^{w_v}t^{-w_v},\qquad w_v=V-\operatorname{rank}(v)+1.$$

内部动量固定，合并 component 的能量是根顶点能量之和，其领头幂由其中最大权重给出。把完整 dlog connection
沿曲线拉回为单变量矩阵 $A(t)$，并取 $t\to 0$ 留数

$$R=\operatorname{Res}_{t=0}A(t)=R_0 .$$

以三顶点 ppp 生产几何为例（`MSInitTree`/`MSDLogDE` 后，首点参数），$R$ 的结构为：
对角元 $\alpha_1,\dots,\alpha_d$（$d$ 为主积分总数）＋ 严格块上三角的小整数 contact 项
（$\pm1,\pm2$），sector 对角块纯对角，块非零元布局与 contact DAG 严格一致。

## 2. indicial 幂次的来源（Hankel/deformed mode 展开）

对 generic Hankel 阶数，乘上 deformed mode 的前因子后，两条小宗量幂序列为

$$h(\nu,0;z)=z^{s_\nu\nu}H_\nu(z)=A_0\bigl(1+O(z^2)\bigr)+B_0z^{2\,\texttt{formulaNu}}\bigl(1+O(z^2)\bigr),\qquad \texttt{formulaNu}=s_\nu\nu.$$

指数从 0 和 $2\,\texttt{formulaNu}$ 开始；不能把未乘前因子的 Hankel 幂次直接用于 h。导数态满足

$$h(\nu,1;z)=\frac{d}{dz}h(\nu,0;z)=O(z)+2\,\texttt{formulaNu}\,B_0z^{2\,\texttt{formulaNu}-1}\bigl(1+O(z^2)\bigr).$$

第一条序列的常数项求导为零，但其高阶项仍产生 O(z)。边界记录中 bit=0 的锚点取常数序列，bit=1 的锚点取指数为 $2\,\texttt{formulaNu}-1$ 的序列；这不表示删除另一条序列，也不是按模长选择唯一占优项。退化参数须取相应极限并检查 log，不直接代入 generic 系数的零分母。

对 sector $s$、state bits $\bm a$、component $C$，所选锚点的时间幂为

$$q_C=p_C+\sum_{\sigma\in H(s,C)}a_\sigma\bigl(2\,\texttt{formulaNu}_\sigma-1\bigr).$$

在收敛阻尼域，

$$\int_0^\infty u^{q_C}e^{-Ku}\,du=\Gamma(q_C+1)K^{-q_C-1},\qquad \Re q_C>-1,\quad \Re K>0.$$

沿既定分支解析延拓后，能量幂仍为 $-q_C-1$；轮廓相位、动量幂及 Hankel 支系数进入物理权重。因此

$$\alpha_C=q_C+1=p_C+1-\sum_{\sigma\in H(s,C)}a_\sigma\bigl(1-2\,\texttt{formulaNu}_\sigma\bigr).$$

其中 $p_C$ 是该 component 的 base time power，$a_\sigma\in\{0,1\}$ 是 slot 的 state bit。
对应代码：`Kernel/Numerics/Boundary.wl` 的 `msGenericSectorLeadingRecord`，
其中 `frobeniusPower` 逐 component 计算。

master $(s,\bm a)$ 的 Frobenius 指数是各 component 幂按 rank weight 的加权求和：

$$\lambda_{s,\bm a}=\sum_C w_C\,\alpha_C,\qquad w_C=V-\operatorname{rank}_C+1,$$

同函数内 `exponent = Total[rankWeight * frobeniusPower]`。物理权重
$W_{s,\bm a}$（normalization×shared 因子×component 系数×anchor 因子）也在同一记录里
解析给出，不经过任何线性求解；本文件只处理首项向量 $C$。

**对角恒等式**：按 globalIndex 次序把各 master 的 $\lambda$ 排成向量 $\bm\lambda$，
则 $\texttt{SameQ}[\mathrm{Diagonal}[R],\bm\lambda]$ 在生产几何上精确成立
（三顶点 ppp 首点生产请求实测 True）。
这正是 $R$ 对角元本身就是各 master 的 Hankel 指数这一事实的矩阵表述。

## 3. 首项系统

分支 $(s,\bm a)$（root 指标 $r$）的 Frobenius ansatz $Y(t)=t^{\lambda_r}(C+O(t))$
代入拉回 DE 的 $t^{-1}$ 层给出首项方程

$$(\lambda_r I-R)\,C=0 .$$

生产语义固定归一化 $C[r]=1$，并把非祖先闭包行置零：

$$C[j]=0\quad\text{当 }j\notin \overline{\mathcal A}(s)\;(\text{sector } s\text{ 的祖先闭包，含自身}).$$

高 sector 中由低 sector source 诱导的分量恰好来自 $R$ 的 contact 块：
非对角块位于 supersector 行、subsector 列。解的源项从低 sector 向高 sector 传播，
与源码中按 supersector→subsector 登记的 contact transition 方向相反。
quotient 的 embedding/projection 必须保留，不能将 raw 状态向量直接当作 master 向量。

## 4. 结构定理与商公式

### 4.1 前提（结构证书，整请求判一次）

代码 `msBoundaryAnalyticStructureCertificate`（`Kernel/Numerics/Boundary.wl:810`）：

1. **对角恒等式**：$\mathrm{Diagonal}[R]=\bm\lambda$（逐项 `SameQ`）；
2. **sector 对角块纯对角**：$R$ 在每个 sector 自身行列块内除对角元外为零；
3. **列支撑**：sector $T$ 各列的非零元只出现在 $T$ 自身的对角位与 $T$ 的祖先闭包行块内；
4. **拓扑次序**：每条 contact transition 在 `sectorOrder` 中由前 sector（super）指向后 sector（sub）。

### 4.2 商公式与证明

设 root 分支 $(s,r)$，$\alpha=\lambda_r$。断言：在前提 1–4 且闭包内无第二行与 root 同指数时，

$$C[r]=1,\qquad C[j]=0\;\;(j\notin\overline{\mathcal A}(s)\ \text{或}\ j\in s,\,j\ne r),$$
$$C[J]=\frac{R[J,P]\,C[P]}{\alpha-\bm\lambda[J]}\quad\text{（只遍历严格祖先，按 Reverse[sectorOrder] 自下而上），}$$

其中 $P$ 是已确定行块的列集合（初值 $\{r\}$）。闭包外可以存在同根，其分量保持零，不进入除法。
这是一条显式块三角回代，省去通用受限线性求解器，但仍完整使用 contact 数据。

**证明**：按行写出首项方程，$(\alpha-\lambda_j)C_j=\sum_{k\ne j}R[j,k]C_k$。

- root 行及同 sector 兄弟行：由前提 2、3，非对角项只来自更低 subsector 列，
  这些列在本分支中全为零。root 行两侧均为零；兄弟行因 $\lambda_j\ne\alpha$ 得到 $C_j=0$。
- 祖先闭包外的行：若右端存在非零 $C_k$，则 $k$ 属于 root 的祖先闭包；若同时
  $R[j,k]\ne0$，前提 3 又要求 $j$ 是 $k$ 的祖先。由祖先关系传递性，$j$ 也在 root 的
  祖先闭包内，与假设矛盾。因此右端为零，$C_j=0$ 满足方程，即使 $\lambda_j=\alpha$ 也无需除法。
- 严格祖先行：按从低到高的顺序，所有可能非零的 subsector 分量都已确定。
  同 sector 非对角项为零，右端完全已知；$\alpha-\lambda_j\ne0$ 唯一给出商公式。

归纳终了全部行方程成立。在相同的支撑、root 归一化与非退化条件下，构造与受限解一致。$\square$

回代还可展开为有限路径和。把每个非零非对角元 $R[j,k]$ 看作由列 $k$ 向行 $j$ 的源项连接，则

$$C_j=\sum_{r=i_0\to i_1\to\cdots\to i_m=j}\prod_{h=1}^{m}\frac{R[i_h,i_{h-1}]}{\lambda_r-\lambda_{i_h}}.$$

root 对应空路径、乘积为 1；闭包外没有路径、和为零。块三角性保证路径有限。
实现采用逐块回代复用公共子路径，无需枚举所有路径。

### 4.3 同根多方向／退化／log 情形的 fail-closed 触发

商公式分母为零或本征空间退化时不给构造、直接回退 rooted 受限求解（`msBoundaryRootedIndicialVector`，
`Kernel/Numerics/Boundary.wl:571`）。触发条件（`msBoundaryAnalyticLeadingVector`，`:865`）：

| 触发 | 标签 | 含义 |
| --- | --- | --- |
| 祖先行同根 | `AnalyticSameRootAncestor` | 闭包内严格祖先的 $\lambda_j=\lambda_r$，商公式分母为零 |
| 同 sector 兄弟同根 | `AnalyticSameRootSector` | 本征空间退化，rooted 受限解不唯一，两路线不必一致 |
| 对角恒等式破坏 | `AnalyticDiagonalIdentityFailed` | 前提 1 失败 |
| 列支撑破坏 | `AnalyticColumnSupportViolation` | 前提 3 失败 |
| 次序非拓扑 | `AnalyticSectorOrderNotTopological` | 前提 4 失败 |

log/Jordan 结构不在本层判定：纯幂首项不保证后续阶数无 log，正整数根差处仍须检查相容性；
log 首现阶由 power-log 证书与后端递推处理。回退求解若不满足边界条件也须报告失败，
不能把回退解释为任意退化输入都可求解。

## 5. 与现行实现的对应

| 环节 | 代码位置 |
| --- | --- |
| 边界曲线、拉回、$R$ | `msGenericBoundaryCurveData` |
| component 幂 $\alpha_C$、master 指数 $\lambda$、物理权重 | `msGenericSectorLeadingRecord` |
| 结构证书（前提 1–4） | `msBoundaryAnalyticStructureCertificate` |
| 商公式构造 + 同根 fail-closed | `msBoundaryAnalyticLeadingVector` |
| 回退受限求解 | `msBoundaryRootedIndicialVector` |
| 祖先闭包 | `msBoundaryAncestorSectorKeys`，含自身 |
| 调用点与路线记录 | `msGenericSectorFrobeniusData`，`leadingVectorRoute` 记 `analyticQuotient`/`rootedIndicialSolve` |

证书在 `msGenericSectorFrobeniusData` 内整请求判一次；成立时逐分支走商公式，
失败（结构级）或触发（分支级）即回退，序列化分支 `{a,b,C}` 合同与归一化
$C[r]=1$ 两路线一致。

## 6. 验证记录

同请求（三顶点 ppp 首点、25 master、80 位）：

- 25/25 分支 `analyticQuotient` 路线，`analyticEqualsRooted` 精确 True、
  Together 残差精确为 0；解析与 rooted 两路线的序列化载荷逐字符一致（201,591 字符）。
- 四类 fail-closed 触发各自返回正确标签，回退 rooted 结果与未触发解析构造一致，
  无关分支不受影响。
