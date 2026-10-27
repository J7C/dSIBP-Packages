# Single-massive sunrise

运行 `wolframscript -file main.wl`。本例固定为两个顶点、三条平行传播子和两个圈动量：第一条传播子是 massive h，另外两条是 massless exponential，并显式给出两个 ISP。

脚本把两个顶点的外腿能量统一为 `kE`，并保留独立圈外模长 `kL`，因此是两标度输入。两个 ISP 选成在两条 massless 平行传播子交换下互换的一对；`symmetryRules` 同时实现顶点交换和 massless 传播子/ISP 交换。

脚本只执行 `DSKinematics`、`DSInit`、`DSSeeds/DSAllSeeds`，构造全部 contact-reachable sectors 的 general IBP 模板；连续指标始终保持符号。它同时读取初始化产生的 `{ss11,kE}` general 参数微分算符，并用一个 general top integral 展示两个公开 `ds` 结果。`ss11` 通过 `2 ss11 d/d sp[kL,kL]` 链式法则实现，`kE` 同时作用于两个顶点的共同能量。

本例不调用 `DSMetaSeedRange`、`DSGenerateIBP`、`DSLinear`、Kira、DE 或 scaling，也不选择数值点、target 或 master，不写任何运行产物。

本例不登记为第三套待运行的 Kira reduction。外部 reduction 当前只包含 pure massive bubble 与 mixed bubble+tree；若将来单独扩展 sunrise 数值约化，必须先另行修改独立检验任务书和版本边界。

本目录是 package 中唯一的 sunrise example；其它 examples 不应再复制或改名形成第二套 sunrise 输入。
