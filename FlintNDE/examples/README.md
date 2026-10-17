# Examples

仓库内示例自动发现程序包根：发布布局直接取 `examples/` 上级目录，发展仓布局取
`versions/` 下唯一的当前版本目录；同一份示例文件在两种布局下字节一致地运行，不复制
程序包源码，也不读取项目外其它代码包。发布用户应先安装 `FlintNDE`，随后在任意调用
目录直接 `import flintnde`；路径段只服务于尚未安装的仓库内示例。

当前保留八个典型示例：

- `qnm_2x2.py`：统一 `u` 的 exact 2x2 一阶系统。分别从字面量 horizon `1` 和
  infinity `"inf"` 的 `{a,b,C}` 边界出发，输运到公共匹配点并检查另一端的禁戒分量；
  同时运行一个轻微偏移频率作为反事实对照。第三块从 horizon 单端输出 QNM 解全曲线：
  tortoise 坐标 `r* = x+ln(x-1)` 均匀网格 `[-6,160]` 共 721 点（含精确 `x=4`）一次规划
  执行，1e-10 相对精度；271 个规划节点覆盖全部请求点，585 点共享段 Taylor 系数（最大桶
  32 点），双链端点差 1.5e-17，与独立规划 horizon 链在 `x=4` 一致到 5.2e-32；样本曲线
  写入 `transport/qnm_2x2_curve_samples.json`。
- `regular_singular_save.py`：标量正则奇点，从 `{a,b,C}` 起点输运到普通点，并保存可复用奇点边界和普通终点。
- `exponential_boundary_save.py`：严格解耦的指数型奇点，从 `{phi,a,b,C}` 起点输运，保存并重新读入同一边界。
- `exponential_polynomial_boundary.py`：三阶极点（本性奇点）的双 sector 对角系统，两条
  `{phi,a,b,C}` 记录各含两项负幂的指数多项式相位；闭式解在终点为 `(e^{1/4}, e^{-11/6})`，
  与输运结果做相对差检查。
- `ep_parallel.py`：Python 入口对不同固定 `ep` 做有界多进程 NDE；`parallel_task_count`
  缺省 12，实际并行数自动取任务数与该值的较小者。
- `ep_series_reconstruction.py`：从上游符号证书给出的 `leading_power=-1` 出发，自适应重构
  `1/ep` 与有限项；演示 7 个显式生产候选点中首轮按内部最高幂只取 4 个、验证失败后复用
  旧值并新增 2 个、最终留下 1 个未使用候选点，以及缺省 12 个外层进程。显式验证点始终不
  进入拟合，FlintNDE 也不从数值样本猜最低阶；同一脚本还用开角域 `(-pi/3,pi/3)` 演示
  最多三条内部复射线，模长仍由目标精度自动决定。
- `ep_parallel_mathematica.wl`：Wolfram 入口使用
  `ParallelTaskCount -> 12`（缺省）运行相同的固定 `ep` 并行合同。
- `mathematica_interface_example.wl`：Wolfram 入口从 exact 有理矩阵构造系统，显式规划并执行
  普通点路径，与闭式解比较，并把正式结果写入 example 自身的 `results/`。

运行：

```powershell
python qnm_2x2.py
# 或显式选择 convention
python qnm_2x2.py --config ../config/qnm_u_unified_it0_3_it1_minus1.json
python regular_singular_save.py
python exponential_boundary_save.py
python exponential_polynomial_boundary.py
python ep_parallel.py
python ep_series_reconstruction.py
wolframscript -file ep_parallel_mathematica.wl
```

两个 Wolfram 示例都把 `"WorkDirectory"` 指向 example 目录下唯一的 `results_temp/`；程序内部
只追加短的 `bridge/` 和 16 位 token 文件名。用户若显式改目录，所给路径同样表示运行根本身。

正式摘要写入 `results/qnm_2x2/summary/qnm_2x2_summary.json`；初始化布局说明写入
`results/qnm_2x2/configuration/output_layout.json`。两条路径均由调用脚本传入
`__file__` 后通过 package 公共接口生成，不依赖启动时的工作目录。
`it0/it1` 只从实际 config 读取并写回 summary。无穷远自动使用 start-only
`formal_exponential_asymptotic`；最近指数根差和缺省四倍最小项阶数规则自动生成首匹配点，
结果固定保留指定 N 阶并记录后五阶/前五阶矢量块比、五阶相对 refinement 和辅助最小项信息。
块比不小于 1 时 warning 但仍保存结果；自动点和附近点另与独立二阶标量递推交叉验证。

后两个脚本抽自当前 `tests/test_save_points.py` 已执行的公开调用配置。它们把输出写到本 example 调用目录下的 `results/`；发布前从本目录 fresh 执行并记录结果。
