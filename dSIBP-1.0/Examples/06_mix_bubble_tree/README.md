# Mix bubble+tree

运行 `wolframscript -file main.wl`。本例固定三条传播子：`v1-v2` 间是一条 massive h cycle 传播子和一条 massless exponential cycle 传播子，`v2-v3` 间是一条 massless exponential bridge 传播子。

两条 cycle propagator 分别取全偶，即对各自完整的 `{b,n1,n2}` 指标施加偶性约束。第三条传播子是 fixed tree/bridge，不属于圈积分，没有连续 `b` 指标，也不施加 parity。

`loopExternalMomenta={k1+k2}` 定义 loop 外方向，`independentExternalMomenta={k1,k2}` 定义两个无圈模长；脚本展示缺省 `ss11/sE1/sE2`、自定义 `loopScale/legScale1/legScale2`、顶点独立相位，以及 cycle/fixed 传播子 pack 和 pinched/pinch 后的 sector 数据。

本例不增加第二条 massive 传播子，避免把 compound topology 的角色、massless convention 和 bridge/cycle contraction 演示混入不必要的 massive function-system 组合。

本例保留已认证 Kira 检验的轻量记录：`kira_input_summary.wl` 给出 `+++` 分支、两条 cycle-even parity、最终 seed 范围 `a1,a2∈[0,5]`、`a3∈[-1,1]`、`b1,b2∈[-4,7]`、精确规则与输入规模（205434 条方程、85306 个积分）；`kira_result_summary.wl` 给出 903 targets、84 个稳定候选主积分及其 sector/总幂次分布；`master_integrals.wl` 保存机器可读 J 清单。target 扩到 a₃=0/1 的补测保持同一 84 个。本例不运行 DE 或 scaling。完整边界台账见版本外 `example-test-ledger/029_example6.md`。

完整 Kira 输入由 Windows Wolfram 在仓库外验证工作区生成，reduction 在 WSL 中运行，并未在本 example 目录执行。发布前已删除 `init/`、`kira/`、database、save、日志、reduction table、cache 和 DE 运行目录。两份摘要只用于把本例输入与正式结果对应起来，不作为程序包或 importer 输入；实际证明保存在独立检验报告中。
