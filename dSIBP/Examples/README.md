# dSIBP examples

本目录只保存可直接阅读和长期维护的 Wolfram 源文件，以及 Kira 成品案例的轻量输入/结果摘要。外部 reduction 不在本目录运行，也不把 `init/`、`kira/`、database、save、日志、reduction table、cache 或 DE 运行目录提交到仓库。

- `01_mixed_bubble_workflow.wl` 只在内存中构造预约化 plan，不导出或运行 Kira。
- `03_single_massive_sunrise/` 只生成 general seeds 和参数微分算符，不进入 Kira、DE 或 scaling。
- `04_pure_massive_bubble_closed_loop/` 给出完整 export/import 脚本；运行时必须通过 `DSIBP_KIRA_WORKSPACE` 指定仓库外工作区。目录内的 `kira_input_summary.wl` 和 `kira_result_summary.wl` 是 fresh WSL Kira 检验的轻量记录。
- `06_mix_bubble_tree/` 是 mixed bubble+tree 正式例子；其 Kira 输入规模、parity、精确点和最终结果也分别保存在两份轻量摘要中。

两套正式 reduction 均由 Windows Wolfram 前端在仓库外生成输入，再在 WSL 中运行 Kira 2.3。发布前已删除完整中间工作树；摘要不作为程序包输入，也不能代替 `000-report/` 中的独立检验报告。
