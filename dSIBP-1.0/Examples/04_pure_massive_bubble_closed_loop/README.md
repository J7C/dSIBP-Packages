# Pure massive bubble closed loop

本例固定 `--` branch、even parity 和等能量约束。package 与 reference 的能量满足 `P_pkg=P0=-P_ref`。两个顶点的零点均加入通用符号 `analyticRegulator`，在符号 seed 生成后才取 `analyticRegulator -> 0`；这使时间端点门禁能够区分已正规化 family 与未定义的发散 IBP。其余流程从初始化到 formal Kira 输入、外部 reduction、结果取回、19 维微分方程与 Eq. (51)/(64) 标度检查。顶点交换 symmetry 只在本例等能量约束下启用；独立 `P1/P2` family 不得复用。

`reference_user_mi_basis.wl` 只保存 reference 的 19 个候选线性组合（全部 active，没有辅助关系）、`activeIndices=Range[19]`、physical `ks=ss11`、导数变量和 19 个统一 scaling degrees。`main.wl` 把这些数据直接交给 package 的 `DSUserMI`；线性秩、可逆 `J/userMI` 映射、backend ID、导数闭包和 manifest 均由 package 生成，example 不实现自己的 basis adapter。`dlog_basis.wl` 只保留 reference-readable 记号，不应单独作为 formal basis。

当前 target 包络为 `a[v1],a[v2] in [0,3]`、`b[1],b[2] in [-3,5]`，直接覆盖 19 项 basis 及导数支撑；生成时使用 `SeedEnvelopeCut->Automatic` 与 `TotalPowerLimit->None`。fresh 结果为 20228 条生成 IBP、2442 条 canonical relations、2460 条 formal equations、1317 个源积分、300 targets。Kira 的 24 个全局 masters 中有 19 个 active masters 和 5 个额外全局项；验收以目标 RHS 为准，19 维 DE 无残留，post-Kira 23/23、三套 reference 矩阵各 361/361。复现旧未截断输入时须显式关闭两个截断选项。

本仓库没有在该 example 目录内运行 Kira。正式检验由 Windows Wolfram 在仓库外工作区生成输入，再由 WSL Kira 2.3 运行。发布副本已删除 `init/`、`kira/`、database、save、日志、reduction table、cache 和 DE 运行目录，只保留 `kira_input_summary.wl` 与 `kira_result_summary.wl` 两份轻量摘要。

1. 在 PowerShell 中把 `DSIBP_KIRA_WORKSPACE` 设为仓库外目录，再运行 `wolframscript -file main.wl`。脚本会在其 `pure_massive_bubble/` 子目录写出 `init/` 和 `kira/`，不会启动 Kira；首次正常状态为 `awaitingExternalKira`。
2. 从 WSL 进入同一外置工作区的 `pure_massive_bubble/kira/`，运行 `kira --parallel=10 jobs.yaml`。完整结果写入外置工作区，不复制回 example 目录。
3. 再运行 `main.wl`，或只运行 `post_kira_check.wl`。后者不重新生成输入，只执行 `DSKiraImport -> DSDE[{ss11,P0}] -> DSScaleCheck[{1,1}]`，结果仍写入外置工作区。

最终 DE 与同序 master list 写入外置工作区的 `results/dlogDE/`，闭环摘要写入其 `results/post_kira_summary.wl`。若修改 family、active basis、package 或 Kira 输入，必须重新生成并运行 reduction；脚本会拒绝读取早于当前 export manifest 的旧结果。仓库内两份摘要只记录正式检验的输入规模和最终结论，不作为 importer 输入。
