// FlintNDE 引擎公开接口：请求校验、prefix 认证输运与固定 ep 批量调度。
#ifndef FLINTNDE_ENGINE_HPP
#define FLINTNDE_ENGINE_HPP

#include <string>

#include "json.hpp"

namespace fnde {

// 极点留数乘法内核选择：auto 按留数矩阵稀疏度自适应。
enum class ResidueMode { Auto, Sparse, Dense };

struct EngineOptions {
    long flint_threads = 1;      // FLINT 内核线程池上限（对标 ctx.threads）
    long worker_override = 0;    // 批量 worker 数覆盖；0 表示用请求 JSON 里的值
    ResidueMode residue_mode = ResidueMode::Auto;
};

ResidueMode parse_residue_mode(const std::string& text);
std::string residue_mode_text(ResidueMode mode);

// 处理单个 evaluate 请求（schema 与 MadStree↔FlintNDE 桥一致）。
fjson::Value run_request(const fjson::Value& request, const EngineOptions& options);

// 处理固定 ep 批量请求（madstree_flintnde_ep_batch_v1）；
// raw 是完整批请求文本，任务列表以流式扫描分片，worker 各自解析，避免整批 DOM。
fjson::Value run_batch(const std::string& raw, const EngineOptions& options);

}  // namespace fnde

#endif  // FLINTNDE_ENGINE_HPP
