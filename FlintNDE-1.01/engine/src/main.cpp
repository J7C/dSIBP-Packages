// FlintNDE C++ 引擎命令行入口：对齐 flintnde_transport.py 的 INPUT OUTPUT 文件协议。
// 未实现的功能（如 ep series control）fail-closed 返回结构化拒绝，不静默降级。
// 批量请求按 schema 流式分片，不整批建 DOM。
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "engine.hpp"
#include "flint/flint.h"

namespace {

const char* kUsage =
    "FlintNDE engine: C++/FLINT twin kernel of FlintNDE (MadStree bridge contract)\n"
    "usage: flintnde_engine INPUT.json OUTPUT.json [--workers N]\n"
    "                  [--flint-threads N] [--residues auto|sparse|dense]\n"
    "\n"
    "  --workers N         批量请求 worker 数覆盖（默认用请求 JSON 里的 parallelTaskCount）\n"
    "  --flint-threads N   FLINT 内核线程池上限（默认 1，对标 python ctx.threads=1）\n"
    "  --residues          留数乘法内核（默认 auto：按留数非零元密度自适应选\n"
    "                      sparse 稀疏跳乘或 dense 稠密 acb_mat_mul）\n";

std::string read_file(const std::string& path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot open input file: " + path);
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    return buffer.str();
}

void write_file(const std::string& path, const std::string& text) {
    std::ofstream stream(path, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("cannot open output file: " + path);
    stream << text;
    stream.flush();
    if (!stream) throw std::runtime_error("cannot write output file: " + path);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc == 1) {
        std::fputs(kUsage, stdout);
        return 0;
    }
    fnde::EngineOptions options;
    bool have_input = false;
    bool have_output = false;
    std::string input_path;
    std::string output_path;
    for (int index = 1; index < argc; ++index) {
        const std::string flag = argv[index];
        auto next_value = [&]() -> std::string {
            if (index + 1 >= argc)
                throw std::runtime_error("missing value for flag: " + flag);
            return argv[++index];
        };
        if (flag == "--residues") {
            options.residue_mode = fnde::parse_residue_mode(next_value());
        } else if (flag == "--workers") {
            options.worker_override = std::stol(next_value());
            if (options.worker_override < 1)
                throw std::runtime_error("--workers must be a positive integer");
        } else if (flag == "--flint-threads") {
            options.flint_threads = std::stol(next_value());
            if (options.flint_threads < 1)
                throw std::runtime_error("--flint-threads must be a positive integer");
        } else if (!flag.empty() && flag[0] == '-') {
            throw std::runtime_error("unknown flag: " + flag);
        } else if (!have_input) {
            input_path = flag;
            have_input = true;
        } else if (!have_output) {
            output_path = flag;
            have_output = true;
        } else {
            throw std::runtime_error("unexpected extra argument: " + flag);
        }
    }
    if (!have_input || !have_output) {
        std::fputs(
            "usage: flintnde_engine INPUT.json OUTPUT.json [options]\n", stderr);
        return 1;
    }
    // FLINT 内核线程池：默认 1（python ctx.threads 同缺省），可显式提高。
    flint_set_num_threads(static_cast<long>(options.flint_threads));
    fjson::Value result;
    int exit_code = 0;
    try {
        const std::string raw = read_file(input_path);
        const std::string schema = fjson::scan_schema_field(raw);
        if (schema.empty()) {
            throw std::runtime_error("request must be a JSON object with a schema string");
        }
        if (schema == "madstree_flintnde_ep_batch_v1") {
            // 批量：只流式扫 tasks 分片，各 worker 解析自己的任务文本。
            result = fnde::run_batch(raw, options);
        } else if (schema == "madstree_flintnde_evaluate_v1") {
            result = fnde::run_request(fjson::parse(raw), options);
        } else if (schema == "madstree_flintnde_ep_series_control_v1") {
            // 结构化拒绝：引擎只复刻数值输运合同。
            fjson::Value refusal;
            refusal.type = fjson::Value::Type::ObjectT;
            refusal.fields["status"] = fjson::Value(std::string("seriesControlNotImplemented"));
            refusal.fields["reason"] = fjson::Value(std::string(
                "the FlintNDE engine covers the numeric transport contract only; "
                "ep series control requests are refused"));
            result = refusal;
        } else {
            throw std::runtime_error("unsupported request schema: " + schema);
        }
    } catch (const std::exception& error) {  // CLI 边界必须序列化所有错误
        fjson::Value failed;
        failed.type = fjson::Value::Type::ObjectT;
        failed.fields["status"] = fjson::Value(std::string("failed"));
        failed.fields["errorType"] = fjson::Value(std::string(typeid(error).name()));
        failed.fields["error"] = fjson::Value(std::string(error.what()));
        result = failed;
        exit_code = 1;
        std::fprintf(stderr, "flintnde_engine: %s: %s\n", typeid(error).name(), error.what());
    }
    write_file(output_path, fjson::dump(result));
    flint_cleanup();
    return exit_code;
}
