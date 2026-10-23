// FlintNDE 引擎极简 JSON：解析 + 输出 + 批量流式扫描，仅覆盖桥接 schema 所需特性。
// 输入为 ASCII 为主的 UTF-8 文本；数字按 double 存储（数值载荷全部走字符串）。
#ifndef FLINTNDE_ENGINE_JSON_HPP
#define FLINTNDE_ENGINE_JSON_HPP

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fjson {

class Value;
using Array = std::vector<Value>;
using Object = std::map<std::string, Value>;

class Value {
public:
    enum class Type { Null, Bool, Number, String, ArrayT, ObjectT };

    Type type = Type::Null;
    bool boolean = false;
    double number = 0.0;
    std::string text;
    Array items;
    Object fields;

    Value() = default;
    explicit Value(bool b) : type(Type::Bool), boolean(b) {}
    explicit Value(double x) : type(Type::Number), number(x) {}
    explicit Value(std::string s) : type(Type::String), text(std::move(s)) {}

    bool is_null() const { return type == Type::Null; }
    bool is_bool() const { return type == Type::Bool; }
    bool is_number() const { return type == Type::Number; }
    bool is_string() const { return type == Type::String; }
    bool is_array() const { return type == Type::ArrayT; }
    bool is_object() const { return type == Type::ObjectT; }

    const Value* find(const std::string& key) const {
        if (type != Type::ObjectT) return nullptr;
        auto it = fields.find(key);
        return it == fields.end() ? nullptr : &it->second;
    }
    const Value& at(const std::string& key) const {
        const Value* item = find(key);
        if (item == nullptr) throw std::runtime_error("missing JSON field: " + key);
        return *item;
    }
    long long as_integer(const char* label) const {
        if (type != Type::Number || !std::isfinite(number) ||
            number != std::floor(number)) {
            throw std::runtime_error(std::string("expected integer for ") + label);
        }
        return static_cast<long long>(number);
    }
};

struct TopLevelScan;

// 基础解析器：作用域为 [data, data+size) 的只读文本，不持有内存。
// position 公开给流式扫描包装使用；skip_value 为不建 DOM 的快速跳过。
class Parser {
public:
    size_t position = 0;

    Parser(const char* data, size_t size) : data(data), size(size) {}
    explicit Parser(const std::string& source) : data(source.data()), size(source.size()) {}

    Value parse() {
        skip();
        Value result = parse_value(0);
        skip();
        if (position != size) throw std::runtime_error("trailing JSON content");
        return result;
    }

    // 顶层对象流式扫描：stream_field 数组的元素只记录 span，其余字段建 DOM。
    void parse_streaming_object(const std::string& stream_field, TopLevelScan& scan);

    // 只读顶层 "schema" 字符串字段，不解析其余内容（分发用）。
    std::string scan_schema();

private:
    const char* data;
    size_t size;

    void skip() {
        while (position < size) {
            char c = data[position];
            if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
                ++position;
            } else {
                break;
            }
        }
    }

    char peek() {
        if (position >= size) throw std::runtime_error("unexpected JSON end");
        return data[position];
    }

    void expect(char c) {
        if (position >= size || data[position] != c) {
            throw std::runtime_error(std::string("expected '") + c + "' in JSON");
        }
        ++position;
    }

    bool literal(const char* word, size_t length) {
        return position + length <= size && std::memcmp(data + position, word, length) == 0;
    }

    // 无 DOM 快速跳过一个完整 JSON 值（流式扫描数组元素用）。
    void skip_value(int depth) {
        if (depth > 128) throw std::runtime_error("JSON nesting too deep");
        skip();
        char c = peek();
        if (c == '{') {
            ++position;
            skip();
            if (peek() == '}') {
                ++position;
                return;
            }
            while (true) {
                skip();
                skip_string_raw();
                skip();
                if (peek() != ':') throw std::runtime_error("expected ':' in JSON object");
                ++position;
                skip();
                skip_value(depth + 1);
                skip();
                char t = peek();
                if (t == ',') {
                    ++position;
                    continue;
                }
                if (t == '}') {
                    ++position;
                    return;
                }
                throw std::runtime_error("expected ',' or '}' in JSON object");
            }
        }
        if (c == '[') {
            ++position;
            skip();
            if (peek() == ']') {
                ++position;
                return;
            }
            while (true) {
                skip();
                skip_value(depth + 1);
                skip();
                char t = peek();
                if (t == ',') {
                    ++position;
                    continue;
                }
                if (t == ']') {
                    ++position;
                    return;
                }
                throw std::runtime_error("expected ',' or ']' in JSON array");
            }
        }
        if (c == '"') {
            skip_string_raw();
            return;
        }
        while (position < size) {
            char t = data[position];
            if (t == ',' || t == '}' || t == ']' || t == ' ' || t == '\t' ||
                t == '\r' || t == '\n') {
                break;
            }
            ++position;
        }
    }

    void skip_string_raw() {
        if (position >= size || data[position] != '"') {
            throw std::runtime_error("expected JSON string");
        }
        ++position;
        while (true) {
            if (position >= size) throw std::runtime_error("unterminated JSON string");
            char c = data[position];
            if (c == '\\') {
                position += 2;
                continue;
            }
            if (c == '"') {
                ++position;
                return;
            }
            ++position;
        }
    }

    Value parse_value(int depth) {
        if (depth > 128) throw std::runtime_error("JSON nesting too deep");
        char c = peek();
        if (c == '{') return parse_object(depth);
        if (c == '[') return parse_array(depth);
        if (c == '"') {
            Value v;
            v.type = Value::Type::String;
            v.text = parse_string();
            return v;
        }
        if (c == 't' || c == 'f') return parse_bool();
        if (c == 'n') {
            if (!literal("null", 4)) throw std::runtime_error("bad JSON literal");
            position += 4;
            Value v;
            v.type = Value::Type::Null;
            return v;
        }
        return parse_number();
    }

    Value parse_object(int depth) {
        expect('{');
        Value v;
        v.type = Value::Type::ObjectT;
        skip();
        if (peek() == '}') {
            ++position;
            return v;
        }
        while (true) {
            skip();
            std::string key = parse_string();
            skip();
            expect(':');
            skip();
            v.fields[key] = parse_value(depth + 1);
            skip();
            char c = peek();
            if (c == ',') {
                ++position;
                continue;
            }
            if (c == '}') {
                ++position;
                break;
            }
            throw std::runtime_error("expected ',' or '}' in JSON object");
        }
        return v;
    }

    Value parse_array(int depth) {
        expect('[');
        Value v;
        v.type = Value::Type::ArrayT;
        skip();
        if (peek() == ']') {
            ++position;
            return v;
        }
        while (true) {
            skip();
            v.items.push_back(parse_value(depth + 1));
            skip();
            char c = peek();
            if (c == ',') {
                ++position;
                continue;
            }
            if (c == ']') {
                ++position;
                break;
            }
            throw std::runtime_error("expected ',' or ']' in JSON array");
        }
        return v;
    }

    std::string parse_string() {
        expect('"');
        std::string out;
        while (true) {
            if (position >= size) throw std::runtime_error("unterminated JSON string");
            unsigned char c = static_cast<unsigned char>(data[position]);
            if (c == '"') {
                ++position;
                break;
            }
            if (c == '\\') {
                ++position;
                if (position >= size) throw std::runtime_error("bad JSON escape");
                char e = data[position++];
                switch (e) {
                    case '"': out.push_back('"'); break;
                    case '\\': out.push_back('\\'); break;
                    case '/': out.push_back('/'); break;
                    case 'b': out.push_back('\b'); break;
                    case 'f': out.push_back('\f'); break;
                    case 'n': out.push_back('\n'); break;
                    case 'r': out.push_back('\r'); break;
                    case 't': out.push_back('\t'); break;
                    case 'u': {
                        unsigned cp = parse_hex4();
                        if (cp >= 0xD800 && cp <= 0xDBFF && position + 1 < size &&
                            data[position] == '\\' && data[position + 1] == 'u') {
                            position += 2;
                            unsigned low = parse_hex4();
                            if (low >= 0xDC00 && low <= 0xDFFF) {
                                cp = 0x10000 + ((cp - 0xD800) << 10) + (low - 0xDC00);
                            }
                        }
                        append_utf8(out, cp);
                        break;
                    }
                    default: throw std::runtime_error("bad JSON escape character");
                }
            } else if (c < 0x20) {
                throw std::runtime_error("control character in JSON string");
            } else {
                out.push_back(static_cast<char>(c));
                ++position;
            }
        }
        return out;
    }

    unsigned parse_hex4() {
        unsigned value = 0;
        for (int i = 0; i < 4; ++i) {
            if (position >= size) throw std::runtime_error("short \\u escape");
            char c = data[position++];
            value <<= 4;
            if (c >= '0' && c <= '9') value |= static_cast<unsigned>(c - '0');
            else if (c >= 'a' && c <= 'f') value |= static_cast<unsigned>(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') value |= static_cast<unsigned>(c - 'A' + 10);
            else throw std::runtime_error("bad \\u escape hex digit");
        }
        return value;
    }

    static void append_utf8(std::string& out, unsigned cp) {
        if (cp < 0x80) {
            out.push_back(static_cast<char>(cp));
        } else if (cp < 0x800) {
            out.push_back(static_cast<char>(0xC0 | (cp >> 6)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        } else if (cp < 0x10000) {
            out.push_back(static_cast<char>(0xE0 | (cp >> 12)));
            out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        } else {
            out.push_back(static_cast<char>(0xF0 | (cp >> 18)));
            out.push_back(static_cast<char>(0x80 | ((cp >> 12) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | ((cp >> 6) & 0x3F)));
            out.push_back(static_cast<char>(0x80 | (cp & 0x3F)));
        }
    }

    Value parse_bool() {
        Value v;
        v.type = Value::Type::Bool;
        if (literal("true", 4)) {
            v.boolean = true;
            position += 4;
        } else if (literal("false", 5)) {
            v.boolean = false;
            position += 5;
        } else {
            throw std::runtime_error("bad JSON literal");
        }
        return v;
    }

    Value parse_number() {
        size_t start = position;
        while (position < size) {
            char c = data[position];
            if ((c >= '0' && c <= '9') || c == '-' || c == '+' || c == '.' ||
                c == 'e' || c == 'E') {
                ++position;
            } else {
                break;
            }
        }
        if (start == position) throw std::runtime_error("bad JSON number");
        std::string token(data + start, position - start);
        Value v;
        v.type = Value::Type::Number;
        v.number = std::strtod(token.c_str(), nullptr);
        return v;
    }
};

struct TopLevelScan {
    Value fields;                                 // 除流式数组外的全部顶层字段
    bool has_streamed_array = false;              // 指定数组字段是否存在
    std::vector<std::pair<size_t, size_t>> spans; // 数组元素文本 span
};

inline void Parser::parse_streaming_object(
    const std::string& stream_field, TopLevelScan& scan) {
    scan.fields.type = Value::Type::ObjectT;  // 缺省构造是 Null，必须显式置为对象
    skip();
    if (peek() != '{') throw std::runtime_error("expected JSON object");
    ++position;
    skip();
    if (peek() == '}') {
        ++position;
    } else {
        while (true) {
            skip();
            std::string key = parse_string();
            skip();
            expect(':');
            skip();
            if (key == stream_field) {
                if (peek() != '[') {
                    throw std::runtime_error("streamed field must be an array");
                }
                ++position;
                skip();
                scan.has_streamed_array = true;
                if (peek() != ']') {
                    while (true) {
                        skip();
                        size_t begin = position;
                        skip_value(0);
                        scan.spans.emplace_back(begin, position);
                        skip();
                        char c = peek();
                        if (c == ',') {
                            ++position;
                            continue;
                        }
                        if (c == ']') {
                            ++position;
                            break;
                        }
                        throw std::runtime_error("expected ',' or ']' in JSON array");
                    }
                } else {
                    ++position;
                }
            } else {
                scan.fields.fields[key] = parse_value(0);
            }
            skip();
            char c = peek();
            if (c == ',') {
                ++position;
                continue;
            }
            if (c == '}') {
                ++position;
                break;
            }
            throw std::runtime_error("expected ',' or '}' in JSON object");
        }
    }
    skip();
    if (position != size) throw std::runtime_error("trailing JSON content");
}

inline std::string Parser::scan_schema() {
    skip();
    if (peek() != '{') throw std::runtime_error("expected JSON object");
    ++position;
    skip();
    if (peek() == '}') return std::string();
    while (true) {
        skip();
        std::string key = parse_string();
        skip();
        expect(':');
        skip();
        if (key == "schema") {
            if (peek() != '"') throw std::runtime_error("schema must be a string");
            return parse_string();
        }
        skip_value(0);
        skip();
        char c = peek();
        if (c == ',') {
            ++position;
            continue;
        }
        if (c == '}') break;
        throw std::runtime_error("expected ',' or '}' in JSON object");
    }
    return std::string();
}

inline Value parse(const std::string& text) { return Parser(text).parse(); }

// 解析原始文本的一个片段（批量任务元素的 span）。
inline Value parse_span(const char* data, size_t begin, size_t end) {
    return Parser(data + begin, end - begin).parse();
}

// 批量请求体积可达数百 MB：顶层对象除 tasks 外的小字段建 DOM，tasks 元素只记
// 录文本 span，由 worker 各自解析、用完即弃，避免整批 DOM 常驻内存。
inline TopLevelScan scan_object_streaming(const std::string& text, const std::string& field) {
    Parser parser(text.data(), text.size());
    TopLevelScan scan;
    parser.parse_streaming_object(field, scan);
    return scan;
}

inline std::string scan_schema_field(const std::string& text) {
    Parser parser(text.data(), text.size());
    return parser.scan_schema();
}

inline void escape_into(const std::string& text, std::string& out) {
    out.push_back('"');
    for (unsigned char c : text) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (c < 0x20) {
                    char buffer[8];
                    std::snprintf(buffer, sizeof(buffer), "\\u%04x", c);
                    out += buffer;
                } else {
                    out.push_back(static_cast<char>(c));
                }
        }
    }
    out.push_back('"');
}

// 紧凑输出；基准响应体积可控，缩进交给阅读器。
inline void dump_into(const Value& v, std::string& out) {
    switch (v.type) {
        case Value::Type::Null: out += "null"; break;
        case Value::Type::Bool: out += v.boolean ? "true" : "false"; break;
        case Value::Type::Number: {
            char buffer[40];
            std::snprintf(buffer, sizeof(buffer), "%.17g", v.number);
            out += buffer;
            break;
        }
        case Value::Type::String: escape_into(v.text, out); break;
        case Value::Type::ArrayT: {
            out.push_back('[');
            for (size_t i = 0; i < v.items.size(); ++i) {
                if (i != 0) out.push_back(',');
                dump_into(v.items[i], out);
            }
            out.push_back(']');
            break;
        }
        case Value::Type::ObjectT: {
            out.push_back('{');
            bool first = true;
            for (const auto& entry : v.fields) {
                if (!first) out.push_back(',');
                first = false;
                escape_into(entry.first, out);
                out.push_back(':');
                dump_into(entry.second, out);
            }
            out.push_back('}');
            break;
        }
    }
}

inline std::string dump(const Value& v) {
    std::string out;
    dump_into(v, out);
    return out;
}

}  // namespace fjson

#endif  // FLINTNDE_ENGINE_JSON_HPP
