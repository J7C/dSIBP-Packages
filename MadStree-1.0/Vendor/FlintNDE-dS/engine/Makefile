# FlintNDE 引擎构建脚本：MSYS2 UCRT64 g++ + FLINT 3.x（Acb 并入）。
# 用法：在 MSYS2 UCRT64 shell 中 `make` 生成 build/flintnde_engine.exe。
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall -Wextra -Isrc
LDFLAGS ?= -flto
LDLIBS ?= -lflint

BUILD := build
TARGET := $(BUILD)/flintnde_engine.exe
SOURCES := src/main.cpp src/engine.cpp
HEADERS := src/engine.hpp src/json.hpp src/util.hpp

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(SOURCES) $(HEADERS) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SOURCES) $(LDLIBS) -o $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)
