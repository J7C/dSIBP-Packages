# MadStree 更新日志

## 1.01

统一传播子术语为 (anti-)time-ordered 与 Wightman，内部类型改为 `massiveT`、`masslessT`、`massiveW`、`masslessW`，同步相关字段、示例和双语手册；公开输入与计算定义不变。

修复自动无穷远边界与 `ParameterRules` 固定外能量的衔接：所有外能量现在都从边界匹配值完整输运到首个用户点，再在扫描中保持给定值；受影响的旧调用需要重新计算。

更新 MadStree 程序包论文引用为 arXiv:2610.28526；加载提示和中英文说明现提供对应论文链接。
