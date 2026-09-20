# 校准计算不确定度传播 API

基于 FastAPI + SymPy + Pint，对测量模型 $y=f(x_1,\dots,x_n)$ 做 GUM 一阶不确定度传播：
服务端解析仅含白名单运算的公式并进行量纲校验，利用符号偏导与相关系数矩阵计算结果值、
单位、合成标准不确定度，并按输入量列出含相关项的方差贡献。

## 启动

```bash
docker compose up --build
```

依赖在**容器启动时**安装；服务随后监听 API 端口 **8000**。

## 配置

- 宿主机端口：通过环境变量 `PORT` 覆盖（默认 `8000`），容器内固定监听 `8000`。
  例：`PORT=9000 docker compose up --build`
- 依赖清单：`requirements.txt`（容器启动时由 `docker/entrypoint.sh` 安装）。

## API（端口 8000）

- `GET /health` — 健康检查
- `GET /api/v1/whitelist` — 公式白名单
- `POST /api/v1/propagate` — 不确定度传播；相关矩阵非对称 / 非半正定 / 单位量纲不兼容 / 公式不合法时返回 **422**，错误体 `error.loc` 精确定位（如 `correlation_matrix.0.1`、`inputs.0.unit`、`formula`）
- `POST /api/v1/allocate-uncertainty` — 不确定度分配：指定一个待分配输入量与带单位的目标合成标准不确定度，反解满足目标的**最大非负标准不确定度**；目标不可达 / 待分配输入灵敏度为零（无有限上界）/ 单位量纲不兼容时返回 **422**
- 交互文档：`http://localhost:8000/docs`

## 请求示例

### 不确定度传播

```bash
curl -s http://localhost:8000/api/v1/propagate \
  -H 'Content-Type: application/json' \
  -d '{
    "inputs": [
      {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
      {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02}
    ],
    "formula": "U*I",
    "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]]
  }'
```

### 不确定度分配

在 `/api/v1/propagate` 同款请求体上增加两个字段：

- `allocate_input`：待分配输入量名称（须出现在 `inputs` 中；该输入自带的
  `standard_uncertainty` 会被忽略，视为未知量）；
- `target_combined_standard_uncertainty`：带单位的目标合成标准不确定度
  （`{value, unit}`，无量纲时可省略 `unit`；单位须与公式结果量纲一致）。

将该输入的不确定度记为 $z$，完整保留相关项后总方差是关于 $z$ 的二次式

$$u_c^2(z)=A z^2+Bz+C,\quad
A=c_k^2,\; B=2c_k\sum_{j\ne k}c_j u_j r_{kj},\;
C=\sum_{i,j\ne k}c_i c_j u_i u_j r_{ij},$$

服务端在 $z\ge 0$ 上解二次不等式 $u_c^2(z)\le T^2$，返回可行区间的上界（最大允许标准
不确定度，以该输入自身单位表示）及代回后的合成标准不确定度。

```bash
curl -s http://localhost:8000/api/v1/allocate-uncertainty \
  -H 'Content-Type: application/json' \
  -d '{
    "inputs": [
      {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
      {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02}
    ],
    "formula": "U*I",
    "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]],
    "allocate_input": "U",
    "target_combined_standard_uncertainty": {"value": 0.5, "unit": "W"}
  }'
```

响应：

```json
{
  "result": {
    "input": "U",
    "standard_uncertainty": 0.1845207879911715,
    "unit": "V",
    "achieved_combined_standard_uncertainty": 0.5,
    "target_combined_standard_uncertainty": 0.5,
    "combined_unit": "W",
    "feasible_range": {"lower": 0.0, "upper": 0.1845207879911715, "unit": "V"}
  },
  "number_of_inputs": 2,
  "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]],
  "sensitivity_coefficients": [
    {"input": "U", "value": 2.0, "unit": "W / V"},
    {"input": "I", "value": 10.0, "unit": "W / A"}
  ],
  "formula": "U*I"
}
```

`feasible_range` 给出标准不确定度的完整可行区间：强负相关时二次式顶点可能落在正域，
下界会大于 0（此时不确定度过小反而因相关抵消不足而超标）；上界即返回的最大允许值。

分配接口的 **422** 错误码：

| code | 含义 | loc |
| --- | --- | --- |
| `target_unreachable` | 目标不可达：即便最优选取该输入的不确定度，其余输入量贡献（含相关项）仍超过目标 | `target_combined_standard_uncertainty` |
| `zero_sensitivity` | 待分配输入的灵敏度系数为零，其不确定度不影响合成结果，不存在有限上界 | `allocate_input` |
| `unknown_input` | `allocate_input` 不在 `inputs` 中 | `allocate_input` |
| `incompatible_unit` | 目标单位与公式结果量纲不一致 | `target_combined_standard_uncertainty.unit` |

公式、单位与相关矩阵的校验错误沿用 `/api/v1/propagate` 的错误码与 loc 规则。
