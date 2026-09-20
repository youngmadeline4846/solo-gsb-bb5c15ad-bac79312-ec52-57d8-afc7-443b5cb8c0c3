"""端到端快速验证（不经过 HTTP，直接调传播函数与 API 客户端）。"""
import json
import sys

sys.path.insert(0, "/workspace")

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def show(title, resp):
    print("=" * 70)
    print(title)
    print(json.dumps(resp, ensure_ascii=False, indent=2))


# 1) 欧姆定律 P = U*I，单位校验 + 不相关
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
            {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02},
        ],
        "formula": "U*I",
    },
)
show("1) P=U*I (expect 20 W, uc=sqrt(0.2^2+0.2^2)=0.2828 W)", r.json())
assert r.status_code == 200, r.text
b = r.json()
assert abs(b["result"]["value"] - 20.0) < 1e-12
assert b["result"]["unit"] == "W"
assert abs(b["result"]["standard_uncertainty"] - (0.08 ** 0.5)) < 1e-12

# 2) 完全正相关，相关项使不确定度为 0.4
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
            {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02},
        ],
        "formula": "U*I",
        "correlation_matrix": [[1, 1], [1, 1]],
    },
)
show("2) P=U*I, r=1 (expect uc=0.4 W, share=0.5 each)", r.json())
assert r.status_code == 200
b = r.json()
assert abs(b["result"]["standard_uncertainty"] - 0.4) < 1e-12
assert abs(b["variance_contributions"][0]["share"] - 0.5) < 1e-9

# 3) 负相关抵消
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
            {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02},
        ],
        "formula": "U*I",
        "correlation_matrix": [[1, -1], [-1, 1]],
    },
)
show("3) P=U*I, r=-1 (expect uc=0)", r.json())
assert r.status_code == 200
assert abs(r.json()["result"]["standard_uncertainty"]) < 1e-9

# 4) 量纲不一致：V + A
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
            {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02},
        ],
        "formula": "U + I",
    },
)
show("4) V + A -> 422 dimension_mismatch", r.json())
assert r.status_code == 422
assert r.json()["error"]["code"] == "dimension_mismatch"
assert r.json()["error"]["loc"] == ["formula"]

# 5) 非对称矩阵，定位到 [1,0]
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "y", "value": 2, "unit": "m", "standard_uncertainty": 0.1},
        ],
        "formula": "x+y",
        "correlation_matrix": [[1, 0.5], [0.3, 1]],
    },
)
show("5) asymmetric -> 422 loc [correlation_matrix,1,1?]", r.json())
assert r.status_code == 422
assert r.json()["error"]["code"] == "asymmetric_matrix"
assert r.json()["error"]["loc"][0] == "correlation_matrix"

# 6) 非半正定（三个相关系数互相矛盾）
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "y", "value": 2, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "z", "value": 3, "unit": "m", "standard_uncertainty": 0.1},
        ],
        "formula": "x+y+z",
        "correlation_matrix": [
            [1, 0.9, 0.9],
            [0.9, 1, -0.9],
            [0.9, -0.9, 1],
        ],
    },
)
show("6) non-PSD -> 422", r.json())
assert r.status_code == 422
assert r.json()["error"]["code"] == "not_positive_semidefinite"

# 7) 白名单外函数 __import__
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [{"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1}],
        "formula": "__import__('os').system('id')",
    },
)
show("7) injection -> 422", r.json())
assert r.status_code == 422

# 8) 三角函数 + 角度单位，ln/log10
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "theta", "value": 30, "unit": "deg", "standard_uncertainty": 1},
            {"name": "L", "value": 2.0, "unit": "m", "standard_uncertainty": 0.01},
        ],
        "formula": "L * cos(theta)",
    },
)
show("8) L*cos(30 deg) (expect ~1.732 m)", r.json())
assert r.status_code == 200, r.text
b = r.json()
assert abs(b["result"]["value"] - 2.0 * (3 ** 0.5) / 2) < 1e-9
assert b["result"]["unit"] == "m"

# 9) 不确定度显式单位（mV 与 V 量纲一致）
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {
                "name": "U",
                "value": 10.0,
                "unit": "V",
                "standard_uncertainty": {"value": 100, "unit": "mV"},
            },
        ],
        "formula": "U*2",
    },
)
show("9) u=100 mV for 10 V (expect uc=0.2 V)", r.json())
assert r.status_code == 200, r.text
assert abs(r.json()["result"]["standard_uncertainty"] - 0.2) < 1e-12

# 10) 不确定度单位量纲不一致 -> 422 loc 到字段
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {
                "name": "U",
                "value": 10.0,
                "unit": "V",
                "standard_uncertainty": {"value": 0.1, "unit": "A"},
            },
        ],
        "formula": "U*2",
    },
)
show("10) u unit A vs V -> 422", r.json())
assert r.status_code == 422
assert r.json()["error"]["loc"][-1] == "unit"

# 11) log 量纲错误：log(m)
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [{"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1}],
        "formula": "log(x)",
    },
)
show("11) log(m) -> 422", r.json())
assert r.status_code == 422

# 12) 未声明符号
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [{"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1}],
        "formula": "x + y",
    },
)
show("12) unknown symbol y -> 422", r.json())
assert r.status_code == 422 and r.json()["error"]["code"] == "unknown_symbol"

# 13) 幂次为变量
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "x", "value": 2, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "n", "value": 2, "unit": "dimensionless", "standard_uncertainty": 0.0},
        ],
        "formula": "x**n",
    },
)
show("13) x**n -> 422 non_constant_exponent", r.json())
assert r.status_code == 422 and r.json()["error"]["code"] == "non_constant_exponent"

# 14) 未知单位
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [{"name": "x", "value": 1, "unit": "frobnicate", "standard_uncertainty": 0.1}],
        "formula": "x",
    },
)
show("14) unknown unit -> 422 loc inputs/0/unit", r.json())
assert r.status_code == 422 and r.json()["error"]["loc"] == ["inputs", 0, "unit"]

# 15) 相关系数越界
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "y", "value": 2, "unit": "m", "standard_uncertainty": 0.1},
        ],
        "formula": "x+y",
        "correlation_matrix": [[1, 1.2], [1.2, 1]],
    },
)
show("15) r=1.2 -> 422 loc element", r.json())
assert r.status_code == 422 and r.json()["error"]["code"] == "correlation_out_of_range"
assert r.json()["error"]["loc"] == ["correlation_matrix", 0, 1]

# 16) 对角线 != 1
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "y", "value": 2, "unit": "m", "standard_uncertainty": 0.1},
        ],
        "formula": "x+y",
        "correlation_matrix": [[0.9, 0], [0, 1]],
    },
)
show("16) diag=0.9 -> 422", r.json())
assert r.status_code == 422 and r.json()["error"]["code"] == "diagonal_not_one"

# 17) 矩阵形状不符
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "y", "value": 2, "unit": "m", "standard_uncertainty": 0.1},
        ],
        "formula": "x+y",
        "correlation_matrix": [[1]],
    },
)
show("17) bad shape -> 422", r.json())
assert r.status_code == 422 and r.json()["error"]["code"] == "bad_matrix_shape"

# 18) sqrt(m^2) 量纲
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "a", "value": 3, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "b", "value": 4, "unit": "m", "standard_uncertainty": 0.1},
        ],
        "formula": "sqrt(a^2 + b^2)",
    },
)
show("18) sqrt(a^2+b^2) (expect 5 m)", r.json())
assert r.status_code == 200, r.text
b = r.json()
assert abs(b["result"]["value"] - 5.0) < 1e-9
assert b["result"]["unit"] == "m"

# 19) 重复名称
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [
            {"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1},
            {"name": "x", "value": 2, "unit": "m", "standard_uncertainty": 0.1},
        ],
        "formula": "x",
    },
)
show("19) duplicate name -> 422", r.json())
assert r.status_code == 422 and r.json()["error"]["code"] == "duplicate_name"

# 20) 常数公式
r = client.post(
    "/api/v1/propagate",
    json={
        "inputs": [{"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1}],
        "formula": "pi*E",
    },
)
show("20) pi*E, no x dependence (uc=0)", r.json())
assert r.status_code == 200, r.text
assert abs(r.json()["result"]["standard_uncertainty"]) < 1e-12

# === POST /api/v1/allocate-uncertainty ====================================

ALLOC = "/api/v1/allocate-uncertainty"


def alloc_body(**overrides):
    body = {
        "inputs": [
            {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
            {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02},
        ],
        "formula": "U*I",
        "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]],
        "allocate_input": "U",
        "target_combined_standard_uncertainty": {"value": 0.5, "unit": "W"},
    }
    body.update(overrides)
    return body


# 21) 基本分配：目标 uc=0.5 W，相关项完整保留
#    V(z)=4z^2+0.4z+0.04 <= 0.25 -> z = (-0.4+sqrt(3.52))/8
r = client.post(ALLOC, json=alloc_body())
show("21) allocate U, target uc=0.5 W (expect z=0.18452 V)", r.json())
assert r.status_code == 200, r.text
b = r.json()["result"]
assert abs(b["standard_uncertainty"] - ((-0.4 + 3.52 ** 0.5) / 8)) < 1e-12
assert b["unit"] == "V"
assert b["combined_unit"] == "W"
assert abs(b["achieved_combined_standard_uncertainty"] - 0.5) < 1e-9
assert b["feasible_range"]["lower"] == 0.0
assert abs(b["feasible_range"]["upper"] - b["standard_uncertainty"]) < 1e-12
# 待分配输入自带的 standard_uncertainty 被忽略
r2 = client.post(ALLOC, json=alloc_body(inputs=[
    {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 999.0},
    {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02},
]))
assert r2.status_code == 200
assert abs(r2.json()["result"]["standard_uncertainty"] - b["standard_uncertainty"]) < 1e-12

# 22) 与 /propagate 交叉验证：代回分配值后合成值达标，且再大即超标
z = b["standard_uncertainty"]
rp = client.post("/api/v1/propagate", json={
    "inputs": [
        {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": z},
        {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02},
    ],
    "formula": "U*I",
    "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]],
})
assert abs(rp.json()["result"]["standard_uncertainty"] - 0.5) < 1e-9
rp2 = client.post("/api/v1/propagate", json={
    "inputs": [
        {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": z * 1.001 + 1e-9},
        {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.02},
    ],
    "formula": "U*I",
    "correlation_matrix": [[1.0, 0.5], [0.5, 1.0]],
})
assert rp2.json()["result"]["standard_uncertainty"] > 0.5

# 23) 完全负相关 + 目标 0：可行区间退化为单点 z=0.1 V
r = client.post(ALLOC, json=alloc_body(
    correlation_matrix=[[1.0, -1.0], [-1.0, 1.0]],
    target_combined_standard_uncertainty={"value": 0.0, "unit": "W"},
))
show("23) r=-1, target 0 (single feasible z=0.1 V)", r.json())
assert r.status_code == 200
b = r.json()["result"]
assert abs(b["standard_uncertainty"] - 0.1) < 1e-12
assert abs(b["feasible_range"]["lower"] - 0.1) < 1e-12
assert abs(b["feasible_range"]["upper"] - 0.1) < 1e-12
assert b["achieved_combined_standard_uncertainty"] < 1e-12

# 24) 完全负相关形成有上下界的区间：f=x+y, u_y=0.5, V=(z-0.5)^2, target=sqrt(0.1)
r = client.post(ALLOC, json={
    "inputs": [
        {"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1},
        {"name": "y", "value": 2, "unit": "m", "standard_uncertainty": 0.5},
    ],
    "formula": "x+y",
    "correlation_matrix": [[1.0, -1.0], [-1.0, 1.0]],
    "allocate_input": "x",
    "target_combined_standard_uncertainty": {"value": 0.1 ** 0.5, "unit": "m"},
})
show("24) bounded feasible interval via r=-1", r.json())
assert r.status_code == 200
b = r.json()["result"]
assert abs(b["standard_uncertainty"] - (0.5 + 0.1 ** 0.5)) < 1e-12
assert abs(b["feasible_range"]["lower"] - (0.5 - 0.1 ** 0.5)) < 1e-12

# 25) 目标不可达：其余输入量自身贡献（uc=1 W）已超目标 0.5 W
r = client.post(ALLOC, json=alloc_body(inputs=[
    {"name": "U", "value": 10.0, "unit": "V", "standard_uncertainty": 0.1},
    {"name": "I", "value": 2.0, "unit": "A", "standard_uncertainty": 0.1},
]))
show("25) target unreachable -> 422 target_unreachable", r.json())
assert r.status_code == 422
e = r.json()["error"]
assert e["code"] == "target_unreachable"
assert e["loc"] == ["target_combined_standard_uncertainty"]

# 26) 灵敏度为零（公式不含待分配输入）：无有限上界 -> 422 zero_sensitivity
r = client.post(ALLOC, json=alloc_body(formula="I*2",
                                       target_combined_standard_uncertainty={"value": 0.1, "unit": "A"}))
show("26) zero sensitivity / unbounded -> 422 zero_sensitivity", r.json())
assert r.status_code == 422
e = r.json()["error"]
assert e["code"] == "zero_sensitivity"
assert e["loc"] == ["allocate_input"]

# 27) 灵敏度为零且其余贡献已超目标，同样 422
r = client.post(ALLOC, json=alloc_body(formula="I*2",
                                       target_combined_standard_uncertainty={"value": 0.01, "unit": "A"}))
assert r.status_code == 422 and r.json()["error"]["code"] == "zero_sensitivity"

# 28) 待分配输入未声明 -> 422 unknown_input
r = client.post(ALLOC, json=alloc_body(allocate_input="P"))
show("28) unknown allocate_input -> 422", r.json())
assert r.status_code == 422
e = r.json()["error"]
assert e["code"] == "unknown_input" and e["loc"] == ["allocate_input"]

# 29) 目标单位量纲与结果不一致 -> 422，定位到 unit
r = client.post(ALLOC, json=alloc_body(
    target_combined_standard_uncertainty={"value": 0.5, "unit": "A"}))
show("29) target unit A vs result W -> 422", r.json())
assert r.status_code == 422
e = r.json()["error"]
assert e["code"] == "incompatible_unit"
assert e["loc"] == ["target_combined_standard_uncertainty", "unit"]

# 30) 目标单位自动换算：f=x（m），目标 2000 mm -> z=2 m
r = client.post(ALLOC, json={
    "inputs": [{"name": "x", "value": 5, "unit": "m", "standard_uncertainty": 0.1}],
    "formula": "x",
    "allocate_input": "x",
    "target_combined_standard_uncertainty": {"value": 2000, "unit": "mm"},
})
show("30) target 2000 mm -> z=2 m", r.json())
assert r.status_code == 200
b = r.json()["result"]
assert abs(b["standard_uncertainty"] - 2.0) < 1e-12 and b["unit"] == "m"
assert abs(b["target_combined_standard_uncertainty"] - 2.0) < 1e-12

# 31) 无量纲目标可省略单位
r = client.post(ALLOC, json={
    "inputs": [{"name": "x", "value": 2.0, "unit": "dimensionless",
                "standard_uncertainty": 0.1}],
    "formula": "x^2",
    "allocate_input": "x",
    "target_combined_standard_uncertainty": {"value": 0.4},
})
show("31) dimensionless target without unit (z=0.1)", r.json())
assert r.status_code == 200
assert abs(r.json()["result"]["standard_uncertainty"] - 0.1) < 1e-12

# 32) 复用相关矩阵校验：非对称矩阵
r = client.post(ALLOC, json=alloc_body(
    inputs=[
        {"name": "x", "value": 1, "unit": "m", "standard_uncertainty": 0.1},
        {"name": "y", "value": 2, "unit": "m", "standard_uncertainty": 0.1},
    ],
    formula="x+y",
    correlation_matrix=[[1.0, 0.5], [0.3, 1.0]],
    target_combined_standard_uncertainty={"value": 1.0, "unit": "m"},
))
show("32) asymmetric matrix -> 422", r.json())
assert r.status_code == 422 and r.json()["error"]["code"] == "asymmetric_matrix"

# 33) 复用公式量纲校验
r = client.post(ALLOC, json=alloc_body(formula="U + I",
                                       target_combined_standard_uncertainty={"value": 1, "unit": "V"}))
assert r.status_code == 422 and r.json()["error"]["code"] == "dimension_mismatch"

# 34) 目标为负数 -> 422（请求体校验，带 loc）
r = client.post(ALLOC, json=alloc_body(
    target_combined_standard_uncertainty={"value": -0.5, "unit": "W"}))
show("34) negative target -> 422", r.json())
assert r.status_code == 422
assert r.json()["error"]["loc"][0] == "target_combined_standard_uncertainty"

print("\nALL TESTS PASSED")
