"""不确定度分配：在目标合成标准不确定度下反解某输入量的最大允许不确定度。

设待分配输入为 x_k，其标准不确定度 u(x_k) = z 为未知量（以该输入自身单位表示），
其余输入量的不确定度保持给定值。完整保留相关项时，GUM 一阶总方差为 z 的二次式：

    u_c^2(z) = A z^2 + B z + C

其中 c_i 为公式在输入估计值处的灵敏度系数（已换算到结果显示单位），
r_ij 为相关系数，u_j（j≠k）为其余输入量的标准不确定度：

    A = c_k^2
    B = 2 c_k Σ_{j≠k} c_j u_j r_kj
    C = Σ_{i≠k} Σ_{j≠k} c_i c_j u_i u_j r_ij
      （即令 z=0 时的总方差）

在 z ≥ 0 上求解二次不等式 A z^2 + B z + C ≤ T^2（T 为目标合成标准不确定度），
返回满足目标的最大非负 z 及代回后的合成值。
"""
from __future__ import annotations

import math

from .errors import UnprocessableError
from .propagation import evaluate_model
from .units import (
    Quantity,
    convert_magnitude,
    parse_unit,
    same_dimensionality,
    unit_str,
)


def _quadratic_coefficients(model: dict, k: int) -> tuple[float, float, float]:
    """计算总方差关于待分配输入不确定度 z 的二次系数 (A, B, C)。"""
    inputs = model["inputs"]
    n = model["n"]
    corr = model["corr"]
    c = model["c_mags"]
    u = [p["u_magnitude"] for p in inputs]

    # A = c_k^2；相关矩阵对角元为 1，z^2 项不含相关结构。
    A = c[k] * c[k]

    # B z = 2 c_k z Σ_{j≠k} c_j u_j r_kj
    cross_sum = 0.0
    for j in range(n):
        if j == k:
            continue
        cross_sum += c[j] * u[j] * float(corr[k, j])
    B = 2.0 * c[k] * cross_sum

    # C：z = 0 时的总方差，仅由其余输入量（含彼此相关项）贡献。
    C = 0.0
    for i in range(n):
        if i == k:
            continue
        for j in range(n):
            if j == k:
                continue
            C += c[i] * u[i] * float(corr[i, j]) * c[j] * u[j]
    if C < 0.0 and C > -1e-12:
        C = 0.0

    return A, B, C


def _max_nonnegative_root(
    A: float, B: float, C: float, target_variance: float
) -> tuple[float, float]:
    """求解 A z^2 + B z + C ≤ target_variance 在 z ≥ 0 上的最大 z。

    返回 (最大不确定度 z_max, 可行区间下界 z_min)。
    目标不可达时抛出 422（code=target_unreachable）。
    调用前需保证 A > 0（灵敏度非零）。
    """
    D_const = C - target_variance  # 等价方程 A z^2 + B z + D_const = 0

    # 顶点 z* = -B/(2A)：抛物线上方差最小处（相对容差，随量级缩放）。
    scale = max(1.0, abs(B) / A if A > 0 else 0.0, math.sqrt(max(C, 0.0) / A if A > 0 else 0.0))
    tol = 1e-12 * scale * scale

    z_vertex = -B / (2.0 * A)
    vertex_in_domain = z_vertex >= -1e-12 * scale
    v = z_vertex if vertex_in_domain else 0.0
    min_variance = A * v * v + B * v + C
    if min_variance > target_variance + tol:
        raise UnprocessableError(
            f"目标不可达：即便最优地选取待分配输入的标准不确定度，"
            f"最小合成方差 {min_variance:.6g} 仍大于目标方差 {target_variance:.6g}"
            "（其余输入量的不确定度及其相关项已超过目标）",
            ["target_combined_standard_uncertainty"],
            code="target_unreachable",
        )

    disc = B * B - 4.0 * A * D_const
    if disc < -tol:
        # 顶点检查已通过时不应发生；作为数值兜底。
        raise UnprocessableError(
            "目标不可达：总方差二次方程无实根，给定目标合成标准不确定度无法满足",
            ["target_combined_standard_uncertainty"],
            code="target_unreachable",
        )
    if disc < 0.0:
        disc = 0.0

    # 相切判定：可行域上的最小方差与目标方差在数值容差内相等时，
    # 判别式受浮点舍入可能略呈正值；此时可行解唯一（顶点或边界 z=0）。
    tangent_gap = target_variance - min_variance
    tangent_tol = 1e-12 * max(
        1.0,
        abs(target_variance),
        abs(C),
        B * B / (4.0 * A),
    )
    if 0.0 <= tangent_gap <= tangent_tol:
        z_tangent = max(0.0, z_vertex) if vertex_in_domain else 0.0
        return z_tangent, z_tangent

    root = math.sqrt(disc)
    z_plus = (-B + root) / (2.0 * A)
    z_minus = (-B - root) / (2.0 * A)

    # 可行域为两根之间（disc>0）；与 z ≥ 0 取交集。
    root_zero_tol = 1e-12 * scale
    if abs(z_plus) <= root_zero_tol:
        z_plus = 0.0
    if abs(z_minus) <= root_zero_tol:
        z_minus = 0.0
    z_max = z_plus
    z_lower = max(0.0, z_minus)
    if disc == 0.0:
        z_lower = max(0.0, z_max)
    if z_max < -1e-12 * scale:
        raise UnprocessableError(
            "目标不可达：总方差二次不等式在非负标准不确定度范围内无解",
            ["target_combined_standard_uncertainty"],
            code="target_unreachable",
        )
    if z_max < 0.0:
        z_max = 0.0
    if z_lower > z_max:
        z_lower = z_max
    return z_max, z_lower


def allocate(req) -> dict:
    """执行不确定度分配，返回可直接序列化的结果字典。"""
    model = evaluate_model(req)
    inputs = model["inputs"]
    n = model["n"]
    corr = model["corr"]
    display_unit = model["display_unit"]
    sensitivities = model["sensitivities"]

    # 定位待分配输入。
    target_name = req.allocate_input
    k = next(
        (idx for idx, p in enumerate(inputs) if p["name"] == target_name),
        None,
    )
    if k is None:
        raise UnprocessableError(
            f"待分配输入量 {target_name!r} 不在 inputs 声明清单中",
            ["allocate_input"],
            code="unknown_input",
        )
    target_input = inputs[k]

    # 目标合成标准不确定度：单位量纲须与公式结果一致，并换算到结果显示单位。
    target = req.target_combined_standard_uncertainty
    target_loc = ["target_combined_standard_uncertainty", "unit"]
    if target.unit is None:
        target_unit = parse_unit("dimensionless", target_loc)
    else:
        target_unit = parse_unit(target.unit, target_loc)
    if not same_dimensionality(target_unit, display_unit):
        raise UnprocessableError(
            f"目标合成标准不确定度单位 {target_unit:~} 与公式结果单位 "
            f"{display_unit:~} 量纲不一致",
            target_loc,
            code="incompatible_unit",
        )
    target_mag = convert_magnitude(
        Quantity(float(target.value), target_unit),
        display_unit,
        ["target_combined_standard_uncertainty", "value"],
    )
    target_variance = target_mag * target_mag

    # 总方差关于未知量 z 的二次式。
    A, B, C = _quadratic_coefficients(model, k)

    # 灵敏度为零：未知量对合成不确定度没有贡献，不存在有限上界。
    if A <= 0.0:
        gap_tol = 1e-12 * max(1.0, abs(C), abs(target_variance))
        if C <= target_variance + gap_tol:
            raise UnprocessableError(
                f"待分配输入 {target_name!r} 的灵敏度系数为零，其标准不确定度"
                "不影响合成结果，满足目标的取值不存在有限上界，无法分配",
                ["allocate_input"],
                code="zero_sensitivity",
            )
        raise UnprocessableError(
            f"目标不可达：待分配输入 {target_name!r} 的灵敏度系数为零，"
            "无法通过调整其不确定度降低合成不确定度，而其余输入量的贡献"
            f"（方差 {C:.6g}）已超过目标方差 {target_variance:.6g}",
            ["allocate_input"],
            code="zero_sensitivity",
        )

    z_max, z_lower = _max_nonnegative_root(A, B, C, target_variance)

    # 代回：以分配值替换该输入的 u_k，重新计算合成标准不确定度。
    cu = list(model["cu_mags"])
    cu[k] = model["c_mags"][k] * z_max
    achieved_variance = 0.0
    for i in range(n):
        for j in range(n):
            achieved_variance += cu[i] * float(corr[i, j]) * cu[j]
    if achieved_variance < 0.0 and achieved_variance > -1e-12:
        achieved_variance = 0.0
    achieved_uc = math.sqrt(max(achieved_variance, 0.0))

    input_unit_str = unit_str(target_input["unit"])
    combined_unit_str = unit_str(display_unit)

    return {
        "result": {
            "input": target_name,
            "standard_uncertainty": z_max,
            "unit": input_unit_str,
            "achieved_combined_standard_uncertainty": achieved_uc,
            "target_combined_standard_uncertainty": target_mag,
            "combined_unit": combined_unit_str,
            "feasible_range": {
                "lower": z_lower,
                "upper": z_max,
                "unit": input_unit_str,
            },
        },
        "number_of_inputs": n,
        "correlation_matrix": [[float(x) for x in row] for row in corr],
        "sensitivity_coefficients": sensitivities,
        "formula": req.formula,
    }
