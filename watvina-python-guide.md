# WATVina Python 使用指南

WATVina 的 Python 模块名为 `watvina_python`。它暴露配置、解析结果、`DockingDriver` 和结构化结果对象，不需要通过 subprocess 解析 CLI 文本。

本指南以 `src/python/watvina_python.cpp` 当前 binding 为准。文中所有示例都经过实际运行验证。源码仓库示例的路径假定当前目录是项目根目录；只下载
`watvina.zip` 的用户应将示例里的 receptor、ligand 和 `vina.in` 换成自己的文件，因为精简发布包不包含源码树的 `demo/` 数据。

## 目录

1. [选择正确的 Python 产物](#1-选择正确的-python-产物)
2. [5 分钟快速上手](#2-5-分钟快速上手)
3. [推荐入口：复用 CLI 解析](#3-推荐入口复用-cli-解析)
4. [`ParseResult`](#4-parseresult)
5. [内存输入](#5-内存输入)
6. [Docking、打分与结果](#6-docking打分与结果)
7. [Cookbook 实战配方](#7-cookbook-实战配方)
8. [输出字符串](#8-输出字符串)
9. [Mini 结果](#9-mini-结果)
10. [FHFT 水预测](#10-fhft-水预测)
11. [离散水 docking](#11-离散水-docking)
12. [性能与并行](#12-性能与并行)
13. [API 总览](#13-api-总览)
14. [错误处理与排错](#14-错误处理与排错)
15. [当前边界](#15-当前边界)

## 1. 选择正确的 Python 产物

当前交付的是 CPython 原生扩展模块，不是与平台和 Python 版本无关的wheel。文件必须同时匹配操作系统、CPU 架构和 CPython 主次版本：

| 平台 | 分发文件 | 要求 |
|---|---|---|
| Linux x86-64 | `bin/watvina_python.cpython-312-x86_64-linux-gnu.so` | CPython 3.12；当前产物最高要求 `GLIBC_2.38`、`GLIBCXX_3.4.32` |
| Windows x86-64 | `bin/watvina_python.pyd` | CPython 3.12；模块依赖 `python312.dll` |

直接使用 `bin/` 中的 Linux 模块：

```bash
PYTHONPATH=/path/to/watvina/bin python3.12 -c \
  'import watvina_python as wv; print(wv.__doc__)'
```

首次导入时建议同时核对版本和模块实际加载路径，避免误用环境中残留的旧模块：

```bash
PYTHONPATH=/path/to/watvina/bin python3.12 -c \
  'import watvina_python as wv; print(wv.__version__); print(wv.__file__)'
```

本发布版应输出 `2026.7.20`，并指向刚解压目录中的 `.so`。

Windows PowerShell：

```powershell
$env:PYTHONPATH = "C:\path\to\watvina\bin"
python -c "import watvina_python as wv; print(wv.__doc__)"
```

不要把 CPython 3.12 产物用于 Python 3.11/3.13，也不要混用 Linux `.so`和 Windows `.pyd`。如果目标 Python ABI 不同，应使用该解释器从源码构建。

## 2. 5 分钟快速上手

以下流程使用**源码仓库**自带的 `demo/1a3l` 示例，在普通桌面机上几十秒内完成。二进制发布包不带演示结构；使用发布包时只需把三个路径替换为自己的全氢 receptor、全氢 ligand 和 box 配置文件。

### 2.1 第一次 docking

```python
import watvina_python as wv

parsed = wv.parse_argv([
    "-c", "demo/1a3l/vina.in",      # 搜索盒子（center/size）
    "-r", "demo/1a3l/rec.pdb",      # 受体（显式全氢）
    "-l", "demo/1a3l/i.sdf",        # 配体（显式全氢）
    "--exhaustiveness", "2",        # 演示用小预算；正式计算用默认 8+
    "--num_modes", "3",
    "--seed", "42",                 # 固定种子，结果可复现
    "--cpu", "4",
])
if not parsed.success:
    raise RuntimeError(parsed.error_message)

driver = wv.DockingDriver.from_parse(parsed)
result = driver.run()
if not result.ok:
    raise RuntimeError(result.error_message)

print("输入构象能量:", round(result.initial_energy, 3))
for rank, pose in enumerate(result.poses, start=1):
    print(f"mode {rank}: energy={pose.energy:.3f} "
          f"torsions={pose.pose.num_torsions}")

with open("docked.sdf", "w") as handle:
    handle.write(driver.poses_to_string(result, "sdf",
                                        num_modes=3, energy_range=3.0))
```

实测输出（`--seed 42` 时可复现到 ~0.01 kcal/mol）：

```text
输入构象能量: -6.778
mode 1: energy=-6.987 torsions=3
mode 2: energy=-6.003 torsions=3
mode 3: energy=-5.653 torsions=3
```

### 2.2 第一次水预测

```python
parsed = wv.parse_argv([
    "-r", "demo/1a3l/rec.pdb",
    "-c", "demo/1a3l/vina.in",
    "--predict_water",
])
driver = wv.DockingDriver.from_parse(parsed)
prediction = driver.predict_water()

print("位点数:", len(prediction.sites))
print("收敛:", prediction.diagnostics.converged)
print("格点:", prediction.grid.dimensions, "@", prediction.grid.spacing, "Å")

with open("waters.pdb", "w") as handle:
    handle.write(wv.format_water_pdb(prediction))
with open("occupancy.dx", "w") as handle:
    handle.write(wv.format_water_open_dx(prediction))
```

1a3l 演示盒（约 18 Å 边长）实测得到 34 个水合位点，`converged=True`。FHFT 外部场构建是耗时大头（本例约 40 s，见 [§12](#12-性能与并行)）。

## 3. 推荐入口：复用 CLI 解析

`parse_argv()` 是功能最完整且最不容易与 CLI 漂移的 Python 入口。传入的列表不包含程序名：

```python
import watvina_python as wv

parsed = wv.parse_argv([
    "-c", "vina.in",
    "-r", "receptor.pdb",
    "-l", "ligand.sdf",
    "--cpu", "4",
    "--seed", "42",
    "--num_modes", "10",
])

if not parsed.success:
    raise RuntimeError(parsed.error_message)

driver = wv.DockingDriver.from_parse(parsed)
result = driver.run()
if not result.ok:
    raise RuntimeError(result.error_message)

sdf = driver.poses_to_string(
    result,
    format="sdf",
    num_modes=10,
    energy_range=3.0,
)
with open("docked.sdf", "w") as handle:
    handle.write(sdf)
```

`parse_argv()` 使用与二进制相同的 cxxopts/config 合并逻辑：配置文件先加载，argv 覆盖同名 key；不要在 argv 自身重复同一参数。

**经验法则**：凡是 `watvina --help` / `--help_water` / `--help_ph4` /
`--help_relax` 里出现的选项，都能以字符串形式传给 `parse_argv()`；
`DockingConfig` attribute 直写只覆盖常用子集（见 §5.2）。

## 4. `ParseResult`

可读写字段：

| 字段 | 类型 | 含义 |
|---|---|---|
| `config` | `DockingConfig` | 解析后的运行配置 |
| `ligands` | `list[str]` | ligand 文件名；多个元素表示联合 docking |
| `out_name` | `str` | 输出名 |
| `scaffold_file` | `str` | scaffold SDF |
| `root_indices` | `str` | 0-based ROOT 索引字符串 |
| `receptor_content` | `str` | 内存 receptor 内容 |
| `receptor_format` | `str` | `pdb` 或 `pdbqt` |
| `flex_content` | `str` | 内存 flex PDBQT；PDB flex 应写在 receptor block 内 |
| `ligand_contents` | `list[str]` | 内存 ligand records |
| `ligand_formats` | `list[str]` | 与 contents 对齐的 `sdf`/`pdbqt` |
| `score_only` | `bool` | score-only 标志 |
| `success` | `bool` | 参数解析是否成功 |
| `error_message` | `str` | 失败原因 |

## 5. 内存输入

### 5.1 `DockingDriver.from_strings`

```python
import watvina_python as wv

config = wv.DockingConfig()
config.center_x = 29.931
config.center_y = 15.142
config.center_z = 12.455
config.size_x = 17.477
config.size_y = 18.646
config.size_z = 17.013
config.exhaustiveness = 4
config.population_size = 8
config.ga_searching = 4
config.num_modes = 10
config.cpu = 4

with open("receptor.pdb") as handle:
    receptor = handle.read()
with open("ligand.sdf") as handle:
    ligand = handle.read()

driver = wv.DockingDriver.from_strings(
    config=config,
    receptor=receptor,
    receptor_format="pdb",
    ligands=[ligand],
    ligand_formats=["sdf"],
)
result = driver.run()
```

`ligands` 与 `ligand_formats` 长度必须一致。全部输入仍必须显式全氢。若 receptor 是 PDB，柔性残基必须通过 receptor 字符串内部的 `REMARK FLEX` 表达；`flex=` 只用于 PDBQT flex 内容。

`from_strings` 适合配体来自数据库/SMILES 在线生成、不落盘的场景；但它只能设置 `DockingConfig` 暴露的字段——`seed`、打分权重等仍需
`parse_argv()`（两条路径可在同一项目里混用）。

### 5.2 `DockingConfig` 暴露字段

| 字段 | 类型 | 默认 | 含义 |
|---|---|---|---|
| `rigid_name` / `flex_name` | str | `""` | 受体/柔性残基文件名（记录用途） |
| `center_x/y/z` | float | — | 盒子中心（Å），无默认，必须设置 |
| `size_x/y/z` | float | 25.0 | 盒子边长（Å） |
| `exhaustiveness` | int | 8 | 全局搜索预算倍率 |
| `population_size` | int | 8 | GA 种群大小 |
| `ga_searching` | int | 4 | GA 搜索步数乘子 |
| `num_modes` | int | 10 | RMSD 聚类后保留的 mode 数 |
| `local_steps` | int | 0 | L-BFGS 步数；0 = 按分子大小启发式 |
| `cpu` | int | 0 | 并行搜索任务数；0 = 用满硬件线程 |
| `verbosity` | int | 2 | 日志级别 0–3 |
| `mini` | str | `"off"` | MM 松弛模式 `off`/`ligand`/`site` |
| `mini_optimizer` | str | `"lbfgs"` | MM 优化器 `lbfgs`/`abnr` |
| `mini_steps` | int | 80 | 每个输出 pose 的最大 MM 迭代 |
| `mini_gradient_tolerance` | float | 0.05 | MM 梯度收敛阈值 |
| `mini_dielectric_k` | float | 4.0 | 距离相关介电系数 |
| `mini_forcefield` | str | `"amber14sb"` | 模板力场 |
| `predict_water` | bool | False | 是否进行 FHFT 水预测 |
| `water_ff` | str | `"amber14sb"` | 水外部场力场 |
| `water_model` | str | `"tip3p"` | 刚性水模型 `tip3p`/`tip4p` |
| `water_grid_spacing` | float | 0.5 | FHFT 格点间距（Å） |
| `water_hb_epsilon` | float | 0.75 | 协同 O–O 壳层强度（kcal/mol） |
| `water_directional_cooperativity` | float | 0.75 | 第一壳层方向协同无量纲乘子 |
| `water_min_excess` | float | 0.25 | 位点最小超额占据 |
| `water_min_distance` | float | 2.5 | NMS 基础距离（Å） |
| `water_max_sites` | int | 0 | 位点数上限；0 = 不限 |
| `water_bias_weight` | float | 1.0 | 离散水置换偏置权重 |
| `water_dx_out` / `water_sites_out` | str | `""` | DX/PDB 输出路径；空 = 默认名 |
| `water_reference` | str | `""` | 参考水 PDB（验证用） |

并非所有 C++/CLI 字段都直接暴露为 Python attribute。例如 `seed`、打分权重、`rmsd`、`energy_range`、`template` 和 water 文件名应通过 `parse_argv()` 设置。需要完整 CLI 功能时优先使用 `parse_argv()`。

水取向求积、介电模型和网络化学势字段也没有直接暴露为 attribute：
`water_orientation_axis_order`、`water_orientation_azimuth_count`、`water_orientation_spin_count`、`water_dielectric_mode`、各介电参数以及`water_network_mu` 必须通过 `parse_argv()` 设置。

## 6. Docking、打分与结果

### 6.1 全局 docking

```python
result = driver.run()
if not result.ok:
    raise RuntimeError(result.error_message)
```

### 6.2 Score only

```python
parsed = wv.parse_argv([
    "-c", "vina.in",
    "-r", "receptor.pdb",
    "-l", "pose.sdf",
    "--score_only",
])
driver = wv.DockingDriver.from_parse(parsed)
result = driver.score_only()
print(result.initial_energy)
```

`score_breakdown` 当前没有暴露给 Python。需要完整逐项表格时使用 CLI。

### 6.3 `DockingResult`

| 字段 | 类型 | 含义 |
|---|---|---|
| `ok` | bool | 是否成功 |
| `poses` | list[`PoseOutput`] | 最终 pose，按能量排序 |
| `initial_energy` | float | 输入构象能量 |
| `error_message` | str | 失败原因 |
| `chem_summary` | str | 化学系统摘要 |
| `mini_relax` | `MiniRelaxSummary` | mini 结构化结果 |

### 6.4 `PoseOutput`

| 字段 | 类型 | 含义 |
|---|---|---|
| `energy` | float | 最终 pose 排序能量，包含启用的搜索偏置 |
| `pose` | `PoseConf` | 第一个 ligand 的构象 |
| `ligand_poses` | list[`PoseConf`] | 联合 docking 中全部 ligand 的构象 |

`PoseOutput` 当前不直接暴露 Cartesian coordinates、raw energy 或 channel breakdown。用 `poses_to_string()` 得到带坐标和元数据的 SDF/PDBQT（SDF 中的 `WV_*` 属性标签含义见 README「输出解读」一节）。

### 6.5 `PoseConf`

```python
pose = result.poses[0].pose
print(pose.position)        # (x, y, z) 三元组
print(pose.num_torsions)    # 可旋转键数

pose.position = (30.0, 15.0, 12.0)   # 可写：刚体平移
energy = driver.eval_at_pose(pose)   # 在该构象处重新求能
```

| 成员 | 类型 | 读写 | 含义 |
|---|---|---|---|
| `position` | tuple[float, float, float] | RW | 配体刚体中心 |
| `num_torsions` | int | RO | 可旋转键数 |

Python 可修改刚体平移并读取 torsion 数量；当前 binding 不提供逐 torsion setter。`eval_at_pose()` 当前只接受恰好包含一个 ligand 的 driver；联合 docking 的单个 `PoseConf` 不足以定义其他 ligand 的坐标，因此会明确抛出异常。

## 7. Cookbook 实战配方

### 7.1 批量虚拟筛选并导出 CSV

```python
import csv
import glob

import watvina_python as wv

rows = []
for ligand_path in sorted(glob.glob("library/*.sdf")):
    parsed = wv.parse_argv([
        "-c", "vina.in",
        "-r", "receptor.pdb",
        "-l", ligand_path,
        "--exhaustiveness", "8",
        "--seed", "42", "--cpu", "8",
    ])
    if not parsed.success:
        print(f"[skip] {ligand_path}: {parsed.error_message}")
        continue

    try:
        result = wv.DockingDriver.from_parse(parsed).run()
    except Exception as error:                # 单个失败不拖垮整批
        print(f"[fail] {ligand_path}: {error}")
        continue
    if not result.ok:
        print(f"[fail] {ligand_path}: {result.error_message}")
        continue

    best = result.poses[0].energy if result.poses else float("nan")
    rows.append((ligand_path, best, len(result.poses)))

with open("ranking.csv", "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["ligand", "best_energy", "n_poses"])
    writer.writerows(sorted(rows, key=lambda row: row[1]))
```

注意：每个配体都新建 `DockingDriver`，受体参数化（surface mask、directional cache）会随之重复——这是当前 API 的固有开销，万级配体库应按 §12 做进程级并行。

### 7.2 批量 score_only（重打分）

```python
def score_pose(receptor_pdb, pose_sdf, config="vina.in"):
    parsed = wv.parse_argv([
        "-c", config, "-r", receptor_pdb, "-l", pose_sdf, "--score_only",
    ])
    if not parsed.success:
        raise ValueError(parsed.error_message)
    result = wv.DockingDriver.from_parse(parsed).score_only()
    if not result.ok:
        raise RuntimeError(result.error_message)
    return result.initial_energy
```

### 7.3 水预测 → 保存 → 定量验证

```python
parsed = wv.parse_argv(["-r", "receptor.pdb", "-c", "vina.in", "--predict_water"])
driver = wv.DockingDriver.from_parse(parsed)
prediction = driver.predict_water()

if not prediction.diagnostics.converged:
    print("警告: FHFT 未收敛 —", prediction.diagnostics.message)

open("waters.pdb", "w").write(wv.format_water_pdb(prediction))
open("occupancy.dx", "w").write(wv.format_water_open_dx(prediction))

# 与晶体参考水做 1.4 Å 匹配（精度/召回/AP）
metrics = wv.validate_water_prediction(
    prediction, "crystal_waters.pdb", cutoff=1.4,
)
print(f"P={metrics.precision:.2f} R={metrics.recall:.2f} "
      f"AP={metrics.average_precision:.2f} "
      f"({metrics.matched}/{metrics.predicted} 命中)")
```

`validate_water_prediction()` 按参考 PDB 中的全部水氧计算，不会自动按docking box 过滤；调用者应确保参考水集合与预测区域一致。

### 7.4 水预测 → 水偏置 docking 闭环

`format_water_pdb()` 的输出可直接作为 `-w` 输入（fixed-column 格式）：

```python
# 第一步：预测
pred = wv.DockingDriver.from_parse(wv.parse_argv([
    "-r", "receptor.pdb", "-c", "vina.in", "--predict_water",
])).predict_water()
with open("waters.pdb", "w") as handle:
    handle.write(wv.format_water_pdb(pred))

# 第二步：带置换偏置 docking
parsed = wv.parse_argv([
    "-r", "receptor.pdb", "-l", "ligand.sdf", "-c", "vina.in",
    "-w", "waters.pdb", "--water_bias_weight", "0.1",
    "--seed", "42", "--cpu", "8",
])
driver = wv.DockingDriver.from_parse(parsed)
result = driver.run()
sdf = driver.poses_to_string(result, "sdf", 10, 3.0)   # 含 WV_WATER_BIAS 标签
```

`PoseOutput.energy` 的搜索目标包含 WaterBias；报告的 `Affinity`（SDF 中的 `WV_ENERGY`）仍排除它。

### 7.5 占据场导入 numpy 做体素分析

`WaterPredictionResult.occupancy`（以及 `mean_field`、`effective_energy`）是长度为 `grid.size` 的 `list[float]`，按 **x 最快** 展平
（`flat = x + nx*(y + ny*z)`）。用 numpy C-order reshape 为 `(nz, ny, nx)`即可还原三维场：

```python
import numpy as np

nx, ny, nz = prediction.grid.dimensions
occupancy = np.asarray(prediction.occupancy).reshape((nz, ny, nx))

# 体素 (ix, iy, iz) 的值：occupancy[iz, iy, ix]
```

位点的权威下标是 `site.grid_index`——`site.position` 是 GCMC 网络精炼后的坐标，**不一定落在体素中心**；对 position 做 round 反推体素可能落在相邻格子上（实测可差一行）。因此按位点查场值应直接用：

```python
site = prediction.sites[0]
assert prediction.occupancy[site.grid_index] == site.occupancy
```

### 7.6 `eval_at_pose` 单点能量

```python
pose = result.poses[0].pose
pose.position = (x0 + 0.5, y0, z0)     # 刚体微扰
energy = driver.eval_at_pose(pose)     # 返回 float
```

适合做结合姿态的刚体扫描/敏感性分析。限制：仅单配体 driver；torsion 不可逐键设置（见 §15）。

## 8. 输出字符串

### 8.1 Ligand poses

```python
sdf = driver.poses_to_string(result, "sdf", 10, 3.0)
pdbqt = driver.poses_to_string(result, "pdbqt", 10, 3.0)
```

参数：

- `format`：`sdf` 或 `pdbqt`；
- `num_modes`：最多输出 mode 数；
- `energy_range`：相对最佳 pose 的能量窗口。

SDF 包含 `WV_ENERGY`、`WV_RAW_ENERGY`、`WV_WATER_BIAS`、可选的`WV_MODE` 和 `WV_LIGAND`。

### 8.2 Flex poses

```python
text = driver.flex_poses_to_string(
    result,
    format="pdb",
    num_modes=10,
    energy_range=3.0,
)
```

format 为 `pdb` 或 `pdbqt`。无 flex 时返回空字符串。

### 8.3 Mini site poses

```python
if result.mini_relax.site_pose_count:
    site_pdb = driver.mini_site_poses_to_string(result, "pdb")
```

## 9. Mini 结果

通过 argv 启用：

```python
parsed = wv.parse_argv([
    "-c", "vina.in",
    "-r", "receptor.pdb",
    "-l", "ligand.sdf",
    "--mini", "site",
    "--mini_optimizer", "lbfgs",
    "--mini_steps", "80",
    "--mini_forcefield", "charmm36",
])
result = wv.DockingDriver.from_parse(parsed).run()
```

`MiniRelaxSummary`：

| 字段 | 类型 | 含义 |
|---|---|---|
| `mode` | str | `ligand`/`site` |
| `optimizer` | str | 实际使用的优化器 |
| `worker_count` | int | 并行 worker 数 |
| `elapsed_seconds` | float | 总耗时 |
| `poses` | list[`MiniPoseStats`] | 逐 pose 统计 |
| `setup_log` | str | 初始化日志 |
| `pose_logs` | list[str] | 逐 pose 日志 |
| `parameterization` | `ParameterizationSummary` | 与 CLI 同源的参数来源、覆盖率和退化摘要 |
| `site_pose_count` | int | site 模式输出的 pose 数（只读属性） |

每个 `MiniPoseStats`（均为 float）：

| 字段 | 含义 |
|---|---|
| `initial_energy` / `final_energy` | MM 前后能量 |
| `rms_displacement` / `max_displacement` | 整体位移 |
| `ligand_rms_displacement` / `ligand_max_displacement` | 配体位移 |
| `site_rms_displacement` / `site_max_displacement` | 位点（柔性残基）位移 |

`--mini full`、batch mini 和 multi-ligand mini 当前未启用。

`ParameterizationSummary` 同时用于 mini 和 FHFT。`provenance` 给出`producer/method/version/source`；`bonded_coverage` 给出 atom、bond、angle、torsion、improper 和 CMAP 的 template/generic/fallback/pruned 计数。摘要自身还暴露`atoms_total`、`atoms_explicit`、`atoms_forcefield`、`atoms_nonstandard`、`atoms_element_safety_net`、`atoms_unresolved`、`diagnostics`、`fallbacks` 和`degraded`。这些字段直接来自 watmol `Result/Diagnostics`，不是 Python 端重新推断。

## 10. FHFT 水预测

使用 CLI 参数构造 driver 后直接预测：

```python
parsed = wv.parse_argv([
    "-r", "receptor.pdb",
    "-c", "vina.in",
    "--predict_water",
])
if not parsed.success:
    raise RuntimeError(parsed.error_message)

driver = wv.DockingDriver.from_parse(parsed)
prediction = driver.predict_water()
print(prediction.grid.dimensions)
print(len(prediction.sites), prediction.diagnostics.converged)
```

默认使用 TIP3P、`4r` 蛋白-水介电、`8 x 12 x 6` C2 约化取向求积、0.5 Å FHFT 格点和自动标定的 TIP3P 网络化学势。通过 `parse_argv()` 可传入与 `watvina --help_water` 完全相同的覆盖参数，例如：

```python
parsed = wv.parse_argv([
    "-r", "receptor.pdb", "-c", "vina.in", "--predict_water",
    "--water_dielectric_mode", "sigmoidal",
    "--water_orientation_axis_order", "10",
    "--water_orientation_azimuth_count", "16",
    "--water_orientation_spin_count", "16",
])
```

TIP4P 当前没有 `auto` 网络化学势；选择 `--water_model tip4p` 时必须同时提供经过独立标定的 `--water_network_mu`。

`WaterPredictionResult` 暴露：

| 字段 | 类型 | 含义 |
|---|---|---|
| `sites` | list[`HydrationSite`] | NMS 提取的水合位点，按网络边际占据排序 |
| `receptor_atoms` | int | 受体原子数 |
| `forcefield_atoms` | int | 成功参数化的原子数 |
| `element_safety_net_atoms` | int | 走元素兜底参数的原子数 |
| `metal_steric_atoms` | int | 加入 FHFT steric 排除的金属数 |
| `metal_charge_fallback_atoms` | int | 电荷回退到 formal/氧化态的金属数 |
| `receptor_net_charge` | float | 受体净电荷 |
| `forcefield` / `charge_model` / `water_model` | str | 实际使用的模型名 |
| `parameterization` | `ParameterizationSummary` | 与 CLI 相同的 provenance/coverage/退化摘要 |
| `timings` | `WaterPredictionTimings` | 分阶段耗时 |
| `grid` | `FhftGrid` | 格点规格（只读属性） |
| `occupancy` | list[float] | 标量占据场，x 最快展平（只读属性） |
| `mean_field` | list[float] | 平均场（只读属性） |
| `effective_energy` | list[float] | 有效能场（只读属性） |
| `diagnostics` | `FhftDiagnostics` | 求解诊断（只读属性） |

`FhftGrid`：

| 字段 | 类型 | 含义 |
|---|---|---|
| `origin` | tuple[float, float, float] | 格点原点（Å） |
| `dimensions` | tuple[int, int, int] | `(nx, ny, nz)` |
| `spacing` | float | 体素边长（Å） |
| `size` | int | 总体素数（只读属性） |

`FhftDiagnostics`：

| 字段 | 类型 | 含义 |
|---|---|---|
| `converged` | bool | 是否收敛 |
| `iterations` | int | Picard 迭代数 |
| `bulk_occupancy` | float | 体相占据 ρ_bulk |
| `coarse_volume` | float | 粗粒化体积 |
| `chemical_potential` | float | 化学势 |
| `discrete_kernel_integral` | float | 离散核积分 |
| `structure_factor_zero` | float | 零波数结构因子 |
| `isothermal_compressibility` | float | 等温压缩率 |
| `calibrated_core_epsilon` | float | 标定的 core ε |
| `max_change` / `max_residual` | float | 收敛判据终值 |
| `grand_potential` | float | 巨势 |
| `message` | str | 诊断消息 |

`WaterPredictionTimings` 暴露 `parameterization_seconds`、`external_field_seconds`、`solve_seconds` 和 `site_extraction_seconds`
（均为 float）。FHFT 核心状态只有 `occupancy` 一个标量占据场；显式刚性水取向和边际占据在随后的 GCMC/GCI 网络步骤中求得。

`metal_steric_atoms` 是按常见氧化态离子半径加入 FHFT steric 排除的金属数；
`metal_charge_fallback_atoms` 只统计原始 partial charge 为零、因而回退到 formal charge 或常见氧化态的金属。已有非零模板电荷不会被覆盖。

`HydrationSite` 字段：

| 字段 | 类型 | 含义 |
|---|---|---|
| `position` | tuple[float, float, float] | 氧坐标（Å） |
| `occupancy` | float | GCMC 网络边际占据 [0,1] |
| `excess_occupancy` | float | 相对体相的超额占据 |
| `basin_excess_waters` | float | 盆地内超额水分子数 |
| `effective_energy` | float | 相对体相 Euler 场，**不是**位点结合自由能 |
| `grid_index` | int | 在扁平占据场中的下标 |

`format_water_pdb()` 输出每个位点的氧和 MAP 取向氢。PDB occupancy 列是GCMC 网络边际占据概率；B-factor 列为`-kBT log(p/(1-p))` 的网络自由能表示，不是绝对结合自由能。距离小于2.3 Å 的高占据候选可以同时列出，但它们在网络中只能互斥占据，表示同一高密度区域的替代水位置。

格式化与验证：

```python
dx_text = wv.format_water_open_dx(prediction)
pdb_text = wv.format_water_pdb(prediction)
metrics = wv.validate_water_prediction(
    prediction, "reference_waters.pdb", cutoff=1.4
)
print(metrics.precision, metrics.recall, metrics.average_precision)
```

`WaterMatchMetrics` 字段：

| 字段 | 类型 | 含义 |
|---|---|---|
| `cutoff` | float | 匹配阈值（Å） |
| `predicted` / `reference` | int | 预测/参考位点数 |
| `matched` | int | 成功匹配数 |
| `precision` / `recall` | float | 精度/召回 |
| `mean_distance` | float | 匹配对的平均距离 |
| `average_precision` | float | AP（排序质量） |

这个 Python helper 按给定 PDB 中的全部参考水计算，不会自动按 docking box过滤；调用者应确保参考水集合与预测区域一致。

## 11. 离散水 docking

Python 与 CLI 使用同一条简化路径：

```python
parsed = wv.parse_argv([
    "-r", "receptor.pdb", "-l", "ligand.sdf",
    "-c", "vina.in", "-w", "waters.pdb",
    "--water_bias_weight", "0.1",
])
if not parsed.success:
    raise RuntimeError(parsed.error_message)

result = wv.DockingDriver.from_parse(parsed).run()
```

该路径只加入可导的水位点置换偏置，没有显式水桥或完整 FHFT pose 重排。`PoseOutput.energy` 包含搜索目标中的 WaterBias；报告的 `Affinity` 仍排除它。

## 12. 性能与并行

**GIL 已释放。** `run()`、`score_only()`、`predict_water()` 和`eval_at_pose()` 在执行期间释放 GIL，因此可以用线程池并行多个独立driver（各自独立持有受体副本，无线程安全问题）：

```python
from concurrent.futures import ThreadPoolExecutor

def dock_one(ligand_path):
    parsed = wv.parse_argv(["-c", "vina.in", "-r", "receptor.pdb",
                            "-l", ligand_path, "--seed", "42", "--cpu", "2"])
    return wv.DockingDriver.from_parse(parsed).run()

with ThreadPoolExecutor(max_workers=6) as pool:     # 6 × cpu=2 ≈ 12 核
    results = list(pool.map(dock_one, ligand_paths))
```

**实测开销分布（demo/1a3l，12 核桌面机）：**

| 阶段 | 实测耗时 | 说明 |
|---|---|---|
| docking（exhaustiveness=2, cpu=4） | 数秒 | 与 `exhaustiveness × population × ga_search` 近似成正比 |
| FHFT 参数化 | 0.5 s | `timings.parameterization_seconds` |
| FHFT 外部场构建 | ~42 s | `timings.external_field_seconds`，取向求积主导，随盒体积和求积阶数增长 |
| FHFT 求解 | ~2 s | `timings.solve_seconds` |
| 位点提取 | <0.1 s | `timings.site_extraction_seconds` |

实践要点：

- `--cpu` 在**搜索任务之间**并行，单个简单配体吃不下所有核；批量筛选时「每 driver 少核 + 多 driver 并行」比「单 driver 满核」吞吐更高。
- `--seed 0`（默认）= 非确定（entropy ^ pid ^ 时间戳）；排名可复现性要求时显式设置种子。
- FHFT 外部场是水预测的成本中心；若同一受体要反复预测（如参数扫描），应缓存 `prediction` 或其文件输出（`format_water_pdb` 可无损回读为
  `-w` 输入，见 §7.4），而不是重复调用 `predict_water()`。
- 每个 `DockingDriver` 构造都会重做受体参数化（surface mask、directional cache），批处理时这部分开销无法跨 driver 共享。

## 13. API 总览

模块函数：

```text
parse_argv(argv)
format_water_open_dx(prediction)
format_water_pdb(prediction)
validate_water_prediction(prediction, reference_pdb, cutoff=1.4)
```

公开类：

```text
DockingConfig ParseResult DockingDriver DockingResult
PoseOutput PoseConf
MiniPoseStats MiniRelaxSummary Provenance ParameterCoverage ParameterizationSummary
FhftGrid FhftDiagnostics HydrationSite
WaterPredictionTimings WaterPredictionResult WaterMatchMetrics
```

`DockingDriver` 方法：

```text
from_parse(parse_result)
from_strings(config, receptor, receptor_format, ligands, ligand_formats, flex="")
run()
score_only()
predict_water()
eval_at_pose(pose)
poses_to_string(result, format="sdf", num_modes=10, energy_range=3.0)
mini_site_poses_to_string(result, format="pdb")
flex_poses_to_string(result, format="pdb", num_modes=10, energy_range=3.0)
```

只读属性：`size`、`chem`。`chem` 是内部 ChemicalSystem 的只读引用，不应作为稳定的 Python 化学编辑 API 使用。

## 14. 错误处理与排错

所有长任务都通过结果对象报告业务失败：

```python
if not parsed.success:
    raise ValueError(parsed.error_message)

result = driver.run()
if not result.ok:
    raise RuntimeError(result.error_message)
```

输入解析、格式输出和 FHFT 不收敛等情况也可能直接抛出 Python exception，应在批处理调用点捕获：

```python
try:
    result = driver.run()
except Exception as error:
    print(f"WATVina failed: {error}")
```

**`ImportError: dynamic module does not define module export function`**或 `undefined symbol: PyUnicode_...`：Python 主次版本与扩展 ABI 不匹配
（如用 3.11 导入 cpython-312 产物）。用与构建一致的解释器，或从源码重建。

**`ImportError: libstdc++.so.6: version 'GLIBCXX_3.4.32' not found`**：系统 libstdc++ 太旧。换用更新的系统/conda 环境，或用`make linux_glibc_release` 思路在目标 ABI 下自行构建扩展。

**`parse success=False, error_message 提到未知选项`**：argv 中传了CLI 不存在的选项或重复传参；对照 `watvina --help` 修正。

**`eval_at_pose` 抛异常**：driver 含多个配体（联合 docking）。该接口当前只支持单配体。

**`predict_water` 返回不收敛**：检查 `prediction.diagnostics.message`；金属/异常化学环境可能需要调整 `--water_ff` 或显式 `--water_network_mu`。

## 15. 当前边界

- 没有老版 `WATVina` class、setter 集合或 Python batch helper。
- Python 没有直接暴露 score channel breakdown、Cartesian pose arrays 或 torsion setter。
- `DockingConfig` 只暴露常用字段；完整功能走 `parse_argv()`。
- `--mini full`、batch/multi-ligand mini 未启用。
- TIP4P 水网络没有自动化学势标定值，必须显式提供`--water_network_mu`。
- 第三方显式水的方向和 flex water-direction refresh 尚未完成。
- 水位移偏置不包含显式水桥或逐 pose FHFT 重排，也不属于报告的`Affinity`。
- 当前未提供 wheel；Linux `.so` 和 Windows `.pyd` 都受平台与 CPythonABI 约束。
