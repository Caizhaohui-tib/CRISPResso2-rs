# CRISPResso2 Rust 重写分析与迁移计划

## Summary

当前仓库是 Python CLI 套件，入口定义在 `CRISPResso2-master/pyproject.toml`，包含 `CRISPResso`、`CRISPRessoBatch`、`CRISPRessoPooled`、`CRISPRessoWGS`、`CRISPRessoCompare`、`CRISPRessoPooledWGSCompare`、`CRISPRessoAggregate` 七个命令。

核心结构：

- `CRISPRessoCORE.py`：单样本核心分析，负责参数校验、FASTQ/BAM 预处理、读段比对、indel/substitution 分类、计数、TSV/NPZ/VCF 输出、绘图和报告调度。
- `CRISPRessoShared.py`：共享参数系统、序列工具、JSON 序列化、guide/quantification window 计算、报告辅助、guardrails。
- `CRISPResso2Align.pyx`：Needleman-Wunsch affine gap 全局比对，带 cut-site gap incentive，是最适合优先迁移的热点。
- `CRISPRessoCOREResources.pyx`：从 read/ref alignment 中提取 insertion/deletion/substitution，是第二个优先迁移热点。
- `CRISPRessoBatch/Pooled/WGS`：主要是编排器，读取批量/区域/amplicon 输入，调用 Core 子任务，汇总结果。
- `plots/` 和 `CRISPRessoReports/`：matplotlib/seaborn/Jinja 报告层，工作量大但不属于第一阶段必改。

推荐采用“Rust 核心先行 + Python 报告暂留 + 输出逐字节优先兼容”的路线。

## Key Changes

- 新建 Rust workspace，建议 crate 拆分为：
  - `crispresso-core`：序列、参考、guide、quantification window、alignment、variant classification、summary aggregation。
  - `crispresso-cli`：兼容现有 CLI 参数、目录名和输出文件名。
  - `crispresso-py`：可选 PyO3 绑定，让现有 Python 编排层能调用 Rust 核心，便于渐进迁移。
- 当前已落地：
  - `crispresso-core` 与 `crispresso-cli` 已实现。
  - `crispresso-py` 尚未开始，仍适合作为 Core API 稳定后的后续项。
- 第一阶段只替换可计算核心：
  - Rust 实现 `read_matrix/make_matrix/global_align`，严格复刻 `CRISPResso2Align.pyx` 的 tie-breaking、gap incentive、末端 gap 行为和分数格式。
  - Rust 实现 `find_indels_substitutions` 与 legacy 版本，输出字段保持与当前 payload 一致。
  - Rust 实现 FASTQ 读取、unique-read cache、variant payload 生成和 count vector 聚合。
- 当前验证结果：
  - 最小 Core 路径 `CRISPResso -r1 ... -a ... -g ...` 已跑通。
  - `CRISPResso_quantification_of_editing_frequency.txt` 与 `Nucleotide_frequency_table.txt` 已在 `FANC.Cas9` 数据集上与 Python golden output 对齐。
- 保留现有 Python 报告层：
  - Rust 先产出 `CRISPResso2_info.json`、`Alleles_frequency_table.zip`、`CRISPResso_quantification_of_editing_frequency.txt`、nucleotide/modification tables。
  - Python `CRISPRessoReport.py` 和 `plots/CRISPRessoPlot.py` 暂时继续消费这些输出。
- 外部工具策略：
  - 初期继续调用 `fastp`、`samtools`、`bowtie2`，保持行为兼容。
  - 后续再评估用 `rust-htslib` 替换 BAM/SAM 处理，用 Rust-native FASTQ/FASTA 库替换 shell 管道。
- 数据模型建议：
  - 明确定义 `ReferenceInfo`、`GuideInfo`、`VariantPayload`、`AlignmentResult`、`AlleleRow`、`RunSummary`。
  - 避免直接复制 Python 的动态 dict；Rust 内部用强类型，序列化层映射成现有 JSON/TSV 形状。

## Migration Order

1. 建立黄金测试基线：运行现有 pytest 和 `tests/Makefile` 中的 CLI 集成样例，保存关键输出作为 Rust parity 标准。
2. 迁移 Cython 热点：先做 alignment 和 indel/substitution extraction，并用现有 unit tests 对齐。
3. 迁移 Core 单样本最小路径：`CRISPResso -r1 ... -a ... -g ...`，先覆盖 `tests/CRISPResso_on_FANC.Cas9`。
4. 扩展 Core 参数面：HDR、base editor、prime editing、paired-end、BAM input、VCF output、legacy insertion quantification。
5. 迁移编排器：Batch 最先，因为它主要生成 Core 命令；Pooled/WGS 后移，因为依赖 Bowtie2/Samtools 和区域拆分。
6. 最后处理报告层：只有在 Rust 核心输出稳定后，再决定是否重写 plotting/report。

### Current Status

- 第 1 步已完成：当前直接使用 `tests/expectedResults/CRISPResso_on_FANC.Cas9/` 作为 parity 基线，无需在 Windows 上重新构建 Python Cython 扩展即可完成对比。
- 第 2 步已大体完成：`align.rs` 与 `edits.rs` 已实现，且已补充一批镜像 Python 行为的测试。
- 第 3 步已完成：Rust CLI 已支持 `--fastq_r1/--amplicon_seq/--guide_seq` 的最小工作流，并在 `FANC.Cas9` 上得到与 Python 一致的关键输出。
- 当前优先级已进入“第 2 步收尾 + 第 4 步准备”阶段，接下来优先补测试、固化 golden 对比、再扩参数面。

### Next Priorities

1. 继续补齐 Python 单测迁移：优先移植 `test_CRISPRessoCOREResources.py` 与 `test_CRISPRessoShared.py` 中与 Core 计算直接相关的测试。
2. 把 `FANC.Cas9` golden 对比固化为自动化测试：将 `CRISPResso_quantification_of_editing_frequency.txt` 与 `Nucleotide_frequency_table.txt` 的对比纳入 Rust 测试或集成测试流程，避免回归。
3. 开始扩展参数面：建议先做 `legacy insertion quantification`、`paired-end`、`BAM input`、`VCF output`，再进入 HDR、base editor、prime editing。

### Lessons Learned

- Python 单元测试参数与 Python CLI 默认参数并不完全相同：
  - Python 对齐单测常使用 `EDNAFULL` 与 `gap_open=-1/gap_extend=-1`。
  - Python CLI 实际默认路径使用 `make_matrix(5, -4, -2, -1)` 与 `gap_open=-20/gap_extend=-2`。
  - 因此 parity 不能只对齐 unit test，还必须对齐真实 CLI 输出。
- `default_min_aln_score=60` 是影响 `FANC.Cas9` 结果一致性的关键行为：Rust 若默认用 `0.0`，会把 250 条 reads 全部计为 aligned；修正为 `60.0` 后，与 Python 的 235 aligned 完全一致。
- 输出兼容不仅是算法兼容，也包括格式兼容：百分比字段需要按 Python 一样保留 8 位小数，否则即使统计值正确，golden file 仍无法逐字节匹配。
- 在当前 Windows 环境下，直接复用已有 golden output 比重新构建 Python Cython 环境更高效，也更适合作为第一阶段 parity 验证手段。

## Test Plan

- 单元测试：
  - 直接移植 `test_CRISPResso2Align.py`、`test_CRISPRessoCOREResources.py`、quantification window、VCF left normalization、guide parsing 测试。
- Golden file 测试：
  - 对比 `tests/expectedResults` 里的 quantification、nucleotide frequency、summary 文件。
  - 浮点输出先按逐字节兼容处理；如发现 Python/numpy 平台差异，再明确允许误差的字段。
- CLI 兼容测试：
  - 覆盖 `CRISPResso`、`CRISPRessoBatch`、`CRISPRessoPooled`、`CRISPRessoWGS`、`CRISPRessoCompare` 示例命令。
- 性能测试：
  - 单独 benchmark global alignment、variant classification、FASTQ unique-read cache。
  - 用相同 FASTQ 比较 Python/Cython 版本与 Rust 版本的运行时间和峰值内存。

### Test Status

- 已完成：
  - `align.rs` 中补充了与 Python 行为对齐的测试，包括 EDNAFULL 场景。
  - `FANC.Cas9` 的 quantification 与 nucleotide frequency 已完成 golden 对比验证。
- 待补齐：
  - 系统移植 `test_CRISPRessoCOREResources.py`。
  - 系统移植 `test_CRISPRessoShared.py` 中与 guide、窗口、共享计算逻辑相关的测试。
  - 将 `FANC.Cas9` golden 对比改造成可重复执行的自动化测试，而不是仅作为手工验证步骤。

## Assumptions

- 第一版 Rust 重写以结果兼容优先，不重新设计用户可见输出。
- HTML 报告和绘图先继续使用 Python 层，避免把算法迁移和视觉重写混在一起。
- `fastp`、`samtools`、`bowtie2` 初期仍作为外部依赖保留。
- 当前仓库不是 git 工作区；实施前建议先初始化新分支或新 Rust workspace，避免和原始 Python 源码混乱。

## Contributors

- 项目负责人：用户本人
- AI 协作贡献者：Claude Sonnet 4.6
- AI 协作贡献者：GPT-5.4
