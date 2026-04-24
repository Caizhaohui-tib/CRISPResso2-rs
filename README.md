# CRISPResso2-rs

`CRISPResso2-rs` 是一个使用 Rust 重写 `CRISPResso2` 核心分析流程的实验性项目，目标是在保持结果兼容的前提下，逐步替代原有 Python/Cython 热点模块。

当前仓库同时包含两部分内容：

- `crates/`：正在开发中的 Rust workspace
- `CRISPResso2-master/`：上游 Python 版本代码与测试数据，用于行为对照和 golden output 校验

## 项目目标

- 用 Rust 重写 CRISPResso2 的核心计算路径
- 优先保证输出结果与 Python 版本兼容
- 用自动化测试持续验证 Rust 与 Python golden output 的一致性
- 逐步扩展参数面，最终支撑更完整的 CRISPResso2 工作流

## 当前状态

当前已实现并验证的内容包括：

- Rust workspace：
  - `crispresso-core`
  - `crispresso-cli`
- 已迁移的核心模块：
  - 全局比对
  - indel/substitution 提取
  - guide 相关计算
  - FASTQ 读取与 unique read 计数
  - quantification / nucleotide frequency 输出
- 已支持的最小 CLI 路径：
  - `CRISPResso --fastq_r1 ... --amplicon_seq ... --guide_seq ...`
- 已完成的兼容性验证：
  - `FANC.Cas9` 数据集的 `CRISPResso_quantification_of_editing_frequency.txt`
  - `FANC.Cas9` 数据集的 `Nucleotide_frequency_table.txt`
  - 以上文件已通过自动化测试与 Python golden output 逐字节比对

## 仓库结构

```text
.
├── .cargo/
├── crates/
│   ├── crispresso-core/
│   └── crispresso-cli/
├── CRISPResso2-master/
├── Cargo.toml
├── Cargo.lock
└── README.md
```

## 本地开发

当前开发环境基于 Windows + Conda，Rust 测试在 `rust-env` 环境中执行。

运行测试：

```bash
conda run -n rust-env cargo test
```

运行当前最小 CLI：

```bash
conda run -n rust-env cargo run -p crispresso-cli -- \
  --fastq_r1 CRISPResso2-master/tests/FANC.Cas9.fastq \
  --amplicon_seq CGGATGTTCCAATCAGTACGCAGAGAGTCGCCGTCTCCAAGGTGAAAGCGGAAGTAGGGCCTTCGCGCACCTCATGGAATCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTCAAGCACTACCTACGTCAGCACCTGGGACCCCGCCACCGTGCGCCGGGCCTTGCAGTGGGCGCGCTACCTGCGCCACATCCATCGGCGCTTTGGTCGG \
  --guide_seq GGAATCCCTTCTGCAGCACC
```

## 已有测试保障

- Rust 单元测试覆盖 alignment、edits、guide、FASTQ、输出格式等核心行为
- CLI 侧包含 `FANC.Cas9` golden output 自动化回归测试
- 目前 `cargo test` 已通过

## 当前限制

以下能力尚未完整实现：

- paired-end 处理
- BAM 输入分析
- VCF 输出
- HDR / base editor / prime editing 全量支持
- Batch / Pooled / WGS / Compare 等编排层功能

其中部分参数当前已经接入 CLI，但会显式报出未实现错误，而不是静默忽略。

## 贡献者

- 项目负责人：用户本人
- AI 协作贡献者：Claude Sonnet 4.6
- AI 协作贡献者：GPT-5.4
