# BMI3 Group 1 前端 Demo 设计与实现 Prompt

你是一名熟悉生物信息学可视化的全栈工程师。请在当前 BMI3 Algorithms ICA Group 1 仓库中加入一个可本地运行的 **Protein Conservation Explorer**。它面向课堂演示，也应当像一个小型生信工具一样能接收蛋白 FASTA、执行实际计算，并且让人从逐位点分数追溯到原始比对。借鉴 UniProt 和 SMART 面向生物研究者的结果组织方式，设计属于本项目的浅色科研工作台。

## 1. Analyze · 上传与计算

- 首页主操作为 **Upload FASTA**，支持文件选择和拖放。旁边提供 **Load example**，内置一个可运行的 Pattern C 已对齐示例，以及一个用于展示比对流程的 globin 未对齐示例。加载示例后仍可点击运行真实计算。初始状态简短说明两种输入各自适用的场景。
- 让用户明确选择 **Aligned FASTA → Score alignment** 或 **Unaligned FASTA → Align & score**。前者复用现有解析器校验对齐矩阵和序列 ID，直接评分；后者先运行 MAFFT，环境安装 ClustalW 时也可选择它。以用户选择与输入校验决定流程。显示文件名、序列数、输入模式和已选方法，主按钮随流程切换。
- 评分默认使用 **SCAW-RVNE**，提供 **Original RVNE**。使用冻结 benchmark 参数时标明 **Frozen benchmark settings**，把参数放进可展开的 **Advanced settings**。运行时区分 Reading FASTA、Aligning 和 Scoring 等实际阶段或显示当前可确认的状态。失败信息区分文件格式、缺少比对工具和评分失败，并允许从当前输入重试。

## 2. Results · 比对、保守性曲线与热图

- 计算后在 **Analyze** 内打开结果工作台。顶部紧凑摘要显示文件、序列数、最终对齐长度、可评分列比例、实际比对工具（如有）和评分方法。提供 **New analysis**，并允许下载对齐 FASTA、逐位点 TSV 和可用的分段 TSV 或参数 JSON。
- 主区是带固定序列 ID 列、残基字符、gap 标记和位置标尺的 **MSA viewer**。初始只显示能读清的局部窗口与有限行数，横向平移对齐列，纵向滚动查看序列。残基配色按化学性质分组，字符和 gap 即使在投影上也可辨。标尺使用从 1 起算的 alignment column；有明确映射时，可另标参考残基号。
- 共享同一列坐标的 **Conservation profile** 用折线显示实际逐位点分数。提供 **Line / Heatmap** 视图切换，并允许用热图作为贴在 MSA 上方的紧凑分数轨。分数热图每列对应一个比对列，每行对应一种已计算的方法，提供明示范围的色标。对同一次 SCAW 计算，支持切换原始 RVNE 与最终分数。方法分数尺度差异明显时使用独立标尺或分面，并明确标出比较范围。RVNE 主线分数越高表示越保守。未评分列在折线中保留断点，在热图中用单独样式标为 Unscored。
- 用全局概览和可拖动选区定位任意区间。缩放、平移、点选或键盘移动都同步更新 MSA、折线和热图。点击某列固定选中，详情显示列号、残基计数或比例、非 gap 覆盖率、原始与最终分数，并在结果提供时显示窗口宽度和分段范围。SCAW 拟合片段与合成真值各有独立图层及开关。合成案例可查看已知 motif 或突变率区间，自定义上传结果按实际可用数据绘制。

## 3. Benchmarks · 测试与消融

- 将 **Benchmarks** 设为独立子页面，支持按数据集、Pattern A–D、方法和指标筛选。用一张横向点图或条形图显示比较，用表格核对精确值、统计口径、来源和实际贡献的家族数 n。优先显示 Spearman、variable-site recall 和 boundary F1，再提供 AUROC、AUPRC、可评分比例，以及已有耗时结果。
- Pattern C 展示短保守 motif 和固定窗口平滑的关系，同时呈现 E7 或 Shannon 对该 motif 的检出结果。Pattern D 展示 **E7 retuned RVNE、E8 retuned SCAW-RVNE、Shannon** 在中等突变率位点上的检出差异，以及 **E5 segmentation only** 在该组的 Spearman。说明 D 的二元评估将 medium 与 high 合并为 variable。指标无定义时显示 **N/A** 和原因，例如 Pattern A 的每个家族只有一个类别，按本项目评估口径，AUROC、AUPRC 显示 N/A。结论应注明对应的合成测试集和指标。
- 已有真实 globin 对齐产物时，可在 **Case study** 区比较 MAFFT 和 ClustalW 对同一组序列的影响。两条曲线按 P02144 的 ungapped reference residue 对齐，并标明来源、共同可评分覆盖率和所用的比对工具。归类为对齐敏感性示例。

## 4. 视觉与布局

- 采用浅色科研工作台。背景用 #F5F7FA，图表背景为白色，正文用 #162536，分隔线用 #D8E1E8。SCAW 主曲线用深青 #0F766E，原始 RVNE 用蓝灰虚线，分段边界用琥珀，真值使用独立轨道或纹理。各图层同时用颜色与线型、标记或标签区分。界面用 Inter 或本地无衬线字体，残基、坐标与数值用 IBM Plex Mono 或本地等宽字体。
- 紧凑顶部导航只保留项目名、**Analyze** 与 **Benchmarks** 两个主要入口。结果页优先留空间给 MSA、折线或热图，位点详情在宽屏放右侧，小屏放下方。按 1280×720 和 1440×900 检查投影可读性。默认界面使用英文。控件提供可读标签与键盘焦点，图表图例、坐标和色标清晰可见。

## 5. 仓库接入与验收

- 先读 `code/rvne/README.md`。已对齐评分复用 `python -m msa_conservation.cli score`，未对齐流程复用 `python -m msa_conservation.cli real-data` 或 `code/rvne/src/msa_conservation/real_data.py` 的现有调用路径。两套代码使用相同的 `msa_conservation` 包名，运行主线时将模块路径指向 `code/rvne/src`。用轻量本地接口处理上传、调用 Python 与返回结果，前端读取结果并渲染。优先复用仓库已有技术栈，若没有前端栈，可用 React + Vite 配轻量 Python API。
- 使用逐列输出中的 `position`、`raw_rvne`、`nongap_coverage`、`score`、`window_size` 和分段字段。冻结配置读取 `code/rvne/results/reference_run/parameters.json`。示例可选 Pattern C seed 50124，配对的 aligned FASTA 和 truth TSV 由 `code/data_generation/generate_pattern_data.py` 生成。实时样例与预计算结果都注明来源。测试数据读取 `code/rvne/results/reference_run/` 和 `code/rvne/results/supplemental_run/` 的现有 TSV。缺失产物通过现有脚本生成，或在界面显示明确状态。
- 交付可运行的 Demo 与简短启动说明。验收时实际上传已对齐 FASTA 获得真实分数，再用未对齐 globin 示例获得比对和评分；核对选中列在矩阵、折线、热图和详情中一致；检查非法输入和缺少 MAFFT 时能收到明确反馈；最后在 Benchmarks 核对 Pattern D 的指标可追溯到 TSV，并确认 Pattern A 的 AUROC/AUPRC 显示 N/A。默认演示与示例资源应支持无外网运行，未对齐分析依赖本地安装的比对工具。
