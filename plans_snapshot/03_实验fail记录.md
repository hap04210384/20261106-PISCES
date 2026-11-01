# 实验 Fail / Blocked 记录（夜间自主执行，2026-10-31）

> 纪律：本文件记录环境受限或失败项，裁决与替代方案一并留痕。所有时间均为本地时间。

## B-01 CUDA 工具链获取 blocked（nvcc 不可用）
- **现象**：pip 版 `nvidia-cuda-nvcc-cu12`（12.9.86 / 12.8.93 均验证）Windows wheel 只含 `ptxas.exe`，**不含 nvcc 前端**（与 Linux wheel 不同）。
- **官方安装包下载被阻断**：`developer.download.nvidia.com` 按地理位置 301 重定向至 `developer.download.nvidia.cn`，后者对 12.4.1/12.6.x 安装包一律 404；目录索引亦不可用。winget 的 winget 源更新失败（网络不可达），仅 msstore 源可用但无 CUDA。
- **影响**：E-A03（TPDS 引擎三配置 route-independence）无法在 GPU 上运行；E-B04 的"单机单卡"判据降级为 CPU 单路。
- **裁决（自主）**：GPU 实验记 blocked，不在本机伪造。CPU-only 继续全部 MVP 主线实验（128GB 内存、i9-12900 足够 2000 基因组 × 5.5 万特征规模）。nvcc 获取待明日用户网络环境（或手机热点）重试，或改用 WSL2/Colab 复现 GPU 部分。
- **佐证**：`code/envs/cuda-pip/` 仅 `nvidia/cuda_nvcc/bin/ptxas.exe`；`wcheck/` 内留 12.8.93 wheel 备查。

## B-02 GTDB 元数据直连弃用（网速不可行）
- **现象**：`data.gtdb.ecogenomic.org` 实测 ~57 B/s（280 s 仅 16 KB/289 MB）。
- **裁决（自主）**：GTDB 元数据非 MVP 关键路径——PanGBank 已提供 K. pneumoniae 预计算泛基因组（GTDB_refseq R232 集合，2001 基因组，54590 基因家族，完整性/污染度均质控），直接满足 E-B01/E-B02 数据需求。GTDB 元数据下载放弃，不改实验设计。

## B-03 PanGBank 前端不可爬，改用官方 REST API
- **现象**：Web 前端为 Reflex SPA，无静态 API 路由；经源码检索找到官方 CLI（labgem/PanGBank-cli），由其源码确认 REST API 基地址 `https://pangbank-api.genoscope.cns.fr/`。
- **结果**：`/pangenomes/?pangenome_name=s__Klebsiella_pneumoniae` 命中 2 条记录；选定 GTDB_refseq 2.0.0 版（id=11756，2001 基因组）为 MVP 数据集，GTDB_all（id=15033）留作扩展。
- **下载限速**：单连接 ~33 KB/s，1.12 GB 文件需 9+ 小时；裁决用 16 路 HTTP Range 并行分块下载（聚合 ~500 KB/s，约 30–45 分钟）。

## B-05 基线环境修补记录（Windows venv 内，不入 git）
- **pyseer 1.3.11**：`pysam`、`glmnet_python/cvglmnet*` 无 Windows 版，已打桩（桩仅在使用 VCF/enet 代码路径时抛错，本研究只用 TSV + FE/LMM，不受影响）。**投稿前必须在 Linux/WSL2 重跑基线复验**。
- **scoary 1.1.4**：两处兼容性修补——`open(...,"rU")` → `"r"`（Python 3.11+ 移除 rU 模式）；`scipy.stats.binom_test` → `binomtest(...).pvalue`（SciPy ≥1.12 移除 binom_test）。修补仅涉及文件打开模式与检验函数名，统计逻辑未动。
- **下载进程中断**：16 路并行 curl 在 18:51 前后全部静默退出（原因未查明，疑为长时间连接被中间设备断开）；已用 `curl -C -` 断点续传重启全部 19 个分块，续传有效（697 MB → 持续增长）。

## B-04 数据方案替换留痕（E-B01/E-B02 设计变更）
- **原计划**：GTDB 清单 → NCBI 下载基因组 → Panaroo 泛基因组。
- **实际执行**：PanGBank 预计算 PPanggolin 泛基因组（persistent/shell/cloud 分区基因家族）替代自注释流程。
- **理由**：① 一晚上跑不动数百基因组下载+注释流水线；② PanGBank 基因家族即 OG 等价特征（L1 层），persistent/shell/cloud 分区提供额外先验；③ 变更只影响特征构建工程路径，不影响方法学主张（输入均为 genome × gene-family presence/absence 矩阵）。
- **风险记录**：审稿人可能质疑 Panaroo vs PPanggolin 差异 → 回复策略：方法对特征构建黑箱不可知（任何 presence/absence 矩阵均可输入），扩展实验可补 Panaroo 对拍。
