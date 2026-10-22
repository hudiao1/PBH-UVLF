# PBH-UVLF

[English](README.md) | 中文

本仓库提供研究原初黑洞（PBH）与星系紫外光度函数（UVLF）所用的分析代码和支持数据，覆盖七个红移节点：z = 6、7、8、9、10、11 和 12.5。

## 文件内容

| 路径 | 内容 |
| --- | --- |
| `mcmc/` | 作者提供的 16 参数模型、MPI/emcee 采样程序、UVLF 似然输入数据和尘埃校准表 |
| `plot/UVLF_full_WithPBH_noRobertson_20261016.ipynb` | 使用已有理论曲线缓存绘制七个红移面板的 UVLF 图 |
| `plot/prl_free_sf_accretion_rate_v22.ipynb` | 绘制作者提供的六个红移箱的爱丁顿比约束及对应吸积率 |
| `plot/statistics and plot_F_216_consistent_MCSamples_shortest68.ipynb` | 作者提供的后验分析 notebook，已调整为读取便于分发的样本文件 |
| `plot/*.txt` | 文献观测数据台账、理论曲线缓存，以及绘图输入和导出数据 |
| `posterior/posterior_samples.npz` | 按原 notebook 的处理规则选出的 6512 个等权重样本 |
| `posterior/metadata.json` | 原始链的校验值、处理设置、剔除的 walker 编号和数组维度 |
| `posterior/parameters.csv` | 后验数组的参数列顺序及定义 |
| `scripts/export_posterior.py` | 从完整 HDF5 链重新导出后验样本文件 |
| `figures/` | 运行所提供 notebook 后生成的 PDF 预览图 |

本次发布保留作者提供的科学计算设置。整理过程中没有重新运行 MCMC，也没有验证论文的统计主张。

## 安装

使用 Python 3.12 或兼容的更新版本，并创建独立环境：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

依赖版本清单记录了本次发布检查使用的环境，不能据此还原历史 MCMC 链生成时的软件环境。本次没有测试从空环境开始的完整安装流程。

## 重新绘图

```bash
cd plot
jupyter lab
```

打开相应绘图 notebook，按顺序运行所有单元。请将 notebook 的工作目录设为 `plot/`，其中的相对文件名均以该目录为基准。UVLF notebook 读取已有模型曲线缓存，不重新拟合。吸积率 notebook 绘制作者提供的数值摘要，不从后验样本中重新推断这些数值。程序指定 Times New Roman 字体；系统缺少该字体时，可能使用替代字体。

后验分析 notebook 读取 `../posterior/posterior_samples.npz`。原处理规则为：丢弃前 30000 步，剔除 burn-in 后发生移动的步数比例低于 0.02 的 walker，再按 emcee 的 `thin=1000` 规则抽稀。最终保留 88 个 walker，每个 walker 有 74 个采样点，每个采样点包含 16 个参数。完整分辨率的轨迹和自相关诊断需要原始链；不能把抽稀后的轨迹当作完整链已收敛的证据。

## 完整 MCMC 采样（可选，计算量较大）

先安装 MPI 实现，例如 Open MPI，再安装可选依赖：

```bash
python -m pip install -r requirements-mpi.txt
mpiexec -n 4 python mcmc/emcee_UVLF_P18_a16_merged_z20.py
```

脚本保留原先验和模型设置，使用 112 个 walker，每次调用运行 60000 步。新建 HDF5 后端时，程序从给定初始化范围生成初始位置；后端已有链时继续采样，不重置原链。输出文件为 `mcmc/emcee_UVLF_P18_a16_merged_z20.h5`。脚本没有固定随机种子，因此新链不会逐字节复现历史样本。本次发布整理没有执行这项耗时采样。

完整原始链约 1.9 GB，本 GitHub 仓库没有包含该文件。`posterior/metadata.json` 记录了它的 SHA-256 校验值和样本筛选细节。仓库提供处理后的后验样本，不包含完整链诊断所需的全部数据，也没有声明完整链的公开下载地址。

## 数据来源与单位

拟合输入文件 `mcmc/UVLF_PBH.txt` 包含六列：红移、M_UV（AB mag）、星等箱宽（mag）、Phi_UV（Mpc^-3 mag^-1）、上侧误差和下侧误差。所提供的似然程序会施加最低 20% 的误差，并执行其中定义的尘埃及宇宙学转换；请保留原输入表。文件表头标注了 HST 数据来源（arXiv:2102.07775），高红移数据块的来源仍需确认。

`mcmc/Beta_parameters.txt` 提供尘埃修正所用的紫外光谱斜率校准，并注明来源为 arXiv:1306.2950 的表 3。`plot/UVLF_literature_data.txt` 和导出的 `*_data.txt` 逐行记录文献、原表位置、红移区间、误差定义、上下限类型和显示标记。仅用于比较的观测点没有额外加入似然。不同文献的观测宇宙学和选择函数尚未全部统一。

理论数据表中，M_UV 的单位为 AB mag，数密度单位为 Mpc^-3 mag^-1；表头区分 `phi_model` 和 `phi_observed`。后验参数的对数以 10 为底；`log_prob` 是 emcee 保存的对数后验，不是独立重新计算的对数似然。参数顺序见 `posterior/parameters.csv`。所提供吸积率图采用的转换关系为 mdot = 9.73 lambda_E。

## 引用与获取

仓库地址：[hudiao1/PBH-UVLF](https://github.com/hudiao1/PBH-UVLF)。

目前没有添加论文引用或归档 DOI，因为作者尚未提供已核验的发表标识或 DOI。复用文献测量时，请同时引用原始观测论文和本仓库。作者确定论文对应的最终版本后，可再创建带标签的发布版本，并存入长期归档仓库。

## 权利与许可证

本仓库暂未指定许可证。公开访问不等于授予通用的复用或再分发许可。作者编写的代码仍由相应权利人保留权利，复用授权请联系权利人。第三方数据和软件保留各自的权利及使用条件。来源引用用于记录出处，不能替代使用许可；不能假定整理后的文献数据表已经获得新的授权。
