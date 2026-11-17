# 宇宙线输运与核反应

编译和运行（需要 C++20）：

```sh
make
./cr_solver model.cfg
./cr_solver
make test
```

所有运行参数统一写在 `model.cfg`。无参数运行读取当前工作目录的 `model.cfg`，
也可传入另一个完整配置文件。核素选择与源权重必须在文件中逐项指定，
列表顺序不限制反应顺序，但不能重复核素。
已有输出目录中的同名结果会被覆盖。

网格、风场、扩散参数、气体密度、源强和核素源丰度可通过 `model.cfg` 设置。
源空间形状由 `initialize_CR_source` 设置，是位于 R=z=10 kpc 的高斯分布，
谱形为 Rg^-2；这里只作为运行示例，源强没有做观测归一化。
初始密度为零，主程序将该形状作为持续源项注入。

## OpenMP

`make` 生成串行版 `cr_solver`。`make openmp` 生成 `cr_solver_omp`，
`make test-openmp` 编译和运行 OpenMP 测试（默认 4 线程）。
macOS 上 Apple Clang 需要 Homebrew 的 `libomp`，Linux 上默认使用 `-fopenmp`。
其他编译器可通过 `CXX`、`OMP_CXXFLAGS` 和 `OMP_LDLIBS` 指定。

```sh
make openmp
OMP_NUM_THREADS=4 ./cr_solver_omp model.cfg
```

源注入、空间对流和空间扩散按刚度平面分配线程。刚度演化和核反应按空间径向行分配线程，
核反应在每行内保持母核先于子核。运行时可通过 `OMP_NUM_THREADS` 限定线程数。

## 算子与单位

每个核素对应一个 `Field3D`，保存单元平均 q=dn/dRg。
空间单位 kpc、时间 yr、刚度 GV；D0 在 1 GV 处定义，配置项 `D0_cm2_per_s` 使用 cm²/s，程序在构造扩散模型前换算为 kpc²/yr；vA 输入 km/s。
现有空间对流、空间扩散、绝热能变与随机再加速模块保持原有接口。
每步顺序执行源注入、各核素的空间扩散/对流/动量演化，以及核反应。
核反应使用 backward Euler，因此整体分裂推进为时间一阶。
主程序在初始化风场后检查一次对流 CFL，使用单元总流出率乘以时间步长，
当前均匀网格、PLM + SSPRK2 采用保守上限 0.5。超限会报错并给出最大允许 dt_yr。
所有步长指标由 `configuration.cpp` 的 `calculate_cfl_numbers` 在初始化时计算，
并写入 `settings.txt`。定义为 dt 乘以实际离散算子的最大对角损失率，
包含空间扩散 R/z/总和、刚度扩散、绝热能变、合并动量算子，以及核反应 Γdt。
均匀笛卡尔内部单元的单轴扩散指标对应 2Ddt/dx²，需注意与 Ddt/dx² 的区别。
空间 CN–ADI 各轴指标 <=2、合并动量 CN 指标 <=2，分别是显式右端非负的充分条件，
这些隐式指标仅记录。核反应 backward Euler 的 Γdt 用于评估时间精度。
固定加性源没有独立 CFL 限制。
CN 算子仍有大步长负密度限制；修改网格或参数后需要重新检查步长与收敛。

核反应从 `cross_section.hpp` 读取常数截面：对角项为总破坏截面，
非对角项为母核到子核的产生截面。反应率为
nH * beta * c * sigma * seconds_per_year，单位 yr^-1。
beta 按母核刚度单元中心、共同核子质量 0.9315 GeV 计算。
默认只考虑氢靶，`helium_enhancement=1`；需要近似氦增强时可以显式修改该参数。
总破坏截面包括未选入网络的通道，不能用已选产生截面的和代替。

## 每核子动能守恒的刚度映射

共同核子质量近似下，每核子动能不变意味着每核子动量不变：

```
alpha = (A_d / Z_d) / (A_p / Z_p)
Rg_d = alpha * Rg_p
Q_d(Rg_d) = reaction_rate * q_p(Rg_d / alpha) / alpha
```

预计算映射后的母核区间与各子核单元的重叠权重。
母核谱按单元内常数重建，速度和截面也在该母核单元内取常数。
保留落在刚度域内的产生量；域外产物视为离开计算域，不重新归一化。
按截面表的母核到子核顺序进行隐式更新，支持任意选择顺序的核素子集。
该处理依赖当前截面表为无反向通道的上三角网络，不支持循环反应网络。

`nuclear_solver.calculate_terms` 提供瞬时产生项和损失项（均为 q/yr）；
`advance` 则直接隐式推进。调用者提供与初始化匹配的核素/网格维度，
`calculate_terms` 的输入、产生场和损失场必须为不同的存储。
核反应只更新物理空间单元，空间 ghost 由下一次空间算子刷新。

## 输出与范围

每个核素输出一个二进制文件，例如 `C12.bin`、`B11.bin`，保存最终时刻的
物理单元平均密度 q=dn/dRg。文件无头，使用本机字节序的原始 `double`
（当前平台为小端 float64），存储顺序为 `[nRg][nR][nz]`，z 最快变化，
不包含 ghost 单元，不乘空间体积或刚度宽度。同名文件每次覆盖。
维度及网格范围从同目录的 `resolved.cfg` 读取，`settings.txt` 记录运行摘要。

```python
q = np.fromfile("output/cr_solver/C12.bin", dtype=np.float64).reshape(nRg, nR, nz)
# q[k, i, j] 对应 C++ 的 density[s](k, i, j)
```

其中 `np` 为 NumPy，各维度使用该次运行 `resolved.cfg` 中的值。
总粒子数由 sum(q * volume * dRg) 得到。

质子已加入核素列表，但原截面表没有质子反应数据，对应行和列为零，
因此 `proton` 模式求解纯输运与已有动量演化。
Be10、Be7 等目前仅按表进行碰撞反应，不包含放射性衰变。
当前能损仅含原有绝热项，不含电离、库仑或强子连续能损。
现有十核素截面数值保持不变，未引入新的实验截面或能量依赖。

## 代码阅读顺序

1. `main.cpp`：按运行选项、物理参数与网格、核素选择、初始化、演化、输出六个阶段阅读。
2. `field.hpp` / `field.cpp`：先确认 `q=dn/dRg`、空间 ghost 和柱坐标体积；
   `Field3D` 的存储顺序是 `[rigidity][R][z]`，z 连续。
3. `initialization.cpp`：查看源、风、气体和扩散系数的具体模型。
   当前风场仅最底层单元散度非零，因此绝热项仅作用于该层。
4. `spatial_solver.cpp`：对流使用共享面通量，扩散使用两个 ADI 阶段。
5. `momentum_solver.cpp`：先构造刚度面通量，再组装 CN 三对角系统，最后用 Thomas 算法求解。
6. `nuclear_solver.cpp`：构造时缓存破坏系数、母子核刚度映射和反应顺序；
   演化时先更新母核，再用更新后的母核计算子核源项。

核反应中的 `s,d` 是所选核素列表的位置，`id,daughter_id` 是截面表的固定索引；
自定义核素子集时两者通常不同。`order_` 只改变更新顺序，不改变密度场或输出的排列。
顺序算子分裂和核反应 backward Euler 使完整推进为时间一阶，
不能由单个 CN 算子的时间二阶推断完整算法也是二阶。

## 配置文件运行

```sh
make
./cr_solver model.cfg
```

编辑 `model.cfg` 即可更改参数，无须重新编译。使用 `key = value`，支持空行、
`#` 注释以及行尾注释；字符串不加引号，数值使用十进制或科学计数法。
配置文件以 `.cfg` 为扩展名，可复制成不同实验文件，例如 `carbon_boron.cfg`。
相对路径均相对于运行命令时的工作目录。

`mode` 支持 `proton`、`network` 或核素列表，例如：

```ini
mode = O16,C12,B11
source_weights = 1,1,0
steps = 200
dt_yr = 10000
output_directory = output/oxygen_carbon_boron
```

`source_weights` 与所选核素逐项对应，零表示无直接注入。改变 `mode` 时需同步修改
权重。全部配置项均须填写；`configuration.hpp` 仅声明参数类型。缺失参数、未知参数、
重复参数、不合法数值以及核素/权重数量不匹配会报错。
单位写在参数名或示例文件注释中；`source_width_kpc` 是高斯指数中的宽度，
不是标准差，`source_index` 为刚度幂律谱的正指数。

每次输出包含 `resolved.cfg`，记录实际生效的全部可配置参数和源权重，可再次作为输入。
运行模式、步数和输出目录统一在配置文件中设置。
边界条件、数值格式和截面表仍沿用现有实现，不属于本次新增配置项。

D0 的换算为 `D0_internal = D0_cm2_per_s * seconds_per_year / cm_per_kpc^2`，其中 seconds_per_year = 31557600、cm_per_kpc = 3.0856775814913673e21。
默认值约为 `3.0171515377e28 cm²/s`，与原来的 `1e-7 kpc²/yr` 等价。
`resolved.cfg` 保存输入单位 cm²/s，`settings.txt` 同时记录输入和内部单位的 D0。
旧配置中的 `D0_kpc2_per_yr` 需改名并换算数值，避免混淆单位。

物理常数、单位换算与核数据统一定义在 `cross_section.hpp`（常数位于 `physical_constants` 命名空间）；所有输运和核反应均采用 Julian year（365.25 天）。
