% =========================================================================
% 逆有限元（Inverse FEM）节点力场重构程序——中文详注版
%
% 功能概述：
%   1. 从 Abaqus 输入文件中读取节点编号和节点坐标；
%   2. 从 Abaqus 导出的 .mtx 文件中组装全局刚度矩阵；
%   3. 删除固定边界自由度以及节点编号跳跃造成的“空自由度”；
%   4. 根据测点位移构造力—位移灵敏度矩阵；
%   5. 使用 Tikhonov 正则化反演全场节点力；
%   6. 输出部分节点的反演结果，并绘制三维节点力云图。
%
% 重要假设：
%   - 当前程序针对三维实体模型，每个节点只考虑 U1、U2、U3 三个平动自由度；
%   - 刚度矩阵、位移、节点坐标和反演力必须采用一致的单位制；
%   - get_model.m 已经位于 MATLAB 路径中，并能读取指定的 Abaqus .inp 文件；
%   - Abaqus 导出的刚度矩阵采用五列格式：
%       [节点I, 自由度I, 节点J, 自由度J, 刚度值]。
% =========================================================================

clear;
clc;

%% 1. 配置输入文件和基本参数

% Abaqus 模型文件：用于读取节点 ID 和节点坐标。
inp_filepath = 'element_nodes_2.inp';

% Abaqus 刚度矩阵文件：要求为五列节点—自由度格式。
stiffness_filepath = 'get_matrix-3_STIF2.mtx';

% 测点位移文件：
% 第 1 列为节点 ID，第 2~4 列依次为 U1、U2、U3。
measurement_filepath = 'Abaqus_Nodal_U_and_Pressure.csv';

% 固定节点文件：可以是单列或单行节点 ID。
fixed_nodes_filepath = 'fixed_nodes.csv';

% 三维实体模型每个节点具有 3 个平动自由度：U1、U2、U3。
% 注意：壳单元通常还包含转动自由度，不能只把这里改成 6；
% 测点数据读取、力向量重组和可视化部分也需要同步修改。
num_dofs_per_node = 3;

% Tikhonov 正则化参数。
% lambda 越大，结果越平滑、数值越稳定，但反演力的幅值可能被压低；
% lambda 越小，越贴合测量位移，但对噪声和病态矩阵更敏感。
lambda = 1e-6;

% 命令行最多打印多少个节点的力结果。
print_limit_max = 20;

%% 2. 读取模型节点信息

% get_model 应返回：
%   node_ids    : [N_nodes x 1] 节点编号
%   node_coords : [N_nodes x 3] 节点坐标
[node_ids, node_coords] = get_model(inp_filepath);

% 基本完整性检查。
if isempty(node_ids) || isempty(node_coords)
    error('未能从 %s 中读取有效节点数据。', inp_filepath);
end

if size(node_coords, 2) < 3
    error('当前程序需要三维节点坐标，但 node_coords 只有 %d 列。', size(node_coords, 2));
end

node_ids = node_ids(:);  % 强制转换为列向量，便于后续集合运算。

%% 3. 读取并组装 Abaqus 全局刚度矩阵

% 读取五列刚度矩阵数据：
%   第 1 列：行节点编号 Node_I
%   第 2 列：行自由度编号 DOF_I
%   第 3 列：列节点编号 Node_J
%   第 4 列：列自由度编号 DOF_J
%   第 5 列：刚度值 Value
mtx_data = load(stiffness_filepath);

if size(mtx_data, 2) < 5
    error('刚度矩阵文件 %s 必须至少包含 5 列。', stiffness_filepath);
end

node_i = mtx_data(:, 1);
dof_i  = mtx_data(:, 2);
node_j = mtx_data(:, 3);
dof_j  = mtx_data(:, 4);
values = mtx_data(:, 5);

% 将“节点编号 + 节点内自由度编号”映射为全局自由度编号。
% 对于节点 n：(节点三个自由度方向的求全局的自由度编号)
%   U1 -> 3*n-2
%   U2 -> 3*n-1
%   U3 -> 3*n
row_idx = num_dofs_per_node * (node_i - 1) + dof_i;
col_idx = num_dofs_per_node * (node_j - 1) + dof_j;

% 检查节点内自由度编号是否与三维实体假设一致。
% 意思是这里面的自由度编号必须在 1~3 之间，否则说明刚度矩阵可能来自壳单元或其他类型单元。
if any(dof_i < 1 | dof_i > num_dofs_per_node) || ...
   any(dof_j < 1 | dof_j > num_dofs_per_node)
    error(['刚度矩阵中出现了超出 1~%d 的节点内自由度编号。' ...
           '当前脚本仅适用于三维实体单元的平动自由度。'], num_dofs_per_node);
end

% 矩阵尺寸由最大全局自由度编号确定。
% 节点 ID 不连续时，中间会产生全零行/列，后续将通过 existing_dofs 删除。
max_dof = max([row_idx; col_idx]);

% 使用 sparse 直接组装稀疏矩阵；若文件中存在重复项，MATLAB 会自动累加。
K_half = sparse(row_idx, col_idx, values, max_dof, max_dof);
% S = sparse(行号, 列号, 数值, 总行数, 总列数);

% Abaqus 常只导出对称矩阵的一个三角部分，因此需要镜像补全。
% 减去对角线是为了避免对角项被重复计算。
K_raw = K_half + K_half' - spdiags(diag(K_half), 0, max_dof, max_dof);

%% 4. 读取测点位移数据

measured_results = readmatrix(measurement_filepath);

if size(measured_results, 2) < 4
    error('测点文件 %s 至少需要 4 列：节点 ID、U1、U2、U3。', measurement_filepath);
end

% 删除全空行，避免 CSV 尾部空白被 readmatrix 读取为 NaN 行。
measured_results = measured_results(~all(isnan(measured_results), 2), :);

id_measured = measured_results(:, 1);
U_measured_mat = measured_results(:, 2:4);

if any(isnan(id_measured)) || any(isnan(U_measured_mat), 'all')
    error('测点文件中包含 NaN，请检查节点 ID 和 U1/U2/U3 数据。');
end

% 按节点顺序展开为列向量：
% [U1_x; U1_y; U1_z; U2_x; U2_y; U2_z; ...]
U_measured = reshape(U_measured_mat', [], 1);

% 已知测点外力向量。
% 当前假设测点没有额外已知集中力，因此初始化为零。
% 若实际存在已知测点力，应按与 U_measured 完全相同的顺序赋值。
% 结构当然受到未知力，否则不会产生位移。只是这些力正是程序准备反演的对象。
F_measured = zeros(size(U_measured));

%% 5. 设置虚拟受力节点与固定边界

% 本程序把所有模型节点都作为候选虚拟受力节点。
id_virtual = node_ids;

% 读取固定节点编号，并强制整理为行向量。
fixed_nodes_data = readmatrix(fixed_nodes_filepath);
fixed_nodes_data = fixed_nodes_data(~isnan(fixed_nodes_data));
fixed_nodes = unique(reshape(fixed_nodes_data, 1, []));

% 将固定节点转换为受约束的全局自由度编号。
% 使用向量化写法代替循环动态扩容。
if isempty(fixed_nodes)
    constrained_dofs = zeros(1, 0);
else
    constrained_dofs = reshape( ...
        num_dofs_per_node * (fixed_nodes - 1) + (1:num_dofs_per_node)', ...
        1, []);
end

% existing_dofs 表示 .mtx 文件中真正出现过的物理自由度。
% 这样可以删除节点 ID 跳号造成的全零“幽灵行/列”。
existing_dofs = unique([row_idx; col_idx]);

% 活动自由度 = 真实存在的自由度 - 固定边界自由度。
active_dofs = setdiff(existing_dofs, constrained_dofs);

if isempty(active_dofs)
    error('删除固定自由度后没有剩余活动自由度，请检查固定节点列表。');
end

% 提取施加边界条件后的活动刚度矩阵。
K_global = K_raw(active_dofs, active_dofs);

%% 6. 过滤无效节点，并建立节点 ID 到活动自由度的映射

% 先删除明确属于固定端的候选虚拟节点。
id_virtual_valid = setdiff(id_virtual, fixed_nodes, 'stable');

% 找到候选虚拟节点在坐标数组中的位置。
[is_virtual_in_model, loc_in_coords] = ismember(id_virtual_valid, node_ids);

% 理论上 id_virtual 来自 node_ids，因此应全部匹配；仍保留检查以增强健壮性。
if any(~is_virtual_in_model)
    warning('有 %d 个虚拟节点未在 node_ids 中找到，已自动删除。', sum(~is_virtual_in_model));
    id_virtual_valid = id_virtual_valid(is_virtual_in_model);
    loc_in_coords = loc_in_coords(is_virtual_in_model);
end

coords_virtual_valid = node_coords(loc_in_coords, 1:3);

% 为每个候选虚拟节点构造其 3 个全局自由度编号。
virtual_global_dof_matrix = ...
    num_dofs_per_node * (id_virtual_valid - 1) + (1:num_dofs_per_node);

% 只有当节点的 U1、U2、U3 三个自由度都存在于 active_dofs 中时，
% 才把该节点保留为可反演的虚拟节点。
valid_mask_v = all(ismember(virtual_global_dof_matrix, active_dofs), 2);

id_virtual_valid = id_virtual_valid(valid_mask_v);
coords_virtual_valid = coords_virtual_valid(valid_mask_v, :);
virtual_global_dof_matrix = virtual_global_dof_matrix(valid_mask_v, :);

% 对测点做同样的合法性检查。
measured_global_dof_matrix = ...
    num_dofs_per_node * (id_measured - 1) + (1:num_dofs_per_node);
valid_mask_m = all(ismember(measured_global_dof_matrix, active_dofs), 2);

if any(~valid_mask_m)
    warning('有 %d 个测点不具备完整活动自由度，已从反演中删除。', sum(~valid_mask_m));
end

id_measured = id_measured(valid_mask_m);
U_measured_mat_valid = U_measured_mat(valid_mask_m, :);
measured_global_dof_matrix = measured_global_dof_matrix(valid_mask_m, :);

% 过滤节点后，重新生成测点位移和已知力向量，保证维度严格一致。
U_measured = reshape(U_measured_mat_valid', [], 1);
F_measured = zeros(size(U_measured));

if isempty(id_measured)
    error('过滤后没有有效测点，无法执行力场反演。');
end

if isempty(id_virtual_valid)
    error('过滤后没有有效虚拟节点，无法构造待求力向量。');
end

% 按节点顺序展开全局自由度：
% [node1-U1; node1-U2; node1-U3; node2-U1; ...]
global_dofs_measured = reshape(measured_global_dof_matrix', [], 1);
global_dofs_virtual  = reshape(virtual_global_dof_matrix', [], 1);

% 将原始全局自由度编号映射到 K_global 的局部行列编号。
[is_measured_active, dof_measured] = ismember(global_dofs_measured, active_dofs);
[is_virtual_active,  dof_virtual]  = ismember(global_dofs_virtual, active_dofs);

if any(~is_measured_active) || any(~is_virtual_active)
    error('自由度映射失败：过滤逻辑与 active_dofs 不一致。');
end

%% 7. 构造虚拟力到测点位移的灵敏度矩阵 H

N_total_dof = size(K_global, 1);
N_virtual_dof = numel(dof_virtual);

% 在每一个虚拟自由度上依次施加单位力。
% F_unit_virtual 的第 j 列表示：在第 j 个虚拟自由度上施加 1 单位力。
F_unit_virtual = sparse( ...
    dof_virtual, ...              % 非零元素所在行
    1:N_virtual_dof, ...          % 每个虚拟自由度对应一列
    1, ...                        % 单位载荷幅值
    N_total_dof, ...              % 总行数
    N_virtual_dof);               % 总列数

% 有限元平衡方程为 K*u=f。
% 对所有单位虚拟载荷同时求解，可得到柔度矩阵中所需的列：
%   K_global * C_virtual_cols = F_unit_virtual
%   C_virtual_cols = K_global \ F_unit_virtual
%
% 这里使用反斜杠求解线性方程，避免显式计算 inv(K_global)。
C_virtual_cols = K_global \ F_unit_virtual;

% 只保留测点自由度对应的响应行，得到灵敏度矩阵 H。
% H 的每一列表示某一虚拟单位力在全部测点上产生的位移响应。
H = C_virtual_cols(dof_measured, :);

%% 8. 计算已知测点力引起的基准位移

if any(F_measured ~= 0)
    % 将已知测点力放入活动自由度空间中的稀疏列向量。
    % 原代码使用 1:length(dof_measured) 作为列索引但矩阵仅有 1 列，
    % 在 F_measured 非零时会产生列索引越界；此处统一使用列索引 1。
    F_measured_global = sparse( dof_measured, ones(size(dof_measured)), F_measured, N_total_dof, 1);

    % 求解已知力产生的全场位移，并提取测点位置的位移。
    U_from_F_measured_full = K_global \ F_measured_global;
    U_from_F_measured = U_from_F_measured_full(dof_measured);
else
    % 当前默认没有已知测点外力，因此基准位移为零。
    U_from_F_measured = zeros(size(U_measured));
end

%% 9. 运行前矩阵健康检查

disp(' ');
disp('====== 矩阵健康状态诊断 ======');
fprintf('1. 固定端节点数量: %d 个\n', numel(fixed_nodes));
fprintf('2. 有效测点数量: %d 个\n', numel(id_measured));
fprintf('3. 有效虚拟节点数量: %d 个\n', numel(id_virtual_valid));
fprintf('4. 活动自由度数量: %d 个\n', numel(active_dofs));

if isempty(fixed_nodes)
    warning('固定节点数量为 0，模型可能保留刚体位移并导致刚度矩阵奇异。');
end

fprintf('5. U_measured 是否包含 NaN: %s\n', ...
    mat2str(any(isnan(U_measured), 'all')));
fprintf('6. K_global 是否包含 NaN: %s\n', ...
    mat2str(any(isnan(K_global), 'all')));
fprintf('7. H 是否包含 NaN: %s\n', ...
    mat2str(any(isnan(H), 'all')));

if any(isnan(K_global), 'all') || any(isinf(K_global), 'all')
    error('K_global 中包含 NaN 或 Inf，无法继续计算。');
end

if any(isnan(H), 'all') || any(isinf(H), 'all')
    error('灵敏度矩阵 H 中包含 NaN 或 Inf，无法继续计算。');
end

% condest 适用于稀疏矩阵，用于估算 1-范数条件数。
% 条件数越大，线性系统越病态，对测量噪声越敏感。
cond_K = condest(K_global);
fprintf('8. K_global 条件数估计: %.6e\n', cond_K);

if ~isfinite(cond_K) || cond_K > 1e15
    warning(['K_global 严重病态或接近奇异。请检查固定边界、连接关系、' ...
             '未约束刚体运动以及刚度矩阵导出设置。']);
end

disp('================================');
disp(' ');

%% 10. 使用 Tikhonov 正则化反演虚拟节点力

% 测量位移中尚未被已知测点力解释的部分。
Delta_U = U_measured - U_from_F_measured;

% 使用稀疏单位矩阵，避免 eye 在虚拟自由度很多时占用大量内存。
I_reg = speye(N_virtual_dof);

% 求解目标：
%   min ||H*F_virtual - Delta_U||_2^2 + lambda*||F_virtual||_2^2
%
% 对应的正规方程：
%   (H'*H + lambda*I)*F_virtual = H'*Delta_U
%
% MATLAB 的反斜杠用于求解线性系统，不显式计算矩阵逆。
F_virtual_reconstructed = ...
    (H' * H + lambda * I_reg) \ (H' * Delta_U);

%% 11. 将自由度力向量重组为节点三向力

% F_virtual_reconstructed 的排列为：
% [F1x; F1y; F1z; F2x; F2y; F2z; ...]
%
% 重组后 F_nodes_3D 的每一行为：
% [Fx, Fy, Fz]
F_nodes_3D = reshape(F_virtual_reconstructed, num_dofs_per_node, [])';

% 计算节点合力幅值，便于后续保存或绘图。
F_nodes_magnitude = vecnorm(F_nodes_3D, 2, 2);

disp('========================================================');
disp('                 虚拟节点力场重构完成');
disp('========================================================');

N_valid_nodes = numel(id_virtual_valid);
print_limit = min(print_limit_max, N_valid_nodes);

for i = 1:print_limit
    node_id = id_virtual_valid(i);
    coords = coords_virtual_valid(i, :);

    fprintf('虚拟节点 ID: %d，坐标: [%.6g, %.6g, %.6g]\n', ...
        node_id, coords(1), coords(2), coords(3));
    fprintf('  X 方向重构力: %12.6g N\n', F_nodes_3D(i, 1));
    fprintf('  Y 方向重构力: %12.6g N\n', F_nodes_3D(i, 2));
    fprintf('  Z 方向重构力: %12.6g N\n', F_nodes_3D(i, 3));
    fprintf('  合力幅值:       %12.6g N\n', F_nodes_magnitude(i));
    disp('--------------------------------------------------------');
end

if N_valid_nodes > print_limit
    fprintf('\n其余 %d 个节点未在命令行中展开。\n', N_valid_nodes - print_limit);
end

% 将完整结果整理成表格，方便在工作区查看或导出。
force_result_table = table( ...
    id_virtual_valid, ...
    coords_virtual_valid(:, 1), ...
    coords_virtual_valid(:, 2), ...
    coords_virtual_valid(:, 3), ...
    F_nodes_3D(:, 1), ...
    F_nodes_3D(:, 2), ...
    F_nodes_3D(:, 3), ...
    F_nodes_magnitude, ...
    'VariableNames', { ...
        'NodeID', 'X', 'Y', 'Z', 'Fx', 'Fy', 'Fz', 'ForceMagnitude'});

% 如需自动写出 CSV，可取消下一行注释。
% writetable(force_result_table, 'reconstructed_nodal_force.csv');

%% 12. 绘制三维节点力云图

disp('正在渲染三维力场云图...');

X_virt = coords_virtual_valid(:, 1);
Y_virt = coords_virtual_valid(:, 2);
Z_virt = coords_virtual_valid(:, 3);

% 查找固定节点坐标。
% ismember 未找到的节点会返回位置 0，必须先过滤，否则 node_coords(0,:) 会报错。
[is_fixed_in_model, loc_fixed] = ismember(fixed_nodes, node_ids);
loc_fixed = loc_fixed(is_fixed_in_model);
coords_fixed = node_coords(loc_fixed, 1:3);

if any(~is_fixed_in_model)
    warning('有 %d 个固定节点未在 node_ids 中找到，绘图时已忽略。', sum(~is_fixed_in_model));
end

% 当前颜色映射使用 Z 方向法向力。
% 若想显示合力大小，可改为：Color_Data = F_nodes_magnitude;
Color_Data = F_nodes_3D(:, 3);

figure( ...
    'Name', '全模型逆向受力云图', ...
    'Color', 'w', ...
    'Position', [100, 100, 850, 650]);

hold on;
grid on;

% 第一层：固定端节点，用灰色小点显示几何边界。
if ~isempty(coords_fixed)
    scatter3( ...
        coords_fixed(:, 1), ...
        coords_fixed(:, 2), ...
        coords_fixed(:, 3), ...
        15, ...
        [0.7, 0.7, 0.7], ...
        'filled', ...
        'MarkerEdgeColor', 'none', ...
        'DisplayName', '固定端边界');
end

% 第二层：有效虚拟节点，用颜色表示对应的重构力。
scatter3( ...
    X_virt, Y_virt, Z_virt, ...
    40, Color_Data, ...
    'filled', ...
    'MarkerEdgeColor', 'none', ...
    'DisplayName', '反演受力节点');

view(45, 30);
axis equal;

xlabel('X 坐标', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y 坐标', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Z 坐标', 'FontSize', 12, 'FontWeight', 'bold');
title('三维逆有限元节点力场重构结果', 'FontSize', 15, 'FontWeight', 'bold');

set(gcf, 'Renderer', 'painters');
colormap('jet');

cb = colorbar;
ylabel(cb, '重构节点 Z 向力', 'FontSize', 12, 'FontWeight', 'bold');

% 对拉力与压力使用关于 0 对称的颜色范围。
c_max = max(abs(Color_Data));
if isfinite(c_max) && c_max > 0
    caxis([-c_max, c_max]);
end

legend('Location', 'best');
hold off;

disp('云图渲染完成。完整结果保存在 force_result_table 和 F_nodes_3D 中。');
