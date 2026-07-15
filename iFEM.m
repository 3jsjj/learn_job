% =========================================================================
% 纯净版逆有限元 (Inverse FEM) + 力场重构框架 (基于真实输入)
% =========================================================================
clear; clc;

%% 1. 用户数据输入区 (请替换以下占位符)

% --- A. 模型基础信息 ---
% 导入所有节点的坐标矩阵，维度为 [N_nodes, 2] (2D) 或 [N_nodes, 3] (3D)
inp_filepath='element_nodes_2.inp';
[node_ids, node_coords] = get_model(inp_filepath);
% 导入处理过边界条件后的全局刚度矩阵 K (建议为 sparse 稀疏矩阵格式)
% 2. 导入并组装 Abaqus 刚度矩阵 (处理 5 列的节点-自由度格式)
mtx_data = load('get_matrix-3_STIF2.mtx'); 

% 拆解 Abaqus 的 5 列数据: [Node_I, DOF_I, Node_J, DOF_J, Value]
node_i = mtx_data(:, 1);
dof_i  = mtx_data(:, 2);
node_j = mtx_data(:, 3);
dof_j  = mtx_data(:, 4);
values = mtx_data(:, 5);

% 将“节点ID + 自由度方向”映射为矩阵的绝对行/列索引
% 这里的 3 代表 3D 实体单元(Solid)每个节点有 3 个平动自由度(U1, U2, U3)
% （如果你在这层薄膜模型里用的是壳单元 Shell，请把这里的 3 改成 6）
num_dofs_per_node = 3; 
row_idx = num_dofs_per_node * (node_i - 1) + dof_i;
col_idx = num_dofs_per_node * (node_j - 1) + dof_j;

% 获取模型最大的自由度编号，以确定矩阵的绝对尺寸
max_dof = max(max(row_idx), max(col_idx));

% 直接使用 sparse 函数构建稀疏矩阵 (它比 spconvert 更灵活且自动累加重复项)
K_half = sparse(row_idx, col_idx, values, max_dof, max_dof);

% Abaqus 默认只输出对称矩阵的下三角/上三角，需补全为完整的对称矩阵
K_raw = K_half + K_half' - diag(diag(K_half));


% --- B. 真实测点输入 ---
% 读取你在 Abaqus 提取的测点位移结果文件
measured_results = readmatrix('Abaqus_Nodal_U_and_Pressure.csv');
id_measured = measured_results(:, 1); % 第一列是测点 ID

% 提取 U1, U2, U3 位移矩阵，并展平为一维列向量 [U1_x; U1_y; U1_z; U2_x; ...]
U_measured_mat = measured_results(:, 2:4); 
U_measured = reshape(U_measured_mat', [], 1); 

% 假设测点本身没有承受已知的集中外力干扰
% 【修正】测点没有已知外力，赋全 0 向量，且维度必须和 U_measured 保持一致
F_measured = zeros(size(U_measured));


% --- C. 虚拟位点输入 ---
% 将全场所有节点作为推演目标
id_virtual = node_ids;          % 【新增这行】全场节点 ID
coords_virtual = node_coords;   % 用于最后 print 打印坐标

% 正则化参数 (根据位移数据的信噪比调整。无噪音填0，有噪音尝试 1e-6 ~ 1e-2)
lambda = 1e-6; 


%%
% 3. 施加边界条件消除奇异性 (这一步不可省略！)
% 【完美自动化替换】
% 直接读取 Python 提取出的固定节点 ID 列表
fixed_nodes_data = readmatrix('fixed_nodes.csv');
fixed_nodes = reshape(fixed_nodes_data, 1, []); % 确保它是行向量

constrained_dofs = [];
for n = fixed_nodes
    % 3D 实体模型每个节点有 3 个自由度 (U1, U2, U3)
    constrained_dofs = [constrained_dofs, 3*n-2, 3*n-1, 3*n]; 
end

% =========================================================
% 【核心修复】：找出矩阵中真正有刚度数值的物理自由度
% 这一步直接剔除了所有因为 Abaqus 节点 ID 跳跃而产生的全 0 幽灵空行
existing_dofs = unique([row_idx; col_idx]); 

% 在真实存在的自由度中，进一步扣除被固定的边界自由度
active_dofs = setdiff(existing_dofs, constrained_dofs);
% =========================================================
% 提取最终可计算的非奇异全局刚度矩阵
K_global = K_raw(active_dofs, active_dofs);


%% 2. 节点 ID 到矩阵自由度 (DOF) 的精确映射与安全过滤

% 1. 剔除全场节点中明确属于固定端的部分
id_virtual_valid = setdiff(id_virtual, fixed_nodes);

% 2. 找出虚拟节点 ID 在 node_ids 列表中的“真实行号”
[~, loc_in_coords] = ismember(id_virtual_valid, node_ids);
coords_virtual_valid = node_coords(loc_in_coords, :);

% =========================================================================
% 【防崩溃核心逻辑】：强制过滤“幽灵节点”
% 模型中可能包含参考点、刚体节点、或是没有输出到 .mtx 的多余节点
% 我们必须进行体检，只保留 3 个方向自由度都完全存在于 K_global 中的节点
% =========================================================================

% --- 过滤虚拟节点 ---
valid_mask_v = ismember(3*id_virtual_valid-2, active_dofs) & ...
               ismember(3*id_virtual_valid-1, active_dofs) & ...
               ismember(3*id_virtual_valid, active_dofs);

% 仅保留合法的虚拟节点及其坐标
id_virtual_valid = id_virtual_valid(valid_mask_v);
coords_virtual_valid = coords_virtual_valid(valid_mask_v, :);

% --- 过滤测点节点 ---
valid_mask_m = ismember(3*id_measured-2, active_dofs) & ...
               ismember(3*id_measured-1, active_dofs) & ...
               ismember(3*id_measured, active_dofs);

% 仅保留合法的测点，并同步更新位移和力矩阵维度
id_measured = id_measured(valid_mask_m);
U_measured_mat_valid = U_measured_mat(valid_mask_m, :);
U_measured = reshape(U_measured_mat_valid', [], 1);
F_measured = zeros(size(U_measured));


% 3. 展开为 3D 全局自由度向量 [U1x; U1y; U1z; U2x; ...]
global_dofs_measured = reshape([3*id_measured-2, 3*id_measured-1, 3*id_measured]', [], 1);
global_dofs_virtual  = reshape([3*id_virtual_valid-2, 3*id_virtual_valid-1, 3*id_virtual_valid]', [], 1);

% 4. 获取在消除边界后的活动刚度矩阵 (active_dofs) 中的绝对行/列索引
% 经过上面的严格过滤，这里的映射绝对 100% 成功，不可能再出现 0！
[~, dof_measured] = ismember(global_dofs_measured, active_dofs);
[~, dof_virtual]  = ismember(global_dofs_virtual, active_dofs);


%% 3. 高效构建灵敏度矩阵 (避免全矩阵求逆)
N_total_dof = size(K_global, 1);
N_virtual_dof = length(dof_virtual);

% 构造虚拟节点位置的单位力载荷矩阵 (稀疏矩阵)
% 列数为虚拟节点的数量，每列代表在一个虚拟节点施加单位力
F_unit_virtual = sparse(dof_virtual, 1:N_virtual_dof, 1, N_total_dof, N_virtual_dof);

% 直接求解全局系统得到对应的位移场 (即柔度矩阵中与虚拟节点相关的列)
% K * C_col = F_unit  =>  C_col = K \ F_unit
C_virtual_cols = K_global \ F_unit_virtual;

% 提取灵敏度观测矩阵 H (仅提取测点自由度所在行的位移响应)
H = C_virtual_cols(dof_measured, :);


%% 4. 计算测点已知力产生的基准位移 (如果有)
if any(F_measured)
    % 与上面同理，计算测点已知力在全场产生的位移
    F_unit_measured = sparse(dof_measured, 1:length(dof_measured), F_measured, N_total_dof, 1);
    U_from_F_measured_full = K_global \ F_unit_measured;
    U_from_F_measured = U_from_F_measured_full(dof_measured);
else
    U_from_F_measured = zeros(size(U_measured));
end

%%
% =========================================================================
% 🚑 救命诊断探针 (运行健康检查)
% =========================================================================
disp(' ');
disp('====== 🚑 矩阵健康状态诊断 ======');
fprintf('1. 识别到的固定端节点数量: %d 个\n', length(fixed_nodes));

if length(fixed_nodes) == 0
    warning('致命错误：固定节点数量为 0！刚度矩阵将无法求逆，必然产生 NaN。');
end

fprintf('2. 测点位移 (U_measured) 是否含 NaN: %s\n', mat2str(any(isnan(U_measured))));
fprintf('3. 全局刚度矩阵 (K_global) 是否含 NaN: %s\n', mat2str(any(isnan(K_global), 'all')));

% condest 用于估算稀疏矩阵的条件数。如果 > 1e15，说明矩阵严重奇异（存在刚体位移）
cond_K = condest(K_global);
fprintf('4. K_global 矩阵条件数估算: %e\n', cond_K);
if cond_K > 1e15
    warning('致命错误：全局刚度矩阵严重奇异 (条件数极大)！边界条件未有效约束刚体位移。');
end

fprintf('5. 灵敏度矩阵 (H) 是否含 NaN: %s\n', mat2str(any(isnan(H), 'all')));
disp('=================================');
disp(' ');




%% 5. Tikhonov 正则化求解虚拟节点力场
% 计算需要由“虚拟力”来弥补的位移残差
Delta_U = U_measured - U_from_F_measured;

% 构建单位矩阵用于正则化惩罚项
I_reg = eye(N_virtual_dof);

% 核心求解公式: F_v = inv(H'*H + lambda*I) * H' * Delta_U
F_virtual_reconstructed = (H' * H + lambda * I_reg) \ (H' * Delta_U);


%% 6. 结果重组与输出 (3D 升级版)
% 此时的 F_virtual_reconstructed 是展开的列向量 [F1x; F1y; F1z; F2x; ...]
% 重组为 [N个有效虚拟节点, 3个方向] 的矩阵
F_nodes_3D = reshape(F_virtual_reconstructed, 3, [])';

disp('========================================================');
disp('                  虚拟位点力场重构完成');
disp('========================================================');

% 限制打印前 20 个节点，防止全场节点数量庞大导致命令行卡死
N_valid_nodes = length(id_virtual_valid);
print_limit = min(20, N_valid_nodes);

for i = 1:print_limit
    node_id = id_virtual_valid(i);
    coords = coords_virtual_valid(i, :);

    fprintf('虚拟节点 ID: %d  (坐标: [%.3f, %.3f, %.3f])\n', node_id, coords(1), coords(2), coords(3));
    fprintf('  -> X方向重构力: %10.4f N\n', F_nodes_3D(i, 1));
    fprintf('  -> Y方向重构力: %10.4f N\n', F_nodes_3D(i, 2));
    fprintf('  -> Z方向重构力: %10.4f N\n', F_nodes_3D(i, 3));
    disp('--------------------------------------------------------');
end

if N_valid_nodes > 20
    fprintf('\n... (省略其余 %d 个节点的输出)\n', N_valid_nodes - 20);
    disp('*** 全部推演受力数据已结构化保存在 F_nodes_3D 矩阵中，可用于后续云图绘制。 ***');
end

%% 7. 高清三维力场云图可视化 (全模型对比版)
disp('正在渲染全景三维力场云图...');

% 1. 提取虚拟节点的坐标数据 (用于彩色映射)
X_virt = coords_virtual_valid(:, 1);
Y_virt = coords_virtual_valid(:, 2);
Z_virt = coords_virtual_valid(:, 3);

% 2. 提取固定节点的坐标数据 (用于绘制灰色基座轮廓，补全模型)
[~, loc_fixed] = ismember(fixed_nodes, node_ids);
coords_fixed = node_coords(loc_fixed, :);

% 3. 选择你要映射为颜色的物理量 (Z方向法向力)
Color_Data = F_nodes_3D(:, 3); 

% 4. 创建绘图窗口
figure('Name', '全模型逆向受力云图', 'Color', 'w', 'Position', [100, 100, 850, 650]);
hold on; grid on;

% =======================================================
% 核心画图逻辑：分层渲染
% =======================================================

% 第一层：画出被固定的节点 (使用灰色小点，作为几何边界对照)
if ~isempty(coords_fixed)
    scatter3(coords_fixed(:,1), coords_fixed(:,2), coords_fixed(:,3), ...
             15, [0.7 0.7 0.7], 'filled', 'MarkerEdgeColor', 'none', ...
             'DisplayName', '固定端边界');
end

% 第二层：画出推演出受力的虚拟节点 (使用彩色大点，展示应力梯度)
scatter3(X_virt, Y_virt, Z_virt, ...
         40, Color_Data, 'filled', 'MarkerEdgeColor', 'none', ...
         'DisplayName', '推演受力点');

% =======================================================

% 5. 视角与坐标轴美化
view(45, 30);            
axis equal;              
xlabel('X 坐标 (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y 坐标 (mm)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Z 坐标 (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title('多层薄膜全模型三维力场云图', 'FontSize', 15, 'FontWeight', 'bold');

% 6. 渲染器与色彩映射配置
set(gcf, 'Renderer', 'painters'); 
colormap('jet');                  
cb = colorbar;                    
ylabel(cb, '重构节点法向力 (N)', 'FontSize', 12, 'FontWeight', 'bold');

% 对称化颜色条 (让受压和受拉的颜色对比更强烈)
c_max = max(abs(Color_Data));
if c_max > 0
    caxis([-c_max, c_max]); 
end

legend('Location', 'best');
hold off;
disp('云图渲染完成！现在你可以看到整个模型的完整轮廓了。');
