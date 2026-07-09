clear;
clc;

%% ============================================================
%  使用说明
%
%  本文件不参与训练和预测，仅用于说明运行顺序。
%
%  文件列表：
%
%  1. train_mlp_model.m
%     功能：
%     - 生成训练数据；
%     - 训练 MLP；
%     - 验证模型；
%     - 保存 trained_pressure_model.mat。
%
%  2. predict_force_field.m
%     功能：
%     - 加载 trained_pressure_model.mat；
%     - 输入新的 9 个传感器值；
%     - 输出 21×21 压力场；
%     - 绘制三维图和二维等高线图。
%
%  3. estimate_contact_params.m
%     功能：
%     - 根据 9 点传感器值估计接触中心；
%     - 估计接触半径；
%     - 估计高斯幅值和背景值。
%
%  4. gaussian_fit_error.m
%     功能：
%     - 计算高斯拟合误差；
%     - 供 estimate_contact_params.m 调用。
%
%% ============================================================

%% 第一次使用

% 第一步：
% 在 MATLAB 当前文件夹中放入全部 .m 文件。
%
% 第二步：
% 运行：
%
% train_mlp_model
%
% 训练结束后会生成：
%
% trained_pressure_model.mat
%
% 第三步：
% 打开 predict_force_field.m，
% 修改 current_sensor_values 中的 9 个传感器值。
%
% 第四步：
% 运行：
%
% predict_force_field

%% ============================================================
%  后续预测
%
%  后续预测不需要重新运行 train_mlp_model.m。
%
%  只需要：
%
%  1. 修改 predict_force_field.m 中的 current_sensor_values；
%  2. 运行 predict_force_field.m。
%
%% ============================================================

disp('请先运行 train_mlp_model.m 完成训练。');
disp('训练完成后，再运行 predict_force_field.m 进行预测。');
