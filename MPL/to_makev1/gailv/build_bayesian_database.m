clear;
clc;
close all;

%% 1. 文件路径
x_file = fullfile('results', 'X_data.csv');
y_file = fullfile('results', 'Y_data.csv');

if ~isfile(x_file)
    error('未找到文件：%s', x_file);
end
if ~isfile(y_file)
    error('未找到文件：%s', y_file);
end

%% 2. 读取 Abaqus 数据
X_table = readtable(x_file, 'VariableNamingRule', 'preserve');
Y_table = readtable(y_file, 'VariableNamingRule', 'preserve');

if height(X_table) ~= height(Y_table)
    error('X_data 和 Y_data 的样本数量不一致。');
end

case_name = string(X_table{:, 1});
side_length_mm = X_table{:, 2};
area_mm2 = X_table{:, 3};
depth_mm = X_table{:, 4};

X_sensor = X_table{:, 5:13};
Y_pressure = Y_table{:, 5:445};

if size(X_sensor, 2) ~= 9
    error('X_data 必须包含9个传感器输入。');
end
if size(Y_pressure, 2) ~= 441
    error('Y_data 必须包含441个压力点。');
end

if any(~isfinite(X_sensor), 'all') || any(~isfinite(Y_pressure), 'all')
    error('数据中存在 NaN 或 Inf。');
end

%% 3. 设置传感器噪声模型
relative_noise = 0.05;

sensor_scale = std(X_sensor, 0, 1);
fallback_scale = max(abs(X_sensor), [], 1);

sensor_scale(sensor_scale < eps) = fallback_scale(sensor_scale < eps);
sensor_scale(sensor_scale < eps) = 1;

sigma_sensor = relative_noise .* sensor_scale;
sigma_floor = max(1e-6, 0.001 .* max(abs(X_sensor), [], 1));
sigma_sensor = max(sigma_sensor, sigma_floor);

%% 4. 设置均匀先验
N = size(X_sensor, 1);
prior_probability = ones(N, 1) / N;

%% 5. 建立21×21网格
x_grid = 0:1:20;
y_grid = 0:1:20;
[X_grid, Y_grid] = meshgrid(x_grid, y_grid);
grid_points = [X_grid(:), Y_grid(:)];

%% 6. 保存数据库
save('bayesian_inverse_database.mat', ...
    'case_name', ...
    'side_length_mm', ...
    'area_mm2', ...
    'depth_mm', ...
    'X_sensor', ...
    'Y_pressure', ...
    'sigma_sensor', ...
    'prior_probability', ...
    'relative_noise', ...
    'x_grid', ...
    'y_grid', ...
    'X_grid', ...
    'Y_grid', ...
    'grid_points');

fprintf('\n概率反演数据库构建完成。\n');
fprintf('样本数量：%d\n', N);
fprintf('输入维度：%d\n', size(X_sensor,2));
fprintf('输出维度：%d\n', size(Y_pressure,2));
fprintf('数据库文件：bayesian_inverse_database.mat\n');
disp('各传感器噪声标准差：');
disp(sigma_sensor);
