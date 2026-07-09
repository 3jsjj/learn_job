%% ============================================================
%  局部函数 1：
%  根据 9 个传感器力值估计接触中心、半径、峰值力
% ============================================================
function [cx, cy, r_est, amp_est, F_base, peak_force] = ...
    estimate_contact_params(F_sensor, sensor_points)

    F_sensor = F_sensor(:);

    %% --------------------------------------------------------
    %  1. 估计基础力
    %
    %  这里取最小传感器力值作为基础力。
    %  这样可以削弱整体偏置对接触中心估计的影响。
    %% --------------------------------------------------------

    F_base = min(F_sensor);

    %% --------------------------------------------------------
    %  2. 构造加权质心用的权重
    %
    %  先减去基础力，再使用平方权重。
    %  平方权重可以让大力点对中心位置的影响更强。
    %% --------------------------------------------------------

    F_weight = F_sensor - F_base;
    F_weight(F_weight < 0) = 0;

    if sum(F_weight) <= eps
        F_weight = F_sensor;
    end

    p = 2;
    W = F_weight .^ p;

    if sum(W) <= eps
        [~, idx_max] = max(F_sensor);
        cx = sensor_points(idx_max, 1);
        cy = sensor_points(idx_max, 2);
    else
        cx = sum(W .* sensor_points(:,1)) / sum(W);
        cy = sum(W .* sensor_points(:,2)) / sum(W);
    end

    %% --------------------------------------------------------
    %  3. 自动拟合接触半径
    %
    %  不人为指定接触半径。
    %  在给定范围内搜索使 9 个传感器点拟合误差最小的半径。
    %% --------------------------------------------------------

    r_min = 1.0;      % 最小搜索半径，单位 mm
    r_max = 12.0;     % 最大搜索半径，单位 mm

    obj_fun = @(r) gaussian_fit_error( ...
        r, ...
        sensor_points, ...
        F_sensor, ...
        cx, ...
        cy, ...
        F_base);

    r_est = fminbnd(obj_fun, r_min, r_max);

    %% --------------------------------------------------------
    %  4. 在最优半径下，用高斯最小二乘拟合峰值幅值
    %% --------------------------------------------------------

    d2 = (sensor_points(:,1) - cx).^2 + ...
         (sensor_points(:,2) - cy).^2;

    g = exp(-d2 / (2 * r_est^2));

    amp_est = (g' * (F_sensor - F_base)) / (g' * g);

    % 防止拟合出负幅值
    amp_est = max(amp_est, 0);

    peak_force = F_base + amp_est;

end