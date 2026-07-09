%% ============================================================
%  局部函数 1：
%  根据 9 个传感器力值估计接触中心、半径、峰值力
% ============================================================
function [cx, cy, r_est, amp_est, F_base, peak_force] = ...
    estimate_contact_params(F_sensor, sensor_points)

    F_sensor = F_sensor(:);

    %% --------------------------------------------------------
    %  1. 估计基础力
    %% --------------------------------------------------------

    F_base = min(F_sensor);

    %% --------------------------------------------------------
    %  2. 构造加权质心用的权重
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
    %% --------------------------------------------------------

    r_min = 1.0;
    r_max = 12.0;

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

    denominator = g' * g;

    if denominator <= eps
        amp_est = 0;
    else
        amp_est = (g' * (F_sensor - F_base)) / denominator;
    end

    amp_est = max(amp_est, 0);

    peak_force = F_base + amp_est;

end
