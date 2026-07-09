function [cx, cy, r_est, amp_est, F_base, peak_force] = ...
    estimate_contact_params(F_sensor, sensor_points)
%ESTIMATE_CONTACT_PARAMS
% 根据 9 个传感器力值估计接触中心、接触半径、
% 高斯幅值、基础值和峰值。

    F_sensor = F_sensor(:);

    %% 1. 估计基础值

    F_base = min(F_sensor);

    %% 2. 使用平方权重计算加权质心

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

        cx = sum(W .* sensor_points(:, 1)) / sum(W);
        cy = sum(W .* sensor_points(:, 2)) / sum(W);

    end

    %% 3. 搜索最优接触半径

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

    %% 4. 最小二乘估计高斯幅值

    d2 = ...
        (sensor_points(:, 1) - cx).^2 + ...
        (sensor_points(:, 2) - cy).^2;

    g = exp(-d2 / (2 * r_est^2));

    denominator = g' * g;

    if denominator <= eps
        amp_est = 0;
    else
        amp_est = ...
            (g' * (F_sensor - F_base)) / denominator;
    end

    amp_est = max(amp_est, 0);

    peak_force = F_base + amp_est;

end
