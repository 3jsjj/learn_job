function err = gaussian_fit_error( ...
    r, sensor_points, F_sensor, cx, cy, F_base)
%GAUSSIAN_FIT_ERROR
% 计算给定高斯半径下，拟合值与传感器实测值之间的均方误差。

    F_sensor = F_sensor(:);

    d2 = ...
        (sensor_points(:, 1) - cx).^2 + ...
        (sensor_points(:, 2) - cy).^2;

    g = exp(-d2 / (2 * r^2));

    denominator = g' * g;

    if denominator <= eps
        amp = 0;
    else
        amp = ...
            (g' * (F_sensor - F_base)) / denominator;
    end

    amp = max(amp, 0);

    F_fit = F_base + amp * g;

    err = mean((F_fit - F_sensor).^2);

end
