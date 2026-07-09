
function err = gaussian_fit_error(r, sensor_points, F_sensor, cx, cy, F_base)

    F_sensor = F_sensor(:);

    d2 = (sensor_points(:,1) - cx).^2 + ...
         (sensor_points(:,2) - cy).^2;

    g = exp(-d2 / (2 * r^2));

    amp = (g' * (F_sensor - F_base)) / (g' * g);

    amp = max(amp, 0);

    F_fit = F_base + amp * g;

    err = mean((F_fit - F_sensor).^2);

end