clear;
clc;
close all;

%% 仅用于测试代码流程，正式使用时请换成 Abaqus 数据。

if ~isfolder('results')
    mkdir('results');
end

side_list = [2,4,6,8,10];
depth_list = [0.5,1.0,1.5,2.0,2.5,3.0];

[X_grid, Y_grid] = meshgrid(0:1:20, 0:1:20);

sensor_points = [
    0,0;
    10,0;
    20,0;
    0,10;
    10,10;
    20,10;
    0,20;
    10,20;
    20,20
];

N = numel(side_list) * numel(depth_list);

case_name = strings(N,1);
side_length_mm = zeros(N,1);
area_mm2 = zeros(N,1);
depth_mm = zeros(N,1);

X_sensor = zeros(N,9);
Y_pressure = zeros(N,441);

k = 0;

for side = side_list
    for depth = depth_list

        k = k + 1;

        sigma = max(side/2.5, 0.8);
        amplitude = 1000 * depth / max(side^2,1);

        F_map = amplitude .* exp( ...
            -((X_grid-10).^2 + (Y_grid-10).^2) ...
            ./ (2*sigma^2));

        sensor_values = zeros(9,1);

        for i = 1:9
            x = sensor_points(i,1);
            y = sensor_points(i,2);
            sensor_values(i) = F_map(y+1, x+1);
        end

        case_name(k) = sprintf('sq%.1f_d%.1f',side,depth);
        side_length_mm(k) = side;
        area_mm2(k) = side^2;
        depth_mm(k) = depth;

        X_sensor(k,:) = sensor_values';
        Y_pressure(k,:) = F_map(:)';

    end
end

X_table = table(case_name,side_length_mm,area_mm2,depth_mm);

for i = 1:9
    X_table.(sprintf('S%d',i)) = X_sensor(:,i);
end

Y_table = table(case_name,side_length_mm,area_mm2,depth_mm);

for i = 1:441
    Y_table.(sprintf('P%d',i)) = Y_pressure(:,i);
end

writetable(X_table, fullfile('results','X_data.csv'));
writetable(Y_table, fullfile('results','Y_data.csv'));

fprintf('测试数据已生成。\n');
fprintf('依次运行：\n');
fprintf('1. build_bayesian_database.m\n');
fprintf('2. predict_by_bayesian_inverse.m\n');
