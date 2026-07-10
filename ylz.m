% 虚拟点与真实力源之间的位移
h = X_virtual - source_position;

r = sqrt(sum(h.^2, 2));

eps_r = 1e-9;
r_hat = h ./ max(r, eps_r);

% source_direction 和 virtual_direction 均为二维单位向量
ni = source_direction;
nj = virtual_direction;

% 二维叉积的 z 分量
cross_i = ni(1) .* r_hat(:,2) - ni(2) .* r_hat(:,1);
cross_j = nj(:,1) .* r_hat(:,2) - nj(:,2) .* r_hat(:,1);

direction_difference = ni - nj;

a = cross_i .* cross_j ...
    + beta .* sum(direction_difference .* r_hat, 2) ...
    - beta^2;

% 建议采用正值方向映射
direction_kernel = exp(gamma .* (a - 1));

% 径向有限支撑核
radial_kernel = zeros(size(r));

near_region = r <= r_min;
transition_region = r > r_min & r < r_c;

radial_kernel(near_region) = 1;

radial_kernel(transition_region) = ...
    cos( ...
        pi .* (r(transition_region) - r_min) ./ ...
        (2 .* (r_c - r_min)) ...
    ).^(2 .* zeta);

% YLZ 启发的各向异性核
K = radial_kernel .* direction_kernel;

% 将核最大值归一化为 1
K = K ./ (max(K) + eps);

% 单点力源预测
F_pred = source_force .* K;