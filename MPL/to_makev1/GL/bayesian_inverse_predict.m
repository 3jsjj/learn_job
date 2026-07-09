function result = bayesian_inverse_predict( ...
    observed_sensor, ...
    X_database, ...
    Y_database, ...
    sigma_sensor, ...
    prior_probability)
%BAYESIAN_INVERSE_PREDICT
%
% 根据当前 9 个传感器观测值，在候选数据库中进行离散贝叶斯反演。
%
% 输入：
% observed_sensor      1×9 当前传感器值
% X_database           N×9 候选样本传感器响应
% Y_database           N×441 候选样本力场
% sigma_sensor         1×9 传感器噪声标准差
% prior_probability    N×1 候选样本先验概率
%
% 输出：
% result.posterior_probability
% result.map_index
% result.map_pressure
% result.posterior_mean_pressure
% result.posterior_std_pressure
% result.posterior_entropy
% result.effective_sample_size

    observed_sensor = observed_sensor(:)';

    [N, n_sensor] = size(X_database);

    if numel(observed_sensor) ~= n_sensor
        error('观测输入维度与数据库传感器维度不一致。');
    end

    sigma_sensor = sigma_sensor(:)';

    if numel(sigma_sensor) ~= n_sensor
        error('sigma_sensor 维度不正确。');
    end

    sigma_sensor = max(sigma_sensor, eps);

    prior_probability = prior_probability(:);

    if numel(prior_probability) ~= N
        error('prior_probability 长度必须等于候选样本数。');
    end

    if any(prior_probability < 0)
        error('先验概率不能为负。');
    end

    prior_probability = prior_probability / sum(prior_probability);

    %% ========================================================
    %  1. 计算高斯似然
    %
    %  假设每个传感器噪声相互独立，且服从高斯分布：
    %
    %  s_i = s_i_hat + noise
    %
    %  noise ~ N(0, sigma_i^2)
    %% ========================================================

    residual = X_database - observed_sensor;

    normalized_residual = residual ./ sigma_sensor;

    squared_mahalanobis = sum(normalized_residual.^2, 2);

    log_likelihood = ...
        -0.5 .* squared_mahalanobis ...
        -sum(log(sigma_sensor)) ...
        -0.5 .* n_sensor .* log(2*pi);

    %% ========================================================
    %  2. 加入先验概率
    %% ========================================================

    log_prior = log(max(prior_probability, realmin));

    log_posterior_unnormalized = log_likelihood + log_prior;

    %% ========================================================
    %  3. log-sum-exp 归一化，避免数值下溢
    %% ========================================================

    max_log = max(log_posterior_unnormalized);

    posterior_unnormalized = exp( ...
        log_posterior_unnormalized - max_log);

    posterior_probability = ...
        posterior_unnormalized / sum(posterior_unnormalized);

    %% ========================================================
    %  4. MAP 最优候选
    %% ========================================================

    [~, map_index] = max(posterior_probability);

    map_pressure = Y_database(map_index, :);

    %% ========================================================
    %  5. 后验均值力场
    %% ========================================================

    posterior_mean_pressure = posterior_probability' * Y_database;

    %% ========================================================
    %  6. 后验标准差力场
    %% ========================================================

    pressure_difference = Y_database - posterior_mean_pressure;

    posterior_variance = ...
        posterior_probability' * (pressure_difference.^2);

    posterior_std_pressure = sqrt(max(posterior_variance, 0));

    %% ========================================================
    %  7. 后验熵和有效候选数
    %% ========================================================

    p_nonzero = posterior_probability(posterior_probability > 0);

    posterior_entropy = -sum(p_nonzero .* log(p_nonzero));

    effective_sample_size = 1 / sum(posterior_probability.^2);

    %% ========================================================
    %  8. 输出结果
    %% ========================================================

    result.posterior_probability = posterior_probability;
    result.map_index = map_index;
    result.map_pressure = map_pressure;
    result.posterior_mean_pressure = posterior_mean_pressure;
    result.posterior_std_pressure = posterior_std_pressure;
    result.posterior_entropy = posterior_entropy;
    result.effective_sample_size = effective_sample_size;
    result.squared_mahalanobis = squared_mahalanobis;

end
