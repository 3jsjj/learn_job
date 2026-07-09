function result = bayesian_inverse_predict(observed_sensor, X_database, Y_database, sigma_sensor, prior_probability)
%BAYESIAN_INVERSE_PREDICT 离散贝叶斯反演核心函数

observed_sensor = observed_sensor(:)';
[N, n_sensor] = size(X_database);

if numel(observed_sensor) ~= n_sensor
    error('观测输入维度与数据库传感器维度不一致。');
end

sigma_sensor = max(sigma_sensor(:)', eps);
prior_probability = prior_probability(:);
prior_probability = prior_probability / sum(prior_probability);

residual = X_database - observed_sensor;
normalized_residual = residual ./ sigma_sensor;
squared_mahalanobis = sum(normalized_residual.^2, 2);

log_likelihood = -0.5 .* squared_mahalanobis ...
    - sum(log(sigma_sensor)) ...
    - 0.5 .* n_sensor .* log(2*pi);

log_prior = log(max(prior_probability, realmin));
log_post = log_likelihood + log_prior;

max_log = max(log_post);
posterior = exp(log_post - max_log);
posterior = posterior / sum(posterior);

[~, map_index] = max(posterior);
map_pressure = Y_database(map_index, :);

posterior_mean_pressure = posterior' * Y_database;

diff_pressure = Y_database - posterior_mean_pressure;
posterior_variance = posterior' * (diff_pressure.^2);
posterior_std_pressure = sqrt(max(posterior_variance, 0));

p_nonzero = posterior(posterior > 0);
posterior_entropy = -sum(p_nonzero .* log(p_nonzero));
effective_sample_size = 1 / sum(posterior.^2);

result.posterior_probability = posterior;
result.map_index = map_index;
result.map_pressure = max(map_pressure, 0);
result.posterior_mean_pressure = max(posterior_mean_pressure, 0);
result.posterior_std_pressure = posterior_std_pressure;
result.posterior_entropy = posterior_entropy;
result.effective_sample_size = effective_sample_size;
result.squared_mahalanobis = squared_mahalanobis;
end
