%% inverse_fem_force_reconstruction_professional.m
% =========================================================================
% Professional inverse finite-element force reconstruction workflow
%
% Purpose
%   Reconstruct equivalent nodal forces from measured nodal displacements
%   using an Abaqus global stiffness matrix and zero-order Tikhonov
%   regularization.
%
% Governing equations
%   K * U = F
%   Delta_U = H * F_unknown
%   F_unknown = argmin ||H*F - Delta_U||_2^2 + lambda*||F||_2^2
%
% Numerical strategy
%   1) Assemble the sparse global stiffness matrix from Abaqus 5-column MTX.
%   2) Remove constrained and non-physical/empty DOFs.
%   3) Build the observation matrix H from the measurement side:
%         H = P_m * K^{-1} * P_v^T
%      For symmetric K, this avoids computing a full flexibility matrix.
%   4) Solve the regularized inverse problem with the dual Tikhonov form:
%         F = H' * (H*H' + lambda*I)^{-1} * Delta_U
%      This avoids forming the potentially huge matrix H'*H.
%   5) Validate the reconstructed force field by predicting measured
%      displacements and reporting residual errors.
%
% Node-reading strategy
%   If get_model.m is available on the MATLAB path, the script uses it.
%   Otherwise, a built-in lightweight Abaqus *Node parser is used.
%
% Expected input files
%   1) Abaqus INP: contains node IDs and coordinates.
%   2) Abaqus stiffness MTX, five columns:
%        [Node_I, DOF_I, Node_J, DOF_J, Value]
%   3) Measured displacement CSV:
%        [Node_ID, U1, U2, U3, ...]
%   4) Fixed-node CSV:
%        one or more columns containing fixed node IDs.
%   5) Optional candidate-force-node CSV:
%        node IDs where unknown force is allowed.
%   6) Optional known-force CSV:
%        [Node_ID, F1, F2, F3]
%
% Assumptions
%   - Linear static finite-element model.
%   - The exported stiffness matrix is symmetric after assembly.
%   - Fixed nodes listed in fixed_nodes.csv are fully constrained in the
%     translational DOFs specified by cfg.fixed_directions.
%   - Units are consistent across K, displacement data, geometry and force.
%
% Important modeling note
%   Using every model node and all three directions as unknown forces usually
%   creates a severely underdetermined inverse problem. For physically
%   meaningful reconstruction, restrict cfg.force_directions and preferably
%   provide cfg.candidate_nodes_file.
% =========================================================================

clear;
clc;
close all;

%% 1. Configuration
cfg = struct();

% --------------------------- Input files -------------------------------
cfg.inp_filepath         = 'element_nodes_2.inp';
cfg.stiffness_mtx_file   = 'get_matrix-3_STIF2.mtx';
cfg.measurement_file     = 'Abaqus_Nodal_U_and_Pressure.csv';
cfg.fixed_nodes_file     = 'fixed_nodes.csv';

% Optional: restrict unknown loads to physically plausible nodes.
% Set to '' to use all valid non-fixed model nodes.
cfg.candidate_nodes_file = '';

% Optional known nodal forces, format: [Node_ID, F1, F2, F3].
% Set to '' when no known nodal forces are present.
cfg.known_force_file     = '';

% ------------------------- Model dimensions ----------------------------
% 3D solid element: 3 translational DOFs per node.
cfg.dofs_per_node        = 3;

% Measured displacement components contained in columns 2:4 of the CSV.
cfg.measurement_directions = [1, 2, 3];

% Unknown force components to reconstruct.
% Use [3] for purely normal/Z-direction force reconstruction when justified.
cfg.force_directions     = [1, 2, 3];

% Fully constrained directions for every node in fixed_nodes.csv.
cfg.fixed_directions     = [1, 2, 3];

% -------------------------- Regularization -----------------------------
% lambda = lambda_relative * ||H||_2^2 when use_relative_lambda = true.
% This makes lambda scale with the forward operator instead of using a
% dimensionless hard-coded number with no relation to the matrix magnitude.
cfg.use_relative_lambda  = true;
cfg.lambda_relative      = 1e-4;

% Used only when cfg.use_relative_lambda = false.
cfg.lambda_absolute      = 1e-6;

% ---------------------------- Diagnostics ------------------------------
cfg.run_condest          = true;
cfg.condest_warning      = 1e15;
cfg.symmetry_tolerance   = 1e-10;

% ----------------------------- Output ----------------------------------
cfg.output_force_csv     = 'reconstructed_nodal_forces.csv';
cfg.output_fit_csv       = 'measured_displacement_fit.csv';
cfg.print_limit          = 20;
cfg.make_plot            = true;
cfg.plot_direction       = 3;  % Color by reconstructed F3 when available.

fprintf('\n============================================================\n');
fprintf('  Inverse FEM force reconstruction - professional workflow\n');
fprintf('============================================================\n');

%% 2. Validate configuration and required dependencies
validateConfiguration(cfg);

assertFileExists(cfg.inp_filepath,       'Abaqus INP file');
assertFileExists(cfg.stiffness_mtx_file, 'stiffness MTX file');
assertFileExists(cfg.measurement_file,   'measurement file');
assertFileExists(cfg.fixed_nodes_file,   'fixed-node file');

if ~isempty(cfg.candidate_nodes_file)
    assertFileExists(cfg.candidate_nodes_file, 'candidate-force-node file');
end

if ~isempty(cfg.known_force_file)
    assertFileExists(cfg.known_force_file, 'known-force file');
end

%% 3. Read model nodes and coordinates
fprintf('\n[1/9] Reading model nodes and coordinates...\n');
if exist('get_model', 'file') == 2
    [node_ids, node_coords] = get_model(cfg.inp_filepath);
    fprintf('      Node reader: external get_model.m\n');
else
    [node_ids, node_coords] = readAbaqusNodeCoordinates(cfg.inp_filepath);
    fprintf('      Node reader: built-in Abaqus *Node parser\n');
end

node_ids = node_ids(:);

if isempty(node_ids) || isempty(node_coords)
    error('get_model returned an empty node list or coordinate matrix.');
end

if size(node_coords, 1) ~= numel(node_ids)
    error('node_ids and node_coords must contain the same number of nodes.');
end

if numel(unique(node_ids)) ~= numel(node_ids)
    error('Duplicate node IDs were detected in node_ids.');
end

if size(node_coords, 2) < 2 || size(node_coords, 2) > 3
    error('node_coords must have either 2 columns (2D) or 3 columns (3D).');
end

% Internally pad 2D coordinates with Z = 0 for robust output and plotting.
if size(node_coords, 2) == 2
    node_coords = [node_coords, zeros(size(node_coords, 1), 1)];
end

fprintf('      Model nodes: %d\n', numel(node_ids));

%% 4. Assemble the Abaqus global stiffness matrix
fprintf('\n[2/9] Assembling global stiffness matrix...\n');
[K_raw, stiffness_info] = assembleAbaqusStiffness( ...
    cfg.stiffness_mtx_file, cfg.dofs_per_node);

fprintf('      Raw matrix size: %d x %d\n', size(K_raw, 1), size(K_raw, 2));
fprintf('      Nonzero entries: %d\n', nnz(K_raw));
fprintf('      MTX storage detected: %s\n', stiffness_info.storage_mode);
fprintf('      Relative symmetry error: %.3e\n', stiffness_info.symmetry_error);

if stiffness_info.symmetry_error > cfg.symmetry_tolerance
    warning(['The assembled stiffness matrix is not sufficiently symmetric. ' ...
             'The matrix has been symmetrized, but the export format should be verified.']);
end

%% 5. Read constraints and construct the reduced active system
fprintf('\n[3/9] Applying constraints and removing empty DOFs...\n');

fixed_nodes = readNodeIdFile(cfg.fixed_nodes_file);
fixed_nodes = intersect(fixed_nodes, node_ids, 'stable');

if isempty(fixed_nodes)
    error(['No valid fixed nodes were found. The reduced stiffness matrix may ' ...
           'contain rigid-body modes and become singular.']);
end

constrained_global_dofs = expandNodeDofs( ...
    fixed_nodes, cfg.fixed_directions, cfg.dofs_per_node);

% Detect DOFs that actually participate in the assembled stiffness matrix.
row_has_stiffness = full(sum(spones(K_raw), 2)) > 0;
col_has_stiffness = full(sum(spones(K_raw), 1)).' > 0;
existing_global_dofs = find(row_has_stiffness | col_has_stiffness);

active_global_dofs = setdiff( ...
    existing_global_dofs, constrained_global_dofs, 'stable');

if isempty(active_global_dofs)
    error('No active DOFs remain after applying constraints.');
end

K_global = K_raw(active_global_dofs, active_global_dofs);
K_global = sparse((K_global + K_global.') / 2);  % Numerical symmetry cleanup.

fprintf('      Valid fixed nodes: %d\n', numel(fixed_nodes));
fprintf('      Constrained DOFs: %d\n', numel(constrained_global_dofs));
fprintf('      Active DOFs: %d\n', numel(active_global_dofs));
fprintf('      Reduced K size: %d x %d\n', size(K_global, 1), size(K_global, 2));

if any(~isfinite(nonzeros(K_global)))
    error('K_global contains NaN or Inf values.');
end

if cfg.run_condest
    fprintf('      Estimating condition number of K_global...\n');
    cond_K = condest(K_global);
    fprintf('      condest(K_global): %.3e\n', cond_K);

    if ~isfinite(cond_K) || cond_K > cfg.condest_warning
        warning(['K_global is singular or severely ill-conditioned. Verify ' ...
                 'boundary conditions, disconnected parts, reference nodes and MPC constraints.']);
    end
else
    cond_K = NaN;
end

% Try a Cholesky factorization check. A linear-elastic constrained stiffness
% matrix is normally symmetric positive definite when no unsupported modes remain.
[~, chol_flag] = chol(K_global);
if chol_flag ~= 0
    warning(['K_global is not strictly positive definite according to CHOL. ' ...
             'The direct backslash solver may still work, but the model should be checked.']);
end

%% 6. Read and validate measured displacements
fprintf('\n[4/9] Reading and validating measured displacements...\n');

measurement_data = readmatrix(cfg.measurement_file);

required_measurement_columns = 1 + max(cfg.measurement_directions);
if size(measurement_data, 2) < required_measurement_columns
    error(['Measurement file does not contain enough columns. Expected Node_ID ' ...
           'plus the requested displacement components.']);
end

id_measured_all = measurement_data(:, 1);
U_measured_all = measurement_data(:, 1 + cfg.measurement_directions);

valid_measurement_rows = all(isfinite([id_measured_all, U_measured_all]), 2);
id_measured_all = id_measured_all(valid_measurement_rows);
U_measured_all = U_measured_all(valid_measurement_rows, :);

if isempty(id_measured_all)
    error('No finite measured displacement records were found.');
end

% Keep only nodes that exist in the model and whose requested measurement
% DOFs all survive in the active reduced system.
[measurement_node_mask, measurement_global_dof_matrix] = validateNodeDirections( ...
    id_measured_all, cfg.measurement_directions, cfg.dofs_per_node, ...
    node_ids, active_global_dofs);

id_measured = id_measured_all(measurement_node_mask);
U_measured_mat = U_measured_all(measurement_node_mask, :);
measurement_global_dof_matrix = measurement_global_dof_matrix(measurement_node_mask, :);

if isempty(id_measured)
    error('No valid measured nodes remain after DOF and model filtering.');
end

if numel(unique(id_measured)) ~= numel(id_measured)
    error(['Duplicate measured node IDs were detected. Combine or remove duplicate ' ...
           'measurement records before reconstruction.']);
end

% Node-major ordering:
% [U(node1,d1); U(node1,d2); ...; U(node2,d1); ...]
U_measured = reshape(U_measured_mat.', [], 1);
measurement_global_dofs = reshape(measurement_global_dof_matrix.', [], 1);

[is_mapped_m, dof_measured] = ismember(measurement_global_dofs, active_global_dofs);
if ~all(is_mapped_m)
    error('Internal error: at least one measured DOF failed to map into active_global_dofs.');
end

fprintf('      Valid measured nodes: %d\n', numel(id_measured));
fprintf('      Measured displacement DOFs: %d\n', numel(dof_measured));

%% 7. Define candidate unknown-force nodes and directions
fprintf('\n[5/9] Building the candidate unknown-force space...\n');

if isempty(cfg.candidate_nodes_file)
    candidate_nodes_all = node_ids;
    fprintf('      Candidate-node file not supplied: using all model nodes initially.\n');
else
    candidate_nodes_all = readNodeIdFile(cfg.candidate_nodes_file);
end

% Remove fixed nodes before DOF validation.
candidate_nodes_all = setdiff(candidate_nodes_all, fixed_nodes, 'stable');

[candidate_node_mask, candidate_global_dof_matrix] = validateNodeDirections( ...
    candidate_nodes_all, cfg.force_directions, cfg.dofs_per_node, ...
    node_ids, active_global_dofs);

id_virtual = candidate_nodes_all(candidate_node_mask);
candidate_global_dof_matrix = candidate_global_dof_matrix(candidate_node_mask, :);

if isempty(id_virtual)
    error('No valid candidate force nodes remain after filtering.');
end

candidate_global_dofs = reshape(candidate_global_dof_matrix.', [], 1);

[is_mapped_v, dof_virtual] = ismember(candidate_global_dofs, active_global_dofs);
if ~all(is_mapped_v)
    error('Internal error: at least one candidate force DOF failed to map into active_global_dofs.');
end

% Coordinates corresponding exactly to the retained candidate node order.
[tf_virtual_coords, loc_virtual_coords] = ismember(id_virtual, node_ids);
if ~all(tf_virtual_coords)
    error('Internal error: a retained candidate node has no coordinate entry.');
end
coords_virtual = node_coords(loc_virtual_coords, :);

fprintf('      Candidate force nodes: %d\n', numel(id_virtual));
fprintf('      Unknown force DOFs: %d\n', numel(dof_virtual));
fprintf('      Measurement/unknown ratio: %.6f\n', ...
    numel(dof_measured) / numel(dof_virtual));

if numel(dof_virtual) > 10 * numel(dof_measured)
    warning(['The inverse problem is strongly underdetermined: unknown force DOFs ' ...
             'greatly exceed measured displacement DOFs. Restrict candidate nodes ' ...
             'and/or force directions for better physical identifiability.']);
end

%% 8. Compute known-load baseline displacement, if supplied
fprintf('\n[6/9] Computing baseline displacement from known loads...\n');

F_known_active = sparse(size(K_global, 1), 1);

if isempty(cfg.known_force_file)
    U_from_known_force = zeros(size(U_measured));
    fprintf('      No known-force file supplied; baseline displacement set to zero.\n');
else
    known_force_data = readmatrix(cfg.known_force_file);

    if size(known_force_data, 2) < 1 + cfg.dofs_per_node
        error('Known-force file must contain [Node_ID, F1, F2, F3].');
    end

    known_node_ids = known_force_data(:, 1);
    known_force_values = known_force_data(:, 2:(1 + cfg.dofs_per_node));

    finite_known_rows = all(isfinite([known_node_ids, known_force_values]), 2);
    known_node_ids = known_node_ids(finite_known_rows);
    known_force_values = known_force_values(finite_known_rows, :);

    % Assemble all valid known nodal force components into the reduced vector.
    for d = 1:cfg.dofs_per_node
        global_dofs_d = cfg.dofs_per_node * (known_node_ids - 1) + d;
        [tf_d, active_idx_d] = ismember(global_dofs_d, active_global_dofs);

        if any(tf_d)
            F_known_active = F_known_active + sparse( ...
                active_idx_d(tf_d), ...
                ones(nnz(tf_d), 1), ...
                known_force_values(tf_d, d), ...
                size(K_global, 1), 1);
        end
    end

    U_known_full = K_global \ F_known_active;
    U_from_known_force = U_known_full(dof_measured);

    fprintf('      Known-load vector nonzeros: %d\n', nnz(F_known_active));
end

Delta_U = U_measured - U_from_known_force;

if any(~isfinite(Delta_U))
    error('Delta_U contains NaN or Inf values.');
end

%% 9. Build the sensitivity matrix H efficiently
fprintf('\n[7/9] Building sensitivity matrix H from the measurement side...\n');

N_total_dof = size(K_global, 1);
N_measured_dof = numel(dof_measured);
N_virtual_dof = numel(dof_virtual);

% Estimate memory required by the dense sensitivity matrix H. MATLAB's
% direct sparse solve usually returns a dense response matrix here.
estimated_H_gb = double(N_measured_dof) * double(N_virtual_dof) * 8 / 1024^3;
fprintf('      Estimated dense H storage: %.3f GB\n', estimated_H_gb);
if estimated_H_gb > 4
    warning(['The sensitivity matrix alone may require more than 4 GB of memory. ' ...
             'Reduce candidate force nodes and/or reconstructed directions.']);
end

% Each column applies one unit load at one measured DOF.
E_measured = sparse( ...
    dof_measured, ...
    (1:N_measured_dof).', ...
    1, ...
    N_total_dof, ...
    N_measured_dof);

% Solve K * X = E_m. For symmetric K:
%   H = P_m*K^{-1}*P_v' = (P_v*K^{-1}*P_m')'
X_measured = K_global \ E_measured;
H = X_measured(dof_virtual, :).';

clear X_measured E_measured;

if any(~isfinite(H(:)))
    error('Sensitivity matrix H contains NaN or Inf values.');
end

fprintf('      H size: %d x %d\n', size(H, 1), size(H, 2));

%% 10. Select regularization strength and solve the inverse problem
fprintf('\n[8/9] Solving the Tikhonov-regularized inverse problem...\n');

if cfg.use_relative_lambda
    H_norm_est = normest(H);
    lambda = cfg.lambda_relative * H_norm_est^2;
else
    H_norm_est = NaN;
    lambda = cfg.lambda_absolute;
end

if ~isfinite(lambda) || lambda <= 0
    error('Regularization parameter lambda must be finite and strictly positive.');
end

fprintf('      lambda: %.6e\n', lambda);
if isfinite(H_norm_est)
    fprintf('      Estimated ||H||_2: %.6e\n', H_norm_est);
end

% Dual zero-order Tikhonov solution:
%   F = H' * (H*H' + lambda*I)^(-1) * Delta_U
% This avoids the large N_virtual_dof x N_virtual_dof matrix H'*H.
A_dual = H * H.' + lambda * speye(N_measured_dof);
dual_variable = A_dual \ Delta_U;
F_virtual_reconstructed = H.' * dual_variable;

if any(~isfinite(F_virtual_reconstructed))
    error('Reconstructed force vector contains NaN or Inf values.');
end

%% 11. Validate reconstruction against measured displacements
U_predicted = H * F_virtual_reconstructed + U_from_known_force;
residual = U_measured - U_predicted;

absolute_residual = norm(residual, 2);
relative_residual = absolute_residual / max(norm(U_measured, 2), eps);
force_norm = norm(F_virtual_reconstructed, 2);

fprintf('\n      Reconstruction diagnostics\n');
fprintf('      --------------------------\n');
fprintf('      ||U_measured - U_predicted||_2 : %.6e\n', absolute_residual);
fprintf('      Relative displacement residual : %.6e\n', relative_residual);
fprintf('      ||F_reconstructed||_2           : %.6e\n', force_norm);

if relative_residual > 0.20
    warning(['The displacement fit residual exceeds 20%%. Possible causes include ' ...
             'an incorrect candidate load region, inadequate boundary conditions, ' ...
             'unit inconsistency, model mismatch, or excessive regularization.']);
end

%% 12. Reconstruct node-wise force components
fprintf('\n[9/9] Organizing and exporting results...\n');

% Use NaN for directions that were not reconstructed, avoiding the false
% implication that an unestimated component is physically zero.
F_nodes = nan(numel(id_virtual), cfg.dofs_per_node);

F_virtual_matrix = reshape( ...
    F_virtual_reconstructed, numel(cfg.force_directions), []).';
F_nodes(:, cfg.force_directions) = F_virtual_matrix;

force_magnitude = sqrt(sum(F_nodes(:, cfg.force_directions).^2, 2));

force_table = table( ...
    id_virtual, ...
    coords_virtual(:, 1), ...
    coords_virtual(:, 2), ...
    coords_virtual(:, 3), ...
    F_nodes(:, 1), ...
    F_nodes(:, 2), ...
    F_nodes(:, 3), ...
    force_magnitude, ...
    'VariableNames', { ...
        'NodeID', 'X', 'Y', 'Z', 'F1', 'F2', 'F3', 'ForceMagnitude'});

writetable(force_table, cfg.output_force_csv);

% Measurement fit output in the same node-major/component-major order.
U_predicted_mat = reshape( ...
    U_predicted, numel(cfg.measurement_directions), []).';
residual_mat = reshape( ...
    residual, numel(cfg.measurement_directions), []).';

fit_table = table(id_measured, 'VariableNames', {'NodeID'});
for k = 1:numel(cfg.measurement_directions)
    d = cfg.measurement_directions(k);
    fit_table.(sprintf('U%d_Measured', d)) = U_measured_mat(:, k);
    fit_table.(sprintf('U%d_Predicted', d)) = U_predicted_mat(:, k);
    fit_table.(sprintf('U%d_Residual', d)) = residual_mat(:, k);
end

writetable(fit_table, cfg.output_fit_csv);

fprintf('      Force results saved to: %s\n', cfg.output_force_csv);
fprintf('      Measurement fit saved to: %s\n', cfg.output_fit_csv);

%% 13. Print a compact result preview
fprintf('\n============================================================\n');
fprintf('                  Reconstruction completed\n');
fprintf('============================================================\n');

print_limit = min(cfg.print_limit, height(force_table));

for i = 1:print_limit
    fprintf('Node %d  Coord [%.6g, %.6g, %.6g]\n', ...
        force_table.NodeID(i), ...
        force_table.X(i), force_table.Y(i), force_table.Z(i));

    for d = 1:cfg.dofs_per_node
        if ismember(d, cfg.force_directions)
            fprintf('   F%d = % .6e\n', d, F_nodes(i, d));
        end
    end
end

if height(force_table) > print_limit
    fprintf('... %d additional nodes omitted from console output.\n', ...
        height(force_table) - print_limit);
end

%% 14. Visualization
if cfg.make_plot
    plot_direction = cfg.plot_direction;

    if ~ismember(plot_direction, cfg.force_directions)
        warning(['cfg.plot_direction is not among cfg.force_directions. ' ...
                 'The force-field plot was skipped.']);
    else
        color_data = F_nodes(:, plot_direction);

        [tf_fixed_coords, loc_fixed_coords] = ismember(fixed_nodes, node_ids);
        coords_fixed = node_coords(loc_fixed_coords(tf_fixed_coords), :);

        figure( ...
            'Name', 'Inverse FEM reconstructed nodal force field', ...
            'Color', 'w', ...
            'Position', [100, 100, 900, 700]);
        hold on;
        grid on;

        if ~isempty(coords_fixed)
            scatter3( ...
                coords_fixed(:, 1), coords_fixed(:, 2), coords_fixed(:, 3), ...
                15, [0.7, 0.7, 0.7], 'filled', ...
                'MarkerEdgeColor', 'none', ...
                'DisplayName', 'Fixed boundary nodes');
        end

        scatter3( ...
            coords_virtual(:, 1), coords_virtual(:, 2), coords_virtual(:, 3), ...
            40, color_data, 'filled', ...
            'MarkerEdgeColor', 'none', ...
            'DisplayName', sprintf('Reconstructed F%d', plot_direction));

        axis equal;
        view(45, 30);
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        title(sprintf('Inverse FEM reconstructed nodal force component F%d', ...
            plot_direction));

        colormap('jet');
        cb = colorbar;
        ylabel(cb, sprintf('Reconstructed nodal force F%d', plot_direction));

        finite_color_data = color_data(isfinite(color_data));
        if isempty(finite_color_data)
            cmax = NaN;
        else
            cmax = max(abs(finite_color_data));
        end

        if isfinite(cmax) && cmax > 0
            caxis([-cmax, cmax]);
        end

        legend('Location', 'best');
        hold off;
    end
end

fprintf('\nDone.\n');

%% ========================================================================
% Local functions
% ========================================================================

function validateConfiguration(cfg)
%VALIDATECONFIGURATION Validate configuration values before any file I/O.

    if ~isscalar(cfg.dofs_per_node) || cfg.dofs_per_node < 1 || ...
            cfg.dofs_per_node ~= floor(cfg.dofs_per_node)
        error('cfg.dofs_per_node must be a positive integer.');
    end

    validateDirectionVector( ...
        cfg.measurement_directions, cfg.dofs_per_node, ...
        'cfg.measurement_directions');

    validateDirectionVector( ...
        cfg.force_directions, cfg.dofs_per_node, ...
        'cfg.force_directions');

    validateDirectionVector( ...
        cfg.fixed_directions, cfg.dofs_per_node, ...
        'cfg.fixed_directions');

    if cfg.use_relative_lambda
        if ~isscalar(cfg.lambda_relative) || ...
                ~isfinite(cfg.lambda_relative) || cfg.lambda_relative <= 0
            error('cfg.lambda_relative must be finite and > 0.');
        end
    else
        if ~isscalar(cfg.lambda_absolute) || ...
                ~isfinite(cfg.lambda_absolute) || cfg.lambda_absolute <= 0
            error('cfg.lambda_absolute must be finite and > 0.');
        end
    end
end

function validateDirectionVector(directions, dofs_per_node, variable_name)
%VALIDATEDIRECTIONVECTOR Validate requested local DOF/component indices.

    if isempty(directions) || ~isvector(directions) || ...
            any(~isfinite(directions)) || ...
            any(directions ~= floor(directions)) || ...
            any(directions < 1) || any(directions > dofs_per_node) || ...
            numel(unique(directions)) ~= numel(directions)
        error('%s contains invalid or duplicate DOF directions.', variable_name);
    end
end

function assertFileExists(filepath, description)
%ASSERTFILEEXISTS Raise a descriptive error if an input file is missing.

    if exist(filepath, 'file') ~= 2
        error('%s not found: %s', description, filepath);
    end
end

function node_ids = readNodeIdFile(filepath)
%READNODEIDFILE Read node IDs from one or more numeric columns.

    data = readmatrix(filepath);
    data = data(isfinite(data));
    node_ids = unique(data(:), 'stable');

    if isempty(node_ids)
        error('No finite node IDs were found in file: %s', filepath);
    end

    if any(node_ids ~= floor(node_ids)) || any(node_ids < 1)
        error('Node IDs in %s must be positive integers.', filepath);
    end
end

function global_dofs = expandNodeDofs(node_ids, directions, dofs_per_node)
%EXPANDNODEDOFS Expand node IDs to global DOF numbers in node-major order.
%
% Example with 3 DOFs/node:
%   node_ids = [10; 20], directions = [1 3]
%   output = [28; 30; 58; 60]

    node_ids = node_ids(:);
    directions = directions(:).';

    dof_matrix = dofs_per_node * (node_ids - 1) + directions;
    global_dofs = reshape(dof_matrix.', [], 1);
end

function [valid_node_mask, global_dof_matrix] = validateNodeDirections( ...
        node_ids_to_check, directions, dofs_per_node, model_node_ids, active_global_dofs)
%VALIDATENODEDIRECTIONS Retain nodes whose requested DOFs all exist and are active.

    node_ids_to_check = node_ids_to_check(:);
    directions = directions(:).';

    node_exists = ismember(node_ids_to_check, model_node_ids);

    global_dof_matrix = ...
        dofs_per_node * (node_ids_to_check - 1) + directions;

    dofs_are_active = ismember(global_dof_matrix, active_global_dofs);
    valid_node_mask = node_exists & all(dofs_are_active, 2);
end

function [K_raw, info] = assembleAbaqusStiffness(filepath, dofs_per_node)
%ASSEMBLEABAQUSSTIFFNESS Assemble a symmetric sparse matrix from Abaqus MTX.
%
% Expected numeric columns:
%   [Node_I, DOF_I, Node_J, DOF_J, Value]
%
% The function detects whether the file contains:
%   - only upper-triangular entries,
%   - only lower-triangular entries, or
%   - both triangles / a complete matrix.
%
% Duplicate entries are automatically accumulated by sparse().

    % .mtx is a plain-text numeric file, but MATLAB does not always infer
    % the file type from the .mtx extension. Explicitly force text import.
    % This avoids: "'.mtx' is an unrecognized file extension".
    try
        mtx_data = readmatrix(filepath, 'FileType', 'text');
    catch ME_readmatrix
        % Compatibility fallback for older MATLAB releases or unusual
        % delimiter settings. Abaqus matrix exports are normally comma- or
        % whitespace-delimited text.
        fid = fopen(filepath, 'r');
        if fid < 0
            error('Unable to open stiffness MTX file: %s', filepath);
        end
        cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

        raw = textscan(fid, '%f%f%f%f%f', ...
            'Delimiter', {',', ' ', '\t'}, ...
            'MultipleDelimsAsOne', true, ...
            'CollectOutput', true, ...
            'CommentStyle', {'**', '#'});

        mtx_data = raw{1};

        if isempty(mtx_data)
            error('Failed to read numeric data from MTX file %s. Original readmatrix error: %s', ...
                filepath, ME_readmatrix.message);
        end
    end

    if size(mtx_data, 2) < 5
        error(['Stiffness MTX must contain at least five numeric columns: ' ...
               '[Node_I, DOF_I, Node_J, DOF_J, Value].']);
    end

    mtx_data = mtx_data(:, 1:5);
    finite_rows = all(isfinite(mtx_data), 2);
    mtx_data = mtx_data(finite_rows, :);

    if isempty(mtx_data)
        error('No valid numeric stiffness entries were found in %s.', filepath);
    end

    node_i = mtx_data(:, 1);
    dof_i  = mtx_data(:, 2);
    node_j = mtx_data(:, 3);
    dof_j  = mtx_data(:, 4);
    values = mtx_data(:, 5);

    if any(node_i < 1 | node_j < 1 | ...
           node_i ~= floor(node_i) | node_j ~= floor(node_j))
        error('MTX node IDs must be positive integers.');
    end

    if any(dof_i < 1 | dof_i > dofs_per_node | ...
           dof_j < 1 | dof_j > dofs_per_node | ...
           dof_i ~= floor(dof_i) | dof_j ~= floor(dof_j))
        error(['MTX DOF indices are inconsistent with cfg.dofs_per_node = %d.'], ...
            dofs_per_node);
    end

    % Remove explicit numerical zeros before detecting participating DOFs.
    nonzero_rows = values ~= 0;
    node_i = node_i(nonzero_rows);
    dof_i  = dof_i(nonzero_rows);
    node_j = node_j(nonzero_rows);
    dof_j  = dof_j(nonzero_rows);
    values = values(nonzero_rows);

    row_idx = dofs_per_node * (node_i - 1) + dof_i;
    col_idx = dofs_per_node * (node_j - 1) + dof_j;

    max_dof = max([row_idx; col_idx]);
    K_entries = sparse(row_idx, col_idx, values, max_dof, max_dof);

    has_upper = nnz(triu(K_entries, 1)) > 0;
    has_lower = nnz(tril(K_entries, -1)) > 0;

    if has_upper && has_lower
        % The file appears to contain both triangles. Do not double the
        % off-diagonal terms; instead enforce symmetry by averaging.
        storage_mode = 'full/both triangles';
        K_raw = (K_entries + K_entries.') / 2;
    elseif has_upper
        storage_mode = 'upper triangular';
        K_raw = K_entries + K_entries.' - ...
            spdiags(diag(K_entries), 0, max_dof, max_dof);
    elseif has_lower
        storage_mode = 'lower triangular';
        K_raw = K_entries + K_entries.' - ...
            spdiags(diag(K_entries), 0, max_dof, max_dof);
    else
        storage_mode = 'diagonal only';
        K_raw = K_entries;
    end

    K_raw = sparse(K_raw);

    denominator = max(norm(K_raw, 'fro'), eps);
    symmetry_error = norm(K_raw - K_raw.', 'fro') / denominator;

    info = struct();
    info.storage_mode = storage_mode;
    info.symmetry_error = symmetry_error;
end

function [node_ids, node_coords] = readAbaqusNodeCoordinates(filepath)
%READABAQUSNODECOORDINATES Read simple Abaqus *Node sections from an INP file.
%
% Supported node records:
%   node_id, x, y
%   node_id, x, y, z
%
% The parser can read multiple *Node sections as long as node IDs are globally
% unique. It does not apply assembly instance transformations. For complex
% assembly models, provide a project-specific get_model.m on the MATLAB path.

    fid = fopen(filepath, 'r');
    if fid < 0
        error('Unable to open Abaqus INP file: %s', filepath);
    end

    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    node_ids = zeros(0, 1);
    node_coords = zeros(0, 3);
    inside_node_block = false;

    while ~feof(fid)
        line = strtrim(fgetl(fid));

        if isempty(line) || startsWith(line, '**')
            continue;
        end

        if startsWith(line, '*')
            inside_node_block = startsWith(lower(line), '*node');
            continue;
        end

        if ~inside_node_block
            continue;
        end

        numeric_line = strrep(line, ',', ' ');
        values = sscanf(numeric_line, '%f').';

        if numel(values) < 3
            warning('Skipping malformed node line: %s', line);
            continue;
        end

        node_id = values(1);
        coords = values(2:end);

        if node_id < 1 || node_id ~= floor(node_id)
            error('Invalid Abaqus node ID encountered: %g', node_id);
        end

        if numel(coords) == 2
            coords = [coords, 0];
        elseif numel(coords) >= 3
            coords = coords(1:3);
        end

        node_ids(end + 1, 1) = node_id; %#ok<AGROW>
        node_coords(end + 1, :) = coords; %#ok<AGROW>
    end

    if isempty(node_ids)
        error('No Abaqus *Node records were found in %s.', filepath);
    end

    if numel(unique(node_ids)) ~= numel(node_ids)
        error(['Duplicate node IDs were found while parsing the INP file. ' ...
               'For assembly models with repeated part-local node IDs, use a ' ...
               'project-specific get_model.m instead of the built-in parser.']);
    end
end

