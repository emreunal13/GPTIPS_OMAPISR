clear all; clc

%% HELPER FUNCTIONS FOR DISTINCT VECTOR COMBINATIONS
% This block is dedicated for generating distinct dimensionless number
% combinations. In the pipe flow study, dimensionless numbers (Re and e/D)
% are directly chosen.

function normalized = normalize_vector(a)
    % Split into two groups
    group1 = a(1:3);
    group2 = a(4:6);

    % Sort each group
    group1 = sort(group1);
    group2 = sort(group2);

    % Sort both groups together to ensure uniqueness
    sorted_groups = sortrows([group1; group2]);

    % Flatten into a row vector
    normalized = sorted_groups(:)';
end

function distinct_combinations = generate_distinct_combinations(num_coefficients)
    if nargin < 1
        num_coefficients = 6;
    end

    % Generate all 0/1 combinations
    all_combinations = dec2bin(1:2^num_coefficients-1) - '0'; % Avoid all-zero case

    % Normalize and deduplicate
    normalized_combinations = zeros(size(all_combinations));
    for i = 1:size(all_combinations, 1)
        normalized_combinations(i, :) = normalize_vector(all_combinations(i, :));
    end

    % Remove duplicates and sort
    distinct_combinations = unique(normalized_combinations, 'rows');
end

% Function to calculate Pi terms
function [pi1, pi2] = calc_pi_terms(a, X, basis1_in, basis2_in, basis3_in)
    % Create mask for zero values
    non_zero_mask = X ~= 0;
    
    % Initialize X_log with zeros
    X_log = zeros(size(X));
    
    % Calculate log for non-zero values
    X_log(non_zero_mask) = log(X(non_zero_mask));
    
    % Calculate dimensionless numbers
    coef_pi1 = a(1) * basis1_in + a(2) * basis2_in + a(3) * basis3_in;
    coef_pi2 = a(4) * basis1_in + a(5) * basis2_in + a(6) * basis3_in;
    
    % Compute pi terms
    pi1_mat = exp(X_log * coef_pi1);
    pi2_mat = exp(X_log * coef_pi2);
    
    % Handle cases where e == 0 (fifth column)
    e_zero_mask = X(:,5) == 0;
    if any(e_zero_mask)
        pi2_mat(e_zero_mask) = 0;
    end
    
    % Ensure output is column vector
    pi1 = pi1_mat;
    pi2 = pi2_mat;
end

%clear all; clc

% [Previous helper functions: normalize_vector, generate_distinct_combinations, calc_pi_terms stay the same]

%% GENETIC PROCESS FUNCTION

% Modified main function to handle multiple runs
function gp = nikuradse_config(gp, combination_index)
    persistent all_data; % Keep data between calls
    persistent all_combinations; % Keep combinations between calls
    persistent results_storage; % Store results from all runs
    
    % Initialize storage and load data on first call
    if isempty(all_data)
        % Load data
        data = readmatrix('NIKURADSE_DATA_NONZERO_SUPERPIPE_INCH_CONVERTED.csv');
        
        % Store all variables
        all_data.rho = data(:, 1);
        all_data.mu = data(:, 2);
        all_data.D = data(:, 3);
        all_data.V = data(:, 4);
        all_data.e = data(:, 5);
        all_data.L = data(:, 6);
        all_data.delta_p = data(:, 10);
        
        % Setup dimensional matrix
        all_data.D_in = [1, 1, 0, 0, 0, 0;
                        -3, -1, 1, 1, 1, 1;
                        0, -1, 0, -1, 0, 0];
        
        % Split matrices
        Din1 = all_data.D_in(:, 1:3);
        Din2 = all_data.D_in(:, 4:end);
        
        % Calculate basis vectors
        x2_basis1 = [1; 0; 0];
        x1_basis1 = -inv(Din1) * Din2 * x2_basis1;
        all_data.basis1_in = [x1_basis1; x2_basis1];
        
        x2_basis2 = [0; 1; 0];
        x1_basis2 = -inv(Din1) * Din2 * x2_basis2;
        all_data.basis2_in = [x1_basis2; x2_basis2];
        
        x2_basis3 = [0; 0; 1];
        x1_basis3 = -inv(Din1) * Din2 * x2_basis3;
        all_data.basis3_in = [x1_basis3; x2_basis3];
        
        % Generate all combinations
        all_combinations = generate_distinct_combinations();
        results_storage = cell(size(all_combinations, 1), 1);
    end
    
    % Use the current combination
    %current_combination = all_combinations(combination_index, :);
    current_combination = [1 0 0 0 1 0]; % Shortcut basis vector

    
    % Calculate Pi terms for the current combination
    X = [all_data.rho all_data.mu all_data.D all_data.V all_data.e all_data.L];
    [Pi1, Pi2] = calc_pi_terms(current_combination, X, all_data.basis1_in, all_data.basis2_in, all_data.basis3_in);
    Pi3 = all_data.delta_p .* 2 .* all_data.D ./ (all_data.rho .* all_data.L .* all_data.V.^2);
    
    % Basic run control parameters
    gp.selection.tournament.p_pareto = 1;   %probability that a pareto tournament will be used for any given selection event.
    gp.treedef.max_mutate_depth = 6;
    gp = gprandom(gp, 2001);
    gp.runcontrol.pop_size = 750;
    gp.runcontrol.num_gen = 1350;
    gp.runcontrol.verbose = 1;
    gp.fitness.complexityMeasure = 0;
    gp.runcontrol.usecache = false;
    if isfield(gp,'fitness') && isfield(gp.fitness,'cache')
        try remove(gp.fitness.cache, keys(gp.fitness.cache)); catch, end
    end
    
    
    % Define function set
    gp.nodes.functions.name = {'times', 'plus', 'minus', 'rdivide', 'power', 'plog', 'exp', 'tanh'};
    
    % Define constants
    gp.nodes.const.num_dec_places = 6;
    gp.nodes.const.range = [-10 10];
    
    % Set maximum tree depth and nodes
    gp.treedef.max_depth = 8;
    gp.treedef.max_nodes = 16;
    gp.genes.max_genes = 2;
    gp.nodes.const.p_ERC = 0.3;
    
    % Store data in gp structure
    gp.userdata.xtrain = [Pi1 Pi2];
    gp.userdata.ytrain = Pi3;

    %–– store raw variables for penalty in fitness function
    gp.userdata.raw.rho = all_data.rho;
    gp.userdata.raw.V   = all_data.V;
    gp.userdata.raw.D   = all_data.D;
    gp.userdata.raw.L   = all_data.L;
    gp.userdata.raw.mu  = all_data.mu;
    

    gp.userdata.inputs = {'Pi1', 'Pi2'};
    gp.userdata.output = 'Pi3';
    gp.userdata.combination_index = combination_index;
    gp.userdata.current_combination = current_combination;
    
    % Fitness function configuration
    gp.fitness.fitfun = @regressmulti_fitfun;
    gp.fitness.minimization = true;
    gp.selection.elite_fraction = 0.15;
    gp.selection.tournament.size = 3;
    gp.operators.mutation.p_mutate = 0.43;
    gp.operators.crossover.p_cross = 0.55;
    gp.operators.directrepro.p_direct = 0.02;
    
    % Allow parallel processing if available
    gp.runcontrol.parallel.enable = true;
    
    % Store results for this combination
    results_storage{combination_index} = struct('combination', current_combination, ...
                                              'Pi1', Pi1, ...
                                              'Pi2', Pi2, ...
                                              'Pi3', Pi3);
    
    % Make results accessible
    gp.userdata.all_results = results_storage;
end
%% MAIN LOOP
% This block is created in order to run the model for each
% distinct dimensionless number combination. In the pipe flow study,
% it is run once for only (Re, e/D) combination.

% Main execution script
%num_combinations = size(generate_distinct_combinations(), 1);
num_combinations = 1; % Re and e/D

all_gp_results = cell(num_combinations, 1);
pareto_data = struct([]);

% Run for each combination
for i = 1:num_combinations
    fprintf('Processing combination %d of %d\n', i, num_combinations);
    gp = rungp(@(gp) nikuradse_config(gp, i));
    gp = fix_best_after_run(gp);
    all_gp_results{i} = gp;

    % % Store pareto front data for this run
    % try
    %     [pareto_front, pareto_indices] = paretofront(gp);
    % 
    %     % Store essential metrics for later comparison
    %     pareto_data(i).combination = gp.userdata.current_combination;
    %     pareto_data(i).fitnesses = gp.results.best.fitness;
    %     pareto_data(i).complexities = arrayfun(@(x) length(x.nodes), gp.pop);
    %     pareto_data(i).r2 = 1 - sum((gp.userdata.ytrain - evalfitness(gp.results.best.individual, gp)).^2) / ...
    %                            sum((gp.userdata.ytrain - mean(gp.userdata.ytrain)).^2);
    % catch ME
    %     fprintf('Warning: Error processing Pareto front for combination %d\n', i);
    %     fprintf('Error message: %s\n', ME.message);
    % end
end

% % Custom Pareto Front Aggregator
% function combined_pareto_report(all_gp_results)
%     num_runs = length(all_gp_results);
%     pareto_metrics = [];
% 
%     for i = 1:num_runs
%         gp = all_gp_results{i};
%         if isempty(gp), continue; end
% 
%         
%         try
%             report = paretoreport(gp);
%             for j = 1:length(report.individuals)
%                 ind = report.individuals(j);
%                 r2 = ind.r2;
%                 complexity = length(ind.nodes);
%                 combination = gp.userdata.current_combination;
%                 pareto_metrics = [pareto_metrics; complexity, r2, i, combination];
%             end
%         catch ME
%             fprintf('Warning: Failed paretoreport for run %d.\nError: %s\n', i, ME.message);
%         end
%     end
% 
%     % Plot combined Pareto fronts
%     figure;
%     scatter(pareto_metrics(:,1), pareto_metrics(:,2), 50, pareto_metrics(:,3), 'filled');
%     xlabel('Model Complexity (Nodes)');
%     ylabel('R² (Goodness of Fit)');
%     title('Combined Pareto Fronts Across All Basis Vector Combinations');
%     colorbar;
%     grid on;
%     legend('Basis Vector Runs', 'Location', 'best');
% 
%     % Find and display the top 5 models
%     [~, sortedIdx] = sort(pareto_metrics(:,2), 'descend');
%     top5 = pareto_metrics(sortedIdx(1:min(5,end)), :);
% 
%     fprintf('\nTop 5 Pareto Solutions:\n');
%     for k = 1:size(top5, 1)
%         fprintf('Run %d | R² = %.4f | Complexity = %d | Combination: [%s]\n', ...
%             top5(k,3), top5(k,2), top5(k,1), num2str(top5(k,4:end)));
%     end
% end

%%
% Convert the best SR model to a symbolic expression
modelSym = gpmodel2sym(gp, 'best');
modelSym = vpa(modelSym, 4); % Optional: Simplify the equation

% Load combined data (Nikuradse + smooth-pipe lab sets)
data = readtable("NIKURADSE_SUPERPIPE_DIMENSONLESS_roughness_corr5.xlsx");

Re_data = data.Re;
f_data  = data.f;
eD_data = data.eD;

nData = height(data);

% Index ranges 
idx_Oregon    = 363:380;
idx_Princeton = 381:406;
all_idx       = 1:nData;
idx_Nikuradse = setdiff(all_idx, [idx_Oregon, idx_Princeton]);

% For info (approx e/D of smooth-pipe datasets)
eD_Oregon    = mean(eD_data(idx_Oregon));
eD_Princeton = mean(eD_data(idx_Princeton));

% Define the Reynolds number range and relative roughness values for plotting
Re_plot = logspace(3, 8, 100).';   % full grid, column
relative_roughness_plot = [ ...
    3.6463e-06, ...
    0.00100603621730382, ...
    0.00202429149797571, ...
    0.00404203718674212, ...
    0.00821692686935086, ...
    0.0164338537387017, ...
    0.0331674958540630 ];

% Symbolic Regression Discovered Relation: SR_f(Re,e/D)
SR_f    = sym('SR_f', [length(relative_roughness_plot) length(Re_plot)]);
SR_Eq_f = sym('SR_Eq_f', size(SR_f));

for i = 1:length(relative_roughness_plot)
    for j = 1:length(Re_plot)
        x1 = Re_plot(j);                  % Re
        x2 = relative_roughness_plot(i);  % e/D

        SR_Eq_f(i, j) = SR_f(i, j) == subs(modelSym, {'x1', 'x2'}, {x1, x2});
        SR_f(i, j)    = vpasolve(SR_Eq_f(i, j), SR_f(i, j));
    end
end
SR_f = double(SR_f); % Convert symbolic results to numeric


% Haaland correlation on the same grid
f_H = @(Re,eD) ( -1.8 .* log10( ((eD)./3.7).^1.11 + 6.9 ./ Re ) ).^(-2);

Haaland_f = zeros(length(relative_roughness_plot), length(Re_plot));
for i = 1:length(relative_roughness_plot)
    eD     = relative_roughness_plot(i);
    Re_vec = Re_plot(:);
    Haaland_f(i,:) = f_H(Re_vec, eD * ones(size(Re_vec))).';
end

% Laminar regime 
Re_laminar = logspace(log10(800), log10(2000), 100); % Laminar range
f_laminar  = 64 ./ Re_laminar;                       % Laminar friction factor

% Moody-style: f(Re) comparison (SR vs Haaland)

% Map SR/Haaland grids to OMA-style naming
f_SR_grid = SR_f;        % size: nED x nRe
f_H_grid  = Haaland_f;

eD_vals = relative_roughness_plot(:).';      % row
eD_plot_vals      = eD_vals;

Re_plot = Re_plot(:);                        % column
nED     = numel(eD_plot_vals);
nRe     = numel(Re_plot);

% For turbulent curves we only draw from Re >= 4e3
mask_turb = Re_plot >= 4e3;
Re_turb   = Re_plot(mask_turb);

% styling helpers
apply_axes_style = @(ax) set(ax, ...
    'FontSize',28,'LineWidth',1.1,'TickDir','out','Layer','top', ...
    'XMinorGrid','on','YMinorGrid','on');

set_grid_colors = @(ax) set(ax, ...
    'GridAlpha',0.25,'MinorGridAlpha',0.40, ...
    'GridColor',0.6*[1 1 1], 'MinorGridColor',0.85*[1 1 1], ...
    'MinorGridLineStyle',':');

% ---- figure & axes ----
figM = figure('Name','f_vs_Re_SR_vs_Haaland', ...
              'NumberTitle','off', ...
              'Units','pixels', ...
              'Position',[100 100 1100 700]);
axM  = axes('Parent',figM);
hold(axM,'on');

% 1) Nikuradse + smooth-pipe datasets 

% Nikuradse data
hNik = loglog(axM, Re_data(idx_Nikuradse), f_data(idx_Nikuradse), 'o', ...
              'MarkerSize',7, ...
              'MarkerEdgeColor',[0.15 0.15 0.15], ...
              'MarkerFaceColor',[0.85 0.85 0.85], ...
              'LineWidth',0.8, ...
              'LineStyle','none', ...
              'DisplayName','Nikuradse (1950)');

% Swanson et al. (Oregon smooth pipe)
hOre = loglog(axM, Re_data(idx_Oregon), f_data(idx_Oregon), 's', ...
              'MarkerSize',8, ...
              'MarkerEdgeColor','k', ...
              'MarkerFaceColor',[0.25 0.45 0.9], ...
              'LineStyle','none', ...
              'DisplayName',sprintf('Swanson et. al. (2002)'));

% Zagarola & Smits (Princeton smooth pipe)
hPri = loglog(axM, Re_data(idx_Princeton), f_data(idx_Princeton), 'd', ...
              'MarkerSize',8, ...
              'MarkerEdgeColor','k', ...
              'MarkerFaceColor',[0.9 0.45 0.1], ...
              'LineStyle','none', ...
              'DisplayName',sprintf('Zagarola \\& Smits (1998)'));

% 2) Haaland (solid) and SR (dashed) curves, one per Nikuradse e/D

idx_lab_vec = round(linspace(0.60*nRe, 0.88*nRe, nED));  % stagger label x-positions

for i = 1:nED
    % Haaland (solid black, turbulent only)
    f_H_i = f_H_grid(i,:).';
    mask_i = mask_turb & isfinite(f_H_i) & (f_H_i > 0);
    if any(mask_i)
        loglog(axM, Re_plot(mask_i), f_H_i(mask_i), 'k-', ...
               'LineWidth',1.4, 'HandleVisibility','off');
    end

    % SR model (dashed black, turbulent only)
    f_SR_i = f_SR_grid(i,:).';
    mask_s = mask_turb & isfinite(f_SR_i) & (f_SR_i > 0);
    if any(mask_s)
        loglog(axM, Re_plot(mask_s), f_SR_i(mask_s), 'k--', ...
               'LineWidth',1.4, 'HandleVisibility','off');
    end

    % e/D LABELS
    idx_lab = idx_lab_vec(i);
    idx_lab = max(1, min(nRe, idx_lab));       % safety
    Re_lab  = Re_plot(idx_lab);
    f_lab   = f_SR_grid(i, idx_lab) * 1.01;  

    
    if eD_plot_vals(i) < 1e-4
        s_eD = char(compose('%.3g', eD_plot_vals(i)));
    else
        s_eD = char(compose('%.3f', eD_plot_vals(i)));
    end
    labelStr = sprintf('$\\mathbf{e/D = %s}$', s_eD);

    text(axM, Re_lab, f_lab, labelStr, ...
        'Interpreter','latex', ...
        'FontSize',18, ...          
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'Color','k', ...
        'Clipping','on');

end

% 4) Dummy entries for theory lines (legend key)

hH  = plot(axM, NaN,NaN,'k-','LineWidth',1.8, ...
           'DisplayName','Haaland (solid)');
hSR = plot(axM, NaN,NaN,'k--','LineWidth',1.8, ...
           'DisplayName','SR model (dashed)');

% 5) Axes, labels, legend

set(axM,'XScale','log','YScale','log');
apply_axes_style(axM);
set_grid_colors(axM);
grid(axM,'on'); box(axM,'on');

xlabel(axM,'Re', ...
       'Interpreter','latex','FontSize',28);
ylabel(axM,'f', ...
       'Interpreter','latex','FontSize',28);


% Legend: datasets + theory
legM = legend(axM, [hNik hOre hPri hH hSR], ...
    'Location','southwest', ...
    'Interpreter','latex', ...
    'FontSize',20, ...
    'Box','on', ...
    'NumColumns',1);
legM.ItemTokenSize = [28 9];

ti  = axM.TightInset;
pad = 0.02;
axM.Position = [ti(1)+pad, ti(2)+pad, ...
                1-(ti(1)+ti(3))-2*pad, ...
                1-(ti(2)+ti(4))-2*pad];

set(figM,'Color','w');
hold(axM,'off');


%% 3D PARETO
chosenIdx = [714 541 171];
[paretoIdx,domMat] = plot_pareto3(gp, 0.3, chosenIdx);

%% PARETO TRADEOFF PLOT
chosenIdx = [714 541 171];
[paretoIdx, domMat, validIdx] = plot_obj_vs_index(gp, 0.3, chosenIdx);

%% CONSTRAINT CHECK PLOTS
modelStr = '0.2544*x2 + 0.1663*x2^((1/(x1*x2^0.9736) + 0.9937)^118.2) + 0.03273*(0.1026*x2 - (1.0*(32840.0*x2 - 2159.0))/x1)^0.4667 - 0.02735/(x2 + (1/(x1*x2^1.059) + 84.17)/x1)^0.02939 + 0.04604';
f_handle = modelSym_to_fhandle(modelStr);
check_C1C2C3(f_handle);   

%%
%  Validation on new rough Superpipe data (Shockling & Langelandsvik)
%  File: NIKURADSE_SUPERPIPE_DIMENSONLESS_roughness_corr6.xlsx
%  Data sets:
%    363:380   Swanson et al. (Oregon smooth)
%    381:406   Zagarola & Smits (Princeton smooth)
%    407:432   Shockling (2006) rough pipe
%    433:449   Langelandsvik (2008) rough pipe
%
%  Here we compare SR vs Haaland on the
%  four corresponding e/D levels.

modelSym = gpmodel2sym(gp, 'best');

% --- Load new validation data (rough Superpipe sets) ---
data_val = readtable("NIKURADSE_SUPERPIPE_DIMENSONLESS_roughness_corr6.xlsx");

Re_val = data_val.Re;
f_val  = data_val.f;
eD_val = data_val.eD;

% Index ranges (given)
idx_Oregon      = 363:380;
idx_Princeton   = 381:406;
idx_Shockling   = 407:432;
idx_Langeland   = 433:449;

% Representative roughness for each experimental set
eD_Oregon    = mean(eD_val(idx_Oregon));
eD_Princeton = mean(eD_val(idx_Princeton));
eD_Shockling = mean(eD_val(idx_Shockling));
eD_Langeland = mean(eD_val(idx_Langeland));

% 4 roughness levels for theory curves
relative_roughness_plot = [ ...
    eD_Oregon, ...
    eD_Princeton, ...
    eD_Shockling, ...
    eD_Langeland ];

%  SR model on validation grid: SR_f(Re,e/D)
Re_plot = logspace(4, 9, 120).';   

SR_f    = sym('SR_f', [numel(relative_roughness_plot), numel(Re_plot)]);
SR_Eq_f = sym('SR_Eq_f', size(SR_f));

for i = 1:numel(relative_roughness_plot)
    for j = 1:numel(Re_plot)
        x1 = Re_plot(j);                       % Re
        x2 = relative_roughness_plot(i);       % e/D
        SR_Eq_f(i,j) = SR_f(i,j) == subs(modelSym, {'x1','x2'}, {x1, x2});
        SR_f(i,j)    = vpasolve(SR_Eq_f(i,j), SR_f(i,j));
    end
end
SR_f = double(SR_f);   % numeric grid for SR model

%  Haaland correlation on the same validation grid
f_H = @(Re,eD) ( -1.8 .* log10( ((eD)./3.7).^1.11 + 6.9 ./ Re ) ).^(-2);

Haaland_f = zeros(numel(relative_roughness_plot), numel(Re_plot));
for i = 1:numel(relative_roughness_plot)
    eD     = relative_roughness_plot(i);
    Re_vec = Re_plot(:);
    Haaland_f(i,:) = f_H(Re_vec, eD * ones(size(Re_vec))).';
end

%  Moody-style validation: only Superpipe data + 4 roughness curves

apply_axes_style = @(ax) set(ax, ...
    'FontSize',28,'LineWidth',1.1,'TickDir','out','Layer','top', ...
    'XMinorGrid','on','YMinorGrid','on');

set_grid_colors = @(ax) set(ax, ...
    'GridAlpha',0.25,'MinorGridAlpha',0.40, ...
    'GridColor',0.6*[1 1 1], 'MinorGridColor',0.85*[1 1 1], ...
    'MinorGridLineStyle',':');

fsz = 28;

figV = figure('Name','Validation_SR_vs_Haaland_Superpipe', ...
              'NumberTitle','off', ...
              'Units','pixels', ...
              'Position',[150 150 1100 700]);
axV  = axes('Parent',figV); hold(axV,'on');

% Colors (match data & curves) 
colOre = [0.25 0.45 0.90];   % blue
colPri = [0.90 0.45 0.10];   % orange
colSho = [0.20 0.70 0.20];   % green
colLan = [0.1 0.10 0.10];   % reddish 

curveColors = [colOre; colPri; colSho; colLan];

% 1) Experimental data markers 

% Swanson et al. (Oregon smooth)
hOre = loglog(axV, Re_val(idx_Oregon), f_val(idx_Oregon), 's', ...
    'MarkerSize',8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor',colOre, ...
    'LineStyle','none', ...
    'DisplayName',sprintf('Swanson et al. (2002), $\\varepsilon/D\\approx %.2g$', eD_Oregon));

% Zagarola & Smits (Princeton smooth)
hPri = loglog(axV, Re_val(idx_Princeton), f_val(idx_Princeton), 'd', ...
    'MarkerSize',8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor',colPri, ...
    'LineStyle','none', ...
    'DisplayName',sprintf('Zagarola \\& Smits (1998), $\\varepsilon/D\\approx %.2g$', eD_Princeton));

% Shockling (2006) rough pipe
hSho = loglog(axV, Re_val(idx_Shockling), f_val(idx_Shockling), '^', ...
    'MarkerSize',8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor',colSho, ...
    'LineStyle','none', ...
    'DisplayName',sprintf('Shockling et al. (2006), $\\varepsilon/D\\approx %.2g$', eD_Shockling));

% Langelandsvik (2008) rough pipe
hLan = loglog(axV, Re_val(idx_Langeland), f_val(idx_Langeland), 'o', ...
    'MarkerSize',8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor',colLan, ...
    'LineStyle','none', ...
    'DisplayName',sprintf('Langelandsvik et al. (2008), $\\varepsilon/D\\approx %.2g$', eD_Langeland));

% 2) Theory curves: Haaland (solid) and SR (dashed), colored by dataset 

Re_plot = Re_plot(:);
nRe     = numel(Re_plot);
nED     = numel(relative_roughness_plot);

mask_turb = Re_plot >= 4e3;

for i = 1:nED
    eD_i   = relative_roughness_plot(i);
    col_i  = curveColors(i,:);

    % Haaland (solid)
    f_H_i = Haaland_f(i,:).';
    mask_i = mask_turb & isfinite(f_H_i) & (f_H_i > 0);
    if any(mask_i)
        loglog(axV, Re_plot(mask_i), f_H_i(mask_i), '-', ...
               'LineWidth',1.4, ...
               'Color',col_i, ...
               'HandleVisibility','off');
    end

    % SR (dashed)
    f_SR_i = SR_f(i,:).';
    mask_s = mask_turb & isfinite(f_SR_i) & (f_SR_i > 0);
    if any(mask_s)
        loglog(axV, Re_plot(mask_s), f_SR_i(mask_s), '--', ...
               'LineWidth',1.4, ...
               'Color',col_i, ...
               'HandleVisibility','off');
    end
end

% 2b) e/D labels on the right, in normalized coordinates 

[~, orderED] = sort(relative_roughness_plot, 'ascend');
y_norm = linspace(0.25, 0.70, nED);   

for k = 1:nED
    i     = orderED(k);
    eD_i  = relative_roughness_plot(i);
    col_i = curveColors(i,:);

    if eD_i < 1e-4
        s_eD = char(compose('%.2g', eD_i));
    else
        s_eD = char(compose('%.3g', eD_i));
    end
    labelStr = sprintf('$\\varepsilon/D = %s$', s_eD);

    text(axV, 0.985, y_norm(k), labelStr, ...
        'Units','normalized', ...           
        'Interpreter','latex', ...
        'FontSize',16, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'Color',col_i);
end

% 3) Dummy entries for theory in legend 
hH  = plot(axV, NaN,NaN,'k-','LineWidth',1.8, ...
           'DisplayName','Haaland (solid)');
hSR = plot(axV, NaN,NaN,'k--','LineWidth',1.8, ...
           'DisplayName','SR model (dashed)');

% 4) Axes / legend styling
set(axV,'XScale','log','YScale','log');
apply_axes_style(axV);
set_grid_colors(axV);
grid(axV,'on'); box(axV,'on');

xlabel(axV,'$\mathrm{Re} = \rho U D / \mu$', ...
       'Interpreter','latex','FontSize',fsz);
ylabel(axV,'$f = \Delta P D /(0.5\,\rho U^2 L)$', ...
       'Interpreter','latex','FontSize',fsz);

xlim(axV,[1e4 1e7]);          
ylim(axV,[2e-3 0.08]);

legV = legend(axV, [hOre hPri hSho hLan hH hSR], ...
    'Location','southwest', ...
    'Interpreter','latex', ...
    'FontSize',22, ...
    'Box','on', ...
    'NumColumns',1);
legV.ItemTokenSize = [28 9];

ti  = axV.TightInset;
pad = 0.02;
axV.Position = [ti(1)+pad, ti(2)+pad, ...
                1-(ti(1)+ti(3))-2*pad, ...
                1-(ti(2)+ti(4))-2*pad];

set(figV,'Color','w');
hold(axV,'off');

%% ======================== LOCAL FUNCTIONS ============================

function phys = get_phys_params_safe(gp)
% Same defaults/logic as in your regressmulti_fitfun
    phys.D      = 0.05;
    phys.L      = 1.0;
    phys.mu_ref = 1e-3;
    phys.rho    = 1000;

    if isfield(gp,'userdata') && isstruct(gp.userdata)
        if isfield(gp.userdata,'D')      && ~isempty(gp.userdata.D),      phys.D      = gp.userdata.D; end
        if isfield(gp.userdata,'L')      && ~isempty(gp.userdata.L),      phys.L      = gp.userdata.L; end
        if isfield(gp.userdata,'mu_ref') && ~isempty(gp.userdata.mu_ref), phys.mu_ref = gp.userdata.mu_ref; end
        if isfield(gp.userdata,'rho')    && ~isempty(gp.userdata.rho),    phys.rho    = gp.userdata.rho; end
    end
    if isfield(gp,'user') && isstruct(gp.user)
        if isfield(gp.user,'D')      && ~isempty(gp.user.D),      phys.D      = gp.user.D; end
        if isfield(gp.user,'L')      && ~isempty(gp.user.L),      phys.L      = gp.user.L; end
        if isfield(gp.user,'mu_ref') && ~isempty(gp.user.mu_ref), phys.mu_ref = gp.user.mu_ref; end
        if isfield(gp.user,'rho')    && ~isempty(gp.user.rho),    phys.rho    = gp.user.rho; end
    end
end

function pen = compute_C1_penalty_param_smallgrid(f_handle, phys, eD_list, nRe)
% Same as compute_C1_penalty_param, but lets you shrink Re_grid for speed.
% (Everything else kept identical to your function structure.)
    rho    = phys.rho;
    mu_ref = phys.mu_ref;
    D      = phys.D;
    L      = phys.L;

    Re_grid = logspace(4, 9, nRe).';
    V_grid  = (Re_grid .* mu_ref) ./ (rho * D);

    ne = numel(eD_list);
    per_eD = nan(ne,1);

    movmean3 = @(v) conv(v,[1 1 1]/3,'same');

    ReA_min = 1e4; ReA_max = 1e6;
    ReB_min = 1e6; ReB_max = 1e9;

    nA_min  = 1.0;
    nA_max  = 2.4;
    nB_ref  = 2.0;
    tolB    = 0.001;

    for ii = 1:ne
        eD = eD_list(ii);

        f = f_handle(Re_grid, eD*ones(size(Re_grid)));

        ok = isfinite(f) & isreal(f) & (f > 0);
        if nnz(ok) < 30
            per_eD(ii) = 1;
            continue;
        end

        Re_ok = Re_grid(ok);
        V_ok  = V_grid(ok);
        f_ok  = f(ok);

        DP_ok = f_ok .* (rho .* (V_ok.^2) .* L) ./ (2 .* D);
        ok2 = isfinite(DP_ok) & isreal(DP_ok) & (DP_ok > 0);
        if nnz(ok2) < 30
            per_eD(ii) = 1;
            continue;
        end

        Re_ok = Re_ok(ok2);
        V_ok  = V_ok(ok2);
        DP_ok = DP_ok(ok2);

        if numel(DP_ok) >= 5
            DP_sm = movmean3(DP_ok);
            V_sm  = movmean3(V_ok);
        else
            DP_sm = DP_ok;
            V_sm  = V_ok;
        end

        lnDP = log(DP_sm);
        lnV  = log(V_sm);
        nloc = gradient(lnDP, lnV);

        maskA = (Re_ok >= ReA_min) & (Re_ok <= ReA_max);
        if nnz(maskA) < 6
            penA = 1;
        else
            nA = nloc(maskA);
            dev_low  = max(0, nA_min - nA);
            dev_high = max(0, nA - nA_max);
            devA = mean(dev_low + dev_high, 'omitnan') / (nA_max - nA_min + eps);
            penA = min(max(real(devA),0),1);
        end

        maskB = (Re_ok >= ReB_min) & (Re_ok <= ReB_max);
        if nnz(maskB) < 6
            penB = 1;
        else
            nB = nloc(maskB);
            devB_raw = mean( max(0, abs(nB - nB_ref) - tolB), 'omitnan' );
            devB = devB_raw / (0.5 + eps);
            penB = min(max(real(devB),0),1);
        end

        wA = 0.3; wB = 0.7;
        per_eD(ii) = min(1, wA*penA + wB*penB);
    end

    C1_core = max(per_eD, [], 'omitnan');
    if ~isfinite(C1_core), C1_core = 1; end
    pen = min(max(real(C1_core),0),1);
end


%%
function [paretoIdx, domMat, validIdx] = plot_obj_vs_index(gp, fitnessMax, chosenIdx)
%PLOT_OBJ_VS_INDEX  Objectives vs population index (3 stacked subplots).
%
%   [paretoIdx, domMat, validIdx] = plot_obj_vs_index(gp, fitnessMax, chosenIdx)
%
% Expects gp.fitness.multiobj as [N x 3]:
%   col1 = fitness (J_err), col2 = complexity (J_comp), col3 = constraint score (J_phys)

    if nargin < 2 || isempty(fitnessMax), fitnessMax = 1; end
    if nargin < 3, chosenIdx = []; end
    chosenIdx = chosenIdx(:);

    if ~isfield(gp,'fitness') || ~isfield(gp.fitness,'multiobj') || isempty(gp.fitness.multiobj)
        error('gp.fitness.multiobj not found or empty.');
    end

    M = gp.fitness.multiobj;
    if size(M,2) ~= 3
        error('gp.fitness.multiobj must be [N x 3] = [fitness, complexity, constraint].');
    end

    % filtering
    good = all(isfinite(M),2) & (M(:,1) <= fitnessMax);
    validIdx = find(good);
    A = M(validIdx,:);    % [Jerr, Jcomp, Jphys] but stored as [fitness, complexity, constraint]

    if isempty(validIdx)
        error('No valid individuals after filtering (fitnessMax=%.3g).', fitnessMax);
    end

    % Pareto front among valid points (minimize all 3) 
    nv = size(A,1);
    domMat = false(nv,nv);
    for i = 1:nv
        Ai = A(i,:);
        for j = 1:nv
            Aj = A(j,:);
            if all(Aj <= Ai) && any(Aj < Ai)   % j dominates i (minimization)
                domMat(i,j) = true;
            end
        end
    end
    isDominated  = any(domMat,2);
    paretoRelIdx = find(~isDominated);
    paretoIdx    = validIdx(paretoRelIdx);

    % map chosen indices into valid subset 
    chosenRelIdx = nan(size(chosenIdx));
    for k = 1:numel(chosenIdx)
        t = find(validIdx == chosenIdx(k), 1, 'first');
        if ~isempty(t), chosenRelIdx(k) = t; end
    end
    keepChosen = isfinite(chosenRelIdx);
    chosenRelIdx = chosenRelIdx(keepChosen);
    chosenIdx    = chosenIdx(keepChosen);

    % plotting
    fsz = 18;
    LW  = 1.2;

    fig = figure('Color','w', ...
        'Name',sprintf('Objectives vs Index (fitness <= %.3g)',fitnessMax), ...
        'NumberTitle','off');
    set(fig,'Renderer','painters');
    set(fig,'WindowStyle','normal');


    set(fig,'Units','pixels');
    try
        fig.WindowState = 'normal';
    catch
        
    end
    fig.Position = [80 80 1500 980];   

    xAll = validIdx;          % x-axis is population index
    Jerr = A(:,1);
    Jcmp = A(:,2);
    Jphy = A(:,3);

    % styles
    popMS   = 14;   popAlpha = 0.6;     
    pfMS    = 45;   pfAlpha  = 0.95;    
    chMS    = 220;  chLW     = 2.0;     

    popColor = [0.10 0.20 0.90];
    pfFace   = [0.85 0.15 0.15];
    pfEdge   = [0 0 0];

    colorList  = [ 0 0 0;
        0.0 0.7 0.0;
        0.8 0.0 0.8;
        0.0 0.75 0.75;
        0.9 0.4 0.0 ];
    markerList = {'p','s','^','v','d','h'};

    
    tlo = tiledlayout(fig, 3, 1, 'TileSpacing','compact', 'Padding','compact');

    ax1 = nexttile(tlo,1); hold(ax1,'on'); grid(ax1,'on');
    ax2 = nexttile(tlo,2); hold(ax2,'on'); grid(ax2,'on');
    ax3 = nexttile(tlo,3); hold(ax3,'on'); grid(ax3,'on');

    set([ax1 ax2 ax3],'FontSize',fsz,'LineWidth',LW,'TickDir','out');
    linkaxes([ax1 ax2 ax3],'x');

    % vertical guides for chosen indices 
    for k = 1:numel(chosenIdx)
        xline(ax1, chosenIdx(k), '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1.0);
        xline(ax2, chosenIdx(k), '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1.0);
        xline(ax3, chosenIdx(k), '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1.0);
    end

    % 1) J_err 
    hPop1 = scatter(ax1, xAll, Jerr, popMS, 'o', ...
        'MarkerFaceColor', popColor, 'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha', popAlpha);

    hPF1  = scatter(ax1, paretoIdx, M(paretoIdx,1), pfMS, 'o', ...
        'MarkerFaceColor', pfFace, 'MarkerEdgeColor', pfEdge, ...
        'LineWidth',0.8, 'MarkerFaceAlpha', pfAlpha);

    hChosen1 = gobjects(0);
    for k = 1:numel(chosenIdx)
        c = colorList(1 + mod(k-1,size(colorList,1)), :);
        m = markerList{1 + mod(k-1,numel(markerList))};
        xx = chosenIdx(k);
        yy = M(xx,1);
        hChosen1(k) = scatter(ax1, xx, yy, chMS, m, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', [0 0 0], 'LineWidth', chLW);
    end

    ylabel(ax1,'J_{err}','FontSize',fsz);
    xlim(ax1,[min(xAll) max(xAll)]);

    % ---- 2) J_phys ----
    scatter(ax2, xAll, Jphy, popMS, 'o', ...
        'MarkerFaceColor', popColor, 'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha', popAlpha);

    scatter(ax2, paretoIdx, M(paretoIdx,3), pfMS, 'o', ...
        'MarkerFaceColor', pfFace, 'MarkerEdgeColor', pfEdge, ...
        'LineWidth',0.8, 'MarkerFaceAlpha', pfAlpha);

    for k = 1:numel(chosenIdx)
        c = colorList(1 + mod(k-1,size(colorList,1)), :);
        m = markerList{1 + mod(k-1,numel(markerList))};
        xx = chosenIdx(k);
        yy = M(xx,3);
        scatter(ax2, xx, yy, chMS, m, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', [0 0 0], 'LineWidth', chLW);
    end

    ylabel(ax2,'J_{phys}','FontSize',fsz);
    xlim(ax2,[min(xAll) max(xAll)]);
    ylim(ax2,[0 0.5]);   

    % 3) J_comp
    scatter(ax3, xAll, Jcmp, popMS, 'o', ...
        'MarkerFaceColor', popColor, 'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha', popAlpha);

    scatter(ax3, paretoIdx, M(paretoIdx,2), pfMS, 'o', ...
        'MarkerFaceColor', pfFace, 'MarkerEdgeColor', pfEdge, ...
        'LineWidth',0.8, 'MarkerFaceAlpha', pfAlpha);

    for k = 1:numel(chosenIdx)
        c = colorList(1 + mod(k-1,size(colorList,1)), :);
        m = markerList{1 + mod(k-1,numel(markerList))};
        xx = chosenIdx(k);
        yy = M(xx,2);
        scatter(ax3, xx, yy, chMS, m, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', [0 0 0], 'LineWidth', chLW);
    end

    ylabel(ax3,'J_{comp}','FontSize',fsz);
    xlabel(ax3,'Population index','FontSize',fsz);
    xlim(ax3,[min(xAll) max(xAll)]);

    % BIG LEGEND (dummy handles so markers are readable) 
    legFS     = fsz + 2;
    legMS_pop = 10;
    legMS_pf  = 12;
    legMS_ch  = 14;

    hLegPop = plot(ax1, NaN, NaN, 'o', ...
        'MarkerFaceColor', popColor, 'MarkerEdgeColor','none', ...
        'MarkerSize', legMS_pop);

    hLegPF  = plot(ax1, NaN, NaN, 'o', ...
        'MarkerFaceColor', pfFace, 'MarkerEdgeColor', pfEdge, ...
        'LineWidth', 0.8, 'MarkerSize', legMS_pf);

    hLegChosen = gobjects(numel(chosenIdx),1);
    for k = 1:numel(chosenIdx)
        c = colorList(1 + mod(k-1,size(colorList,1)), :);
        m = markerList{1 + mod(k-1,numel(markerList))};
        hLegChosen(k) = plot(ax1, NaN, NaN, m, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', [0 0 0], ...
            'LineWidth', chLW, 'MarkerSize', legMS_ch);
    end

    legLabels = [{'Population'}, {'Pareto front'}, ...
        arrayfun(@(v)sprintf('index %d',v), chosenIdx','UniformOutput',false)];

    
    lg = legend(ax1, [hLegPop hLegPF hLegChosen.'], legLabels, ...
        'Location','northeast', 'Box','on');
    lg.FontSize = legFS;

    try
        lg.ItemTokenSize = [26, 18];
    catch
    end
end
%%
% -------------------------- Candidate models --------------------------

function f = f_SR_cand4(Re, eD, p)

x1 = Re; x2 = eD;

c1  = p(1);
c2  = p(2);
c3  = p(3);
c4  = p(4);
c5  = p(5);
c6  = p(6);
c7  = p(7);
c8  = p(8);
c9  = p(9);
c10 = p(10);
c11 = p(11);
c12 = p(12);
c13 = p(13);
c14 = p(14);
c15 = p(15);
c16 = p(16);

term1 = c1 .* x2;
term2 = -c2 .* tanh(c3./x1 + c4);

base3 = x2.^2 - ( (c6.*x2 - c7./x1) ./ x1 );
term3 = c5 .* (base3).^c8;

term4 = -c9  ./ ( (c10./x2).^(c11./x1) );
term5 =  c12 ./ ( (c13./x2).^(c14./x1) );

term6 = -c15 ./ x1;
term7 = c16;

f = term1 + term2 + term3 + term4 + term5 + term6 + term7;
end

%%

