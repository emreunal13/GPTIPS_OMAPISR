function OMA_VS_HAALAND_MOODY()
% OMA_VS_HAALAND_MOODY
%   Moody-style f(Re) comparison:
%   - Datasets: Nikuradse + Oregon + Princeton
%   - Theory:   Haaland (solid) vs OMA model (dashed)
%
%   Requires:
%       f_OMA_model.m
%       NIKURADSE_SUPERPIPE_DIMENSONLESS_roughness_corr5.xlsx
%
%   Call from command window with:
%       OMA_VS_HAALAND_MOODY

close all;

%% Load combined data (Nikuradse + smooth-pipe lab sets)

data = readtable("NIKURADSE_SUPERPIPE_DIMENSONLESS_roughness_corr5.xlsx");

Re_data = data.Re;
f_data  = data.f;
eD_data = data.eD;

nData = height(data);

% Index ranges (same convention as before)
idx_Oregon    = 363:380;
idx_Princeton = 381:406;
all_idx       = 1:nData;
idx_Nikuradse = setdiff(all_idx, [idx_Oregon, idx_Princeton]);

% Approximate e/D of smooth-pipe datasets
eD_Oregon    = mean(eD_data(idx_Oregon));
eD_Princeton = mean(eD_data(idx_Princeton));  %#ok<NASGU>

%% Define theoretical grids: Re & e/D

% Reynolds number grid (same as SR vs Colebrook version)
Re_plot = logspace(3, 8, 260).';    % column

% Nikuradse-style roughnesses
eD_vals = [1.87569231e-06, 0.00100603621730382, 0.00202429149797571, ...
           0.00404203718674212, 0.00821692686935086, ...
           0.0164338537387017, 0.0331674958540630];


eD_plot_vals    = eD_vals;
eD_plot_vals(1) = 3.6463e-06;   

Re_plot = Re_plot(:);           % ensure column
nED     = numel(eD_plot_vals);
nRe     = numel(Re_plot);

% Only draw turbulent part from Re >= 4e3
mask_turb = Re_plot >= 4e3;
Re_turb   = Re_plot(mask_turb);

%% Define friction correlations: OMA & Haaland

f_OMA = @f_OMA_model;   

f_H  = @(Re,eD) ( -1.8 .* log10( ((eD)./3.7).^1.11 + 6.9 ./ Re ) ).^(-2);

% Allocate theory grids
f_OMA_grid = nan(nED, nRe);
f_H_grid   = nan(nED, nRe);

for i = 1:nED
    eD = eD_vals(i);               % for theory we keep original eD_vals
    Re_vec = Re_plot(:);           % column
    eD_vec = eD * ones(size(Re_vec));

    % OMA (explicit)
    f_OMA_i = f_OMA(Re_vec, eD_vec);    % column
    f_OMA_grid(i,:) = f_OMA_i(:).';

    % Haaland (explicit)
    f_H_i = f_H(Re_vec, eD_vec);
    f_H_grid(i,:) = f_H_i(:).';
end

%% Styling helpers

apply_axes_style = @(ax) set(ax, ...
    'FontSize',28,'LineWidth',1.1,'TickDir','out','Layer','top', ...
    'XMinorGrid','on','YMinorGrid','on');

set_grid_colors = @(ax) set(ax, ...
    'GridAlpha',0.25,'MinorGridAlpha',0.40, ...
    'GridColor',0.6*[1 1 1], 'MinorGridColor',0.85*[1 1 1], ...
    'MinorGridLineStyle',':');

%% Figure & axis

figM = figure('Name','f_vs_Re_OMA_vs_Haala nd', ...
              'NumberTitle','off', ...
              'Units','pixels', ...
              'Position',[100 100 1100 700]);
axM  = axes('Parent',figM);
hold(axM,'on');

%% Plot experimental datasets (same as SR plot)

% Nikuradse data
hNik = loglog(axM, Re_data(idx_Nikuradse), f_data(idx_Nikuradse), 'o', ...
              'MarkerSize',7, ...
              'MarkerEdgeColor',[0.15 0.15 0.15], ...
              'MarkerFaceColor',[0.85 0.85 0.85], ...
              'LineWidth',0.8, ...
              'LineStyle','none', ...
              'DisplayName','Nikuradse data');

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

%% Plot Haaland (solid) and OMA (dashed) curves per e/D

% Staggered label x-positions
idx_lab_vec = round(linspace(0.60*nRe, 0.88*nRe, nED));

for i = 1:nED
    % Haaland (solid black, turbulent only) 
    f_H_i = f_H_grid(i,:).';
    mask_i = mask_turb & isfinite(f_H_i) & (f_H_i > 0);
    if any(mask_i)
        loglog(axM, Re_plot(mask_i), f_H_i(mask_i), 'k-', ...
               'LineWidth',1.4, 'HandleVisibility','off');
    end

    % OMA (dashed black, turbulent only)
    f_OMA_i = f_OMA_grid(i,:).';
    mask_s = mask_turb & isfinite(f_OMA_i) & (f_OMA_i > 0);
    if any(mask_s)
        loglog(axM, Re_plot(mask_s), f_OMA_i(mask_s), 'k--', ...
               'LineWidth',1.4, 'HandleVisibility','off');
    end

    % BIGGER, BOLD e/D LABELS 
    idx_lab = idx_lab_vec(i);
    idx_lab = max(1, min(nRe, idx_lab));       % safety
    Re_lab  = Re_plot(idx_lab);
    f_lab   = f_OMA_grid(i, idx_lab) * 1.01;   % a bit above dashed line

    s_eD     = short_eD_lbl(eD_plot_vals(i));
    labelStr = sprintf('{\\boldmath$\\varepsilon/D = %s$}', s_eD);

    text(axM, Re_lab, f_lab, labelStr, ...
        'Interpreter','latex', ...
        'FontSize',18, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'Color','k', ...
        'Clipping','on');
end

%% Dummy entries for theory lines in legend

hHaaland = plot(axM, NaN,NaN,'k-','LineWidth',1.8, ...
                'DisplayName','Haaland (solid)');
hOMA     = plot(axM, NaN,NaN,'k--','LineWidth',1.8, ...
                'DisplayName','OMA Model(dashed)');

%% Axes, labels, legend

set(axM,'XScale','log','YScale','log');
apply_axes_style(axM);
set_grid_colors(axM);
grid(axM,'on'); box(axM,'on');

xlabel(axM,'$\mathrm{Re}$','Interpreter','latex','FontSize',28);
ylabel(axM,'$f$','Interpreter','latex','FontSize',28);

xlim(axM,[1e3 1e8]);
ylim(axM,[3e-3 0.1]);

% Legend: datasets + theory
legM = legend(axM, [hNik hOre hPri hHaaland hOMA], ...
    'Location','southwest', ...
    'Interpreter','latex', ...
    'FontSize',24, ...
    'Box','on', ...
    'NumColumns',1);
legM.ItemTokenSize = [28 9];

% Tight layout
ti  = axM.TightInset;
pad = 0.02;
axM.Position = [ti(1)+pad, ti(2)+pad, ...
                1-(ti(1)+ti(3))-2*pad, ...
                1-(ti(2)+ti(4))-2*pad];

set(figM,'Color','w');
hold(axM,'off');

end  

% local helper: short label for e/D 
function s = short_eD_lbl(x)
if x < 1e-4
    s = char(compose('%.3g', x));
else
    s = char(compose('%.3f', x));
end
end
