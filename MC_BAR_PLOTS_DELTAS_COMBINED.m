function MC_BAR_PLOTS_DELTAS_COMBINED()
%% BAR CHARTS : DELTA J_err + DELTA J_phys
% Combines the two plots into a single figure (two panels) with ONE legend.
% X-axes: delta_J_err and delta_J_phys (both computed relative to nominal p0).
%
% Run:
%   MC_BAR_PLOTS_DELTAS_COMBINED
%
%
% Notes:
%   - Here we define pointwise J_err as dimensionless signed relative error:
%       J_err_point = (f_data - f_model)/f_data
%     then delta_J_err = J_err_point - J_err_point_nominal.
%   - We define pointwise J_phys as your local C1 signed deviation:
%       J_phys_point = C1dev_signed
%     then delta_J_phys = C1dev - C1dev_nominal.
%
% If you prefer J_err to be absolute (e.g. |...|) or squared, tell me and I’ll swap it.

% -------------------- Settings --------------------
nSamples  = 100000;   % Monte Carlo samples per model
relSigma  = 0.01;     % ±1% coefficient noise (uniform)
rng(1);

% C1 zone parameters (match your paper narrative)
Re_switch = 1e6;
nA_min    = 1.0;
nA_max    = 2.4;
nB_ref    = 2.0;

% finite-difference step for local n(V)
hV = 0.02; % 2% velocity perturbation

% -------------------- Load data & pick 3 points from CSV --------------------
csvFile = 'NIKURADSE_DATA_NONZERO_SUPERPIPE_INCH_CONVERTED.csv';
T = readtable(csvFile);

targets = [1e6 1.9e-6; 1e6 0.008; 1e4 0.008];
idx = zeros(3,1);
for i = 1:3
    d = (log10(T.Re) - log10(targets(i,1))).^2 + ...
        (log10(T.eD) - log10(targets(i,2))).^2;
    [~, idx(i)] = min(d);
end

Re_pts = T.Re(idx);
eD_pts = T.eD(idx);
f_data = T.f(idx);

Re_pts = [1e6; 1e6; 1e4];
eD_pts = [1.9e-6; 0.008; 0.008];

nPts = numel(Re_pts);
ptLabel = cell(nPts,1);
for i = 1:nPts
    ptLabel{i} = sprintf('Point %d: Re=%.1e, $\\epsilon/D=%.3g$', i, Re_pts(i), eD_pts(i));
end

% Phys params (used only for local n(V) definition) 
phys.rho    = 1000;    % kg/m^3
phys.mu_ref = 1.0e-3;  % Pa*s
phys.D      = 0.050;   % m
phys.L      = 1.0;     % m

% Nominal parameter vectors 
p_nom1 = [0.2544, 0.1663, 0.9736, 0.9937, 118.2, ...
          0.03273, 0.1026, 32840.0, 2159.0, 0.4667, ...
          0.02735, 1.059, 84.17, 0.02939, 0.04604];

p_nom2 = [0.2332, 0.1197, 0.8878, 0.9958, 222.8, ...
          0.06218, 0.2163, 435.8, 0.2801, 0.004616];

p_nom3 = [0.2547, 0.03494, 1.064, 107.4, 0.02899, ...
          0.1706, 0.9766, 0.9937, 114.3, ...
          0.03829, 0.001307, 20940.0, 1450.0, 0.481, 0.05518];

p_nom4 = [ ...
    0.319, 0.0345, 32950.0, 1.38, ...
    0.1316, 178.2, 22260.0, 0.1584, ...
    0.0345, 0.3391, 11940.0, ...
    0.03535, 0.1135, 13930.0, ...
    27.38, 0.03416 ];

% Haaland “params”
p_haal = [1.8; 3.7; 1.11; 6.9];

% Model handles (must accept (Re,eD,p)) 
f_eval = cell(5,1);
p_nom  = cell(5,1);

f_eval{1} = @(Re,eD,p) f_SR_cand1(Re,eD,p); p_nom{1} = p_nom1(:);
f_eval{2} = @(Re,eD,p) f_SR_cand2(Re,eD,p); p_nom{2} = p_nom2(:);
f_eval{3} = @(Re,eD,p) f_SR_cand3(Re,eD,p); p_nom{3} = p_nom3(:);
f_eval{4} = @(Re,eD,p) f_SR_cand4(Re,eD,p); p_nom{4} = p_nom4(:);
f_eval{5} = @(Re,eD,p) f_Haaland_new(Re,eD,p); p_nom{5} = p_haal(:);

modelNames = {'1','2','3','4','Haaland'};
nModels = numel(f_eval);

% Nominal (baseline) J_err and J_phys at points 
Jerr_nom  = nan(nModels, nPts);
Jphys_nom = nan(nModels, nPts);

for m = 1:nModels
    p0 = p_nom{m};

    f0 = f_eval{m}(Re_pts, eD_pts, p0);
    if any(~isfinite(f0)) || any(~isreal(f0)) || any(f0<=0)
        warning('Model %s has invalid nominal prediction at points.', modelNames{m});
        continue;
    end

    % pointwise "J_err" definition (dimensionless signed)
    Jerr_nom(m,:) = ((f_data(:) - f0(:)) ./ f_data(:)).';

    % pointwise "J_phys" definition
    for i = 1:nPts
        Jphys_nom(m,i) = local_C1_signed_dev_at_point( ...
            f_eval{m}, p0, Re_pts(i), eD_pts(i), phys, hV, ...
            Re_switch, nA_min, nA_max, nB_ref);
    end
end

% Storage for deltas 
dJerr_pts  = cell(nModels,1);  % [nPts x nSamples]
dJphys_pts = cell(nModels,1);  % [nPts x nSamples]
fracBadPred = zeros(nModels,1);

for m = 1:nModels
    dJerr_pts{m}  = nan(nPts, nSamples);
    dJphys_pts{m} = nan(nPts, nSamples);
end

% Monte Carlo 
for m = 1:nModels
    p0 = p_nom{m};
    badCount = 0;

    for k = 1:nSamples
        noise = relSigma * (-1 + 2*rand(size(p0)));
        p_k = p0 .* (1 + noise);

        f_pred = f_eval{m}(Re_pts, eD_pts, p_k);

        if any(~isfinite(f_pred)) || any(~isreal(f_pred)) || any(f_pred <= 0)
            badCount = badCount + 1;
            continue;
        end

        % J_err at points, then delta vs nominal
        Jerr_k = ((f_data(:) - f_pred(:)) ./ f_data(:)).';
        dJerr_pts{m}(:,k) = (Jerr_k(:) - Jerr_nom(m,:).').';

        % J_phys at points, then delta vs nominal
        for i = 1:nPts
            Jphys_k = local_C1_signed_dev_at_point( ...
                f_eval{m}, p_k, Re_pts(i), eD_pts(i), phys, hV, ...
                Re_switch, nA_min, nA_max, nB_ref);

            dJphys_pts{m}(i,k) = Jphys_k - Jphys_nom(m,i);
        end
    end

    fracBadPred(m) = badCount / nSamples;
end

% Combined plotting
titleStr = sprintf('Central 95\\%% confidence interval for %.1f\\%% coefficient noise', 100*relSigma);

fig = figure('Color','w','Position',[50 50 1600 750]);

% Manual layout 
mL = 0.08;   % left margin
mR = 0.03;   % right margin
mT = 0.12;   % top margin (room for sgtitle)
mB = 0.06;   % bottom margin
gapX = 0.06; % gap between the two panels
gapY = 0.03; % gap between plots and legend strip
legH = 0.12; % legend strip height

axW = (1 - mL - mR - gapX)/2;
axH = 1 - mT - mB - legH - gapY;

ax1 = axes('Parent',fig,'Units','normalized', ...
    'Position',[mL, mB+legH+gapY, axW, axH]);
ax2 = axes('Parent',fig,'Units','normalized', ...
    'Position',[mL+axW+gapX, mB+legH+gapY, axW, axH]);

% Legend strip axis 
axLeg = axes('Parent',fig,'Units','normalized', ...
    'Position',[mL, mB, 1-mL-mR, legH], ...
    'Visible','off');
axis(axLeg,'off');
hold(axLeg,'on');

% Plot both panels 
hLegend = plot_ci_bars_points_oldstyle_ax(ax1, dJerr_pts, modelNames, ptLabel, ...
    '$\Delta J_{\mathrm{err}}$', '', true);

plot_ci_bars_points_oldstyle_ax(ax2, dJphys_pts, modelNames, ptLabel, ...
    '$\Delta J_{\mathrm{phys}}$', '', false);

% Only left subplot shows y labels (cleaner)
ax2.YLabel.String = '';
ax2.YTickLabel = {};

% -------- Panel labels (a), (b) --------
text(ax1, 0.02, 0.98, '(a)', 'Units','normalized', ...
    'FontSize',22, 'FontWeight','bold', 'Interpreter','none', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');

text(ax2, 0.02, 0.98, '(b)', 'Units','normalized', ...
    'FontSize',22, 'FontWeight','bold', 'Interpreter','none', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');

% Dummy legend items drawn in legend-strip axes
C = [ ...
    0.0000 0.4470 0.7410; % blue
    0.8500 0.3250 0.0980; % red/orange
    0.4660 0.6740 0.1880  % green
];

hP = gobjects(3,1);
for i = 1:3
    hP(i) = plot(axLeg, nan,nan,'LineWidth',14,'Color',C(i,:));
end
hM = plot(axLeg, nan,nan,'k-','LineWidth',5);

labels = [ptLabel(:); {'Median (50th pct.)'}];

leg = legend(axLeg, [hP; hM], labels, ...
    'Interpreter','latex','Box','on','Orientation','horizontal');
leg.FontSize   = 20;
leg.LineWidth  = 1.5;
leg.NumColumns = 2;
leg.Location   = 'north';    % inside the legend strip

exportgraphics(fig,'BAR_DELTAS_ERR_PHYS_POINTS.png','Resolution',300);

end

%% bar plot helper (AX version) 
function hLegend = plot_ci_bars_points_oldstyle_ax(ax, Vcell, yTickLabels, ptLabel, xlabStr, titleStr, wantLegend)
% Vcell{m} = [nPts x nSamples]

axes(ax); %#ok<LAXES>
cla(ax);
hold(ax,'on'); grid(ax,'on');

nModels = numel(Vcell);
nPts    = size(Vcell{1},1);

% Colors for the 3 points 
C = [ ...
    0.0000 0.4470 0.7410; % blue
    0.8500 0.3250 0.0980; % red/orange
    0.4660 0.6740 0.1880  % green
];

off  = linspace(-0.22, 0.22, nPts);
barH = 0.20;

ax.YDir = 'normal';
ax.Box = 'on';
ax.LineWidth = 1.5;
ax.FontSize  = 20;
ax.TickLabelInterpreter = 'latex';

allVals = [];

for m = 1:nModels
    y0 = m;
    for i = 1:nPts
        vals = Vcell{m}(i,:);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end

        q = quantile(vals, [0.025 0.5 0.975]);
        y = y0 + off(i);

        rectangle(ax, 'Position',[q(1), y-barH/2, (q(3)-q(1)), barH], ...
            'FaceColor', C(i,:), 'EdgeColor','none', 'FaceAlpha',0.95);

        plot(ax, [q(2) q(2)], [y-barH/2 y+barH/2], 'k-', 'LineWidth',5);

        allVals = [allVals; q(:)]; %#ok<AGROW>
    end
end

% Zero line
xline(ax, 0,'k-','LineWidth',2,'Alpha',0.60);

% Labels
yticks(ax, 1:nModels);
yticklabels(ax, yTickLabels);
yl = ylabel(ax, 'Candidate equation','Interpreter','latex','FontSize',24);

yl.Units = 'normalized';
yl.Position(1) = -0.08;   
xlabel(ax, xlabStr,'Interpreter','latex','FontSize',24);
if ~isempty(titleStr)
    title(ax, titleStr,'Interpreter','latex','FontWeight','bold','FontSize',22);
end


% Symmetric x-limits with padding
if ~isempty(allVals)
    xmin = min(allVals); xmax = max(allVals);
    M = max(abs([xmin xmax]));
    pad = 0.08 * max(1e-12, M);
    xlim(ax, [-(M+pad), (M+pad)]);
end

if wantLegend
    hP = gobjects(nPts,1);
    for i = 1:nPts
        hP(i) = plot(ax, nan,nan,'LineWidth',14,'Color',C(i,:));
    end
    hM = plot(ax, nan,nan,'k-','LineWidth',5);

    hLegend.handles = [hP; hM];
    hLegend.labels  = [ptLabel(:); {'Median (50th pct.)'}];
else
    hLegend = struct('handles',[],'labels',[]);
end

hold(ax,'off');
end

%% Signed C1 deviation at a point 
function dev = local_C1_signed_dev_at_point(f_handle, p, Re0, eD0, phys, hV, Re_switch, nA_min, nA_max, nB_ref)

rho = phys.rho; mu = phys.mu_ref; D = phys.D; L = phys.L;

% Convert Re0 -> V0 using Re = rho*V*D/mu
V0 = (Re0 * mu) / (rho * D);
V1 = V0 * (1 - hV);
V2 = V0 * (1 + hV);

Re1 = (rho * V1 * D) / mu;
Re2 = (rho * V2 * D) / mu;

f1 = f_handle(Re1, eD0, p);
f2 = f_handle(Re2, eD0, p);

if ~isfinite(f1) || ~isfinite(f2) || ~isreal(f1) || ~isreal(f2) || f1<=0 || f2<=0 || V1<=0 || V2<=0
    dev = NaN; return;
end

DP1 = f1 * (rho * V1^2) * L / (2*D);
DP2 = f2 * (rho * V2^2) * L / (2*D);

if DP1<=0 || DP2<=0 || ~isfinite(DP1) || ~isfinite(DP2)
    dev = NaN; return;
end

nloc = (log(DP2) - log(DP1)) / (log(V2) - log(V1));

if Re0 < Re_switch
    if nloc < nA_min
        dev = nloc - nA_min;      % negative
    elseif nloc > nA_max
        dev = nloc - nA_max;      % positive
    else
        dev = 0;
    end
else
    dev = nloc - nB_ref;
end
end

%% Candidate models 
function f = f_SR_cand1(Re, eD, p)
a1=p(1); a2=p(2); b1=p(3); b2=p(4); p_exp=p(5);
a3=p(6); c1=p(7); d1=p(8); e1=p(9); p2=p(10);
a4=p(11); b3=p(12); c2=p(13); p3=p(14); a5=p(15);

x1=Re; x2=eD;

term1 = a1.*x2;
term2 = a2.*x2.^((1./(x1.*x2.^b1) + b2).^p_exp);
base3 = c1.*x2 - ((d1.*x2 - e1)./x1);
term3 = a3.*base3.^p2;
denom4 = x2 + (1./(x1.*x2.^b3) + c2)./x1;
term4 = -a4./(denom4.^p3);
term5 = a5;

f = term1 + term2 + term3 + term4 + term5;
end

function f = f_SR_cand2(Re, eD, p)
a1=p(1); a2=p(2); b1=p(3); b2=p(4); pexp=p(5);
a3=p(6); c1=p(7); d1=p(8); p2=p(9); a4=p(10);

x1=Re; x2=eD;

term1 = a1 .* x2;
term2 = a2 .* x2 .^ ((1 ./ (x1 .* x2.^b1) + b2) .^ pexp);
base3 = c1 .* x2 + (x2 + d1) ./ x1;
term3 = a3 .* base3.^p2;
term4 = a4;

f = term1 + term2 + term3 + term4;
end

function f = f_SR_cand3(Re, eD, p)
a1=p(1); a2=p(2); b1=p(3); c1=p(4); p1=p(5);
a3=p(6); d1=p(7); d2=p(8); pexp=p(9);
a4=p(10); e1=p(11); e2=p(12); e3=p(13); p2=p(14); a5=p(15);

x1=Re; x2=eD;

term1 = a1 .* x2;
denom2 = x2 + (1 ./ (x1 .* x2.^b1) + c1) ./ x1;
term2 = -a2 ./ (denom2.^p1);
term3 = a3 .* x2 .^ ((1 ./ (x1 .* x2.^d1) + d2) .^ pexp);
base4 = e1 - ((e2 .* x2 - e3) ./ x1);
term4 = a4 .* base4.^p2;
term5 = a5;

f = term1 + term2 + term3 + term4 + term5;
end

function f = f_SR_cand4(Re, eD, p)
x1 = Re; x2 = eD;

a1=p(1); a2=p(2); A=p(3); B=p(4);
a3=p(5); c=p(6); d=p(7); pexp=p(8);
a4=p(9); q1=p(10); r1=p(11);
a5=p(12); q2=p(13); r2=p(14);
a6=p(15); a7=p(16);

term1 = a1 .* x2;
term2 = -a2 .* tanh(A./x1 + B);
base3 = x2.^2 - ( (c.*x2 - d./x1) ./ x1 );
term3 = a3 .* (base3).^pexp;
term4 = -a4 ./ ( (q1./x2) .^ (r1./x1) );
term5 = a5 ./ ( (q2./x2) .^ (r2./x1) );
term6 = -a6 ./ x1;
term7 = a7;

f = term1 + term2 + term3 + term4 + term5 + term6 + term7;
end

function f = f_Haaland_new(Re, eD, p)
A = p(1); B = p(2); C = p(3); D = p(4);
f = ( -A .* log10( (eD./B).^C + D ./ Re ) ).^(-2);
end
