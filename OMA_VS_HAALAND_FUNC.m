function [C1, C2, C3, C4] = OMA_VS_HAALAND_FUNC()
%OMA_VS_HAALAND_FUNC  Visual C1–C4 diagnostics for the OMA model.
%
%   [C1,C2,C3,C4] = OMA_VS_HAALAND_FUNC()
%
% Uses your OMA correlation:
%   f = f_OMA_model(Re, eD)
%
% Makes the same plots as check_C1C2C3.m:
%   C1: n(U)      = d log(ΔP) / d log(U)         (U sweep)
%   C2: s(e)      = d log(ΔP) / d log(e)         (e sweep)
%   C3: α(μ)      = d log(ΔP) / d log(μ)         (μ sweep)
%   C4: γ(ρ)      = d log(ΔP) / d log(ρ)         (ρ sweep)
%
% PLUS it returns scalar “violation magnitudes”:
%   C1 = max_eD,max_U |n(U) - 2|
%   C2 = max_Re,max_e |s(e)|          (target 0)
%   C3 = max_eD,max_μ |α(μ) - 1|      (target 1)
%   C4 = max_eD,max_ρ |γ(ρ) - 1|      (target 1)

f_SR = @f_OMA_model; 

% --------- defaults (copied from check_C1C2C3) ---------
eD_vals = [3.65e-06, 0.00100603621730382, ...
           0.00202429149797571, 0.00404203718674212, ...
           0.00821692686935086, 0.0164338537387017, ...
           0.0331674958540630];                 % Nikuradse-like
D      = 0.05;
Lpipe  = 1.0;
rho0   = 1000;
mu_ref = 1e-3;

% Grid helpers  
Re_Moody = logspace(log10(2500), 7, 100).';   
Re_C1    = Re_Moody;                         % U-sweep Re grid 
Re_C2    = [1e4, 1e5, 3e5];                  % C2 slices (still row, fine)

e_grid   = logspace(log10(min(eD_vals)*D), ...
                   log10(max(eD_vals)*D), numel(Re_Moody)).';  % 50x1

% Haaland reference
f_H = @(Re,eD) (-1.8 .* log10( ((eD)./3.7).^1.11 + 6.9 ./ Re ) ).^(-2);

% helper: ΔP from f
dp_from_f = @(f,Re,rho,mu) ...
    f .* ( rho .* ((Re .* mu ./ (rho .* D)).^2) .* Lpipe ) ./ (2 .* D);

col_eD = lines(numel(eD_vals));
col_Re = lines(numel(Re_C2));

%  C1 : n = d log ΔP / d log U     (U sweep at fixed μ,ρ,D,L)
Ugrid   = Re_C1 .* mu_ref ./ (rho0 * D);  % U corresponding to Re_C1

figC1 = figure('Name','C1: n','NumberTitle','off');
set(figC1,'Position',[100 100 900 650]);
ax1 = axes('Parent',figC1); hold(ax1,'on'); apply_axes_style(ax1);

% Header text + Re range
txt1 = sprintf(['Constants:  D = %.2f m,  L = %.1f m,  \\rho = %.0f kg/m^3, ', ...
                '\\mu = %.3g Pa·s'], D, Lpipe, rho0, mu_ref);
txt2 = sprintf('Re range: %.1e – %.1e', min(Re_C1), max(Re_C1));
text(ax1,0.01,0.96,txt1,'Units','normalized','Interpreter','tex', ...
    'FontSize',20,'VerticalAlignment','top','BackgroundColor','w');
text(ax1,0.01,0.88,txt2,'Units','normalized','Interpreter','tex', ...
    'FontSize',20,'VerticalAlignment','top','BackgroundColor','w');

yline(ax1,1,'k:','LineWidth',1.2);
yline(ax1,2,'k:','LineWidth',1.2);

legLines = []; legLabs = {};
C1_eD    = nan(numel(eD_vals),1);   % for scalar C1

for i = 1:numel(eD_vals)
    eD   = eD_vals(i);
    col  = col_eD(i,:);

    % Haaland
    fH   = f_H(Re_C1, eD*ones(size(Re_C1)));
    DPH  = dp_from_f(fH, Re_C1, rho0, mu_ref);
    maskH= isfinite(DPH) & DPH>0;
    [U_mid_H, nH] = local_loglog_slope(Ugrid(maskH), DPH(maskH), 7);
    plot(ax1, U_mid_H, nH, '-', 'Color', col, 'LineWidth', 1.7, ...
         'HandleVisibility','off');

    % OMA (plays the SR role, dashed)
    fS   = f_SR(Re_C1, eD*ones(size(Re_C1)));
    DPS  = dp_from_f(fS, Re_C1, rho0, mu_ref);
    maskS= isfinite(DPS) & DPS>0;
    [U_mid_S, nS] = local_loglog_slope(Ugrid(maskS), DPS(maskS), 7);
    plot(ax1, U_mid_S, nS, '--', 'Color', col, 'LineWidth', 1.9, ...
         'HandleVisibility','off');

    legLines(end+1) = plot(ax1, NaN,NaN, '-', 'Color', col, 'LineWidth',2.1); %#ok<AGROW>
    legLabs{end+1}  = eD_label(eD);

    if ~isempty(nS)
        C1_eD(i) = max(abs(nS - 2));   % target n ≈ 2
        [~, idx] = max(abs(nS - 2));
        Re_bad = (rho0*U_mid_S(idx)*D)/mu_ref;
        fprintf('C1 (OMA): %s, max |n-2| at U = %.3g m/s, Re ≈ %.3g, n = %.3f\n', ...
                eD_label(eD), U_mid_S(idx), Re_bad, nS(idx));
    end
end

set(ax1,'XScale','log'); xlim(ax1,[min(Ugrid) max(Ugrid)]);
xlabel(ax1,'$\overline{U}_m\,[\mathrm{m\,s^{-1}}]$','Interpreter','latex','FontSize',34);
ylabel(ax1,'$\chi=\mathrm{\partial}\log(\Delta P)/\mathrm{\partial}\log(\overline{U}_m)$', ...
       'Interpreter','latex','FontSize',34);


% Cross-hair at a reference Re (for orientation)
Re_ref_C1 = 1e5;
U_ref_C1  = Re_ref_C1 * mu_ref / (rho0*D);
xline(ax1,U_ref_C1,'k:','LineWidth',1.2);

hKey1 = plot(ax1,NaN,NaN,'k-','LineWidth',2.1);
hKey2 = plot(ax1,NaN,NaN,'k--','LineWidth',2.1);
labels = [{'Haaland (solid)','OMA Model (dashed)'}  legLabs(:)'];
leg = legend(ax1,[hKey1 hKey2 legLines], labels, ...
    'Location','southeast','Box','on','Interpreter','latex');
leg.ItemTokenSize = [24 10];
leg.FontSize = 20;

C1 = max(C1_eD);   % scalar

%  C2 : s = d log ΔP / d log e      (e sweep, Re = 1e4,1e5,3e5)
figC2 = figure('Name','C2: s','NumberTitle','off');
set(figC2,'Position',[150 90 900 650]);
ax1 = axes('Parent',figC2); hold(ax1,'on'); apply_axes_style(ax1);

txt1 = sprintf(['Constants:  D = %.2f m,  L = %.1f m,  \\rho = %.0f kg/m^3, ', ...
                '\\mu = %.3g Pa·s'], D, Lpipe, rho0, mu_ref);
txt2 = sprintf('Re slices: %.0e, %.0e, %.0e', Re_C2(1), Re_C2(2), Re_C2(3));
text(ax1,0.01,0.96,txt1,'Units','normalized','Interpreter','tex', ...
    'FontSize',20,'VerticalAlignment','top','BackgroundColor','w');
text(ax1,0.01,0.88,txt2,'Units','normalized','Interpreter','tex', ...
    'FontSize',20,'VerticalAlignment','top','BackgroundColor','w');

yline(ax1,0,'k:','LineWidth',1.2);

legLines = []; legLabs = {};
C2_Re    = nan(numel(Re_C2),1);

for k = 1:numel(Re_C2)
    Re0 = Re_C2(k);
    col = col_Re(k,:);
    eD_arg = e_grid./D;

    % Haaland
    fH  = f_H(Re0*ones(size(e_grid)), eD_arg);
    DPH = dp_from_f(fH, Re0, rho0, mu_ref);
    maskH = isfinite(DPH) & DPH>0;
    [e_mid_H, sH] = local_loglog_slope(e_grid(maskH), DPH(maskH), 7);
    plot(ax1, e_mid_H, sH, '-', 'Color', col, 'LineWidth', 1.7, ...
         'HandleVisibility','off');

    % OMA
    fS  = f_SR(Re0*ones(size(e_grid)), eD_arg);
    DPS = dp_from_f(fS, Re0, rho0, mu_ref);
    maskS = isfinite(DPS) & DPS>0;
    [e_mid_S, sS] = local_loglog_slope(e_grid(maskS), DPS(maskS), 7);
    plot(ax1, e_mid_S, sS, '--', 'Color', col, 'LineWidth', 1.9, ...
         'HandleVisibility','off');

    legLines(end+1) = plot(ax1,NaN,NaN,'-','Color',col,'LineWidth',2.1); %#ok<AGROW>
    legLabs{end+1}  = sprintf('Re = %.0e', Re0); %#ok<AGROW>

    if ~isempty(sS)
        C2_Re(k) = max(abs(sS));   % target ≈ 0
        [~, idx] = max(abs(sS));
        fprintf('C2 (OMA): Re = %.0e, worst s at e = %.3g m, s = %.3f\n', ...
                Re0, e_mid_S(idx), sS(idx));
    end
end

set(ax1,'XScale','log');
xlabel(ax1,'$\epsilon\ [\mathrm{m}]$','Interpreter','latex','FontSize',34);
ylabel(ax1,'$s=\mathrm{\partial}\log(\Delta P)/\mathrm{\partial}\log(\epsilon)$', ...
       'Interpreter','latex','FontSize',34);

hKey1 = plot(ax1,NaN,NaN,'k-','LineWidth',2.1);
hKey2 = plot(ax1,NaN,NaN,'k--','LineWidth',2.1);
labels = [{'Haaland (solid)','OMA Model (dashed)'}  legLabs(:)'];
leg = legend(ax1,[hKey1 hKey2 legLines], labels, ...
    'Location','southeast','Box','on','Interpreter','latex');
leg.ItemTokenSize = [24 10];
leg.FontSize = 20;

C2 = max(C2_Re);

%  C3 : alpha(μ) = d log ΔP / d log μ     (μ sweep at fixed U)
U_mu = 0.2;
mu_grid = (rho0 * U_mu * D) ./ Re_Moody(:);   % column vector
Re_mu   = Re_Moody(:);

figC3 = figure('Name','C3: alpha','NumberTitle','off');
set(figC3,'Position',[200 80 900 650]);
ax1 = axes('Parent',figC3); hold(ax1,'on'); apply_axes_style(ax1);

txt1 = sprintf(['Constants:  D = %.2f m,  L = %.1f m,  \\rho = %.0f kg/m^3, ', ...
               'U = %.2f m/s'], D, Lpipe, rho0, U_mu);
txt2 = sprintf('Re range: %.1e – %.1e', min(Re_mu), max(Re_mu));
text(ax1,0.01,0.96,txt1,'Units','normalized','Interpreter','tex', ...
    'FontSize',20,'VerticalAlignment','top','BackgroundColor','w');
text(ax1,0.01,0.88,txt2,'Units','normalized','Interpreter','tex', ...
    'FontSize',20,'VerticalAlignment','top','BackgroundColor','w');

yline(ax1,0,'k:','LineWidth',1.2);
yline(ax1,1,'k:','LineWidth',1.2);

legLines = []; legLabs = {};
C3_eD    = nan(numel(eD_vals),1);

for i = 1:numel(eD_vals)
    eD = eD_vals(i);
    col = col_eD(i,:);

    % Haaland
    fH  = f_H(Re_mu, eD*ones(size(Re_mu)));
    DPH = dp_from_f(fH, Re_mu, rho0, mu_grid);
    maskH = isfinite(DPH) & DPH>0 & isfinite(mu_grid);
    [mu_mid_H, aH] = local_loglog_slope(mu_grid(maskH), DPH(maskH), 7);
    plot(ax1, mu_mid_H, aH, '-', 'Color', col, 'LineWidth',1.7, ...
         'HandleVisibility','off');

    % OMA
    fS  = f_SR(Re_mu, eD*ones(size(Re_mu)));
    DPS = dp_from_f(fS, Re_mu, rho0, mu_grid);
    maskS = isfinite(DPS) & DPS>0 & isfinite(mu_grid);
    [mu_mid_S, aS] = local_loglog_slope(mu_grid(maskS), DPS(maskS), 7);
    plot(ax1, mu_mid_S, aS, '--', 'Color', col, 'LineWidth',1.9, ...
         'HandleVisibility','off');

    legLines(end+1) = plot(ax1,NaN,NaN,'-','Color',col,'LineWidth',2.1); %#ok<AGROW>
    legLabs{end+1}  = eD_label(eD);

    if ~isempty(aS)
        [~,idx] = max(abs(aS-1));
        C3_eD(i) = abs(aS(idx)-1);    % target ≈ 1
        Re_bad = rho0*U_mu*D/mu_mid_S(idx);
        fprintf('C3 (OMA): %s, worst alpha(mu) at mu = %.3g Pa·s, Re ≈ %.3g, alpha = %.3f\n', ...
                eD_label(eD), mu_mid_S(idx), Re_bad, aS(idx));
    end
end

set(ax1,'XScale','log','XDir','reverse');  % rightwards = decreasing mu = higher Re
xlabel(ax1,'$\mu\,[\mathrm{Pa\,s}]$','Interpreter','latex','FontSize',34);
ylabel(ax1,'$\alpha=\mathrm{\partial}\log(\Delta P)/\mathrm{\partial}\log(\mu)$', ...
       'Interpreter','latex','FontSize',34);

% cross-hair at mu_ref
mu_ref_line = mu_ref;
xline(ax1,mu_ref_line,'k:','LineWidth',1.2);

hKey1 = plot(ax1,NaN,NaN,'k-','LineWidth',2.1);
hKey2 = plot(ax1,NaN,NaN,'k--','LineWidth',2.1);
labels = [{'Haaland (solid)','OMA Model (dashed)'}  legLabs(:)'];
leg = legend(ax1,[hKey1 hKey2 legLines], labels, ...
    'Location','southeast','Box','on','Interpreter','latex');
leg.ItemTokenSize = [24 10];
leg.FontSize = 20;

C3 = max(C3_eD);

%  C4 : gamma(ρ) = d log ΔP / d log ρ   (ρ sweep at fixed U)
U_rho = 0.2;
rho_grid = Re_Moody(:) * (mu_ref / (U_rho * D));
Re_rho   = Re_Moody(:);

figC4 = figure('Name','C4: gamma','NumberTitle','off');
set(figC4,'Position',[220 70 900 650]);
ax1 = axes('Parent',figC4); hold(ax1,'on'); apply_axes_style(ax1);

txt1 = sprintf(['Constants:  D = %.2f m,  L = %.1f m,  \\mu = %.3g Pa·s, ', ...
               'U = %.2f m/s'], D, Lpipe, mu_ref, U_rho);
txt2 = sprintf('Re range: %.1e – %.1e', min(Re_rho), max(Re_rho));
text(ax1,0.01,0.96,txt1,'Units','normalized','Interpreter','tex', ...
    'FontSize',20,'VerticalAlignment','top','BackgroundColor','w');
text(ax1,0.01,0.88,txt2,'Units','normalized','Interpreter','tex', ...
    'FontSize',20,'VerticalAlignment','top','BackgroundColor','w');

yline(ax1,0,'k:','LineWidth',1.2);
yline(ax1,1,'k:','LineWidth',1.2);

legLines = []; legLabs = {};
C4_eD    = nan(numel(eD_vals),1);

for i = 1:numel(eD_vals)
    eD = eD_vals(i);
    col = col_eD(i,:);

    % Haaland
    fH  = f_H(Re_rho, eD*ones(size(Re_rho)));
    DPH = dp_from_f(fH, Re_rho, rho_grid, mu_ref);
    maskH = isfinite(DPH) & DPH>0;
    [rho_mid_H, gH] = local_loglog_slope(rho_grid(maskH), DPH(maskH), 7);
    plot(ax1, rho_mid_H, gH, '-', 'Color', col, 'LineWidth',1.7, ...
         'HandleVisibility','off');

    % OMA
    fS  = f_SR(Re_rho, eD*ones(size(Re_rho)));
    DPS = dp_from_f(fS, Re_rho, rho_grid, mu_ref);
    maskS = isfinite(DPS) & DPS>0;
    [rho_mid_S, gS] = local_loglog_slope(rho_grid(maskS), DPS(maskS), 7);
    plot(ax1, rho_mid_S, gS, '--', 'Color', col, 'LineWidth',1.9, ...
         'HandleVisibility','off');

    legLines(end+1) = plot(ax1,NaN,NaN,'-','Color',col,'LineWidth',2.1); %#ok<AGROW>
    legLabs{end+1}  = eD_label(eD);

    if ~isempty(gS)
        [~,idx] = max(abs(gS-1));
        C4_eD(i) = abs(gS(idx)-1);   % target ≈ 1
        Re_bad = rho_mid_S(idx)*U_rho*D/mu_ref;
        fprintf('C4 (OMA): %s, worst gamma at rho = %.3g kg/m^3, Re ≈ %.3g, gamma = %.3f\n', ...
                eD_label(eD), rho_mid_S(idx), Re_bad, gS(idx));
    end
end

set(ax1,'XScale','log');
xlabel(ax1,'$\rho\,[\mathrm{kg\,m^{-3}}]$','Interpreter','latex','FontSize',34);
ylabel(ax1,'$\gamma=\mathrm{\partial}\log(\Delta P)/\mathrm{\partial}\log(\rho)$', ...
       'Interpreter','latex','FontSize',34);

rho_ref_line = rho0;
xline(ax1,rho_ref_line,'k:','LineWidth',1.2);

hKey1 = plot(ax1,NaN,NaN,'k-','LineWidth',2.1);
hKey2 = plot(ax1,NaN,NaN,'k--','LineWidth',2.1);
labels = [{'Haaland (solid)','OMA Model (dashed)'}  legLabs(:)'];
leg = legend(ax1,[hKey1 hKey2 legLines], labels, ...
    'Location','southeast','Box','on','Interpreter','latex');
leg.ItemTokenSize = [24 10];
leg.FontSize = 20;

C4 = max(C4_eD);

fprintf('\nOMA constraint scalars (vs Haaland):\n');
fprintf('   C1 = max_eD |n(U)-2|      ≈ %.3g\n', C1);
fprintf('   C2 = max_{Re,e} |s|    ≈ %.3g\n', C2);
fprintf('   C3 = max_eD |alpha-1|     ≈ %.3g\n', C3);
fprintf('   C4 = max_eD |gamma-1|     ≈ %.3g\n', C4);

end   


% helper: robust local log–log slope 
function [x_mid, slope] = local_loglog_slope(x, y, win)
x   = x(:); y = y(:);
mask= isfinite(x) & isfinite(y) & x>0 & y>0;
x   = x(mask); y = y(mask);
n   = numel(x);
if n < 3
    x_mid = []; slope = [];
    return;
end
if nargin < 3 || isempty(win), win = 7; end
win = max(3, 2*floor(win/2)+1);  
h   = (win-1)/2;

lx = log(x); ly = log(y);
x_mid = nan(n,1); slope = nan(n,1);
for i = 1:n
    i1 = max(1, i-h);
    i2 = min(n, i+h);
    idx = i1:i2;
    if numel(idx) >= 3
        p = polyfit(lx(idx), ly(idx), 1);
        slope(i) = p(1);
        x_mid(i) = exp(mean(lx(idx)));
    end
end
mask2 = isfinite(slope) & isfinite(x_mid);
x_mid = x_mid(mask2);
slope = slope(mask2);
end


% helper: legend label for e/D 
function str = eD_label(eD)
if eD < 1e-4
    str = sprintf('$\\varepsilon/D = %.3g$', eD);
else
    str = sprintf('$\\varepsilon/D = %.3f$', eD);
end
end



function apply_axes_style(ax)
set(ax, ...
    'FontSize', 30, ...
    'FontName', 'Times New Roman', ...
    'LineWidth', 1.5, ...
    'TickDir', 'out', ...
    'Layer', 'top', ...
    'XMinorGrid', 'on', ...
    'YMinorGrid', 'on', ...
    'TickLabelInterpreter','latex');

ax.GridAlpha        = 0.25;
ax.MinorGridAlpha   = 0.40;
ax.GridColor        = 0.6*[1 1 1];
ax.MinorGridColor   = 0.85*[1 1 1];
ax.MinorGridLineStyle = ':';

grid(ax,'on');
box(ax,'on');
end

