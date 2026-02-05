function [paretoIdx,domMat] = plot_pareto3(gp, fitnessMax, chosenIdx)
% PLOT_PARETO3  3-D Pareto scatter for GPTIPS population (minimize all 3)
%   [paretoIdx, domMat] = plot_pareto3(gp)
%   [paretoIdx, domMat] = plot_pareto3(gp, fitnessMax)
%   [paretoIdx, domMat] = plot_pareto3(gp, fitnessMax, chosenIdx)
%
% Expects gp.fitness.multiobj as [N x 3]:
%   col1 = fitness, col2 = complexity, col3 = constraintScore (all minimized)



% input / filtering 
if nargin < 2 || isempty(fitnessMax), fitnessMax = 1; end
if nargin < 3, chosenIdx = []; end

if ~isfield(gp,'fitness') || ~isfield(gp.fitness,'multiobj') || isempty(gp.fitness.multiobj)
    error('gp.fitness.multiobj not found or empty.');
end
M = gp.fitness.multiobj;
if size(M,2) ~= 3
    error('multiobj must be [N x 3] = [fitness, complexity, constraint].');
end

% remove NaN/Inf rows
badFinite = any(~isfinite(M),2);
validIdx  = find(~badFinite);
A         = M(validIdx,:);

% cut off very bad fitness
badFit   = A(:,1) > fitnessMax;
validIdx = validIdx(~badFit);
A        = A(~badFit,:);

nv = size(A,1);
if nv == 0
    error('No valid individuals to plot after filtering (check fitnessMax).');
end

% dominance / Pareto set
domMat = false(nv,nv);
for i = 1:nv
    Ai = A(i,:);
    for j = 1:nv
        Aj = A(j,:);
        if all(Aj <= Ai) && any(Aj < Ai)    % j dominates i  (minimization)
            domMat(i,j) = true;
        end
    end
end
isDominated  = any(domMat,2);
paretoRelIdx = find(~isDominated);        % indices in A
paretoIdx    = validIdx(paretoRelIdx);    % population indices

% map chosen indices (can be many) 
if isempty(chosenIdx)
    chosenIdx = [];
else
    chosenIdx = chosenIdx(:);             % column vector
end

nChosen = numel(chosenIdx);
chosenRelIdx = nan(nChosen,1);

for k = 1:nChosen
    iRel = find(validIdx == chosenIdx(k), 1, 'first');
    if isempty(iRel)
        warning('plot_pareto3:chosenFiltered', ...
            'Chosen index %d is not in the finite/fitness-filtered set. It will not be highlighted.', ...
            chosenIdx(k));
    else
        chosenRelIdx(k) = iRel;
    end
end

isChosenValid = isfinite(chosenRelIdx);
chosenRelIdx  = chosenRelIdx(isChosenValid);
chosenIdx     = chosenIdx(isChosenValid);
nChosenValid  = numel(chosenRelIdx);

%  plot 
fsz = 28;
popColor   = [0.15 0.15 0.85];  % population
pfFace     = [0.85 0.15 0.15];  % Pareto
pfEdge     = [0.10 0.10 0.10];

fig = figure('Name',sprintf('3-D Pareto (fitness \\le %.3g)',fitnessMax), ...
             'NumberTitle','off', 'Color','w');
set(fig,'Renderer','opengl');            % better 3-D

ax = axes('Parent',fig); hold(ax,'on'); set(ax,'FontSize',fsz);
ax.Box = 'on';
ax.LineWidth = 1.25;
ax.TickDir   = 'out';
grid(ax,'on');

% population
hAll = scatter3(ax, A(:,1), A(:,2), A(:,3), 28, ...
    'Marker','o','MarkerFaceColor',popColor,'MarkerEdgeColor','none');

% Pareto front
hPF  = scatter3(ax, A(paretoRelIdx,1), A(paretoRelIdx,2), A(paretoRelIdx,3), ...
    90, 'Marker','o','MarkerFaceColor',pfFace,'MarkerEdgeColor',pfEdge,'LineWidth',0.8);

% Optional transparency (guarded for compatibility)
if isprop(hAll,'MarkerFaceAlpha')
    hAll.MarkerFaceAlpha = 0.75;  % fade population
end
if isprop(hPF,'MarkerFaceAlpha')
    hPF.MarkerFaceAlpha  = 0.90;  % keep PF strong
end

% axis labels / view
xlabel(ax,'Fitness','FontSize',fsz);
ylabel(ax,'Complexity','FontSize',fsz);
zlabel(ax,'Constraint score','FontSize',fsz);
view(ax, 40, 20);

% tight limits: central cloud + chosen points
pr = [5 95];                   % 5–95% central region
fitLim  = prctile(A(:,1), pr);
compLim = prctile(A(:,2), pr);
consLim = prctile(A(:,3), pr);

% small padding
padFit  = 0.05*(fitLim(2)  - fitLim(1));
padComp = 0.05*(compLim(2) - compLim(1));
padCons = 0.05*(consLim(2) - consLim(1));

fitLim  = fitLim  + [-padFit padFit];
compLim = compLim + [-padComp padComp];
consLim = consLim + [-padCons padCons];

% expand limits so all highlighted models are visible
if nChosenValid > 0
    margin = 0.03;
    for k = 1:nChosenValid
        p = A(chosenRelIdx(k),:);
        fitLim(1)  = min(fitLim(1),  p(1)*(1 - margin));
        fitLim(2)  = max(fitLim(2),  p(1)*(1 + margin));
        compLim(1) = min(compLim(1), p(2)*(1 - margin));
        compLim(2) = max(compLim(2), p(2)*(1 + margin));
        consLim(1) = min(consLim(1), p(3)*(1 - margin));
        consLim(2) = max(consLim(2), p(3)*(1 + margin));
    end
end

xlim(ax, fitLim);
ylim(ax, compLim);
zlim(ax, consLim);

% chosen models: HALO + boxed labels + leader lines 
hChosen = gobjects(0);
hHalo   = gobjects(0);
% hLine   = gobjects(0);
% hTxt    = gobjects(0);

if nChosenValid > 0
    % colour & marker cycles for highlighted models
    colorList = [ 1.0 0.9 0.0;   
                  0.0 0.7 0.0;   
                  0.0 0.75 0.75; 
                  0.8 0.0 0.8;   
                  0.9 0.4 0.0 ]; 
    markerList = {'p','s','^','v','d','h'};
    nColors  = size(colorList,1);
    nMarks   = numel(markerList);

    % label offset controls
    rx = range(xlim(ax)); if rx==0, rx = 1; end
    ry = range(ylim(ax)); if ry==0, ry = 1; end
    rz = range(zlim(ax)); if rz==0, rz = 1; end

    baseXY = 0.045;   
    baseZ  = 0.050;   

    % spread labels around the point to avoid overlap 
    thetaList = deg2rad([25 155 285 70 200 330 115 245]); % enough unique directions
    nTheta = numel(thetaList);

    hChosen = gobjects(nChosenValid,1);
    hHalo   = gobjects(nChosenValid,1);
    hLine   = gobjects(nChosenValid,1);
    hTxt    = gobjects(nChosenValid,1);

    for k = 1:nChosenValid
        p = A(chosenRelIdx(k),:);    % [fitness, complexity, constraint]
        c = colorList( 1 + mod(k-1,nColors), : );
        if chosenIdx(k) == 171
            c = [1.0 0.0 1.0];   
        end
        m = markerList{ 1 + mod(k-1,nMarks) };

        
        hHalo(k) = scatter3(ax, p(1), p(2), p(3), ...
            520, m, 'MarkerFaceColor',[1 1 1], ...
            'MarkerEdgeColor',[1 1 1], 'LineWidth',1.0);

        
        hChosen(k) = scatter3(ax, p(1), p(2), p(3), ...
            300, m, 'MarkerFaceColor', c, ...
            'MarkerEdgeColor',[0 0 0], 'LineWidth',2.0);

        uistack(hHalo(k),'top');
        uistack(hChosen(k),'top');

    end
end

if isprop(hAll,'DataTipTemplate')
    hAll.DataTipTemplate.DataTipRows(1).Label = 'Fitness';
    hAll.DataTipTemplate.DataTipRows(2).Label = 'Complexity';
    hAll.DataTipTemplate.DataTipRows(3).Label = 'Constraint';
    rowIdxAll = dataTipTextRow('Index', validIdx);
    hAll.DataTipTemplate.DataTipRows(end+1) = rowIdxAll;
end

if isprop(hPF,'DataTipTemplate')
    hPF.DataTipTemplate.DataTipRows(1).Label = 'Fitness';
    hPF.DataTipTemplate.DataTipRows(2).Label = 'Complexity';
    hPF.DataTipTemplate.DataTipRows(3).Label = 'Constraint';
    rowIdxPF = dataTipTextRow('Index', paretoIdx);
    hPF.DataTipTemplate.DataTipRows(end+1) = rowIdxPF;
end

if nChosenValid > 0
    for k = 1:nChosenValid
        if isprop(hChosen(k),'DataTipTemplate')
            hChosen(k).DataTipTemplate.DataTipRows(1).Label = 'Fitness';
            hChosen(k).DataTipTemplate.DataTipRows(2).Label = 'Complexity';
            hChosen(k).DataTipTemplate.DataTipRows(3).Label = 'Constraint';
            rowIdxChosen = dataTipTextRow('Index', chosenIdx(k));
            hChosen(k).DataTipTemplate.DataTipRows(end+1) = rowIdxChosen;
        end
    end
end

% --------- legend ----------
legHandles = [hAll hPF];
legStrings = {'population','Pareto front'};

for k = 1:nChosenValid
    legHandles(end+1) = hChosen(k); %#ok<AGROW>
    legStrings{end+1} = sprintf('index %d', chosenIdx(k)); %#ok<AGROW>
end

leg = legend(ax, legHandles, legStrings, ...
    'Location','northeast','Box','on','Interpreter','none','FontSize',30);
leg.ItemTokenSize = [28 9];

hold(ax,'off');

% summary
fprintf('Pareto front: %d points out of %d valid individuals (fitness <= %.3g).\n', ...
        numel(paretoIdx), nv, fitnessMax);
fprintf('Pareto indices (population indices):\n');
fprintf(' %d', paretoIdx(:));
fprintf('\n');

if nChosenValid > 0
    for k = 1:nChosenValid
        p = A(chosenRelIdx(k),:);
        fprintf('Chosen model index %d at [fitness=%.4g, complexity=%.4g, constraint=%.4g].\n', ...
            chosenIdx(k), p(1), p(2), p(3));
    end
end
end
