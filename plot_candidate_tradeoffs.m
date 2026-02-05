function plot_candidate_tradeoffs()
% PLOT_CANDIDATE_TRADEOFFS
% Visualise J_err, J_phys, J_comp for the three chosen candidates.

% data 
candID  = [1 2 3];                   % Candidate number in the paper
eqIdx   = [714 541 171];             % Equation indices
J_err   = [0.0120 0.0151 0.0122];
J_comp  = [45     26     41   ];
J_phys  = [0.168  0.154  0.178];

% For x–axis we use positions 1..3, and label them with the equation IDs
x = 1:numel(candID);
xtickLabels = arrayfun(@num2str, eqIdx, 'UniformOutput', false);

% figure + common style 
fig = figure('Color','w','Position',[100 100 900 600]);

fs_axis   = 14;   % axis tick font size
fs_label  = 16;   % axis label font size
fs_title  = 16;   % subplot title font size
ms        = 7;    % marker size
lw        = 1.5;  % line width

% Small helper for repeated formatting
    function style_axes()
        set(gca,'XTick',x,'XTickLabel',xtickLabels, ...
                'FontSize',fs_axis,'LineWidth',1);
        grid on;
        xlim([0.8 3.2]);
    end

% Fitness vs candidate
subplot(3,1,1);
plot(x, J_err,'-o','MarkerFaceColor','w', ...
     'MarkerSize',ms,'LineWidth',lw);
style_axes();
ylabel('fitness $J_{\mathrm{err}}$','Interpreter','latex', ...
       'FontSize',fs_label);
title('Fitness trade-off','FontSize',fs_title);

for k = 1:numel(x)
    text(x(k)+0.05, J_err(k), sprintf('Cand %d',candID(k)), ...
        'FontSize',11,'VerticalAlignment','bottom');
end

% Physics score vs candidate
subplot(3,1,2);
plot(x, J_phys,'-o','MarkerFaceColor','w', ...
     'MarkerSize',ms,'LineWidth',lw);
style_axes();
ylabel('constraint score $J_{\mathrm{phys}}$','Interpreter','latex', ...
       'FontSize',fs_label);
title('Physics-compatibility trade-off','FontSize',fs_title);

for k = 1:numel(x)
    text(x(k)+0.05, J_phys(k), sprintf('Cand %d',candID(k)), ...
        'FontSize',11,'VerticalAlignment','bottom');
end

% Complexity vs candidate
subplot(3,1,3);
plot(x, J_comp,'-o','MarkerFaceColor','w', ...
     'MarkerSize',ms,'LineWidth',lw);
style_axes();
ylabel('complexity $J_{\mathrm{comp}}$','Interpreter','latex', ...
       'FontSize',fs_label);
xlabel('Equation ID','FontSize',fs_label);
title('Structural complexity trade-off','FontSize',fs_title);

for k = 1:numel(x)
    text(x(k)+0.05, J_comp(k), sprintf('Cand %d',candID(k)), ...
        'FontSize',11,'VerticalAlignment','bottom');
end

% tighten layout a bit
subplot(3,1,1); pos1 = get(gca,'Position');
subplot(3,1,2); pos2 = get(gca,'Position');
subplot(3,1,3); pos3 = get(gca,'Position');

% manually reduce vertical gaps
pos1(2) = 0.70; pos1(4) = 0.22;
pos2(2) = 0.40; pos2(4) = 0.22;
pos3(2) = 0.10; pos3(4) = 0.22;
subplot(3,1,1); set(gca,'Position',pos1);
subplot(3,1,2); set(gca,'Position',pos2);
subplot(3,1,3); set(gca,'Position',pos3);

end
