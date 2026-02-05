function gp = spatial_evol_config(gp)

%SPATIAL_EVO_CONFIG Multigene config file for y = 1/(1+x1^-4) + 1/(1+x2^-4)

 

%run control

gp.runcontrol.pop_size = 1000; %population of 250 models                    

gp.runcontrol.runs = 1; %perform 2 runs that are merged at the end

gp.runcontrol.timeout = 600; %each run terminates after 60 seconds

gp.runcontrol.parallel.auto = true; %enable Parallel Computing if installed

gp.runcontrol.num_gen = 200; % Number of generations
 

%selection

gp.selection.tournament.size = 100;  % i.e. number og generations

gp.selection.tournament.p_pareto = 0.3; %encourages less complex models 

gp.selection.elite_fraction = 0.3; % approx. 1/3 models copied to next gen

 

%genes

gp.genes.max_genes = 5;  



%2d training grid in range -5 to 5 (676 points)

[x1,x2] = meshgrid(-5:0.4:5,-5:0.4:5);

x1 = x1(:); x2 = x2(:); %convert grid to vectors

y = (1./(1+x1.^-4)) + (1./(1+x2.^-4));

gp.userdata.ytrain = y;

gp.userdata.xtrain = [x1 x2];

  

%test grid in range -5 to 5 (2601 points)

[x1,x2] = meshgrid(-5:0.2:5,-5:0.2:5);

x1 = x1(:); x2 = x2(:); 

y = (1./(1+x1.^-4)) + (1./(1+x2.^-4));

gp.userdata.ytest = y;

gp.userdata.xtest = [x1 x2];



%function nodes

% gp.nodes.functions.name = {'times','minus','plus','rdivide','square',...
%     'sin','cos','exp','mult3','add3','sqrt','cube','power','negexp',...
%     'neg','abs','log'};

gp.nodes.functions.name = {'times','minus','plus','power','neg','rdivide'};