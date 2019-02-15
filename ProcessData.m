%% PHASE AVERAGE ALL MACH 2 DATA ORION
clear
cmq = [];
cma = [];
emq = [];

% Square
in.diam = 5;
in.base_dir   ='/home/fagleyc/Projects/Capsule/Exp/2018Fall/Orion/8_28_Square';
in.aoff       ='Run001';    %% Aoff file indec
in.fprefix    ='Run_m2_';   %% file name prefix
in.files      =2:1:91;        %% Array of files
in.dcfile     ='Run000';    %% DC File index
[data o] = run_multiple_datasets(in);

a0 = o.sum.a0_resamp;
cmq = [cmq; o.sum.Cmq_resamp];
cma = [cma; o.sum.Cma_resamp];
emq = [emq; o.sum.error_freq_resamp];

% Chirp
in.base_dir   ='/home/fagleyc/Projects/Capsule/Exp/2018Fall/Orion/8_30_Chirp/m020';
in.aoff       ='Run001';    %% Aoff file indec
in.fprefix    ='Run';   %% file name prefix
in.files      =2:1:87;        %% Array of files
in.dcfile     ='Run000';    %% DC File index
[data o] = run_multiple_datasets(in);

% a0 = o.sum.a0_resamp;
cmq = [cmq; o.sum.Cmq_resamp];
cma = [cma; o.sum.Cma_resamp];
emq = [emq; o.sum.error_freq_resamp];

% 

% Smaller diameter rod larger capsule
in.diam = 6;
in.base_dir   ='/home/fagleyc/Projects/Capsule/Exp/2018Fall/Orion/10_30_Square/';
in.aoff       ='Run001';    %% Aoff file indec
in.fprefix    ='Run_m2_';   %% file name prefix
in.files      =1:1:91;        %% Array of files
in.dcfile     ='Run000';    %% DC File index
[data o] = run_multiple_datasets(in);

cmq = [cmq; o.sum.Cmq_resamp];
cma = [cma; o.sum.Cma_resamp];
emq = [emq; o.sum.error_freq_resamp];

in.base_dir   ='/home/fagleyc/Projects/Capsule/Exp/2018Fall/Orion/10_30_Square/';
in.aoff       ='Run001';    %% Aoff file indec
in.fprefix    ='Run_m2_A3_';   %% file name prefix
in.files      =1:1:91;        %% Array of files
in.dcfile     ='Run000';    %% DC File index
[data o] = run_multiple_datasets(in);

cmq = [cmq; o.sum.Cmq_resamp];
cma = [cma; o.sum.Cma_resamp];
emq = [emq; o.sum.error_freq_resamp];

%
Figure CMQ
% errorbar(o.sum.a0,o.sum.Cmq,o.sum.error,'.')
plot(a0,mean(cmq,1))
hold all
hline(0,':')
xlabel('$\alpha$')
ylabel('$C_{mq}$')

Figure CMA
plot(a0,mean(cma,1),'o')
