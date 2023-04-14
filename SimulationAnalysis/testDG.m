clear;clc;
close all;
%% current parallel pool
if isempty(gcp('nocreate'))
    parpool(32);
end
enhancer_idx = [51,55,60,65,70,75,80,85,90,95];
promoter_idx = [49,45,40,35,30,25,20,15,10,5];
for idx = 1:length(enhancer_idx)
    %% Parameters setting
    % Upstream: Chromatin conformation
    input_options.simulated_on      = true;
    input_options.EP_flag           = true;
    input_options.attraction_coef   = 0.1; % K, E-P communication coefficient
    input_options.enhancer_index    = enhancer_idx(idx); % enhancer index number
    input_options.promoter_index    = promoter_idx(idx); % promoter index number
    
    % Downstream: Transcriptional bursting
    input_options.k_on_max       = 1; % off to on max
    input_options.k_on           = 0.06;  % off to on min
    input_options.k_off          = 0.2; % on to off
    input_options.mu_max         = 3; % generate mRNA
    input_options.mu             = 0.5; % generate mRNA
    input_options.delta          = 0.1; % mRNA degradation
    input_options.friction_coef  = 50; % friction coefficient
    input_options.simulation_num    = 32; % simulation number
    input_options.simulation_reaction_step = 1000000; %
    input_options.simulation_time   = 600000; % [s]
    
    input_options.result_base_folder = fullfile(pwd, 'ResultDG');
    input_options.filename = sprintf('//E_%d_P_%d_%f_%f_%f_%f_%f_%f_%f',[input_options.enhancer_index,...
        input_options.promoter_index,input_options.k_on_max,input_options.k_on,input_options.k_off,...
        input_options.mu_max,input_options.mu,input_options.delta,input_options.friction_coef]);
    if exist([input_options.result_base_folder,input_options.filename],'file') == 0
        mkdir(input_options.result_base_folder,input_options.filename);
    end
    
    %% Parallel computing
    tic;
    if input_options.simulated_on
        parfor s_idx = 1:input_options.simulation_num
            % Load parameters
            params = ParametersBurst(input_options);
            % simulation
            result = SimulateBurst(params,s_idx);
        end
    end
    timerVal = toc;
    X = ['Total simulation time:',num2str(timerVal)];
    disp(X)
    %% Analysing Burst
    % Only parameters are needed, and data are loaded from the .mat files
    tic;
    params = ParametersBurst(input_options);
    results = AnalyseBurst(params);
    filename = sprintf('//E_%d_P_%d_%f_%f_%f_%f_%f_%f_%f.mat',[input_options.enhancer_index,...
        input_options.promoter_index,input_options.k_on_max,input_options.k_on,input_options.k_off,...
        input_options.mu_max,input_options.mu,input_options.delta,input_options.friction_coef]);
    save([input_options.result_base_folder,filename],'results');
    timerVal = toc;
    X = ['Analysing time:',num2str(timerVal)];
    disp(X)
    hold on
    if results.params.simulated_on == true
        h = histogram(results.mRNA_total,'BinEdges',0:max(results.mRNA_total),"Normalization",'probability','DisplayName','simulation');
    end
    plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.mRNA_Prob,'LineWidth',1,'DisplayName','mixed')
%      plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.PDF.mRNA_Prob_v1,'LineWidth',1,'DisplayName','fast')
%      plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.PDF.mRNA_Prob_v2,'LineWidth',1,'DisplayName','slow')
    set(figure1,'position',[300 400 280 190]);
    set(gca,'TickLength',[0.02,0.025]);
    legend1 = legend(gca,'show');
    results.weight
    box on
    xlim([0 200])
   title(['E-P ',num2str(abs(enhancer_idx(idx)-promoter_idx(idx))),...
       '  weight ',num2str(results.weight(1)) ]);
   data.mRNAdistri{1,idx} = results.mRNA_Prob;
   data.EPgenomedistance(idx) = abs(enhancer_idx(idx)-promoter_idx(idx));
   data.mRNAdistribin{1,idx} = 0:1:length(results.mRNA_Prob)-1;
end

