clear;clc;
close all;
%% current parallel pool
if isempty(gcp('nocreate'))
    parpool(32);
end
K_coefficient =  0.05:0.1:0.95;
for idx =  1:length(K_coefficient)
    %% Parameters setting
    % Upstream: Chromatin conformation
    input_options.simulated_on      = true;
    input_options.EP_flag           = true;
    input_options.attraction_coef   = K_coefficient(idx); % K, E-P communication coefficient
    input_options.enhancer_index    = 25; % enhancer index number
    input_options.promoter_index    = 75; % promoter index number
    
    % Downstream: Transcriptional bursting
    input_options.k_on_max       = 1; % off to on max
    input_options.k_on           = 0.06;  % off to on min
    input_options.k_off          = 0.20; % on to off
    input_options.mu_max         = 3; % generate mRNA
    input_options.mu             = 0.5; % generate mRNA
    input_options.delta          = 0.1; % mRNA degradation
    input_options.friction_coef  = 50; % friction coefficient
    input_options.simulation_num    = 32; % simulation number
    input_options.simulation_reaction_step = 1000000; %
    input_options.simulation_time   = 600000; % [s]
    
    input_options.result_base_folder = fullfile(pwd, 'ResultKEP');
    input_options.filename = sprintf('//%f_%f_%f_%f_%f_%f_%f_%f',[input_options.attraction_coef,...
        input_options.k_on_max,input_options.k_on,input_options.k_off,input_options.mu_max...
        input_options.mu,input_options.delta,input_options.friction_coef]);
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
    filename = sprintf('//%f_%f_%f_%f_%f_%f_%f_%f.mat',[input_options.attraction_coef,...
        input_options.k_on_max,input_options.k_on,input_options.k_off,input_options.mu_max...
        input_options.mu,input_options.delta,input_options.friction_coef]);
    save([input_options.result_base_folder,filename],'results');
    timerVal = toc;
    X = ['Analysing time:',num2str(timerVal)];
    hold on
    if results.params.simulated_on == true
        h = histogram(results.mRNA_total,'BinEdges',0:max(results.mRNA_total),"Normalization",'probability','DisplayName','simulation');
        scatter(0.5+(0:3:length(h.BinEdges)-2),h.Values(1:3:end),'MarkerEdgeColor',[0.0392156876623631 0.658823549747467 0.250980406999588],...
    'LineWidth',1,...
    'MarkerFaceColor',[0.894117653369904 0.941176474094391 0.901960790157318]);
    end
    plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.mRNA_Prob,'LineWidth',1,'DisplayName','mixed')
%     plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.PDF.mRNA_Prob_v1,'LineWidth',1,'DisplayName','fast')
%     plot(0.5+(0:1:length(results.mRNA_Prob)-1),results.PDF.mRNA_Prob_v2,'LineWidth',1,'DisplayName','slow')
    set(figure1,'position',[300 400 280 190]);
    set(gca,'TickLength',[0.02,0.025]);
    legend1 = legend(gca,'show');
    results.weight
    box on
%     title(['kep ',num2str(K_coefficient(idx)),'  weight ',num2str(results.weight) ]);
end








