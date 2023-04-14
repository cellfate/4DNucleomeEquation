% clear;clc;
% close all;
%% current parallel pool
format long
K_coefficient =  10.^(-2.5:0.02:-0);
dG = 1:4:81;
for idx =  38
    TheoryCalculateDistribution(idx,K_coefficient,dG,true)
end
for idx = 38
    TheoryCalculateDistribution(idx,K_coefficient,dG,false)
end

peak_num_all = zeros(length(dG),1);
peak_index = zeros(length(dG),3);
peak_prob_value = zeros(length(dG),3);
for idx =  38
    load(sprintf('PeakNum//peak_num_%d.mat',idx));
    peak_num_all = peak_numbers;
    load(sprintf('PeakInformation_kep//peak_infor_%d.mat',idx));
    for jdx = 1:length(dG)
        peak_index(jdx,1:peak_numbers(jdx)) = peakInfor.dg(jdx).peak_index;
        peak_prob_value(jdx,1:peak_numbers(jdx)) = peakInfor.dg(jdx).peak_prob_value;
    end
end




figure1 = figure;
X = [dG,dG,dG]';
Y = [peak_index(:,1); peak_index(:,2); peak_index(:,3)];
Z = [peak_prob_value(:,1); peak_prob_value(:,2); peak_prob_value(:,3)];
scatter(X,Y,[],Z,'LineWidth',0.7,'Marker','^')
hold on
patch(X(1:length(dG)),Y(1:length(dG)),Z(1:length(dG)),'EdgeColor','flat','MarkerFaceColor','flat','FaceColor','none')
patch(X(length(dG) + 1:2*length(dG)),Y(length(dG) + 1:2*length(dG)),Z(length(dG) + 1:2*length(dG)),'EdgeColor','flat','MarkerFaceColor','flat','FaceColor','none')
patch(X(2*length(dG) + 1:3*length(dG)),Y(2*length(dG) + 1:3*length(dG)),Z(2*length(dG) + 1:3*length(dG)),'EdgeColor','flat','MarkerFaceColor','flat','FaceColor','none')
set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);
box on


figure1 = figure;
set(figure1,'position',[300 400 280 190]);
imagesc(dG,log10(K_coefficient),BC_all)
set(gca,'TickLength',[0.02,0.025]);
axis square
axis xy
hold on
% scatter(dG(points_y),log10(K_coefficient(points_x)))

figure1 = figure;
set(figure1,'position',[300 400 280 190]);
imagesc(peak_num_all)
set(gca,'TickLength',[0.02,0.025]);
axis square
axis xy



function TheoryCalculateDistribution(idx,K_coefficient,dG,flag)
BC = zeros(1,length(dG));
peak_numbers = zeros(1,length(dG));
fig = figure('Visible','off','units','normalized','position',[0.1,0.1,0.9,0.75]);
for jdx = 1:length(dG)
    %% Parameters setting
    % Upstream: Chromatin conformation
    input_options.simulated_on      = false;
    input_options.EP_flag           = true;
    input_options.attraction_coef   = K_coefficient(idx); % K, E-P communication coefficient
    input_options.enhancer_index    = 1; % enhancer index number
    input_options.promoter_index    = input_options.enhancer_index + dG(jdx); % promoter index number
    
    % Downstream: Transcriptional bursting
    input_options.k_on_max       = 2; % off to on max
    input_options.k_on           = 0.06;  % off to on min
    input_options.k_off          = 0.1; % on to off
    input_options.mu_max         = 9; % generate mRNA
    input_options.mu             = 3; % generate mRNA
    input_options.delta          = 0.1; % mRNA degradation
    input_options.friction_coef  = 50; % friction coefficient
    
    
    input_options.result_base_folder = fullfile(pwd, 'Result');
    filename = sprintf('//%f_%d_%f_%f_%f_%f_%f_%f_%f.mat',[input_options.attraction_coef,dG(jdx),...
        input_options.k_on_max,input_options.k_on,input_options.k_off,input_options.mu_max...
        input_options.mu,input_options.delta,input_options.friction_coef]);
    if flag == true
        params = ParametersBurst(input_options);
        results = AnalyseBurst(params);
        save([input_options.result_base_folder,filename],'results');
    else
        load([input_options.result_base_folder,filename]);
        index = find(cumsum(results.mRNA_Prob) > 0.999,1,'first');
        mRNA_Prob = results.mRNA_Prob(1:index);
        Prob_normal = mRNA_Prob(2:end-1);
        Prob_forward = mRNA_Prob(1:end-2);
        Prob_backward = mRNA_Prob(3:end);
        % find the extremun points
        extremum_flag = ((Prob_forward - Prob_normal).*(Prob_normal - Prob_backward) < 0);
        extremum_index = find(extremum_flag ==1) + 1;
        extremum_value = Prob_normal(extremum_flag);
        if mRNA_Prob(2) - mRNA_Prob(1) < 0
            extremum_value = [mRNA_Prob(1), extremum_value];
            extremum_index = [1, extremum_index];
        end
        %
        if length(extremum_value) == 1
            peak_num = 1;
            peak_index = extremum_index;
            peak_prob_value = mRNA_Prob(peak_index);
        else
            extremum_value_forward = [0,0,extremum_value];
            extremum_value_normal = [0,extremum_value,0];
            extremum_value_backward = [extremum_value,0,0];
            peak_flag = (extremum_value_forward - extremum_value_normal) < 0 & (extremum_value_forward - extremum_value_normal) < -2e-4 ...
                & (extremum_value_normal - extremum_value_backward) > 0 & (extremum_value_normal - extremum_value_backward) > 2e-4;
            peak_num = sum(peak_flag);
            temp =  find(peak_flag ==1) - 1;
            peak_index = extremum_index(temp);
            peak_prob_value = mRNA_Prob(peak_index);
        end
        peak_numbers(1,jdx) = peak_num;
        BC(1,jdx) = results.BC;
        peakInfor.dg(jdx).dgindex = jdx;
        peakInfor.dg(jdx).peak_index = peak_index;
        peakInfor.dg(jdx).peak_prob_value = peak_prob_value;
    end
end
if flag == false
    save(sprintf('PeakNum//peak_num_%d.mat',idx),"peak_numbers");
    save(sprintf('PeakInformation_kep//peak_infor_%d.mat',idx),"peakInfor")
end

end


