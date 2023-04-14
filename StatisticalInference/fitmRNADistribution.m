clear;clc;
 if isempty(gcp('nocreate'))
     parpool(24);
 end
%
load('dataToFit.mat')
params.simulatetime     = 1000;
params.delta            = 0.01;
params.distance_T       = 0.12;
params.distance_05      = 0.25;
params.H                = 3;

% data collection
data.EPgenomedistance   = [0 22 8 5 3 1];
lb = [1e-3;1e-3;1e-4;1e-4;1e-4;1e-3;1e-3;1];
ub = [1e1;1e1;1e0;1e0;5e0;1e2;1e2;2e2];
tic;
% initial point: KNN; KEP; alpha_max; alpha_min; beta; mu_max; mu_min; gamma
parfor iter = 1:params.simulatetime
    theta0 = log10(lb) + (log10(ub) - log10(lb)).*rand(length(lb),1);
    if theta0(1) < theta0(2)
        temp = theta0(1);
        theta0(1) = theta0(2);
        theta0(2) = temp;
    end
    if theta0(3) < theta0(4)
        temp = theta0(3);
        theta0(3) = theta0(4);
        theta0(4) = temp;
    end
    if theta0(6) < theta0(7)
        temp = theta0(6);
        theta0(6) = theta0(7);
        theta0(7) = temp;
    end
    theta0(1:end) = 10.^theta0(1:end)
    CEInferenceBurst(params,data,theta0,lb,ub,iter);
end
toc;

% collect data
result_all = zeros(params.simulatetime,length(lb)+2);
for idx = 1:params.simulatetime
    filename = sprintf("result_CE//%d.mat",idx);
    load(filename)
    result_all(idx,:) = result_CE;
end
result_all = sortrows(result_all,1);

% figure
theta_est = result_all(select_idx,3:end);
idx = 1;
figure1 = figure;
[mRNA_P_mix,mRNA_P_v1,mRNA_P_v2,weight] = calculateTheoryProb(theta_est,data,params,idx);
bar(data.mRNAdistribin{1,idx},data.mRNAdistri{1,idx});
hold on
plot(data.mRNAdistribin{1,idx},mRNA_P_mix,'DisplayName','mix','LineWidth',1);
set(figure1,'position',[300 400 170 190]);
set(gca,'TickLength',[0.02,0.025]);
xlim(gca,[-0.77 10]);

idx = 2;
figure1 = figure;
[mRNA_P_mix,mRNA_P_v1,mRNA_P_v2,weight] = calculateTheoryProb(theta_est,data,params,idx);
bar(data.mRNAdistribin{1,idx},data.mRNAdistri{1,idx});
hold on
plot(data.mRNAdistribin{1,idx},mRNA_P_mix,'DisplayName','mix','LineWidth',1);
set(figure1,'position',[300 400 170 190]);
set(gca,'TickLength',[0.02,0.025]);
set(gca,'TickLength',[0.02,0.025],'XTick',[5 55 105],'XTickLabel',...
    {'5','55','105'});
xlim(gca,[-5 115]); %2

idx = 3;
figure1 = figure;
[mRNA_P_mix,mRNA_P_v1,mRNA_P_v2,weight] = calculateTheoryProb(theta_est,data,params,idx);
bar(data.mRNAdistribin{1,idx},data.mRNAdistri{1,idx});
hold on
plot(data.mRNAdistribin{1,idx},mRNA_P_mix,'DisplayName','mix','LineWidth',1);
% plot(5:10:195,mRNA_P_mix,'DisplayName','mix','LineWidth',1);
set(figure1,'position',[300 400  170 190]);
set(gca,'TickLength',[0.02,0.025]);
set(gca,'TickLength',[0.02,0.025],'XTick',[5 55 105],'XTickLabel',...
    {'5','55','105'});
ylim(gca,[0 0.3]);% 3
xlim(gca,[-5 110]);% 3

idx = 4;
figure1 = figure;
[mRNA_P_mix,mRNA_P_v1,mRNA_P_v2,weight] = calculateTheoryProb(theta_est,data,params,idx);
bar(data.mRNAdistribin{1,idx},data.mRNAdistri{1,idx});
hold on
plot(data.mRNAdistribin{1,idx},mRNA_P_mix,'DisplayName','mix','LineWidth',1);
% plot(5:10:195,mRNA_P_mix,'DisplayName','mix','LineWidth',1);
set(figure1,'position',[300 400  170 190]);
set(gca,'TickLength',[0.02,0.025]);
set(gca,'TickLength',[0.02,0.025],'XTick',[5 55 105],'XTickLabel',...
    {'5','55','105'});
ylim(gca,[0 0.3]);% 3
xlim(gca,[-5 135]);% 4

idx = 5;
figure1 = figure;
[mRNA_P_mix,mRNA_P_v1,mRNA_P_v2,weight] = calculateTheoryProb(theta_est,data,params,idx);
bar(data.mRNAdistribin{1,idx},data.mRNAdistri{1,idx});
hold on
plot(data.mRNAdistribin{1,idx},mRNA_P_mix,'DisplayName','mix','LineWidth',1);
% plot(5:10:195,mRNA_P_mix,'DisplayName','mix','LineWidth',1);
set(figure1,'position',[300 400  170 190]);
set(gca,'TickLength',[0.02,0.025]);
set(gca,'TickLength',[0.02,0.025],'XTick',[5 55 105],'XTickLabel',...
    {'5','55','105'});
ylim(gca,[0 0.3]);% 3
xlim(gca,[-5 135]);% 5

idx = 6;
figure1 = figure;
[mRNA_P_mix,mRNA_P_v1,mRNA_P_v2,weight] = calculateTheoryProb(theta_est,data,params,idx);
bar(data.mRNAdistribin{1,idx},data.mRNAdistri{1,idx});
hold on
plot(data.mRNAdistribin{1,idx},mRNA_P_mix,'DisplayName','mix','LineWidth',1);
% plot(5:10:195,mRNA_P_mix,'DisplayName','mix','LineWidth',1);
set(figure1,'position',[300 400  170 190]);
set(gca,'TickLength',[0.02,0.025]);
set(gca,'TickLength',[0.02,0.025],'XTick',[5 55 105],'XTickLabel',...
    {'5','55','105'});
ylim(gca,[0 0.3]);% 3
xlim(gca,[-5 175]);% 6






