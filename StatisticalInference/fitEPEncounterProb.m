% clear;clc;
%% ===================================================================
load('dGvsEncounterProb.mat');
load('dataToFit.mat');
K_NN = result_all(1,3);
K_EP = result_all(1,4);
dG_r = distdata(13:22)./5; % 一个单体代表5kb
dG_l =  (abs(distdata(2:12)))./5;
dG = [eps, 3:3:60];
distance_T = params.distance_T ;
distance_05      = params.distance_05;
EncounterProb_slove = erf(distance_T./sqrt(2)./(sqrt(4e-3.*(K_NN./dG + K_EP).^(-1)))) ...
    - sqrt(2/pi)*(sqrt(4e-3.*(K_NN./dG + K_EP).^(-1))).^(-1).*...
    distance_T.*exp(-distance_T.^2./(2.*(sqrt(4e-3.*(K_NN./dG + K_EP).^(-1))).^2));
figure1 = figure;
hold on
scatter([dG_l;dG_r],cpdata(2:end),...
    'MarkerEdgeColor',[0 0.447058826684952 0.721568644046783],...
    'LineWidth',1,...
    'MarkerFaceColor',[0.87058824300766 0.921568632125854 0.980392158031464]);
hold on
plot([0, 3:3:60],[EncounterProb_slove])
set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);
xlim(gca,[-3 63])
box on
axis square

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter([dG_l;dG_r],cpdata(2:end),...
    'MarkerEdgeColor',[0 0.447058826684952 0.721568644046783],...
    'LineWidth',1,...
    'SizeData',20,...
    'MarkerFaceColor',[0.87058824300766 0.921568632125854 0.980392158031464]);
% 创建 loglog
loglog([0, 3:3:60],[EncounterProb_slove]);
set(figure1,'position',[300 400 280 190]);
xlim(axes1,[3 70]);
ylim(axes1,[0.04 0.4]);
box(axes1,'on');
axis(axes1,'square');
hold(axes1,'off');
set(axes1,'TickLength',[0.02 0.025],'XMinorTick','on','XScale','log',...
    'YMinorTick','on','YScale','log');

%% ========================================================================
theta_est = result_all(1,5:end);
dG_all = dG_l(1:end-1);
Contact_model_all = erf(distance_T./sqrt(2)./(sqrt(4e-3.*(K_NN./dG_all + K_EP).^(-1)))) ...
    - sqrt(2/pi)*(sqrt(4e-3.*(K_NN./dG_all + K_EP).^(-1))).^(-1).*...
    distance_T.*exp(-distance_T.^2./(2.*(sqrt(4e-3.*(K_NN./dG_all + K_EP).^(-1))).^2));
Contact_model_all = [0;Contact_model_all; 0.99];
Contact_data_all = [0,0.0563830638888889,0.0687399518518519,0.0507527666666667,0.0734207305555556,0.0653178592592592,0.0795144666666666,0.0883390481481481,0.101503166666667,0.139436488888889,0.355659407407407,0.999000000000000];
dG_all = [0; dG_all; 1];
mRNAmean_model_all = zeros(1,length(dG_all));
errormean_data_all = [0.542065665268867;3.35523654952613;4.03434668989636;2.56555648542236;3.77221404806007;4.77647391273676;3.10330094866007;5.21077383692754;3.39215283743600;5.94583727321590;8.04662237282550;8.72962229449797];
mRNAmean_data_all = [0.566265060240964;9.22397240708034;14.6240510267009;16.0414346402208;14.2027982266248;16.0956512493138;22.5781904102926;24.6041465929351;32.8068532004963;42.4832649585687;50.3095265562443;57.0435609925083];
CV_model_all = zeros(1,length(dG_all));
errorCV_data_all = [0.009566718449048;0.364573146778925;0.0981325356165627;0.0274426038321208;0.110308902996131;0.111446570833234;0.0647340680279242;0.0442243163273898;0;0.0156542272337047;0.0387517641419734;0.0214517021802704];
CV_data_all = [1.035271116564503;1.16885796531064;0.952398372675985;0.865395332507849;0.891855512171577;0.834733132075547;0.714191936844776;0.686095429331910;0.572992805733194;0.527635437094363;0.484752160611272;0.485488194725970];

for idx = 1:length(mRNAmean_model_all)
    k_on_max = theta_est(1);
    k_on_min = theta_est(2);
    k_off    = theta_est(3);
    mu_max   = theta_est(4);
    mu_min   = theta_est(5);
    friction_coef = theta_est(6);
    H        = 3;
    delta    = 0.01;
    b        = 0.12;
    
    if dG_all(idx) == 0
        Mmax = 100;
        x = 0:1:Mmax;
        mRNA_P_mix = poissbeta(k_on_min/delta,k_off/delta,mu_min/delta,x);
        mRNAmean_model_all(idx) = (x*mRNA_P_mix);
        mRNA_var = sum(mRNA_P_mix.*(x'.^2)) - mRNAmean_model_all(idx)^2;
        CV_model_all(idx) = sqrt(mRNA_var/(mRNAmean_model_all(idx)^2));
    else
        Kaa = sqrt(4e-3*(K_NN/dG_all(idx) + K_EP)^(-1));
        d_EP = 0.001:0.001:3;
        P_EP = sqrt(2./pi).*Kaa.^(-3).* d_EP.^2.*exp(- d_EP.^2./(2.*Kaa.^(2)));
        dt = d_EP(2)-d_EP(1);
        d_T = distance_T;
        d_05 = distance_05;
        Hill = (d_EP > d_T).* (1./(1+((d_EP-d_T)./(d_05-d_T)).^H));
        Hill(1:d_T/dt-1) = 0;
        %%
        lambda_on = (d_EP <= d_T).*(k_on_max./delta) + (d_EP > d_T).*...
            ((k_on_min./delta) + (k_on_max -k_on_min)./delta.*Hill);
        lambda_off = k_off/delta;
        lambda_mu = (d_EP <= d_T).*(mu_max./delta) + (d_EP > d_T).*...
            ((mu_min./delta) + (mu_max - mu_min)./delta.*Hill);
        kon = sum(lambda_on.*P_EP.*dt);
        koff = lambda_off;
        ksyn = sum(lambda_mu.*P_EP.*dt);
        
        Mmax = 500;
        x = 0:1:Mmax;
        % fast
        mRNA_P_v1 = poissbeta(kon,koff,ksyn,x);
        % slow
        mRNA_P = zeros(length(lambda_on),size(x,2));
        for i = 1:length(lambda_on)
            mRNA_P(i,:) = poissbeta(lambda_on(i),lambda_off,lambda_mu(i),x)';
        end
        mRNA_P_v2 = mRNA_P'*(P_EP.*dt)';
        
        %% mu
        minlambda = min([k_on_min, mu_min]);
        maxvelocity = (K_NN/dG(idx) +K_EP)*...
            d_EP(find(cumsum(P_EP) > 99.9,1, 'first'))/(b*friction_coef);
        omega = minlambda/maxvelocity;
        mRNA_P_mix = 1./(1+omega).*mRNA_P_v1+omega./(1+omega).*mRNA_P_v2;
        mRNAmean_model_all(idx) = (x*mRNA_P_mix);
        mRNA_var = sum(mRNA_P_mix.*(x'.^2)) - mRNAmean_model_all(idx)^2;
        CV_model_all(idx) = sqrt(mRNA_var/(mRNAmean_model_all(idx)^2));
        
    end
end


%% ==================================================================
theta_est = result_all(1,5:end);
dG = [22 8 5 3 1];
Contact_model = erf(distance_T./sqrt(2)./(sqrt(4e-3.*(K_NN./dG + K_EP).^(-1)))) ...
    - sqrt(2/pi)*(sqrt(4e-3.*(K_NN./dG + K_EP).^(-1))).^(-1).*...
    distance_T.*exp(-distance_T.^2./(2.*(sqrt(4e-3.*(K_NN./dG + K_EP).^(-1))).^2));
Contact_model = [0,Contact_model(1:end-1) 0.99];
Contact_data = [0,0.061506791736824,0.110107864565724,0.183087987508813,0.414545766095449,0.990000000000000];
dG = [0 22 8 5 3 1];
mRNAmean_model = zeros(1,length(dG));
mRNAmean_data = [0.566265060240964,20.7595,27.3038,43.660083160083160,51.1621,63.066666666666670];
CV_model = zeros(1,length(dG));
CV_data = [3.283607601333832,0.807528717929076,0.691406688488989,0.459397978717772,0.467006433151523,0.454659181709376];
for idx = 1:length(mRNAmean_model)
    k_on_max = theta_est(1);
    k_on_min = theta_est(2);
    k_off    = theta_est(3);
    mu_max   = theta_est(4);
    mu_min   = theta_est(5);
    friction_coef = theta_est(6);
    H        = 3;
    delta    = 0.01;
    b        = 0.12;
    
    if dG(idx) == 0
        Mmax = 100;
        x = 0:1:Mmax;
        mRNA_P_mix = poissbeta(k_on_min/delta,k_off/delta,mu_min/delta,x);
        mRNAmean_model(idx) = (x*mRNA_P_mix);
        mRNA_var = sum(mRNA_P_mix.*(x'.^2)) - mRNAmean_model(idx)^2;
        CV_model(idx) = sqrt(mRNA_var/(mRNAmean_model(idx)^2));
    else
        Kaa = sqrt(4e-3*(K_NN/dG(idx) + K_EP)^(-1));
        d_EP = 0.001:0.001:3;
        P_EP = sqrt(2./pi).*Kaa.^(-3).* d_EP.^2.*exp(- d_EP.^2./(2.*Kaa.^(2)));
        dt = d_EP(2)-d_EP(1);
        d_T = distance_T;
        d_05 = distance_05;
        Hill = (d_EP > d_T).* (1./(1+((d_EP-d_T)./(d_05-d_T)).^H));
        Hill(1:d_T/dt-1) = 0;
        %%
        lambda_on = (d_EP <= d_T).*(k_on_max./delta) + (d_EP > d_T).*...
            ((k_on_min./delta) + (k_on_max -k_on_min)./delta.*Hill);
        lambda_off = k_off/delta;
        lambda_mu = (d_EP <= d_T).*(mu_max./delta) + (d_EP > d_T).*...
            ((mu_min./delta) + (mu_max - mu_min)./delta.*Hill);
        kon = sum(lambda_on.*P_EP.*dt);
        koff = lambda_off;
        ksyn = sum(lambda_mu.*P_EP.*dt);
        
        Mmax = 500;
        x = 0:1:Mmax;
        % fast
        mRNA_P_v1 = poissbeta(kon,koff,ksyn,x);
        % slow
        mRNA_P = zeros(length(lambda_on),size(x,2));
        for i = 1:length(lambda_on)
            mRNA_P(i,:) = poissbeta(lambda_on(i),lambda_off,lambda_mu(i),x)';
        end
        mRNA_P_v2 = mRNA_P'*(P_EP.*dt)';
        
        %% mu
        minlambda = min([k_on_min, mu_min]);
        maxvelocity = (K_NN/dG(idx) +K_EP)*...
            d_EP(find(cumsum(P_EP) > 99.9,1, 'first'))/(b*friction_coef);
        omega = minlambda/maxvelocity;
        mRNA_P_mix = 1./(1+omega).*mRNA_P_v1+omega./(1+omega).*mRNA_P_v2;
        mRNAmean_model(idx) = (x*mRNA_P_mix);
        mRNA_var = sum(mRNA_P_mix.*(x'.^2)) - mRNAmean_model(idx)^2;
        CV_model(idx) = sqrt(mRNA_var/(mRNAmean_model(idx)^2));
    end
end
%% figure mean
figure1 = figure;
axes1 = axes('Parent',figure1);
errorbar(Contact_data_all,mRNAmean_data_all,errormean_data_all,'MarkerSize',3,...
    'MarkerFaceColor',[0.501960813999176 0.501960813999176 0.501960813999176],...
    'MarkerEdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0.501960813999176 0.501960813999176 0.501960813999176]);
hold on
% 创建 scatter
scatter(Contact_model_all,mRNAmean_model_all,'MarkerEdgeColor',[1 0 0],'LineWidth',0.2,...
    'MarkerFaceColor',[0.937254905700684 0.866666674613953 0.866666674613953]);

box(axes1,'on');
axis(axes1,'square');

% 设置其余坐标区属性
set(axes1,'TickLength',[0.02 0.025]);
set(figure1,'position',[300 400 280 190]);
axis square

scatter(Contact_data,mRNAmean_data,'Marker','^')
hold on
scatter(Contact_model,mRNAmean_model,'Marker','v')

set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);
box on
axis square
hold(axes1,'off');


%% figure CV
figure1 = figure;
axes1 = axes('Parent',figure1);
errorbar(Contact_data_all(2:end),CV_data_all(2:end),errorCV_data_all(2:end),'MarkerSize',3,...
    'MarkerFaceColor',[0.501960813999176 0.501960813999176 0.501960813999176],...
    'MarkerEdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0.501960813999176 0.501960813999176 0.501960813999176]);
hold on
% 创建 scatter
scatter(Contact_model_all(2:end),CV_model_all(2:end),'MarkerEdgeColor',[1 0 0],'LineWidth',0.2,...
    'MarkerFaceColor',[0.937254905700684 0.866666674613953 0.866666674613953]);

box(axes1,'on');
axis(axes1,'square');

% 设置其余坐标区属性
set(axes1,'TickLength',[0.02 0.025]);

set(figure1,'position',[300 400 280 190]);
axis square
scatter(Contact_data,CV_data,'Marker','^')
hold on
scatter(Contact_model,CV_model,'Marker','v')

set(figure1,'position',[300 400 280 190]);
set(gca,'TickLength',[0.02,0.025]);
box on
axis square
hold(axes1,'off');


% figure1 = figure;
% set(figure1,'position',[300 400 280 190]);
% % 创建 axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
%
% % 创建 scatter
% scatter(Contact_data(1),mRNAmean_data(1),'MarkerEdgeColor',[1 0 0],'LineWidth',0.7,...
%     'MarkerFaceColor',[0.937254905700684 0.866666674613953 0.866666674613953],...
%     'Marker','o');
%
% % 创建 scatter
% scatter(Contact_model(1),mRNAmean_model(1),'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
%     'LineWidth',0.7,...
%     'MarkerFaceColor',[0.803921580314636 0.878431379795074 0.968627452850342],...
%     'Marker','o');
%
% % 创建 scatter
% scatter(Contact_data(2),mRNAmean_data(2),'MarkerEdgeColor',[1 0 0],'LineWidth',0.7,...
%     'MarkerFaceColor',[0.937254905700684 0.866666674613953 0.866666674613953],...
%     'Marker','hexagram');
%
% % 创建 scatter
% scatter(Contact_model(2),mRNAmean_model(2),'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
%     'LineWidth',0.7,...
%     'MarkerFaceColor',[0.803921580314636 0.878431379795074 0.968627452850342],...
%     'Marker','hexagram');
%
% % 创建 scatter
% scatter(Contact_data(3),mRNAmean_data(3),'MarkerEdgeColor',[1 0 0],'LineWidth',0.7,...
%     'MarkerFaceColor',[0.937254905700684 0.866666674613953 0.866666674613953],...
%     'Marker','^');
%
% % 创建 scatter
% scatter(Contact_model(3),mRNAmean_model(3),'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
%     'LineWidth',0.7,...
%     'MarkerFaceColor',[0.803921580314636 0.878431379795074 0.968627452850342],...
%     'Marker','^');
%
% % 创建 scatter
% scatter(Contact_data(4),mRNAmean_data(4),'MarkerEdgeColor',[1 0 0],'LineWidth',0.7,...
%     'MarkerFaceColor',[0.937254905700684 0.866666674613953 0.866666674613953],...
%     'Marker','v');
%
% % 创建 scatter
% scatter(Contact_model(4),mRNAmean_model(4),'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
%     'LineWidth',0.7,...
%     'MarkerFaceColor',[0.803921580314636 0.878431379795074 0.968627452850342],...
%     'Marker','v');
%
% % 创建 scatter
% scatter(Contact_data(5),mRNAmean_data(5),'MarkerEdgeColor',[1 0 0],'LineWidth',0.7,...
%     'MarkerFaceColor',[0.937254905700684 0.866666674613953 0.866666674613953],...
%     'Marker','square');
%
% % 创建 scatter
% scatter(Contact_model(5),mRNAmean_model(5),'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
%     'LineWidth',0.7,...
%     'MarkerFaceColor',[0.803921580314636 0.878431379795074 0.968627452850342],...
%     'Marker','square');
%
% % 创建 scatter
% scatter(Contact_data(6),mRNAmean_data(6),'MarkerEdgeColor',[1 0 0],'LineWidth',0.7,...
%     'MarkerFaceColor',[0.937254905700684 0.866666674613953 0.866666674613953],...
%     'Marker','diamond');
%
% % 创建 scatter
% scatter(Contact_model(6),mRNAmean_model(6),'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
%     'LineWidth',0.7,...
%     'MarkerFaceColor',[0.803921580314636 0.878431379795074 0.968627452850342],...
%     'Marker','diamond');
%
% box(axes1,'on');
% axis(axes1,'square');
% hold(axes1,'off');
% % 设置其余坐标区属性
% set(axes1,'TickLength',[0.02 0.025]);
%
%
%
%
%
