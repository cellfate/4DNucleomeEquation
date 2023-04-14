enhancer_idx = [51,55,60,65,70,75,80,85,90,95];
promoter_idx = [49,45,40,35,30,25,20,15,10,5];

k_on_max       = [1,5,1]; % off to on max
k_on           = 0.06;  % off to on min
k_off          = 0.2; % on to off
mu_max         = 3; % generate mRNA
mu             = [1,1,0.5]; % generate mRNA
delta          = 0.1; % mRNA degradation
friction_coef  = 50; % friction coefficient
simulated_on = true;
variable_size = size(mu,2);
output_size = size(enhancer_idx,2);
mRNA_mean = zeros(variable_size,output_size);
mRNA_CV = zeros(variable_size,output_size);
mRNA_theor_mean = zeros(variable_size,output_size);
mRNA_theor_CV = zeros(variable_size,output_size);
weight_theor = zeros(variable_size,output_size);

%% loading data
for j = 1:1:variable_size
    for i = 1:1:output_size
        filename = sprintf('ResultDG//E_%d_P_%d_%f_%f_%f_%f_%f_%f_%f.mat',[enhancer_idx(i),...
            promoter_idx(i),k_on_max(j),k_on,k_off,mu_max,mu(j),delta,friction_coef]);
        load(filename);
        if results.params.simulated_on == true
            mRNA_mean(j,i) = results.mRNA_mean;
            mRNA_CV(j,i) = sqrt(results.mRNA_noise);
        end
        mRNA_theor_mean(j,i) = results.mRNA_theor_mean;
        mRNA_theor_CV(j,i) = results.mRNA_theor_CV;
        weight_theor(j,i) = results.weight(1);
    end
end

%% figure  % kep变化使用
%% mRNA_mean
figure1 = figure;
set(figure1,'position',[300 400 280 190]);
% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(enhancer_idx-promoter_idx,mRNA_theor_mean,'LineWidth',1,'Parent',axes1);
set(plot1(2),...
    'Color',[0.0392156876623631 0.658823549747467 0.250980406999588]);
set(plot1(1),...
    'Color',[0.952941179275513 0.149019613862038 0.258823543787003]);
set(plot1(3),'Color',[0 0.447058826684952 0.721568644046783]);
if simulated_on == true
    % 创建 scatter
    scatter(enhancer_idx-promoter_idx,mRNA_mean(2,:),...
        'MarkerEdgeColor',[0.0392156876623631 0.658823549747467 0.250980406999588],...
        'LineWidth',1,...
        'MarkerFaceColor',[0.894117653369904 0.941176474094391 0.901960790157318]);
    
    % 创建 scatter
    scatter(enhancer_idx-promoter_idx,mRNA_mean(1,:),...
        'MarkerEdgeColor',[0.952941179275513 0.149019613862038 0.258823543787003],...
        'LineWidth',1,...
        'MarkerFaceColor',[0.960784316062927 0.921568632125854 0.921568632125854]);
    
    % 创建 scatter
    scatter(enhancer_idx-promoter_idx,mRNA_mean(3,:),...
        'MarkerEdgeColor',[0 0.447058826684952 0.721568644046783],...
        'LineWidth',1,...
        'MarkerFaceColor',[0.87058824300766 0.921568632125854 0.980392158031464]);
end
box(axes1,'on');
hold(axes1,'off');
% % 设置其余坐标区属性
set(axes1,'TickLength',[0.02 0.025]);
% loglog plot
xlim(axes1,[1.5 110]);
% 取消以下行的注释以保留坐标区的 Y 范围
ylim(axes1,[6 30]);
box(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');

%% CV
figure1 = figure;
set(figure1,'position',[300 400 280 190]);
% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(enhancer_idx-promoter_idx,mRNA_theor_CV,'LineWidth',1,'Parent',axes1);
set(plot1(2),...
    'Color',[0.0392156876623631 0.658823549747467 0.250980406999588]);
set(plot1(1),...
    'Color',[0.952941179275513 0.149019613862038 0.258823543787003]);
set(plot1(3),'Color',[0 0.447058826684952 0.721568644046783]);
if simulated_on == true
    % 创建 scatter
    scatter(enhancer_idx-promoter_idx,mRNA_CV(2,:),...
        'MarkerEdgeColor',[0.0392156876623631 0.658823549747467 0.250980406999588],...
        'LineWidth',1,...
        'MarkerFaceColor',[0.894117653369904 0.941176474094391 0.901960790157318]);
    
    % 创建 scatter
    scatter(enhancer_idx-promoter_idx,mRNA_CV(1,:),...
        'MarkerEdgeColor',[0.952941179275513 0.149019613862038 0.258823543787003],...
        'LineWidth',1,...
        'MarkerFaceColor',[0.960784316062927 0.921568632125854 0.921568632125854]);
    
    % 创建 scatter
    scatter(enhancer_idx-promoter_idx,mRNA_CV(3,:),...
        'MarkerEdgeColor',[0 0.447058826684952 0.721568644046783],...
        'LineWidth',1,...
        'MarkerFaceColor',[0.87058824300766 0.921568632125854 0.980392158031464]);
end
box(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'TickLength',[0.02 0.025]);
% 对数
xlim(axes1,[1.5 110]);
% 取消以下行的注释以保留坐标区的 Y 范围
ylim(axes1,[6 30]);
box(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log');


