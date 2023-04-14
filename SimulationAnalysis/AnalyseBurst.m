function result = AnalyseBurst(params)
if params.simulated_on
    % This function analyzes data (numerical & theroetical) after parallel computing
    %% =====================numerical =======================================
    %% load data
    mRNA_level_total = [];
    for i = 1:1:params.simulation_num
        filename = sprintf('//%d.mat',i);
        temp = load([params.result_base_folder,params.filename,filename]);
        mRNA_level_total = [mRNA_level_total,temp.result.mRNA_level(1:end)];
    end
end
%% =====================theroteical =======================================
if params.EP_flag == false
    kon = params.k_on/params.delta; koff = params.k_off/params.delta; mu = params.mu/params.delta;
    % mRNA prob
    Mmax = 200;
    if params.simulated_on; Mmax = max(mRNA_level_total)+ 50; end
    x = 0:1:Mmax;

    mRNA_P = Poissbeta(kon,koff,mu,x)';
    mRNA_theor_mean = sum(mRNA_P.*x);
    mRNA_theor_var = sum(mRNA_P.*(x.^2)) - mRNA_theor_mean^2;
    mRNA_theor_noise = mRNA_theor_var/mRNA_theor_mean^2;
else
    Kaa = sqrt(params.diffusion_const*params.friction_coef*(params.spring_const(1,1)/...
        (abs(params.promoter_index - params.enhancer_index))+params.attraction_coef)^(-1));
    dt = 0.001;
    d_EP = 0.001:dt:5;
    P = sqrt(2./pi).*Kaa.^(-3).*d_EP.^2.*exp(-d_EP.^2./(2.*Kaa.^(2)));
    d_T = params.distance_T;
    d_05 = params.distance_05;
    Hill = (d_EP > d_T).* (1./(1+((d_EP-d_T)./(d_05-d_T)).^params.H));
    Hill(1:d_T/dt-1) = 0;
    %% V1
    lambda_on = (d_EP <= d_T).*(params.k_on_max./params.delta) + (d_EP > d_T).*...
        ((params.k_on./params.delta) + (params.k_on_max -params.k_on)./params.delta.*Hill);
    lambda_off = params.k_off/params.delta;
    lambda_mu = (d_EP <= d_T).*(params.mu_max./params.delta) + (d_EP > d_T).*...
        ((params.mu./params.delta) + (params.mu_max -params.mu)./params.delta.*Hill);
    on = sum(lambda_on.*P.*dt);
    off = sum(lambda_off.*P.*dt);
    mu = sum(lambda_mu.*P.*dt);
    % mRNA
    Mmax = 200;
    if params.simulated_on; Mmax = max(mRNA_level_total)+ 50; end
    x = 0:1:Mmax;
    mRNA_P_v1 = Poissbeta(on,off,mu,x)';
    mRNA_theor_v1 =  sum(mRNA_P_v1.*x);
    mRNA_theor_var_v1 = sum(mRNA_P_v1.*(x.^2)) - mRNA_theor_v1^2;
    %% V2
    % mRNA
    x = 0:1:Mmax;
    mRNA_P = zeros(length(lambda_on),size(x,2));
    for i = 1:length(lambda_on)
        mRNA_P(i,:) = Poissbeta(lambda_on(i),lambda_off,lambda_mu(i),x)';
    end
    mRNA_P_v2 = (mRNA_P'*(P.*dt)')';
    mRNA_theor_v2 = sum(mRNA_P_v2.*x);
    mRNA_theor_var_v2 = sum(mRNA_P_v2.*(x.^2)) - mRNA_theor_v2^2;
    %% mu
    minlambda = min([params.k_on,params.mu]);
    maxvelocity = (params.spring_const(1,1)/(abs(params.promoter_index - ...
        params.enhancer_index))+params.attraction_coef)*d_EP(find(cumsum(P)>...
        99.9,1, 'first'))/(params.b*params.friction_coef);
    omega = minlambda/maxvelocity;
    weight = 1./(1+omega);
    mRNA_P_mix = weight.*mRNA_P_v1+(1-weight).*mRNA_P_v2;
    mRNA_theor_mean = sum(mRNA_P_mix.*x);
    mRNA_theor_var = sum(mRNA_P_mix.*(x.^2)) - mRNA_theor_mean^2;
    mRNA_theor_noise = mRNA_theor_var/mRNA_theor_mean^2;
    mRNA_theor_CV = sqrt(mRNA_theor_noise);
end


%% ====================Saving data=========================
result.params = params;
if params.simulated_on
    result.mRNA_total = mRNA_level_total;
    result.mRNA_mean = mean(mRNA_level_total);
    result.mRNA_var = var(mRNA_level_total);
    result.mRNA_noise = var(mRNA_level_total)/result.mRNA_mean^2;
end
if params.EP_flag == false
    result.mRNA_Prob = mRNA_P;
    result.mRNA_theor = mRNA_theor_mean;
    result.mRNA_theor_var = mRNA_theor_var;
    result.mRNA_theor_noise = mRNA_theor_noise;
else
    result.weight = [1/(1+omega), omega/(1+omega)];
    result.omega = omega;
    result.PDF.mRNA_Prob_v1 = mRNA_P_v1;
    result.PDF.mRNA_Prob_v2 = mRNA_P_v2;
    result.mRNA_Prob = mRNA_P_mix;
    result.mRNA_theor_mean = mRNA_theor_mean;
    result.mRNA_theor_var = mRNA_theor_var;
    result.mRNA_theor_noise = mRNA_theor_noise;
    result.mRNA_theor_CV = mRNA_theor_CV;
end

end





