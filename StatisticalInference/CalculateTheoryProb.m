function [mRNA_P_mix,mRNA_P_v1,mRNA_P_v2,weight] = CalculateTheoryProb(theta,data,params,idx)
K_NN     = theta(1);
K_EP     = theta(2);
k_on_max = theta(3);
k_on_min = theta(4);
k_off    = theta(5);
mu_max   = theta(6);
mu_min   = theta(7);
friction_coef = theta(8);
H     = params.H; 
delta    = params.delta;
if idx == 1
    vals = data.mRNAdistribin{1,idx};
    Mmax = floor(vals(end)*1.2);
    x = 0:1:Mmax;
    mRNA_P_mix = Poissbeta(k_on_min/delta,k_off/delta,mu_min/delta,x);
    weight = nan;
    mRNA_P_v1 = nan;
    mRNA_P_v2 = nan;
else
    Kaa = sqrt(4e-3*(K_NN/data.EPgenomedistance(idx) + K_EP)^(-1));
    d_EP_max = sqrt(2)*Kaa*4;
    d_EP = 0.001:0.001:d_EP_max;
    P_EP = sqrt(2./pi).*Kaa.^(-3).* d_EP.^2.*exp(- d_EP.^2./(2.*Kaa.^(2)));
    dt = d_EP(2)-d_EP(1);
    d_T = params.distance_T;
    d_05 = params.distance_05;
    Hill = (d_EP > d_T).* (1./(1+((d_EP-d_T)./(d_05-d_T)).^H));
    if length(Hill)>= d_T/dt
        Hill(1:d_T/dt-1) = 0;
    end
    %%
    lambda_on = (d_EP <= d_T).*(k_on_max./delta) + (d_EP > d_T).*...
        ((k_on_min./delta) + (k_on_max -k_on_min)./delta.*Hill);
    lambda_off = k_off/delta;
    lambda_mu = (d_EP <= d_T).*(mu_max./delta) + (d_EP > d_T).*...
        ((mu_min./delta) + (mu_max - mu_min)./delta.*Hill);
    kon = sum(lambda_on.*P_EP.*dt);
    koff = lambda_off;
    ksyn = sum(lambda_mu.*P_EP.*dt);
    
    vals = data.mRNAdistribin{1,idx};
    Mmax = floor(vals(end)*2);
    x = 0:1:Mmax;
    % fast
    mRNA_P_v1 = Poissbeta(kon,koff,ksyn,x);
    % slow
    mRNA_P = zeros(length(lambda_on),size(x,2));
    for i = 1:length(lambda_on)
        mRNA_P(i,:) = poissbeta(lambda_on(i),lambda_off,lambda_mu(i),x)';
    end
    mRNA_P_v2 = mRNA_P'*(P_EP.*dt)';

    %% mu
    b = params.distance_T;
    minlambda = min([k_on_min, mu_min]);
    maxvelocity = (K_NN/data.EPgenomedistance(idx) + K_EP)*...
        d_EP(find(cumsum(P_EP) > 99,1, 'first'))/(b*friction_coef);
    omega = minlambda/maxvelocity;
    weight = [1./(1+omega) omega./(1+omega)];
    mRNA_P_mix = 1./(1+omega).*mRNA_P_v1+omega./(1+omega).*mRNA_P_v2;
end
%% bins
% binnedmat = [];
if data.mRNAdistribinsize(idx) > 1
    binnedmat=kron(eye(floor(length(mRNA_P_mix)/data.mRNAdistribinsize(idx))),...
        ones(1,data.mRNAdistribinsize(idx)));
    if rem(length(mRNA_P_mix),data.mRNAdistribinsize(idx))>0
        binnedmat(end+1,end+1:end+rem(length(mRNA_P_mix),data.mRNAdistribinsize(idx)))=...
            ones(1,rem(length(mRNA_P_mix),data.mRNAdistribinsize(idx)));
    end
    mRNA_P_mix = binnedmat*mRNA_P_mix;
    mRNA_P_v1 = binnedmat*mRNA_P_v1;
    mRNA_P_v2 = binnedmat*mRNA_P_v2;
end
vals_prob = data.mRNAdistri{1, idx};
mRNA_P_mix = mRNA_P_mix(1:length(vals_prob));
if data.mRNAdistribinsize(idx) > 1
    mRNA_P_v1 = mRNA_P_v1(1:length(vals_prob));
    mRNA_P_v2 = mRNA_P_v2(1:length(vals_prob));
end

end


