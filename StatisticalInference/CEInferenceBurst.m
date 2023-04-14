function CEInferenceBurst(params,data,theta0,lb,ub,iter)
A = [-1 1 0 0 0 0 0 0 
    0 0 -1 1 0 0 0 0 
    0 0 0 0 0 -1 1 0 
    0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0];
b = [0 0 0 0 0 0 0 0];
Aeq = [];
beq = [];
options = optimoptions('fmincon','Algorithm','interior-point','HessianApproximation','lbfgs','MaxFunctionEvaluations',6000);
[theta_est,fval,exitflag,output] = fmincon(@(theta) LogLikelihood(theta,data,params),theta0,A,b,Aeq,beq,lb,ub,[],options);
result_CE = [fval,exitflag,theta_est'];
filename = sprintf("result_CE//%d.mat",iter);
save(filename,"result_CE");
end

function LogLik = LogLikelihood(theta,data,params)
LogLik = 0;
for idx = 1:length(data.mRNAdistribinsize)
    [mRNA_P_mix,~] = CalculateTheoryProb(theta,data,params,idx);
    vals_prob = data.mRNAdistri{1, idx};
    lddata = (vals_prob) * log(mRNA_P_mix + eps);
    if isnan(lddata)==1 || isinf(lddata)==1
        LogLik = 10000;
    else
        LogLik = LogLik - sum(lddata);
    end
end
end




