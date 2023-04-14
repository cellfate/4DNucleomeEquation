function prob = Poissbeta(kon,koff,ksyn,x)
x = x(:);
[bp,wf] = gaussJacob(50,kon-1,koff-1);
A = 1/beta(kon,koff) * 2^(1-kon-koff);
p = exp( repmat(x,1,50).*log(repmat(ksyn*(1+bp')/2,length(x),1)+eps)...
    - repmat(gammaln(x+1),1,50)...
    - repmat(ksyn*(1+bp')/2,length(x),1) );
prob = sum(A*p*wf,2);
end


function [bp,wf] = gaussJacob(n,alph,bet)
%GaussJacob compute abscissas and weight factors for Gaussian quadratures
%CALL:  [bp,wf] = GaussJacob(n,alpha,beta)
%  bp = base points (abscissas)
%  n  = number of base points (abscissas) (integrates a (2n-1)th order polynomial exactly)
%  p(x) = (1-x)^alpha*(1+x)^beta  (alpha>-1, beta>-1)
a = zeros(n,1);
b = zeros(n-1,1);

ab = bet + alph;
abi = 2 + ab;
muzero = 2^(ab + 1) * beta(alph + 1, bet + 1);
a(1) = (alph - bet)/abi;
b(1) = sqrt(4*(1 + bet)*(1 + alph)/((abi + 1)*abi^2));
a2b2 = alph^2 - bet^2;

i = (2:n-1)';
abi = 2*i + ab;
a(i) = a2b2./((abi - 2).*abi);
a(n) =a2b2./((2*n - 2+ab).*(2*n+ab));
b(i) = sqrt(4*i.*(i + bet).*(i + alph).*(i + ab)./((abi.^2 - 1).*abi.^2));

[v,d] = eig( diag(a) + diag(b,1) + diag(b,-1) );
% wf = v(1,:);
wf = v(1,:)';

[bp,i] = sort( diag(d) );
wf = wf(i);

% bp = bp';
wf = muzero.* wf.^2;
end