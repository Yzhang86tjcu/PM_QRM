function  RSinv = SigComp(X,disttype, distparams, Q)
ql = Q(1); qc = Q(2); qu = Q(3);
RY = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 -1 0; 0 0 0 1 0 -1];  
n = length(X);
D = n^(-1)*(X'*X);  
Dinv = D^(-1);
% normal distribution
if disttype == 1
    dp1 = distparams(1); dp2 = distparams(2);
    [~,V]=normstat(dp1, dp2);
    Fql=icdf('norm', ql, dp1, dp2); pfql=pdf('norm', Fql, dp1, dp2); ssql = (sqrt(V)*pfql)^(-1); 
    Fqc=icdf('norm', qc, dp1, dp2); pfqc=pdf('norm', Fqc, dp1, dp2); ssqc = (sqrt(V)*pfqc)^(-1); 
    Fqu=icdf('norm', qu, dp1, dp2); pfqu=pdf('norm', Fqu, dp1, dp2); ssqu = (sqrt(V)*pfqu)^(-1); 
end

% t-distribution
if disttype == 2
    dp1 = distparams(1); 
    [~,V]=tstat(dp1);
    Fql=icdf('t', ql, dp1); pfql=pdf('t', Fql, dp1); ssql = (sqrt(V)*pfql)^(-1); 
    Fqc=icdf('t', qc, dp1); pfqc=pdf('t', Fqc, dp1); ssqc = (sqrt(V)*pfqc)^(-1); 
    Fqu=icdf('t', qu, dp1); pfqu=pdf('t', Fqu, dp1); ssqu = (sqrt(V)*pfqu)^(-1); 
end

% gamma distribution
if disttype == 3
    dp1 = distparams(1); dp2 = distparams(2);
    [~,V]=gamstat(dp1, dp2);
    Fql=icdf('gam', ql, dp1, dp2); pfql=pdf('gam', Fql, dp1, dp2); ssql = (sqrt(V)*pfql)^(-1); 
    Fqc=icdf('gam', qc, dp1, dp2); pfqc=pdf('gam', Fqc, dp1, dp2); ssqc = (sqrt(V)*pfqc)^(-1); 
    Fqu=icdf('gam', qu, dp1, dp2); pfqu=pdf('gam', Fqu, dp1, dp2); ssqu = (sqrt(V)*pfqu)^(-1); 
end

SIGMA = [qc*(1-qc)*ssqc^2*Dinv, qc*(1-qu)*ssqu*ssqc*Dinv, ql*(1-qc)*ssqc*ssql*Dinv;
    qc*(1-qu)*ssqu*ssqc*Dinv, qu*(1-qu)*ssqu^2*Dinv, ql*(1-qu)*ssql*ssqu*Dinv; 
    ql*(1-qc)*ssqc*ssql*Dinv, ql*(1-qu)*ssql*ssqu*Dinv, ql*(1-ql)*ssql^2*Dinv];
 
RS = RY*SIGMA*RY'; 
RSinv = RS^(-1);
end