function ERSinv = EstSQ(m, X, beta, Q, hntype, disttype, distparams)
n = length(X);
Xm=mean(X);
RY = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 -1 0; 0 0 0 1 0 -1];  
ql = Q(1); qc = Q(2); qu = Q(3);
D = n^(-1)*(X'*X);  
Dinv = D^(-1);
Fql=icdf('norm', ql, 0, 1); pfql=pdf('norm', Fql, 0, 1);
Fqc=icdf('norm', qc, 0, 1); pfqc=pdf('norm', Fqc, 0, 1);
Fqu=icdf('norm', qu, 0, 1); pfqu=pdf('norm', Fqu, 0, 1);
alpha=0.05; Za=icdf('norm', 1-alpha/2, 0, 1); 
if hntype == 1
    hnl = n^(-0.2)*(4.5*pfql^4/(2*Fql^2+1)^2)^0.2;  
    hnc = n^(-0.2)*(4.5*pfqc^4/(2*Fqc^2+1)^2)^0.2; 
    hnu = n^(-0.2)*(4.5*pfqu^4/(2*Fqu^2+1)^2)^0.2; 
end

if hntype == 2
    hnl = n^(-1/3)*Za^(2/3)*(1.5*pfql^2/(2*Fql^2+1))^(1/3);
    hnc = n^(-1/3)*Za^(2/3)*(1.5*pfqc^2/(2*Fqc^2+1))^(1/3);
    hnu = n^(-1/3)*Za^(2/3)*(1.5*pfqu^2/(2*Fqu^2+1))^(1/3);
end

if hntype == 3
    hnl = Za*sqrt(ql*(1-ql)/n);
    hnc = Za*sqrt(qc*(1-qc)/n);
    hnu = Za*sqrt(qu*(1-qu)/n);
end

Sq = zeros(m,3);
ysample = RandomSample(m, X, beta, disttype, distparams);   %% % Generate profile samples
for i = 1:m
    qmlu = rq_fnm(X, ysample(:,i), ql+hnl);
    qmll = rq_fnm(X, ysample(:,i), ql-hnl);
    qmcu = rq_fnm(X, ysample(:,i), qc+hnc);
    qmcl = rq_fnm(X, ysample(:,i), qc-hnc);
    qmuu = rq_fnm(X, ysample(:,i), qu+hnu); 
    qmul = rq_fnm(X, ysample(:,i), qu-hnu); 
    sql = Xm*(qmlu-qmll)/(2*hnl);
    sqc = Xm*(qmcu-qmcl)/(2*hnc);
    squ = Xm*(qmuu-qmul)/(2*hnu);
    Sq(i,:) = [sql, sqc, squ];
end
SQ = mean(Sq);
ssql = SQ(1);
ssqc = SQ(2);
ssqu = SQ(3);

SIGMA = [qc*(1-qc)*ssqc^2*Dinv, qc*(1-qu)*ssqu*ssqc*Dinv, ql*(1-qc)*ssqc*ssql*Dinv;
qc*(1-qu)*ssqu*ssqc*Dinv, qu*(1-qu)*ssqu^2*Dinv, ql*(1-qu)*ssql*ssqu*Dinv; 
ql*(1-qc)*ssqc*ssql*Dinv, ql*(1-qu)*ssql*ssqu*Dinv, ql*(1-ql)*ssql^2*Dinv];

RS = RY*SIGMA*RY'; 
ERSinv = RS^(-1);
end