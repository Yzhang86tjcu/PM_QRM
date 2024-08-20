function  Beta0 = EstICPara(m, X, beta, disttype, distparams, Q)
%%% % Estimate IC parameters based on the generated m profile samples
ql = Q(1); qc = Q(2); qu = Q(3);
ysample = RandomSample(m, X, beta, disttype, distparams);   %% % Generate profile samples
Betae = zeros(m,6);
for i = 1:m
    qml = rq_fnm(X, ysample(:,i), ql); %% estimate the QRM parameters
    qm = rq_fnm(X, ysample(:,i), qc); %% estimate the QRM parameters
    qmu = rq_fnm(X, ysample(:,i), qu); %% estimate the QRM parameters
    Betae(i,:) = [qml', qm', qmu']; 
end
Beta0 = mean(Betae);
end