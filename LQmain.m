%%%%%% **** Finding LQ in the UCL of EWMAQ chart *******
clear
clc
%% %% Initial IC parameters ****
A = 3; B = 2;
beta = [B A];
ql = 0.25; qc = 0.5; qu = 0.75;   %%%% qth quantile
Q = [ql, qc, qu];
disttype = 1;  % 1 - normal; 2 - t; 3 - gamma;
distparams = [0,1]; % parameters of the error distribution
lamda = 0.2;
%% %%%%*** Finding final LQ for different n *****
N = [30 50 100 300 500 1000];
LQT = 13.8; %% Initial LQ
ARL = [];
LQF = [];
for rn = 1:length(N)       
    n =  N(rn);    
    x = zeros(n,1);
    for i=1:n
        x(i) = 2+ (i-1)*(8-2)/n;
    end
    u=ones(n,1);
    X=[x, u];
    m = 5000;
    Beta0 = EstICPara(m, X, beta, disttype, distparams, Q); %% % Estimate IC parameters
    RSinv = SigComp(X,disttype, distparams, Q); %% Compute the SIGMA in Statistic
    [LQ,ARL0,SDRL0]=LQFSearch(X, beta, Beta0, RSinv, lamda, LQT, disttype, distparams, Q);
    ARL = [ARL; ARL0 SDRL0];
    LQF = [LQF LQ];
end
