%%%%%% **** Estimating the Steady-state IC AARL of EWMAQ chart with estimated IC parameters*******
clear
clc
%% %% Initial IC parameters
A = 3; B = 2;
beta = [B A];
ql = 0.25; qc = 0.5; qu = 0.75;   %%%% qth quantile
Q = [ql, qc, qu];
disttype = 1;  % 1 - normal; 2 - t; 3 - gamma;
distparams = [0,1]; % parameters of the error distribution
lamda = 0.05;
n = 30;
M = [30 50 100 300 500 1000 3000];
LQ = 11.3;
%% %% Initial the explanatory variables & Determine related parameters
x = zeros(n,1);
for i=1:n
    x(i) = 2+ (i-1)*(8-2)/n;
end
u=ones(n,1);
X=[x, u];
RSinv = SigComp(X,disttype, distparams, Q); %% Compute the SIGMA in Statistic
%% %% *** Compute the SS IC AARL & SDARL
ASDARL = [];
for rm = 1:length(M)
    m = M(rm);
    SA=100;
    rnm = 1;
    ARL = [];
    for rsa = 1:SA  
        Beta0 = EstICPara(m, X, beta, disttype, distparams, Q); %% % Estimate IC parameters based on m historical IC profile samples
        [ARL0,SDRL0]=SSICARLEwmaQ(X, beta, Beta0, RSinv, lamda, LQ, disttype, distparams, Q);  
        ARL = [ARL; ARL0];
    end
    AARL = mean(ARL);
    SDARL = std(ARL);
    ASDARL = [ASDARL; AARL SDARL];
end 
 
