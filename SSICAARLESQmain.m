%%%%%% **** Estimating the Steady-state IC ARL of EWMAQ chart while the SQ is estimated *******
clear
clc
%% %% Initial IC parameters
A = 3; B = 2;
beta = [B A];
ql = 0.25; qc = 0.5; qu = 0.75;   %%%% qth quantile
Q = [ql, qc, qu];
disttype = 1;  % 1 - normal; 2 - t; 3 - gamma;
distparams = [0,1]; % parameters of the error distribution
hntype = 2; % 1 - hn1; 2 - hn2; 3 - hn3;
lamda = 0.05;
M = [50 100 300 500 1000];
LQ = 11.54;
%% %% Initial the explanatory variables & Determine related parameters
n = 100;
x = zeros(n,1);
for i=1:n
    x(i) = 2+ (i-1)*(8-2)/n;
end
u=ones(n,1);
X=[x, u];
m0 = 5000;
Beta0 = EstICPara(m0, X, beta, disttype, distparams, Q); %% % Estimate IC parameters
%% %% *** Compute the SS IC AARL & SDARL
ASDARL = [];
for rm = 1:length(M)
    m = M(rm);
    SA=100;
    rnm = 1;
    ARL = [];
    for rsa = 1:SA  
        ERSinv = EstSQ(m, X, beta, Q, hntype, disttype, distparams); %% Estimate the SQ and the SIGMA in Statistic
        [ARL0,SDRL0]=SSICARLEwmaQ(X, beta, Beta0, ERSinv, lamda, LQ, disttype, distparams, Q);  
        ARL = [ARL; ARL0];
    end
    AARL = mean(ARL);
    SDARL = std(ARL);
    ASDARL = [ASDARL; AARL SDARL];
end 
 
