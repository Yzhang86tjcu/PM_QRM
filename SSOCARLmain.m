%%%%%% **** Estimating the Steady-state OC ARL of EWMAQ chart *******
clear
clc
%% %% Initial IC parameters
A = 3; B = 2;
beta0 = [B A];
ql = 0.25; qc = 0.5; qu = 0.75;   %%%% qth quantile
Q = [ql, qc, qu];
disttype = 1;  % 1 - normal; 2 - t; 3 - gamma;
distparams = [0,1]; % parameters of the error distribution
lamda = 0.05;
LQT = 11.3;
%% %% Initial the explanatory variables & Determine related parameters
n = 30;
x = zeros(n,1);
for i=1:n
    x(i) = 2+ (i-1)*(8-2)/n;
end
u=ones(n,1);
X=[x, u];
RSinv = SigComp(X,disttype, distparams, Q); %% Compute the SIGMA in Statistic
m = 5000;
Beta0 = EstICPara(m, X, beta0, disttype, distparams, Q); %% % Estimate IC parameters
[LQF,ARL0,SDRL0]=LQFSearch(X, beta0, Beta0, RSinv, lamda, LQT, disttype, distparams, Q);
%%  %%  ***  Intitial OC variables
Delta0 = [0 0.25 0.5 1 1.5 2 3 5];
Delta1 = [0.25 0.5 1 1.5 2 3 5];
Delta2 = [0.4 0.6 0.8 0.9 1.2 1.4 1.6 1.8];
R0=length(Delta0); R1=length(Delta1); R2=length(Delta2); R=R0+R1+R2;
%% %% **** Compute the SS OC ARL & SDRL
ARL = [];  % saving the SS OC ARL
 for r=1:R                     
    if r > R0+R1
        delta0 = 0; delta1 = 0; delta2 = Delta2(r-R0-R1);
    elseif r > R0
        delta0 = 0; delta1 = Delta1(r-R0)*sqrt((x'*x)^(-1)); delta2 = 1;
    else
        delta0 = Delta0(r)*sqrt(n^(-1)); delta1 = 0;  delta2 = 1;
    end     
    Delta = [delta0 delta1 delta2];  
    [ARL1,SDRL1]=SSOCARLEwmaQ(X, beta0, Beta0, Delta, RSinv, lamda, LQF, disttype, distparams,Q);
    ARL = [ARL; ARL1 SDRL1];
 end
 