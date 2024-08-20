function [ARL1,SDRL1]=SSOCARLEwmaQ(X, beta0, Beta0, Delta, RSinv, lamda, LQ, disttype, distparams,Q)
%%%%%%%%---------- Steady-state ARL of the EWMAQ chart
Betal0 = Beta0(:, 1:2)'; Betac0 = Beta0(:, 3:4)'; Betau0 = Beta0(:, 5:6)';
UCLR = LQ*lamda/(2-lamda); 
ql = Q(1); qc = Q(2); qu = Q(3);
RY = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 -1 0; 0 0 0 1 0 -1];  
n = length(X);
delta0 = Delta(1); delta1 = Delta(2); delta2 = Delta(3);
beta = [beta0(1)+delta1, beta0(2)+delta0];  
m = 1;
S=10000;
RL_1 = zeros(1,S);
for s=1:S
    t = 1;  %
    jump = -1;  %
    EWrt = [0; 0; 0; 0];
    while jump < 0   
%% %%%% Generate profile samples
        if t<=50
            ysample = RandomSample(m, X, beta0, disttype, distparams);
        else
            ysample = RandomOCSample(m, X, beta, delta2, disttype, distparams); 
        end
        qml = rq_fnm(X, ysample, ql); %% estimate the QRM parameters
        qmc = rq_fnm(X, ysample, qc); %% estimate the QRM parameters
        qmu = rq_fnm(X, ysample, qu); %% estimate the QRM parameters
%% %%%%%%%%%%%*************EWMAQ Charting Statistic*****************************
        YX = [qmc-Betac0; qmu-Betau0; qml-Betal0]; RX = RY*YX;
        Ewrt = EWrt;
        EWrt = lamda*sqrt(n)*RX + (1-lamda)*Ewrt;              
        ER = EWrt'*RSinv*EWrt;            
%% %%**** Identify whether the statistic exceeds the control limits *****
        if ER > UCLR
            if (1<=t && t<=50)
                EWrt = Ewrt;
                t = t - 1;
            else
                jump = 1;
            end
        end        
        t = t + 1;        
    end
    RL_1(s) = t - 50;
end 
ARL1 = mean(RL_1);
SDRL1 = std(RL_1);
end