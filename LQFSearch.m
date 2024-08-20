function  [LQF,ARL,SDRL]=LQFSearch(X, beta, Beta0, RSinv, lamda, LQ, disttype, distparams, Q) 
%%% %% Searching the final LR to meet the desired ARL0 
JumpLR = -1;
while JumpLR < 0
   %%% %% Compute the SSICARL & SSICSDRL 
    [ARL0,SDRL0]=SSICARLEwmaQ(X, beta, Beta0, RSinv, lamda, LQ, disttype, distparams, Q);       
    if abs(ARL0-200) < 2            
        LQF = LQ;
        ARL = ARL0;
        SDRL = SDRL0;
        JumpLR = 1;
    elseif abs(ARL0-200) < 10
        if ARL0 > 200
            LQ = LQ - 0.01;
        else
            LQ = LQ + 0.01;
        end
    else
        if ARL0 > 200
            LQ = LQ - 0.05;
        else
            LQ = LQ + 0.05;
        end
    end  
end