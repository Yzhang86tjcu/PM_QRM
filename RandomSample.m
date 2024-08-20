function ysample = RandomSample(m, x, beta, type, params)
% generating random profile samples specified by the type and parameters
% m is the number of samples required
% x is the fixed design points with n points in a n-by-p regressor matrix
% type is the distribution type
% params is the cell array containing corresponding paramemters
% output: sample is the matrix n-by-m

[n,p] = size(x);

% normal distribution
if type == 1
    sample1 = normrnd(params(1), params(2), n, m);
    [MNx,Vx]=normstat(params(1), params(2));
    sample = (sample1-MNx)./sqrt(Vx);
end

% t-distribution
if type == 2
    sample1 = trnd(params(1), n, m); 
%     Fqx=icdf('t', 0.5, params(1));  
    [MNx,Vx]=tstat(params(1));
    sample = (sample1-MNx)./sqrt(Vx);
end

% gamma distribution
if type == 3
    sample1 = gamrnd(params(1), params(2), n, m);
    [MNx,Vx]=gamstat(params(1), params(2));
    sample = (sample1-MNx)./sqrt(Vx);
end

ysample = x*beta' + sample;
end

