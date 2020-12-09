function CCCs = calcCCC(Data, W)

% calculate Canonical Correlation Coefficient

Xu = Data.X{1} * W{1};
Yv = Data.X{2} * W{2};
Zw = Data.X{3} * W{3};

CCCs = [corr(Xu, Yv), corr(Yv, Zw), corr(Xu, Zw)];
