function [W, u, v, w] = unAdaSMCCA(Data, opts)

% Uncertainty-aware adaptive sparse multi-view canonical correlation analysis

% data
X = Data.X{1}; p = size(X, 2);
Y = Data.X{2}; d = size(Y, 2);
Z = Data.X{3}; q = size(Z, 2);

% pre-calculate
XX = X' * X; YY = Y' * Y; ZZ = Z' * Z;
XY = X' * Y; YZ = Y' * Z; XZ = X' * Z;

% set parameters
lambda_1 = opts.lambda_1;
beta = opts.beta;
lambda_2 = opts.lambda_2;
lambda_3 = opts.lambda_3;

% initialize u, v and w
u = ones(p, 1); u = u / norm(X * u);
v = ones(d, 1); v = v / norm(Y * v);
w = ones(q, 1); w = w / norm(Z * w);

% initialize loss weights
Xu = X * u; Yv = Y * v; Zw = Z * w;
omega12 = 1 / (norm(Xu - Yv) ^ 2);
omega23 = 1 / (norm(Yv - Zw) ^ 2);
omega13 = 1 / (norm(Xu - Zw) ^ 2);

% iteration
iter = 0; maxIter = 100;
tall = inf; tol = 1e-5;

while (iter < maxIter && tall > tol)
    iter = iter + 1;
    
    % update u
    u_old = u;
    % ----------------------------------------
    DFGL = updateD(u, 'FGL');
    D1 = updateD(u);
    % solve u
    F1 = (omega12 + omega13) * XX + lambda_1 * beta * DFGL + lambda_1 * (1 - beta) * D1;
    b1 = omega12 * XY * v + omega13 * XZ * w;
    u = F1 \ b1;
    % scale u
    u = u / norm(X * u);
    
    % update v
    v_old = v;
    % ----------------------------------------
    D2 = updateD(v);
    % solve v
    F2 = (omega12 + omega23) * YY + lambda_2 * D2;
    b2 = omega12 * XY' * u + omega23 * YZ * w;
    v = F2 \ b2;
    % scale v
    v = v / norm(Y * v);
    
    % update w
    w_old = w;
    % ----------------------------------------
    D3 = updateD(w);
    % solve w
    F3 = (omega13 + omega23) * ZZ + lambda_3 * D3;
    b3 = omega13 * XZ' * u + omega23 * YZ' * v;
    w = F3 \ b3;
    % scale w
    w = w / norm(Z * w);
    
    % update loss weights
    Xu = X * u; Yv = Y * v; Zw = Z * w;
    omega12 = 1 / (norm(Xu - Yv) ^ 2);
    omega23 = 1 / (norm(Yv - Zw) ^ 2);
    omega13 = 1 / (norm(Xu - Zw) ^ 2);
    
    % stopping condition
    % ----------------------------------------
    tu = max(abs(u - u_old));
    tv = max(abs(v - v_old));
    tw = max(abs(w - w_old));
    tall = max([tu, tv, tw]);
end

W{1} = u; W{2} = v; W{3} = w;
