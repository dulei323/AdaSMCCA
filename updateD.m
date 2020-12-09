function D = updateD(w, type)

if nargin == 1
    % for L1-norm
    d = 1 ./ sqrt(w .^ 2 + eps);
elseif strcmpi(type, 'FGL')
    % for FGL-norm
    [n_features, ~] = size(w);
    structure = updateGraph(n_features, 'FGL');
    Gp = 1 ./ sqrt(structure * (w .^ 2) + eps);
    d = [Gp(1), sum(reshape(Gp(2 : end - 1), 2, [])), Gp(end)];
elseif strcmpi(type, 'GGL')
    % for GGL-norm
    [n_features, ~] = size(w);
    structure = updateGraph(n_features, 'GGL');
    Gp = 1 ./ sqrt(structure * (w .^ 2) + eps);
    d = sum(reshape(Gp, n_features - 1, []));
else
    error('Error type.');
end

D = diag(d);

function E = updateGraph(n, type)

if strcmpi(type, 'FGL')
    E = zeros(2 * (n - 1), n);
    for i = 1 : n - 1
        j = i + 1;
        E(2 * i - 1, i) = 1;
        E(2 * i - 1, j) = 1;
        E(2 * i, i) = 1;
        E(2 * i, j) = 1;
    end
elseif strcmpi(type, 'GGL')
    num = 0;
    E = zeros(n * (n - 1), n);
    for i = 1 : n
        for j = 1 : n
            if i ~= j
                num = num + 1;
                E(num, i) = 1;
                E(num, j) = 1;
            end
        end
    end
else
    error('Error type.');
end
