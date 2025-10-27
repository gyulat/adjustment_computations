function [M, e, w] = MFV_matrix(x, nit)
% MFV_MATRIX  Vectorized computation of Most Frequent Value (MFV) and dihesion for each row of a matrix
%
%   [M, e, w] = MFV_MATRIX(x, nit)
%
%   Inputs:
%       x   : data matrix, size (nRows x nCols)
%       nit : max. number of iterations (default = 100)
%
%   Outputs:
%       M   : column vector of MFV values per row (nRows x 1)
%       e   : column vector of dihesion values per row (nRows x 1)
%       w   : weight matrix (same size as x)

    if nargin < 2
        nit = 100;
    end

    % Initialize parameters
    xmax = max(x, [], 2);
    xmin = min(x, [], 2);
    e = sqrt(3.0)/2.0 .* (xmax - xmin);  % initial dihesion per row
    M = median(x, 2);                    % initial most frequent value per row

    [nRows, nCols] = size(x);

    % Expand for broadcasting across columns
    M = M(:); e = e(:);
    Mmat = repmat(M, 1, nCols);
    emat = repmat(e, 1, nCols);

    tol = 1e-6;
    for it = 1:nit
        % compute weights and update
        se = 1.0 ./ (emat.^2 + (x - Mmat).^2).^2;
        num = (x - Mmat).^2 .* se;
        den = se;

        ej1 = sqrt(3.0 * sum(num, 2) ./ sum(den, 2));  % updated dihesion per row
        ej1 = max(ej1, eps); % avoid divide by zero

        ej1mat = repmat(ej1, 1, nCols);
        w = 1.0 ./ (ej1mat.^2 + (x - Mmat).^2);
        Mj1 = sum(w .* x, 2) ./ sum(w, 2);

        % check convergence (vectorized)
        if all(abs(ej1 - e) ./ e < tol)
            break;
        end

        % update
        e = ej1;
        M = Mj1;
        Mmat = repmat(M, 1, nCols);
        emat = repmat(e, 1, nCols);
    end
end

