function [M, e, w] = MFV1(x, nit)
% MFV  Iterative computation of Most Frequent Value and dihesion
%
%   [M, e, w] = MFV(x, nit)
%
%   Inputs:
%       x   : vector of data
%       nit : max. number of iterations (default = 100)
%
%   Outputs:
%       M   : most frequent value
%       e   : dihesion
%       w   : weight vector

    if nargin < 2
        nit = 100;
    end

    e0 = sqrt(3.0)/2.0 * (max(x) - min(x));
    M0 = median(x);

    ej = e0;
    Mj = M0;
    i = 0;

    while i < nit
        se = 1.0 ./ (ej^2 + (x - Mj).^2).^2;
        num = (x - Mj).^2 .* se;
        den = se;
        ej1 = sqrt(3.0 * sum(num) / sum(den));

        % w = ej1^2 ./ (ej1^2 + (x - Mj).^2);
        w = 1.0 ./ (ej1^2 + (x - Mj).^2);
        Mj1 = sum(w .* x) / sum(w);

        if abs(ej1 - ej) / ej < 1e-6
            break;
        end

        ej = ej1;
        Mj = Mj1;
        i = i + 1;
    end

    M = Mj1;
    e = ej1;
end

