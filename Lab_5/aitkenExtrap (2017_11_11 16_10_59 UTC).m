% Given a fixed point iterand parentIter, extrapolate nrExtrap times yielding a
% higher order fixed point iteration.
% INPUT
% parentIter    parent iteration function
% tol           desired tolerance (update based error estimate)
% maxIt         maximum number of iterations
% nrExtrap      number of extrapolations (zero means we just use parentIter)
% OUTPUT
% alpha         approximate root of f
% flag          0: attained tolerance, 1: maxIt reached
% convHist      convergence history
function [alpha, flag, convHist] = aitkenExtrap(parentIter, x0, tol, maxIt,...
    nrExtrap)
    
    % Check the input
    if (nargin < 5), nrExtrap = 1; end
    if (nargin < 4), maxIt = 15; end
    if (nargin < 3), tol = 1E-8; end
    if (nargin < 2), error('Aitken:BadInput', 'Not enough input arguments.'); end

    % Initialise variables
    flag = 1;
    convHist = zeros(maxIt, 1);
    x = x0;

    % Main loop
    for iter = 1 : maxIt

        % Calculate update
        xNew = extrap(parentIter, x, nrExtrap);

        % Update solution
        update = xNew - x;
        x = xNew;

        % Compute error estimate
        convHist(iter) = abs(update);

        % Check convergence
        if (convHist(iter) < tol)
            flag = 0;
            break;
        end
    end

    % Clean up variables
    convHist(iter + 1 : end) = [];
    alpha = x;

end

% Recursively performs extrapolation untill depth = nrExtrap
function xNew = extrap(parentIter, x, depth)

    if (depth == 0)                         % Call parent function
        xNew = parentIter(x);
    else                                    % Extrapolate and increase depth
        depth = depth - 1;

        % Combine two recursive calls
        phix = extrap(parentIter, x, depth);
        phiphix = extrap(parentIter, phix, depth);

        denominator = phiphix - 2 * phix + x;

        % Update of solution
        if (abs(denominator) == 0)          % Prevent division by zero
            xNew = x;
        else
            xNew = x - ((phix - x) ^ 2) / denominator;
        end
    end

end