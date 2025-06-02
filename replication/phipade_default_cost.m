function cost = phipade_default_cost(A, p)
% COST_PHIPADE Cost for the phipade algorithm in default degree deg = 7.
%
% COST returns the costs of phipade measured by the number of matrix
% multiplications.
%
% Reference:
% 
% [1] Bard Skaflestad and Will M. Wright. The scaling and modified squaring 
% method for matrix functions related to the exponential, Appl. Numer.
% MATH., 59 (2009), pp. 783--799.

if p <= 0 || mod(p, 1)~=0
    error('Index p must be a positive integer.')
end

normA = norm(A, inf);

theta = 9.50417899616293194e-01; % cf. [Tab. 4, sect. 4, 1] 
if normA<=theta % no scaling
    s = 0;
else 
    s = ceil( log2(normA/theta) );
end

m = 7; % the default

% number of matrix products needed for Pade approximants to phi_{0:p}
threshold = ceil(m/4) - 1;
if p >= threshold
    M = m - 1 + s*(p+1);
else
    M = floor(m/2) + 2*p + 1 + s*(p+1);
end

D = p + 1; % number of matrix division

cost = M + 4/3 * D;

end