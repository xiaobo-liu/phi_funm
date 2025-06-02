function [m, cost] = select_deg_phipade(A, p)
% SELECT_DEG_PHIPADE Degree and the associated cost for the phipade 
% algorithm.
%
% The degree M takes a value from [3:13] if no scaling is needed; 
% otherwise, M = 13 is in use since the paper [1] appears to suggest so.
% COST retunrs the costs of phipade measured by the number of matrix
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


% try finding the minimal m such that norm(A,inf)<=theta_{m,0} (exponential)
th_vec = get_theta();
theta = min( th_vec(normA <= th_vec) );

if ~isempty(theta) % no scaling
    m = find(th_vec==theta, 1) + 2; 
    s = 0;
else 
    m = 13;
    s = ceil( log2(normA/th_vec(end)) );
end

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

%--------------------------------
% Subfunctions
%--------------------------------

function th = get_theta()
% theta_{m,0} for m=3:13

th = [
    1.495585217958292e-002  % m_vals = 3
    8.536352760102745e-002
    2.539398330063230e-001  
    5.414660951208968e-001
    9.504178996162932e-001  
    1.473163964234804e+000
    2.097847961257068e+000  
    2.811644121620263e+000
    3.602330066265032e+000
    4.458935413036850e+000
    5.371920351148152e+000];% m_vals = 13
end