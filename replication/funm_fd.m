function L = funm_fd(fun,A,E)
% FUNM_FD computes the Frechet derivative of the matrix function FUN at A
% in the direction E.

n = size(A,1);
X = fun( [A E;zeros(n) A]);
L = X(1:n,n+1:end);