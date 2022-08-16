function out = BM_objfunc(pk)
global C d
sz = size(pk);
n = length(pk)-1;
func = @himmelblaufunc;
if sz(2) == 1
    % pk = [x;t0]; % column, x0 is n-dimensional, while t0 is a scalar
    out = pk(n+1)*func(pk(1:n))-sum(log(d-C*pk(1:n))); % equiv to @(x,t) t*func(x) -sum(log(d-C*x));
else % pk contains perturbations for numerical computation of gradient, (n) by (n) since t = pk(n+1,1) is fixed
    out = pk(n+1,1).*func(pk(1:n,1:n))-sum(log(d-C*pk(1:n,1:n)),1); % since t is not a variable, only x is the variable
end
    
end
