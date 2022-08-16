function out = Project_BM_objfunc(pk)
global C d obj1 n_var n N tf
%pk = [xk;t0] (n+1) dimensional
sz = size(pk);
n = length(pk)-1;
%func = @himmelblaufunc;

if obj1
    func = @(x) sum(x(n_var:3:n,:).^2,1)*tf/(2*N); 
else
    func = @(x) (sum(x(n_var:3:n-n_var,:).^2,1) + sum(x(n_var:3:n-n_var,:).*x(2*n_var:3:n,:),1) + sum(x(2*n_var:3:n,:).^2,1) )*tf/(3*N)/2;
end
if sz(2) == 1
    % pk = [x;t0]; % column, x0 is n-dimensional, while t0 is a scalar
    x = pk(1:n); % column
    out = pk(n+1)*func(x)-sum(log(d-C*x)); % equiv to @(x,t) t*func(x) -sum(log(d-C*x));
else % pk contains perturbations for numerical computation of gradient, (n) by (n) since t = pk(n+1,1) is fixed
    x = pk(1:n,1:n); % matrix, row: each (perturbed) point is n-dimensional, column: perturbation along each dimension
    out = pk(n+1,1).*func(x)-sum(log(d-C*x),1); % since t is not a variable, only x is the variable
end
    
end
