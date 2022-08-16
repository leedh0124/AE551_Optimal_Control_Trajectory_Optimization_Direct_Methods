function Hess_f = Hess_Him(X)
% computes (exact) Hessian for Himmelblau function at current X = [x,y] 
x = X(1); y = X(2);
Hess_f = [12*x^2+4*y-42, 4*y+4*x ; 4*y+4*x, 12*y^2+4*x-26];
end
