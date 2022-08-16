function del_f = grad_Him(X)
% computes (exact)gradient for Himmelblau function at current X = [x,y] 
x = X(1); y = X(2);
del_f = [4*x^3+4*x*y-42*x+2*y^2-14;4*y^3+4*x*y-26*y+2*x^2-22];
end
