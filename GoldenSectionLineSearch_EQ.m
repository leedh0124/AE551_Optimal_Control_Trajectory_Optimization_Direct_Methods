function [eta, interval] = GoldenSectionLineSearch_EQ(func, x0, s_k, n, m)
% Evaluates the step size eta_k following GoldenSectionLineSearch 
% func: function that maps R^n to R
% x0  : starting point [x_k;nu_k;g_k] / n,m,n
% s_k : search direction = -residual for equality constrained case 
% eta : (minimizing) search step 
% n   : dimension of x
% m   : number of equality constraints

% 1. Determine Unimodal Interval from current point 
objfunc = @himmelblaufunc;
h = 1e-4;
eps0 = 1e-6;           % initial search step size
x_k = x0(1:n+m);       % x_k now includes primal x_k and dual nu_k / n+m    
g_k = x0(n+m+1:n+m+n);
f_k = norm(func([x_k;g_k]),2); % norm for equality constrained case
x_kp1 = x_k + eps0*s_k; 
g_kp1 = num_grad(objfunc,h,x_kp1(1:n))';
f_kp1 = norm(func([x_kp1;g_kp1]),2);

eps = eps0;
%assert(f_k>f_kp1,'df = f_kp1-f_k < 0 initially')
while f_k > f_kp1 % df = f_kp1-f_k < 0, which is guranteed for feasible search direction s_k
    eps = 2*eps;
    x_k = x_kp1;
    f_k = f_kp1;
    x_kp1 = x_k + eps*s_k;
    g_kp1 = num_grad(objfunc,h,x_kp1(1:n))';
    f_kp1 = norm(func([x_kp1;g_kp1]),2);
end
interval = [x0(1:n+m), x_kp1];

% 2. Perform Golden Section Line Search
a = x0(1:n+m); b = x_kp1;
GR = (sqrt(5)-1)/2;
tol = 1e-8;
err = inf;
while err > tol
    d = GR*(b-a);
    x1 = a + d; % x1 > x2
    g1 = num_grad(objfunc,h,x1(1:n))';
    x2 = b - d;
    g2 = num_grad(objfunc,h,x2(1:n))';
    if norm(func([x1;g1]),2) < norm(func([x2;g2]),2)
        a = x2;
    else
        % func(x1) > func(x2)
        b = x1;
    end
    % check domain
    if norm(b-a,2) < tol
        err = norm(b-a,2);
    end
end

eta_temp = (b-x0(1:n+m))./s_k;
eta = mean(eta_temp);
end
