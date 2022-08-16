function [eta, interval, fmin] = GoldenSectionLineSearch(func, x0, s_k)
% Evaluates the step size eta_k following GoldenSectionLineSearch 
% func: function that maps R^n to R
% x0  : starting point
% s_k : search direction 
% eta : (minimizing) search step 

% 1. Determine Unimodal Interval from current point 
eps0 = 1e-6;           % initial search step size
x_k = x0;             % Column vector
f_k = func(x_k); 
x_kp1 = x0 + eps0*s_k; % initial perturbation in search direction
f_kp1 = func(x_kp1);

eps = eps0;
%assert(f_k>f_kp1,'df = f_kp1-f_k < 0 initially')
while f_k > f_kp1 % df = f_kp1-f_k < 0, which is guranteed for feasible search direction s_k
    eps = 2*eps;
    x_k = x_kp1;
    f_k = f_kp1;
    x_kp1 = x_k + eps*s_k;
    f_kp1 = func(x_kp1);
end
interval = [x0, x_kp1];

% 2. Perform Golden Section Line Search
a = x0; b = x_kp1;
GR = (sqrt(5)-1)/2;
tol = 1e-8;
err = inf;
while err > tol
    d = GR*(b-a);
    x1 = a + d; % x1 > x2
    x2 = b - d;
    if func(x1)< func(x2)
        a = x2;
    else
        % func(x1) > func(x2)
        b = x1;
    end
    % check domain
    if (b-a) < tol
        err = b-a;
        fmin = func(b);
    end
end

eta_temp = (b-x0)./s_k;
eta = mean(eta_temp);
end
