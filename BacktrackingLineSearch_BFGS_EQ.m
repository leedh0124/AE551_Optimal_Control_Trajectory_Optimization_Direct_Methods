function t = BacktrackingLineSearch_BFGS_EQ(func, x0, d_k, n, m)
% Refer: https://www.youtube.com/watch?v=5upFcYJqSwo,
% https://www.ams.jhu.edu/~abasu9/AMS_553-761/lecture05_handout2.pdf
% Evaluates the step size t following BacktrackingLineSearch on norm(residual,2) with Strong Wolfe conditions for termination
% used for Equality Constrained Problem 
%
% Output arguments:
%        t : final stepsize
%
% Input arguments:
%       func: function handle that maps R^n to R = norm(residual,2) function
%       x0 : starting point [x_k;nu_k;g_k] / n,m,n
%       g_k : gradient of func (column vector)
%       d_k  : -residual (n+m dimensional)
%         n : dimension of x
%         m : number of equality constraints

% Initialization of default parameters
global obj1 n_var n N tf
if obj1
    objfunc = @(x) sum(x(n_var:3:n,:).^2,1)*tf/(2*N); % x can be column or matrix
else
    objfunc = @(x) (sum(x(n_var:3:n-n_var,:).^2,1) + sum(x(n_var:3:n-n_var,:).*x(2*n_var:3:n,:),1) + sum(x(2*n_var:3:n,:).^2,1))*tf/(3*N)/2;  
end

%objfunc = @himmelblaufunc;
c1 = 1e-4;
c2 = 0.9; 
t0 = 0;      
t1 = 1.0;
tmax = 5*t1;
h = 1e-4;
maxit = 100;
i = 1;
x_k = x0(1:n+m); % x_k now includes primal x_k and dual nu_k / n+m
g_k = x0(n+m+1:n+m+n);

while true
    f_old = norm(func([x_k + t0*d_k; g_k]),2);
    x_kp1 = x_k + t1*d_k;
    g_kp1 = num_grad(objfunc,h,x_kp1(1:n))'; % g_kp1 is n-dimensional
    if norm(func([x_kp1;g_kp1]),2) > norm(func([x_k;g_k]),2) - c1*t1*norm(d_k,2) || ((i>1) && norm(func([x_kp1;g_kp1]),2) > f_old) % sufficient decrease condition
        t = Zoom_StrongWolfe_EQ(objfunc, func, x_k, g_k, d_k, t0, t1, c1, c2, h, n, m);
        return;
    end
    %if abs(num_grad(func,h,[x_kp1;g_kp1])*d_k) <= c2*norm(d_k,2)
    d_kp1 =  -func([x_kp1;g_kp1]); % residual at x_kp1
    if norm(d_kp1,2) <= c2*norm(d_k,2)
        t = t1;
        return;
    end
    %if num_grad(func,h,x_kp1)*d_k >= 0
    %    t = Zoom_StrongWolfe(func, x_k, g_k, d_k, t1, t0, c1, c2, h);
    %    return;
    %end
    
    if i == maxit
        disp('Maximum number of iteration for Line Search reached');
        t = t1;
        return;
    end
    
    % Update for next loop
    i = i+1;
    t0 = t1;
    t1 = 0.8*t0 + 0.2*tmax;
    
end
end
