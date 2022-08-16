function t = BacktrackingLineSearch_BFGS(func, x_k, g_k, d_k)
% Refer: https://www.youtube.com/watch?v=5upFcYJqSwo,
% https://www.ams.jhu.edu/~abasu9/AMS_553-761/lecture05_handout2.pdf, https://gist.github.com/Heliosmaster/1043132
% Evaluates the step size t following BacktrackingLineSearch for BFGS algorithm with Strong Wolfe conditions for termination
% 
% Output arguments:
%        t : final stepsize
%
% Input arguments:
%       func: function handle that maps R^n to R
%       x_k : starting point where the line search is executed
%       g_k : gradient of func (column vector)
%       d_k : search direction (would be newton_step for Newton's method)/descent direction

% Initialization of default parameters
c1 = 1e-4;
c2 = 0.9; 
t0 = 0;      % initialize search step 
t1 = 1.0;
tmax = 5*t1;
h = 1e-4;
maxit = 100;
i = 1;
while true
    f_old = func(x_k + t0*d_k);
    x_kp1 = x_k + t1*d_k;
    %f_now = func(x_kp1);
    if func(x_kp1) > func(x_k) + c1*t1*g_k'*d_k || ((i>1) && func(x_kp1) > f_old) % sufficient decrease condition
        t = Zoom_StrongWolfe(func, x_k, g_k, d_k, t0, t1, c1, c2, h);
        return;
    end
    if abs(num_grad(func,h,x_kp1)*d_k) <= -c2*g_k'*d_k
        t = t1;
        return;
    end
    if num_grad(func,h,x_kp1)*d_k >= 0
        t = Zoom_StrongWolfe(func, x_k, g_k, d_k, t1, t0, c1, c2, h);
        return;
    end
    
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
