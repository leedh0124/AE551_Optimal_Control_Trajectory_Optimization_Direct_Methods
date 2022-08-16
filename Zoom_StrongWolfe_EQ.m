function t_star = Zoom_StrongWolfe_EQ(objfunc, func, x_k, g_k, d_k, t_lo, t_hi, c1, c2, h, n, m)
% Output arguments:
%       t_star  : final stepsize
%
% Input arguments:
%       fun  : function handle 
%       x_k  : starting point [x_k;nu_k;g_k] / n,m,n
%       d_k  : -residual (n+m dimensional)
%       t_lo : lowest bound of interval
%       t_hi : highest bound of interval
%       c1   : constant c1
%       c2   : constant c2
%       h    : for numerical differentiation
%         n : dimension of x
%         m : number of equality constraints
maxit=20;
j=0;

while true
    t = (t_lo+t_hi)/2; % trial step length as mean of t_low and t_high
    x_kp1 = x_k + t*d_k;
    g_kp1 = num_grad(objfunc,h,x_kp1(1:n))';
    x_klow = x_k + t_lo*d_k;
    g_klow = num_grad(objfunc,h,x_klow(1:n))';
    if norm(func([x_kp1;g_kp1]),2) > norm(func([x_k;g_k]),2) - c1*t*norm(d_k,2) || (norm(func([x_kp1;g_kp1]),2) >= norm(func([x_klow;g_klow]),2))
        t_hi = t; % narrows the interval between [t_low, t_high]
    else
        d_klow = -func([x_klow;g_klow]); % residual at x_klow
        if norm(d_klow,2) <= c2*norm(d_k,2)
            t_star = t;
            return;
        end
        %if num_grad(func,h,x_kp1)*d_k * (t_hi-t_lo) >= 0
        %    t_hi = t_lo;
        %end
        t_lo = t;
    end
    if j==maxit
        t_star = t;           % escape condition
        if t_star < 1e-5
            t_star = 0.5;
        end
        return;
    end
    
    j = j+1;
end




end