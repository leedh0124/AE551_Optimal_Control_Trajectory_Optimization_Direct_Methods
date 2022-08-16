function t_star = Zoom_StrongWolfe(func, x_k, g_k, d_k, t_lo, t_hi, c1, c2, h)
% Output arguments:
%       t_star  : final stepsize
%
% Input arguments:
%       fun  : function handle 
%       x_k  : starting point where the line search is executed
%       g_k : gradient of func (column vector)
%       d_k  : search direction
%       t_lo : lowest bound of interval
%       t_hi : highest bound of interval
%       c1   : constant c1
%       c2   : constant c2
%       h    : for numerical differentiation
maxit=20;
j=0;
while true
    t = (t_lo+t_hi)/2; % trial step length as mean of t_low and t_high
    x_kp1 = x_k + t*d_k;
    if func(x_kp1) > func(x_k) + c1*t*g_k'*d_k || (func(x_kp1) >= func(x_k + t_lo*d_k))
        t_hi = t; % narrows the interval between [t_low, t_high]
    else
        if abs(num_grad(func,h,x_kp1)*d_k) <= -c2*g_k'*d_k
            t_star = t;
            return;
        end
        if num_grad(func,h,x_kp1)*d_k * (t_hi-t_lo) >= 0
            t_hi = t_lo;
        end
        t_lo = t;
    end
    if j==maxit
        t_star = t;           % escape condition
        return;
    end
    
    j = j+1;
end




end