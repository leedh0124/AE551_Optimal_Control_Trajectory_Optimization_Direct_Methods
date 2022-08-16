function t = BacktrackingLineSearch_Himm(func, x_k, g_k, newton_step)
% Evaluates the step size eta_k following BacktrackingLineSearch for Himmelblau function
% func: function that maps R^n to R
% x_k : starting point
% g_k : gradient of func (column vector)
% t : (minimizing) search step 

alpha = 0.3; % between 0 and 0.5 (rand*0.5)
beta = 0.6; % between 0 and 1 (rand)
t = 1;      % initialize search step 1

x_kp1 = x_k + t*newton_step; 
while func(x_kp1) > func(x_k) + alpha*t*g_k'*newton_step 
    t = beta*t;
    x_kp1 = x_k + t*newton_step; 
end
end
