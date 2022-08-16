function [x0, s] = Phase1(A,b)
%% Phase I: Generate strictly feasible solution via Newton's Method
%  Input  : A, b
%  Output : x0, s
[m,n] = size(A);

% Reformulate A into A1 and set c1, x = [x0 ; s] (for Phase 1 problem)
A1 = [A -ones(m,1)];
c1 = [zeros(n,1);1];
x = [zeros(n,1);-min(b)+1]; % s is initially positive

% Parameters for Barrier method. Note: objective function for Phase 1 f(x) = ts + phi(x)
t = 1;
mu = 10;
epsB = 1e-4;
x0 = x(1:end-1);
s = 0;
max_iter = 10^2;

% Newton's method parameters
tol = 1e-8;
alpha = 0.1; %between 0.0 and 0.5
beta = 0.5; %between 0.0 and 1.0

%disp('Initiate basic Phase 1 method');
while 1
    for ii=1:max_iter
        if x(end) < 0
            break; % s<0 is obtained 
        end
        % Step 1: Compute the Newton step and decrement
        grad_f = gradf(A1,b,x,t,c1);
        Hessian_f = Hessf(A1,b,x);
        del_x_nt  = -linsolve(Hessian_f,grad_f);
        nt_dec = -grad_f'*del_x_nt;
        
        % Step 2: Stopping criterion
        error = nt_dec/2;
        if error <= tol
            break;
        end
        
        % Step 3: Choose step size t by backtracking line search
        t_BT = 1;
        while myfunc1(A1,b,x + t_BT*del_x_nt,t,c1) > myfunc1(A1,b,x,t,c1) + alpha*t_BT*grad_f'*del_x_nt
            t_BT = beta*t_BT;
        end
        
        % Step 4: Update x
        x = x + t_BT*del_x_nt;
    end
    
    % Check for s<0
    if (m/t < epsB) || (x(end) < 0)
        x0 = x(1:end-1);
        s = x(end);
        break;
    else
        t = mu*t;
    end
end

% Display minimization result
%disp('basic Phase 1 method complete')
%disp('s is')
%s
end


