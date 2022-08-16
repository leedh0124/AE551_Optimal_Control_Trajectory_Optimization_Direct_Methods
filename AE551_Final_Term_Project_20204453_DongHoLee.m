%% AE551 Final Term Project: Solving Optimal Control Problem via Customized Optimization Solver 
% Method: Barrier Method 
% Collocation: Hermite-Simpson 

close all; clear all; 
global A b C d obj1 n_var n N tf 
g = 9.81; h = 1e-4; conv_tol = 1e-5; BM_tol = 1e-6;
warning('off','all');
%% Define Parametrized Optimization Problem 
N = 200; % Number of collocation nodes 
tf = 10; % final time fixed problem
h_step = tf/N; 
Waypoint = false;
WP = [2,150;5,50]; % for zk

% Define Decision Variables, X (3*N+3)
n_var = 3; % [z,v,u]
n = 3*N+3; % dimension of X
x = zeros(n,1); % column [z0,v0,u0,...,zN,vN,uN], zk = k*n_var+1; vk = k*n_var+2; uk = k*n_var+3, k=0,1.,...,N

% Define Equality Constraints (AX=b or AX-b=0, (N+N+2+2)-by-n)
% 1. Dynamic Constraints (N+N)
A_dyn = zeros(2*N,n);
for k=1:N
    A_dyn(k,[(k-1)*n_var+1, k*n_var+1, (k-1)*n_var+2, k*n_var+2, (k-1)*n_var+3, k*n_var+3]) = [-1, 1, -h_step/2, -h_step/2, h_step^2/12, -h_step^2/12]; % k=1일때, z0,z1,v0,v1,u0,u1 
    A_dyn(k+N,[(k-1)*n_var+2, k*n_var+2, (k-1)*n_var+3, k*n_var+3]) = [-1, 1, h_step/2, h_step/2]; % k=1일때, v0,v1,u0,u1 
end

% 2. Initial and Final Conditions
A_cond = zeros(2+2,n);
k=1; A_cond(1,(k-1)*n_var+1) = 1; A_cond(2,(k-1)*n_var+2) = 1; % z0,v0
k=N; A_cond(3,    k*n_var+1) = 1; A_cond(4,k*n_var+2) = 1;     % zN,vN

% 3. Form matrix A
A = [A_dyn;A_cond];

% 4. Form RHS matrix b, (N+N+2+2)-by-1
b = zeros(2*N+4,1); % column
for k=1:N
    b(k+N) = g*h_step;
end
b([2*N+1,2*N+2,2*N+3,2*N+4]) = [100;10;0;0];

% Define Inequality Constraints (Cx<=d or Cx-d<=0, 2*N+2)
% 1. Upper and Lower bounds on variables 
C = [eye(n);-eye(n)];
d = 1e3*ones(2*n,1); % column

% Define Objective Function 
% 1.objfunc = (u*u')*tf/(2*N); 
% 2.objfunc = (u(1:end-1)*u(1:end-1)' + u(1:end-1)*u(2:end)' + u(2:end)*u(2:end)')*tf/(3*N)/2; 
%obj_func = @(x) sum(x(n_var:3:n,:).^2,1)*tf/(2*N); % u = x(n_var:3:n); % column
obj_func = @(x) (sum(x(n_var:3:n-n_var,:).^2,1) + sum(x(n_var:3:n-n_var,:).*x(2*n_var:3:n,:),1) + sum(x(2*n_var:3:n,:).^2,1) )*tf/(3*N)/2; 
obj1 = false; %true if obj_func1 is used

% Define Waypoint (if any) and add constraints
if Waypoint
    sz = size(WP);
    A_wp = zeros(sz(2),n);
    b_wp = zeros(sz(2),1);
    for i=1:sz(2)
        ti = WP(i,1); % time 
        idx = ti/(tf/N); % index
        A_wp(i,idx*n_var+1) = 1;
        b_wp(i) = WP(i,2);
    end
    A = [A;A_wp];
    b = [b;b_wp];
    % for logging
    trace_data_wp = zeros(100,4);
else
    trace_data = zeros(100,4);
end

%% Solve by Barrier Method 

% First, try with strictly feasible solution x0 = ones(n,1)
x0   = ones(n,1);

x_k  = x0;
p    = length(b); % number of equality constraints 
m    = length(d); % number of inequality constraints
H_k  = eye(n); I = eye(n);
nu_k = zeros(p,1);
t0   = 10; % Initialize centering step 
pk   = [x_k;t0];
BM_obj = @Project_BM_objfunc;

g_k = num_grad(BM_obj, h, pk)';
% Check with analytic gradient
d_prime = 1./(d-C*x_k);
g_k_check = t0*num_grad(obj_func,h,x_k)' + C'*d_prime;
assert(norm(g_k-g_k_check,2)<1e-4,"numerical gradient very different from analytical gradient for BM")

yk = [x_k; nu_k;g_k]; % (n)+(p)+(n)
residual_func = @(yk) [yk(n+p+1:n+p+n)+A'*yk(n+1:n+p); A*yk(1:n)-b]; % equiv to @(x_k,nu_k,g_k) [g_k+A'*nu_k;A*x_k-b]
residual = residual_func(yk);

mu = 10; maxiter=1e2; j=1; done = false;
while m/pk(n+1) > BM_tol && ~done %pk(n+1) = t0
    %pk
    k=1;
    d_prime = 1./(d-C*pk(1:n)); %pk(1:n) = x
    H_k = pk(n+1)*I + C'*diag(d_prime.^2)*C; % initialize Hessian with analytic results
    while norm(residual,2) > conv_tol
        if  k > maxiter 
            break;
        end
        if j>7 && norm(A*pk(1:n)-b,2) < h*1e-1 % j counts the magnitude of pk(n+1) = t or tau
           done = true;
           break;
        end
        % 1. Compute KKT newton_step using current x_k, H_k by LDL' factorisation
        S = [H_k, A';A, zeros(p,p)]; 
        [L,D,P] = ldl(S);
        newton_step = -P*(L^-1)'*(D^-1)*(L^-1)*P' * residual;
        pn_step = newton_step(1:n);
        dn_step = newton_step(n+1:n+p);
        % 2. Backtracking Line Search on Residual
        t = BacktrackingLineSearch_BFGS_EQ(residual_func, yk, -residual, n, p);
        % 3. Update x_k, nu_k, g_k, yk, residual using pk
        pk(1:n) = pk(1:n) + t*pn_step; 
        nu_k = nu_k + t*dn_step; 
        g_k = num_grad(BM_obj,h,pk)';
        yk = [pk(1:n);nu_k;g_k];
        residual = residual_func([pk(1:n);nu_k;g_k]);
        % 4. Update Hessian H_k with BFGS (not inverse)
        y_k = num_grad(BM_obj,h,[pk(1:n)+t*pn_step;pk(n+1)])' - num_grad(BM_obj,h,pk)';
        s_k = t*pn_step;
        H_k = H_k - ((H_k*s_k)*s_k'*H_k')/(s_k'*H_k*s_k) + (y_k*y_k')/(y_k'*s_k);        
        %z = pk(1:n_var:n)
        %v = pk(2:n_var:n)
        %u = pk(3:n_var:n)        
        % Save log
        trace_data((j-1)*100+k,1) = pk(n+1);
        trace_data((j-1)*100+k,2) = norm(residual,2);
        trace_data((j-1)*100+k,3) = norm(A*pk(1:n)-b,2);
        trace_data((j-1)*100+k,4) = obj_func(pk(1:n));   
        text = ['t=', num2str(pk(n+1)), ', norm(residual,2)=', num2str(norm(residual,2)), ', norm(A*x_k-b,2)=', num2str(norm(A*pk(1:n)-b,2)), ', obj_func=', num2str(obj_func(pk(1:n)))]; 
        disp(text)
        % check if BM_objfunc blows up
        assert(isreal(BM_obj(pk)),'BM_objfunc blows up! Increase initial t0')
        k=k+1;
    end
    
    % check norm(residual,2)
    if norm(residual,2) < conv_tol
        disp('norm(residual,2) is smaller than convergence tolerance');
        break;
    else
        % Update t0
        disp('Updating t0 and starting again...')
        pk(n+1) = mu*pk(n+1);
        j=j+1;
    end
    
    
end

% z = pk(1:n_var:n)
% v = pk(2:n_var:n)
% u = pk(3:n_var:n)

%% Save Solution
if ~Waypoint
    Barrier_Sol = pk(1:n);
    text = ['Barrier_Solver_N',num2str(N),'.mat'];
    filename = sprintf(text,k);
    save(filename, 'Barrier_Sol')
    text = ['Barrier_Log_N',num2str(N),'.mat'];
    filename = sprintf(text,k);
    save(filename, 'trace_data')
else
    Barrier_Sol_wp = pk(1:n);
    text = ['Barrier_Solver_N',num2str(N),'_wp','.mat'];
    filename = sprintf(text,k);
    save(filename, 'Barrier_Sol_wp')
    text = ['Barrier_Log_N',num2str(N),'_wp','.mat'];
    filename = sprintf(text,k);
    save(filename, 'trace_data')
end




