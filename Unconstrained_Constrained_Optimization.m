%% Practice Code for Unconstrained Optimization
% 1. Himmelblau Function as benchmark function for minimization problem
% 2. Implement steepest descent method with golden section search
% 3. Implement Newton's method with backtracking line search (later compare with raw Newton method)
% 4. Implement BFGS method with both goldsen section search and backtracking line search

%% Practice Code for Constrained Optimization
% 1. Implement Equality-constrained Newton's method (BFGS is used)
% 2. Implement barrier method for inequality constrained problem

clear all; close all;
h = 1e-4;
global C d % inequality constraints
%% 1. Himmelblau function  
% % https://pythonhosted.org/algopy/examples/minimization/himmelblau_minimization.html
% [x,y] = meshgrid(-7:0.01:7,-7:0.01:7); % easy_init = [3.1; 2.1],hard_init = [-0.27; -0.9]
% scores = himmelblaufunc(x,y);
% figure(1)
% hp = surf(x,y,scores); hold on;                  hp.Annotation.LegendInformation.IconDisplayStyle='off';
% hp = plot3(3,2,0,'rx','LineWidth', 2);           hp.Annotation.LegendInformation.IconDisplayStyle='off';
% hp = plot3(-2.80,3.13,0,'rx','LineWidth', 2);    hp.Annotation.LegendInformation.IconDisplayStyle='off';
% hp = plot3(-3.78,-3.28,0,'rx','LineWidth', 2);   hp.Annotation.LegendInformation.IconDisplayStyle='off';
% hp = plot3(3.58,-1.85,0,'rx','LineWidth', 2);    hp.Annotation.LegendInformation.IconDisplayStyle='off';
% colorbar
% shading interp
% title('Himmelblau Function'); xlabel('x');ylabel('y');zlabel('f(x,y)')
% save Himmelblau.mat scores
% saveas(gcf,'Himmelblau_3D.fig')
% figure(2)
% contour(x,y,scores, 100)
% title('f(x,y)'); xlabel('x');ylabel('y'); hold on;
% hp = plot3(3,2,0,'rx','LineWidth', 2);           hp.Annotation.LegendInformation.IconDisplayStyle='off';
% hp = plot3(-2.80,3.13,0,'rx','LineWidth', 2);    hp.Annotation.LegendInformation.IconDisplayStyle='off';
% hp = plot3(-3.78,-3.28,0,'rx','LineWidth', 2);   hp.Annotation.LegendInformation.IconDisplayStyle='off';
% hp = plot3(3.58,-1.85,0,'rx','LineWidth', 2);    hp.Annotation.LegendInformation.IconDisplayStyle='off';
% saveas(gcf,'Himmelblau_Contour.fig')

trace_data = zeros(2,20000,9);

%% 2. Steepest Descent with Golded Section Search (exact and numerical gradient)
% Search direction sk = -grad_fk
x0 = [-0.7,0.2]';  % column vector
x0 = [-0.27; -0.9];
conv_tol = 1e-6;
func = @himmelblaufunc; %fval = func(x0);

disp('Steepest Descent with Golden Section Search using Exact Gradients')
x_k = x0;
s_k = -grad_Him(x_k); % exact gradient
text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' norm(s_k,2)=',num2str(norm(s_k,2))];
disp(text)
trace_data(:,1,1) = x_k;
k=1;
while norm(s_k,2) > conv_tol
    [eta, interval, fmin] = GoldenSectionLineSearch(func, x_k, s_k);
    x_k = x_k+eta*s_k;
    % use exact gradient
    s_k = -grad_Him(x_k);
    
    text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' norm(s_k,2)=',num2str(norm(s_k,2))];
    disp(text)
    k = k+1;
    trace_data(:,k,1) = x_k;
end
disp('-------------------------------------------------------------------------------------------------------------------------')
disp('Steepest Descent with Golden Section Search using Numerical Gradients')
x_k = x0;
s_k_num = -num_grad(func,h,x_k)'; % numerical gradient (column vector)
text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' norm(s_k,2)=',num2str(norm(s_k_num,2))];
disp(text)
trace_data(:,1,2) = x_k;
k=1;
while norm(s_k_num,2) > conv_tol
    [eta, interval, fmin] = GoldenSectionLineSearch(func, x_k, s_k_num);
    x_k = x_k + eta*s_k_num;
    % use numerical gradient
    s_k_num = -num_grad(func,h,x_k)';
    
    text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' norm(s_k,2)=',num2str(norm(s_k_num,2))];
    disp(text)
    k = k+1;
    trace_data(:,k,2) = x_k;
end

%% 3. (Damped) Newton's method with backtracking line search
%x0=rand(2,1)*5; % initial point MUST be near critical point. Local convexity is required. Otherwise, newton_step may not be descent for nonconvex function
x0 = [4.5787;-3.9610];
Hess_k = Hess_Him(x0);
g_k = grad_Him(x0); % or g_k = num_grad(func,h,x0)
newton_step = -Hess_k\g_k;
assert(g_k'*newton_step < 0, 'Newton step is not descending. Choose another point!')
newton_dec = -g_k'*newton_step; % this is actually newton_dec^2
disp('-------------------------------------------------------------------------------------------------------------------------')
disp('Damped Newton''s method with Backtracking Line Search using Exact Gradients')
x_k=x0;
text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' newton_dec=',num2str(newton_dec)];
disp(text)
trace_data(:,1,3) = x_k;
k=1;
while abs(newton_dec)/2 > conv_tol
    t = BacktrackingLineSearch_Himm(func, x_k, g_k, newton_step);
    x_k = x_k+t*newton_step;
    g_k = grad_Him(x_k);
    Hess_k = Hess_Him(x_k);
    newton_step = -Hess_k\g_k;
    newton_dec = -g_k'*newton_step;
    
    text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' newton_dec=',num2str(newton_dec)];
    disp(text)
    k = k+1;
    trace_data(:,k,3) = x_k;
    assert(g_k'*newton_step < 0, 'Newton step is not descending. Choose another point!')    
end

%% 4. Quasi-Newton Method with BFGS (Golden Section Search and Backtracking LineSearch)
x0=[-0.27; -0.9]; % hard_init -> GoldenSection better than Backtracking for this hard_init example
disp('-------------------------------------------------------------------------------------------------------------------------')
disp('BFGS method with Golden Section Line Search using Numerical Gradients')
x_k = x0; H_k = eye(length(x0)); I = eye(length(x0));
g_k = num_grad(func,h,x_k)';
text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' norm(g_k,2)=',num2str(norm(g_k,2))];
disp(text)
trace_data(:,1,4) = x_k;
k=1;
while norm(g_k,2) > conv_tol
    d_k = -H_k*g_k;
    [eta, interval, fmin] = GoldenSectionLineSearch(func, x_k, d_k); % d_k : search direction
    y_k = num_grad(func,h,x_k+eta*d_k)' - num_grad(func,h,x_k)'; 
    s_k = eta*d_k;
    rho_k = 1/(y_k'*s_k);
    H_k = (I-rho_k*s_k*y_k')*H_k*(I-rho_k*y_k*s_k') + (rho_k*s_k)*s_k';
    x_k = x_k + eta*d_k;
    g_k = num_grad(func,h,x_k)';
    text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' norm(g_k,2)=',num2str(norm(g_k,2))];
    disp(text)
    k = k+1;
    trace_data(:,k,4) = x_k;

end


disp('-------------------------------------------------------------------------------------------------------------------------')
disp('BFGS method with Backtracking LineSearch using Numerical Gradients')
x_k = x0; H_k = eye(length(x0)); I = eye(length(x0));
g_k = num_grad(func,h,x_k)';
text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' norm(g_k,2)=',num2str(norm(g_k,2))];
disp(text)
trace_data(:,1,5) = x_k;
k=1;
while norm(g_k,2) > conv_tol
    d_k = -H_k*g_k;
    t = BacktrackingLineSearch_BFGS(func, x_k, g_k, d_k);
    y_k = num_grad(func,h,x_k+t*d_k)' - num_grad(func,h,x_k)'; 
    s_k = t*d_k;
    rho_k = 1/(y_k'*s_k);
    H_k = (I-rho_k*s_k*y_k')*H_k*(I-rho_k*y_k*s_k') + (rho_k*s_k)*s_k';
    x_k = x_k + t*d_k;
    g_k = num_grad(func,h,x_k)';
    text = ['x=',num2str(x_k(1)), ' y=',num2str(x_k(2)), ' norm(g_k,2)=',num2str(norm(g_k,2))];
    disp(text)
    k = k+1;
    trace_data(:,k,5) = x_k;
end

%% Plot 
%load('Himmelblau.mat');
sz = size(trace_data);
options = [1 0 0;0 1 0;0 0 1;0 1 1;1 0 1;1 0 0;0 1 0;0 0 1]; % red, green, blue, cyan, magenta

fig1 = openfig('Himmelblau_3D.fig');
for i=1:5
    nz = ceil(nnz(trace_data(:,:,i))/2); % number of nonzero entries
    x = trace_data(:,1:nz,i);
    p=plot3(x(1,:),x(2,:),func(x),'-','Color',options(i,:));
    p(1).Marker = '*';
end
I=legend('SteepestDescent(GSLS,ExactGrad)','SteepestDescent(GSLS,NumGrad)','Damped Newton''s Method(BTLS,ExactGrad)','BFGS(GSLS,NumGrad)','BFGS(BTLS,NumGrad)');
I.FontSize = 12;

fig2 = openfig('Himmelblau_Contour.fig');
for i=1:5
    nz = ceil(nnz(trace_data(:,:,i))/2); % number of nonzero entries
    x = trace_data(:,1:nz,i);
    p=plot(x(1,:),x(2,:),'-','Color',options(i,:));
    p(1).Marker = '*';
end
I=legend('Contour','SteepestDescent(GSLS,ExactGrad)','SteepestDescent(GSLS,NumGrad)','Damped Newton''s Method(BTLS,ExactGrad)','BFGS(GSLS,NumGrad)','BFGS(BTLS,NumGrad)');
I.FontSize = 12;


%% (Linear) Equality Constrained Optimization 
% Add an equality constraint y = -x+1
x_constr = linspace(-6,7,100);
y_constr = -x_constr+1; % constr: x+y = 1 -> A = [1,1], x = [x,y]', b = 1

fig1 = openfig('Himmelblau_3D.fig'); hold on;
plot3(x_constr,y_constr, func(x_constr,y_constr) ,'-k','LineWidth',2); 
fig2 = openfig('Himmelblau_Contour.fig');  hold on;
plot(x_constr,y_constr,'-k','LineWidth',2); 


%% Infeasible start Newton Method with BFGS for Hessian update (Golden Section Line Search and Backtracking Line Search) 
conv_tol = 1e-5;
x0 = [3.1; 2.1]; % infeasible point to start
%x0 = [0.1; -4.1]; % infeasible point to start (GoldenSection Line Search fails for this starting point)
%x0 = [-0.27; -0.9]; % infeasible point to start (hard init/GoldenSection Line Search fails)
%x0 = randn(2,1)*2
disp('-------------------------------------------------------------------------------------------------------------------------')
disp('Infeasible start Newton Method with Golden Section Line Search using Numerical Gradients')
x_k = x0; H_k = eye(length(x0)); I = eye(length(x0)); trace_data(:,1,6) = x_k;
g_k = num_grad(func,h,x_k)';
% Equality constraints (Ax=b or Ax-b=0)
A = [1,1];b=1;
n = length(x0); p = length(b); 
nu_k = zeros(p,1); % nu_k is p-dimensional = number of equality constraints
yk = [x_k;nu_k;g_k]; % (n)+(p)+(n)
residual_func = @(yk) [yk(n+p+1:n+p+n)+A'*yk(n+1:n+p); A*yk(1:n)-b]; % equiv to @(x_k,nu_k,g_k) [g_k+A'*nu_k;A*x_k-b]
residual = residual_func(yk);
% Start iteration to solve KKT system with BFGS Hessian update. 
% Let S denote KKT matrix
k=1;  %H_k = [80,20;20,39]; 
while norm(residual,2) > conv_tol
    % 1. Compute KKT newton_step using current x_k, H_k by LDL' factorisation
    S = [H_k, A';A, zeros(p,p)];
    [L,D,P] = ldl(S); % P'SP = LDL' => S = PLDL'P', P orthogonaL (P' = P^-1)
    newton_step = -P*(L^-1)'*(D^-1)*(L^-1)*P' * residual; % S*y = -r(y) => y = -(S^-1)*r(y) where S^-1 = P*(L'^-1)*(D^-1)*(L^-1)*P'
    pn_step = newton_step(1:n); % primal_newton_step
    dn_step = newton_step(n+1:n+p); % dual_newton_step
    % 2. GoldenSection Search on Residual
    [t, interval] = GoldenSectionLineSearch_EQ(residual_func, yk, -residual, n, p);
    % 3. Update x_k, nu_k, g_k, yk, residual  
    x_k = x_k + t*pn_step;
    nu_k = nu_k + t*dn_step;
    g_k = num_grad(func,h,x_k)';
    yk = [x_k;nu_k;g_k];
    residual = residual_func([x_k;nu_k;g_k]);
    % 4. Update Hessian H_k with BFGS (not inverse)
    y_k = num_grad(func,h,x_k+t*pn_step)' - num_grad(func,h,x_k)'; % y_k = g_(k+1) - g_k (g_k:gradient)
    s_k = t*pn_step;
    H_k = H_k - ((H_k*s_k)*s_k'*H_k')/(s_k'*H_k*s_k) + (y_k*y_k')/(y_k'*s_k);
    text = ['x=',num2str(x_k(1)), ', y=',num2str(x_k(2)), ', norm(residual,2)=',num2str(norm(residual,2)), ', A*x_k-b=',num2str(A*x_k-b),', obj_func=',num2str(func(x_k))];
    disp(text)
    plot(x_k(1),x_k(2),'rx');
    k = k+1;
    trace_data(:,k,6) = x_k;
end

%x0 = [0.1; -4.1];
x0 = [-0.27; -0.9];
disp('-------------------------------------------------------------------------------------------------------------------------')
disp('Infeasible start Newton Method with Backtracking LineSearch using Numerical Gradients')
x_k = x0; H_k = eye(length(x0)); I = eye(length(x0)); trace_data(:,1,7) = x_k;
g_k = num_grad(func,h,x_k)';
% Equality constraints (Ax=b or Ax-b=0)
A = [1,1];b=1; 
n = length(x0); p = length(b); 
nu_k = zeros(p,1); % nu_k is p-dimensional = number of equality constraints
yk = [x_k;nu_k;g_k]; % (n)+(m)+(n)
residual_func = @(yk) [yk(n+p+1:n+p+n)+A'*yk(n+1:n+p); A*yk(1:n)-b]; % equiv to @(x_k,nu_k,g_k) [g_k+A'*nu_k;A*x_k-b]
residual = residual_func(yk);
% Start iteration to solve KKT system with BFGS Hessian update. 
% Let S denote KKT matrix
k=1; alpha = 1e-2; beta = 0.8; 

%H_k = [80,20;20,39]; 
while norm(residual,2) > conv_tol
    % 1. Compute KKT newton_step using current x_k, H_k by LDL' factorisation
    S = [H_k, A';A, zeros(p,p)];
    [L,D,P] = ldl(S); % P'SP = LDL' => S = PLDL'P', P orthogonaL (P' = P^-1)
    newton_step = -P*(L^-1)'*(D^-1)*(L^-1)*P' * residual; % S*y = -r(y) => y = -(S^-1)*r(y) where S^-1 = P*(L'^-1)*(D^-1)*(L^-1)*P'
    pn_step = newton_step(1:n); % primal_newton_step
    dn_step = newton_step(n+1:n+p); % dual_newton_step
    % 2. Backtracking Line Search on Residual
    %t=1;
    %while norm(residual_func([x_k+t*pn_step;nu_k+t*dn_step;num_grad(func,h,x_k+t*pn_step)']),2) > (1-alpha*t)*norm(residual,2)
    %    t = beta*t;
    %end
    t = BacktrackingLineSearch_BFGS_EQ(residual_func, yk, -residual, n, p);
    % 3. Update x_k, nu_k, g_k, yk, residual  
    x_k = x_k + t*pn_step;
    nu_k = nu_k + t*dn_step;
    g_k = num_grad(func,h,x_k)';
    yk = [x_k;nu_k;g_k];
    residual = residual_func([x_k;nu_k;g_k]);
    % 4. Update Hessian H_k with BFGS (not inverse)
    y_k = num_grad(func,h,x_k+t*pn_step)' - num_grad(func,h,x_k)'; % y_k = g_(k+1) - g_k (g_k:gradient)
    s_k = t*pn_step;
    H_k = H_k - ((H_k*s_k)*s_k'*H_k')/(s_k'*H_k*s_k) + (y_k*y_k')/(y_k'*s_k);
    text = ['x=',num2str(x_k(1)), ', y=',num2str(x_k(2)), ', norm(residual,2)=',num2str(norm(residual,2)), ', A*x_k-b=',num2str(A*x_k-b),', obj_func=',num2str(func(x_k))];
    disp(text)
    plot(x_k(1),x_k(2),'gx');
    k = k+1;
    trace_data(:,k,7) = x_k;
end

%% Plot 
sz = size(trace_data);
fig1 = openfig('Himmelblau_3D.fig'); hold on;
plot3(x_constr,y_constr, func(x_constr,y_constr) ,'-k','LineWidth',2); 
for i=6:7
    nz = ceil(nnz(trace_data(:,:,i))/2); % number of nonzero entries
    x = trace_data(:,1:nz,i);
    p=plot3(x(1,:),x(2,:),func(x),'-','Color',options(i,:));
    p(1).Marker = '*';
end
I=legend('Equality Constraint (x+y=1)','Constrained Newton''s Method(GSLS)','Constrained Newton''s Method(BTLS)');
I.FontSize = 12;

fig2 = openfig('Himmelblau_Contour.fig');  hold on;
plot(x_constr,y_constr,'-k','LineWidth',2); 
for i=6:7
    nz = ceil(nnz(trace_data(:,:,i))/2); % number of nonzero entries
    x = trace_data(:,1:nz,i);
    p=plot(x(1,:),x(2,:),'-','Color',options(i,:));
    p(1).Marker = '*';
end
I=legend('Contour','Equality Constraint (x+y=1)', 'Constrained Newton''s Method(GSLS)','Constrained Newton''s Method(BTLS)');
I.FontSize = 12;

%% (Linear) Inequality Constrained Optimization
% Add an inequality constraint y<x+3

openfig('Himmelblau_Contour.fig');  hold on;
plot(x_constr,y_constr,'-k','LineWidth',2); % 기존 equality constraint
ineqplot('y>x+3',[-7,7,-7,7]); % x-y<-3 -> C = [1 -1], d = -3

%% Barrier Method for Inequality Constrained Problem
conv_tol = 1e-5; h = 1e-4;
x0 = [0; -4]; % infeasible point to start
trace_data(:,1,8) = x0;
disp('-------------------------------------------------------------------------------------------------------------------------')
disp('(Phase I) Barrier method with Basic Backtracking Line Search using Numerical Gradients')

% Equality Constraints (Ax=b or Ax-b=0)
A = [1,1]; b=1; % b in column
n = length(x0); p = length(b); 
A = [A zeros(p,1)]; % Rewrite A to include s variable, A now p by (n+1)
H_k = eye(length(x0)+1); I = eye(length(x0)+1);
nu_k = zeros(p,1); % nu_k is p-dimensional = number of equality constraints
% Inequality Constraints (Cx<d or Cx-d<0)
C = [1, -1]; d=-3; % d in column 
m = length(d);
% 1. Solve Phase I(Feasibility) Problem to obtain a strictly feasible solution using Barrier Method
t0 = 1;
% Pick s0 larger than max{f(x0)} = max{C*x0-d}
s0 = max(C*x0-d) + 7;
assert(log(s0-C*x0+d)>=0,"s0 is not strictly feasible for phase 1 problem")
pk = [x0;s0;t0]; % column, x0 is n-dimensional, while s0 and t0 are scalars
phase1_obj = @phase1_objfunc;
g_k = num_grad(phase1_obj,h,pk)';
d_tilda = 1./(s0+d-C*x0);
g_k_check = [C'*d_tilda;t0-ones(m,1)'*d_tilda];
assert(norm(g_k-g_k_check,2)<1e-7,"numerical gradient very different from analytical gradient for Phase 1")

x_k = x0; 
yk = [x_k;s0;nu_k;g_k]; % (n+1)+(p)+(n+1) addition of s scalar, now x_k contains point and s
residual_func = @(yk) [yk(n+p+2:n+p+n+2)+A'*yk(n+2:n+p+1); A*yk(1:n+1)-b]; % equiv to @(x_k,nu_k,g_k) [g_k+A'*nu_k;A*x_k-b]
residual = residual_func(yk);

mu = 10; alpha = 1e-2;
while pk(n+1) > 0 % x_k previous is now pk(1:n+1)! pk(1:n+1) = [x_k;s]
    pk 
    while norm(residual,2) > conv_tol
        % 1. Compute KKT newton_step using current x_k, H_k by LDL' factorisation
        S = [H_k, A';A, zeros(p,p)]; % A is now p by (n+1)
        [L,D,P] = ldl(S); 
        newton_step = -P*(L^-1)'*(D^-1)*(L^-1)*P' * residual;
        pn_step = newton_step(1:n+1);
        dn_step = newton_step(n+2:n+p+1);
        % 2. Basic Backtracking Line Search on Residual
        %t_BT=1;
        %while norm(residual_func([pk(1:n+1)+t_BT*pn_step;nu_k+t_BT*dn_step;num_grad(phase1_obj,h,[pk(1:n+1)+t_BT*pn_step;pk(n+2)])']),2) > (1-alpha*t_BT)*norm(residual,2)
        %    t_BT = beta*t_BT;
        %end
        % 2. Backtracking Line Search on Residual
        t = BacktrackingLineSearch_BFGS_EQ(residual_func, yk, -residual, n+1, p);
        % 3. Update x_k, s, nu_k, g_k, yk, residual using pk
        pk(1:n+1) = pk(1:n+1) + t*pn_step; % x_k, s updated together
        nu_k = nu_k + t*dn_step;
        g_k = num_grad(phase1_obj,h,pk)'; % now objective is phase1_obj
        yk = [pk(1:n+1);nu_k;g_k]; 
        residual = residual_func([pk(1:n+1);nu_k;g_k]);
        % 4. Update Hessian H_k with BFGS (not inverse)
        y_k = num_grad(phase1_obj,h,[pk(1:n+1)+t*pn_step;pk(n+2)])' - num_grad(phase1_obj,h,pk)'; % y_k = g_(k+1) - g_k (g_k:gradient)
        s_k = t*pn_step;
        H_k = H_k - ((H_k*s_k)*s_k'*H_k')/(s_k'*H_k*s_k) + (y_k*y_k')/(y_k'*s_k);
        text = ['t=', num2str(pk(n+2)),', x=',num2str(pk(1)), ', y=',num2str(pk(2)), ', s=',num2str(pk(n+1)), ', A*x_k-b=',num2str(A*pk(1:n+1)-b),', phase1_objfunc=',num2str(phase1_objfunc(pk))];
        disp(text)
        %plot(pk(1),pk(2),'mo');
        % check if s < 0
        if pk(n+1) < 0
            text = ['Strictly feasible solution is found, s=',num2str(pk(n+1))];
            disp(text);
            save p1_strictly_feasible_sol.mat pk
            break;
        end
    end
    if pk(n+1) < 0
        break;
    else
        disp('s is still positive at t=',num2str(pk(n+2)));
    end
    % Update t0
    disp('Updating t0 and starting again...')
    pk(n+2) = mu*pk(n+2);
    
end

% load strictly feasible solution
strictly_feasible=load('p1_strictly_feasible_sol.mat'); 
x0 = strictly_feasible.pk(1:n);
trace_data(:,2,8) = x0;

%% 2. Solve the approximated Barrier problem using Barrier method given a strcitly feasible x_k
%x0 = [-4.22;1.86]; % strictly feasible but does not satisfy equality constraints 
disp('-------------------------------------------------------------------------------------------------------------------------')
disp('Barrier method with Backtracking Line Search using Numerical Gradients')
x_k = x0; trace_data(:,1,9) = x_k; %plot(x_k(1),x_k(2),'mx');

% Equality Constraints (Ax=b or Ax-b=0)
A = [1,1]; b=1; % b in column
n = length(x0); p = length(b); 
H_k = eye(length(x0)); I = eye(length(x0));
nu_k = zeros(p,1); % nu_k is p-dimensional = number of equality constraints
% Inequality Constraints (Cx<d or Cx-d<0)
C = [1, -1]; d=-3; % d in column 
m = length(d);
% Initialize t0 for Centering Step
t0 = 5;
pk = [x_k;t0];
BM_obj = @BM_objfunc;
g_k = num_grad(BM_obj,h,pk)';
% check with analytical gradient
d_prime = 1./(d-C*x_k);
g_k_check = t0*num_grad(func,h,x_k)' + C'*d_prime;
assert(norm(g_k-g_k_check,2)<1e-7,"numerical gradient very different from analytical gradient for BM")

yk = [x_k;nu_k;g_k]; % (n)+(p)+(n)
residual_func = @(yk) [yk(n+p+1:n+p+n)+A'*yk(n+1:n+p); A*yk(1:n)-b]; % equiv to @(x_k,nu_k,g_k) [g_k+A'*nu_k;A*x_k-b]
residual = residual_func(yk);

mu = 10; maxiter=1e2; BM_tol = 1e-10; j=1; done = false;
while m/pk(n+1) > BM_tol && ~done %pk(n+1) = t0
    pk
    k=1;
    d_prime = 1./(d-C*pk(1:n));
    H_k = pk(n+1)*I + C'*diag(d_prime.^2)*C; % initialize Hessian with analytic results
    while norm(residual,2) > conv_tol
        if  k > maxiter 
            break;
        end
        if norm(A*pk(1:n)-b,2) < h*1e-1 
            done = true;
            break;
        end
        % 1. Compute KKT newton_step using current x_k, H_k by LDL' factorisation
        S = [H_k, A';A, zeros(p,p)]; 
        [L,D,P] = ldl(S);
        newton_step = -P*(L^-1)'*(D^-1)*(L^-1)*P' * residual;
        pn_step = newton_step(1:n);
        dn_step = newton_step(n+1:n+p);
%        % 2. Basic Backtracking Line Search on Residual
%         t_BT=1;
%         while norm(residual_func([pk(1:n+1)+t_BT*pn_step;nu_k+t_BT*dn_step;num_grad(phase1_obj,h,[pk(1:n+1)+t_BT*pn_step;pk(n+2)])']),2) > (1-alpha*t_BT)*norm(residual,2)
%             t_BT = beta*t_BT;
%         end
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
        text = ['t=', num2str(pk(n+1)),', x=',num2str(pk(1)), ', y=',num2str(pk(2)), ', A*x_k-b=',num2str(A*pk(1:n)-b), ', obj_func=',num2str(func(pk(1:n)))];
        disp(text)
        %plot(pk(1),pk(2),'rx');
        % check if BM_objfunc blows up
        assert(isreal(BM_objfunc(pk)),'BM_objfunc blows up! Increase initial t0')
        k=k+1;
    end
    j=j+1;
    trace_data(:,j,9) = pk(1:n);
    % Update t0
    disp('Updating t0 and starting again...')
    pk(n+1) = mu*pk(n+1);
    
end


%% Plot 
%load('Himmelblau.mat');
sz = size(trace_data);
options = [1 0 0;0 1 0;0 0 1;0 1 1;1 0 1;1 0 0;0 1 0;0 0 1;1 0 0]; % red, green, blue, cyan, magenta

for i=8:9
    nz = ceil(nnz(trace_data(:,:,i))/2); % number of nonzero entries
    x = trace_data(:,1:nz,i);
    p=plot(x(1,:),x(2,:),'-','Color',options(i,:));
    p(1).Marker = '*';
end
I=legend('Contour','Equality Constraint (x+y=1)','Inequality Constraint (x-y<-3)','Phase I by Newton''s method (BFGS,BTLS,NumGrad)','Barrier Method (BFGS,BTLS,NumGrad)');
I.FontSize = 12;









