%% Plot Barrier Solution
close all; clear all;
tf = 10; n_var = 3;
%% 1. Vanilla problem without Waypoint Constraints
Sol_N20 = load('Barrier_Solver_N20.mat'); % struc, Sol_N20.Barrier_Sol
Sol_N50 = load('Barrier_Solver_N50.mat');
Sol_N100 = load('Barrier_Solver_N100.mat');
Sol_N200 = load('Barrier_Solver_N200.mat');

figure(1) 
tiledlayout(1,2)
nexttile
% Plot z(t) 
n = length(Sol_N20.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N20.Barrier_Sol(1:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N50.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N50.Barrier_Sol(1:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N100.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N100.Barrier_Sol(1:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N200.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N200.Barrier_Sol(1:n_var:n), 'LineWidth',2); hold on; grid on;
% Plot v(t)
n = length(Sol_N20.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N20.Barrier_Sol(2:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N50.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N50.Barrier_Sol(2:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N100.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N100.Barrier_Sol(2:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N200.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N200.Barrier_Sol(2:n_var:n), 'LineWidth',2); hold on; grid on;
l = legend('z(t) with N=20','z(t) with N=50','z(t) with N=100','z(t) with N=200','v(t) with N=20','v(t) with N=50','v(t) with N=100','v(t) with N=200');  
set(l, 'FontSize',10);
t = title('Optimal Trajectory by Hermite Simpson Method using Barrier Method'); t.FontSize = 13; t.FontWeight = 'bold';
xl = xlabel('Time [sec]'); xl.FontSize = 12; xl.FontWeight = 'bold';
yl = ylabel('[m] and [m/s]'); yl.FontSize = 12; yl.FontWeight = 'bold';
set(gcf, 'position', [50 450 1800 500]);

nexttile
% Plot u(t)
n = length(Sol_N20.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N20.Barrier_Sol(3:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N50.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N50.Barrier_Sol(3:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N100.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N100.Barrier_Sol(3:n_var:n), 'LineWidth',2); hold on; grid on;
n = length(Sol_N200.Barrier_Sol);plot(linspace(0,tf,n/n_var), Sol_N200.Barrier_Sol(3:n_var:n), 'LineWidth',2); hold on; grid on;
l = legend('N=20','N=50','N=100','N=200');  
set(l, 'FontSize',10);
t = title('Optimal Control Input u(t) by Hermite Simpson Method using Barrier Method'); t.FontSize = 13; t.FontWeight = 'bold';
xl = xlabel('Time [sec]'); xl.FontSize = 12; xl.FontWeight = 'bold';
yl = ylabel('[m/s^{2}]'); yl.FontSize = 12; yl.FontWeight = 'bold';
set(gcf, 'position', [50 450 1800 500]);

% Plot log
log_N20 = load('Barrier_Log_N20.mat');
log_N50 = load('Barrier_Log_N50.mat');
log_N100 = load('Barrier_Log_N100.mat');
log_N200 = load('Barrier_Log_N200.mat');

figure(2)
tiledlayout(1,2)
nexttile
% Plot norm(residual,2) 
nz=min(nnz(log_N20.trace_data(:,1)),100); r=log_N20.trace_data(:,2); plot(1:nz,r(1:nz),'LineWidth',2); hold on; grid on;
nz=min(nnz(log_N50.trace_data(:,1)),100); r=log_N50.trace_data(:,2); plot(1:nz,r(1:nz),'LineWidth',2); 
nz=min(nnz(log_N100.trace_data(:,1)),100); r=log_N100.trace_data(:,2); plot(1:nz,r(1:nz),'LineWidth',2); 
nz=min(nnz(log_N200.trace_data(:,1)),100); r=log_N200.trace_data(:,2); plot(1:nz,r(1:nz),'LineWidth',2);
l = legend('N=20','N=50','N=100','N=200');  
set(l, 'FontSize',10);
t = title('Residual Plot $|\!|\bf{r}|\!|_{2}$ for $\tau$=10', 'interpreter', 'latex');
set(t, 'FontSize',20);
xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
set(gcf, 'position', [50 450 1800 500]);

nexttile
nz=nnz(log_N100.trace_data(:,1)); r=log_N100.trace_data(:,2); plot(101:nz,r(101:nz),'LineWidth',2); hold on; grid on;
nz=nnz(log_N200.trace_data(:,1)); r=log_N200.trace_data(:,2); plot(101:nz,r(101:nz),'LineWidth',2);
l = legend('N=100','N=200');  
set(l, 'FontSize',10);
t = title('Residual Plot $|\!|\bf{r}|\!|_{2}$ for $\tau$=100', 'interpreter', 'latex');
set(t, 'FontSize',20);
xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
set(gcf, 'position', [50 450 1800 500]);

figure(3)
tiledlayout(1,2)
nexttile
% Plot norm(Ax-b,2)
nz=min(nnz(log_N20.trace_data(:,1)),100); r=log_N20.trace_data(:,3); plot(1:nz,r(1:nz),'LineWidth',2); hold on; grid on;
nz=min(nnz(log_N50.trace_data(:,1)),100); r=log_N50.trace_data(:,3); plot(1:nz,r(1:nz),'LineWidth',2); 
nz=min(nnz(log_N100.trace_data(:,1)),100); r=log_N100.trace_data(:,3); plot(1:nz,r(1:nz),'LineWidth',2); 
nz=min(nnz(log_N200.trace_data(:,1)),100); r=log_N200.trace_data(:,3); plot(1:nz,r(1:nz),'LineWidth',2);
l = legend('N=20','N=50','N=100','N=200');  
set(l, 'FontSize',10);
t = title('Norm of Equality Constraints $|\!|\bf{Ax-b}|\!|_{2}$ for $\tau$=10', 'interpreter', 'latex');
set(t, 'FontSize',20);
xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
set(gcf, 'position', [50 450 1800 500]);

nexttile
nz=nnz(log_N100.trace_data(:,1)); r=log_N100.trace_data(:,3); plot(101:nz,r(101:nz),'LineWidth',2); hold on; grid on;
nz=nnz(log_N200.trace_data(:,1)); r=log_N200.trace_data(:,3); plot(101:nz,r(101:nz),'LineWidth',2);
l = legend('N=100','N=200');  
set(l, 'FontSize',10);
t = title('Norm of Equality Constraints $|\!|\bf{Ax-b}|\!|_{2}$ for $\tau$=100', 'interpreter', 'latex');
set(t, 'FontSize',20);
xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
set(gcf, 'position', [50 450 1800 500]);

figure(4)
tiledlayout(1,2)
nexttile
% Plot objective function
nz=min(nnz(log_N20.trace_data(:,1)),100); r=log_N20.trace_data(:,4); plot(1:nz,r(1:nz),'LineWidth',2); hold on; grid on;
nz=min(nnz(log_N50.trace_data(:,1)),100); r=log_N50.trace_data(:,4); plot(1:nz,r(1:nz),'LineWidth',2); 
nz=min(nnz(log_N100.trace_data(:,1)),100); r=log_N100.trace_data(:,4); plot(1:nz,r(1:nz),'LineWidth',2); 
nz=min(nnz(log_N200.trace_data(:,1)),100); r=log_N200.trace_data(:,4); plot(1:nz,r(1:nz),'LineWidth',2);
l = legend('N=20','N=50','N=100','N=200');  
set(l, 'FontSize',10);
t = title('Objective Function $\bf{f_o(x)}$ for $\tau$=10', 'interpreter', 'latex');
set(t, 'FontSize',20);
xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
set(gcf, 'position', [50 450 1800 500]);

nexttile
nz=nnz(log_N100.trace_data(:,1)); r=log_N100.trace_data(:,4); plot(101:nz,r(101:nz),'LineWidth',2); hold on; grid on;
nz=nnz(log_N200.trace_data(:,1)); r=log_N200.trace_data(:,4); plot(101:nz,r(101:nz),'LineWidth',2);
l = legend('N=100','N=200');  
set(l, 'FontSize',10);
t = title('Objective Function $\bf{f_o(x)}$ for $\tau$=100', 'interpreter', 'latex');
set(t, 'FontSize',20);
xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
set(gcf, 'position', [50 450 1800 500]);

%% 2. Problem with Waypoint Constraints
Waypoint = true;
if Waypoint
    WP = [2,150;5,50]; % for zk
    Sol_N20 = load('Barrier_Solver_N20_wp.mat');
    Sol_N50 = load('Barrier_Solver_N50_wp.mat');
    Sol_N100 = load('Barrier_Solver_N100_wp.mat');
    Sol_N200 = load('Barrier_Solver_N200_wp.mat');
    
    figure(5)
    tiledlayout(1,2)
    nexttile
    % Plot z(t)
    n = length(Sol_N20.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N20.Barrier_Sol_wp(1:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N50.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N50.Barrier_Sol_wp(1:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N100.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N100.Barrier_Sol_wp(1:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N200.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N200.Barrier_Sol_wp(1:n_var:n), 'LineWidth',2); hold on; grid on;
    % Plot v(t)
    n = length(Sol_N20.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N20.Barrier_Sol_wp(2:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N50.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N50.Barrier_Sol_wp(2:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N100.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N100.Barrier_Sol_wp(2:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N200.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N200.Barrier_Sol_wp(2:n_var:n), 'LineWidth',2); hold on; grid on;
    l = legend('z(t) with N=20','z(t) with N=50','z(t) with N=100','z(t) with N=200','v(t) with N=20','v(t) with N=50','v(t) with N=100','v(t) with N=200');
    set(l, 'FontSize',10);
    t = title('Optimal Trajectory by Hermite Simpson Method using Barrier Method'); t.FontSize = 13; t.FontWeight = 'bold';
    xl = xlabel('Time [sec]'); xl.FontSize = 12; xl.FontWeight = 'bold';
    yl = ylabel('[m] and [m/s]'); yl.FontSize = 12; yl.FontWeight = 'bold';
    for i=1:length(WP)
        plot(WP(i,1),WP(i,2),'rx','MarkerSize',10, 'LineWidth',2)
    end
    set(gcf, 'position', [50 450 1800 500]);
    
    nexttile
    % Plot u(t)
    n = length(Sol_N20.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N20.Barrier_Sol_wp(3:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N50.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N50.Barrier_Sol_wp(3:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N100.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N100.Barrier_Sol_wp(3:n_var:n), 'LineWidth',2); hold on; grid on;
    n = length(Sol_N200.Barrier_Sol_wp);plot(linspace(0,tf,n/n_var), Sol_N200.Barrier_Sol_wp(3:n_var:n), 'LineWidth',2); hold on; grid on;
    l = legend('N=20','N=50','N=100','N=200');
    set(l, 'FontSize',10);
    t = title('Optimal Control Input u(t) by Hermite Simpson Method using Barrier Method'); t.FontSize = 13; t.FontWeight = 'bold';
    xl = xlabel('Time [sec]'); xl.FontSize = 12; xl.FontWeight = 'bold';
    yl = ylabel('[m/s^{2}]'); yl.FontSize = 12; yl.FontWeight = 'bold';
    set(gcf, 'position', [50 450 1800 500]);
    
    % Plot log
    log_N20 = load('Barrier_Log_N20_wp.mat');
    log_N50 = load('Barrier_Log_N50_wp.mat');
    log_N100 = load('Barrier_Log_N100_wp.mat');
    log_N200 = load('Barrier_Log_N200_wp.mat');
    
    figure(6)
    tiledlayout(1,2)
    nexttile
    % Plot norm(residual,2)
    nz=min(nnz(log_N20.trace_data(:,1)),100); r=log_N20.trace_data(:,2); plot(1:nz,r(1:nz),'LineWidth',2); hold on; grid on;
    nz=min(nnz(log_N50.trace_data(:,1)),100); r=log_N50.trace_data(:,2); plot(1:nz,r(1:nz),'LineWidth',2);
    nz=min(nnz(log_N100.trace_data(:,1)),100); r=log_N100.trace_data(:,2); plot(1:nz,r(1:nz),'LineWidth',2);
    nz=min(nnz(log_N200.trace_data(:,1)),100); r=log_N200.trace_data(:,2); plot(1:nz,r(1:nz),'LineWidth',2);
    l = legend('N=20','N=50','N=100','N=200');
    set(l, 'FontSize',10);
    t = title('Residual Plot $|\!|\bf{r}|\!|_{2}$ for $\tau$=10', 'interpreter', 'latex');
    set(t, 'FontSize',20);
    xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
    set(gcf, 'position', [50 450 1800 500]);
    
    nexttile
    nz=nnz(log_N100.trace_data(:,1)); r=log_N100.trace_data(:,2); plot(101:nz,r(101:nz),'LineWidth',2); hold on; grid on;
    nz=nnz(log_N200.trace_data(:,1)); r=log_N200.trace_data(:,2); plot(101:nz,r(101:nz),'LineWidth',2);
    l = legend('N=100','N=200');
    set(l, 'FontSize',10);
    t = title('Residual Plot $|\!|\bf{r}|\!|_{2}$ for $\tau$=100', 'interpreter', 'latex');
    set(t, 'FontSize',20);
    xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
    set(gcf, 'position', [50 450 1800 500]);
    
    figure(7)
    tiledlayout(1,2)
    nexttile
    % Plot norm(Ax-b,2)
    nz=min(nnz(log_N20.trace_data(:,1)),100); r=log_N20.trace_data(:,3); plot(1:nz,r(1:nz),'LineWidth',2); hold on; grid on;
    nz=min(nnz(log_N50.trace_data(:,1)),100); r=log_N50.trace_data(:,3); plot(1:nz,r(1:nz),'LineWidth',2);
    nz=min(nnz(log_N100.trace_data(:,1)),100); r=log_N100.trace_data(:,3); plot(1:nz,r(1:nz),'LineWidth',2);
    nz=min(nnz(log_N200.trace_data(:,1)),100); r=log_N200.trace_data(:,3); plot(1:nz,r(1:nz),'LineWidth',2);
    l = legend('N=20','N=50','N=100','N=200');
    set(l, 'FontSize',10);
    t = title('Norm of Equality Constraints $|\!|\bf{Ax-b}|\!|_{2}$ for $\tau$=10', 'interpreter', 'latex');
    set(t, 'FontSize',20);
    xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
    set(gcf, 'position', [50 450 1800 500]);
    
    nexttile
    nz=nnz(log_N200.trace_data(:,1)); r=log_N200.trace_data(:,3); plot(101:nz,r(101:nz),'LineWidth',2); grid on;
    l = legend('N=200');
    set(l, 'FontSize',10);
    t = title('Norm of Equality Constraints $|\!|\bf{Ax-b}|\!|_{2}$ for $\tau$=100', 'interpreter', 'latex');
    set(t, 'FontSize',20);
    xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
    set(gcf, 'position', [50 450 1800 500]);
    
    figure(8)
    tiledlayout(1,2)
    nexttile
    % Plot objective function
    nz=min(nnz(log_N20.trace_data(:,1)),100); r=log_N20.trace_data(:,4); plot(1:nz,r(1:nz),'LineWidth',2); hold on; grid on;
    nz=min(nnz(log_N50.trace_data(:,1)),100); r=log_N50.trace_data(:,4); plot(1:nz,r(1:nz),'LineWidth',2);
    nz=min(nnz(log_N100.trace_data(:,1)),100); r=log_N100.trace_data(:,4); plot(1:nz,r(1:nz),'LineWidth',2);
    nz=min(nnz(log_N200.trace_data(:,1)),100); r=log_N200.trace_data(:,4); plot(1:nz,r(1:nz),'LineWidth',2);
    l = legend('N=20','N=50','N=100','N=200');
    set(l, 'FontSize',10);
    t = title('Objective Function $\bf{f_o(x)}$ for $\tau$=10', 'interpreter', 'latex');
    set(t, 'FontSize',20);
    xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
    set(gcf, 'position', [50 450 1800 500]);
    
    nexttile
    nz=nnz(log_N200.trace_data(:,1)); r=log_N200.trace_data(:,4); plot(101:nz,r(101:nz),'LineWidth',2); grid on;
    l = legend('N=200');
    set(l, 'FontSize',10);
    t = title('Objective Function $\bf{f_o(x)}$ for $\tau$=100', 'interpreter', 'latex');
    set(t, 'FontSize',20);
    xl = xlabel('Iteration'); xl.FontSize = 12; xl.FontWeight = 'bold';
    set(gcf, 'position', [50 450 1800 500]);
    
end
