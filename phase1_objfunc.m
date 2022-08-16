function out = phase1_objfunc(pk)
global C d
sz = size(pk);
n = length(pk)-2;
if sz(2) == 1
    % pk = [x;s0;t0]; % column, x0 is n-dimensional, while s0 and t0 are scalars
    out = pk(n+2)*pk(n+1)-sum(log(pk(n+1)-C*pk(1:n)+d)); % equiv to @(x,s,t) t*s -sum(log(s-C*x+d));
else % pk contains perturbations for numerical computation of gradient, (n+1) by (n+1) since t = pk(n+2,1) is fixed
    out = pk(n+2,1).*pk(n+1,1:n+1)-sum(log(pk(n+1,1:n+1)-C*pk(1:n,1:n+1)+d),1); % since t is not a variable, only x and s are variables
end
    
end
