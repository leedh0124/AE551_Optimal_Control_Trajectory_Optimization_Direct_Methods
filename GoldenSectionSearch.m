%% Golden Section Search Method (Zeroth order method)

% Find the minimum of f(x) = x^2-6*x+15 in [0,10]

func = @(x) x^2-6*x+15;
a = 0; b = 10;

% Construct Unimodal Interval 
alp = 0;


GR = (sqrt(5)-1)/2;
tol = 1e-8;
err = inf;
k = 0;
while err > tol
    d = GR*(b-a);
    x1 = a + d; % x1 > x2
    x2 = b - d;
    if func(x1) < func(x2)
        a = x2;
    else
        % func(x1) > func(x2)
        b = x1;
    end
    % check domain
    if abs(b-a) < tol
        err = abs(b-a);
    end
end