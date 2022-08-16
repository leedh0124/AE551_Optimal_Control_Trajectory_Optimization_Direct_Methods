function del_f = num_grad(func,h,x_k)
% computes numerical gradient for any function at current X = [x1,x2,...,xN]'
% using Central Difference Method (2N extra evaluations required)
N = length(x_k);
I = eye(N);
del_f = (func(x_k+h*I)-func(x_k-h*I))/(2*h); % note this is row vector. Convert to column afterwards
end