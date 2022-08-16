function scores = himmelblaufunc(X,Y)
switch nargin
    case 1
        sz = size(X); % 2-by-n
        if sz(2) > 1 % many points
            scores = (X(1,:).^2+X(2,:)-11).^2 + (X(1,:)+X(2,:).^2-7).^2;
        else % one point
            scores = (X(1).^2+X(2)-11).^2 + (X(1)+X(2).^2 - 7).^2;
        end  
    case 2
        scores = (X.^2+Y-11).^2 + (X+Y.^2 - 7).^2;
        
end