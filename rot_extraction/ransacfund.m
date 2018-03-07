function [f, inliers,e1,e2] = ransacfund(x1,x2,t,feedback)

%x1 and x2 are 3byn homogenous 2d points.

%F       - The 3x3 fundamental matrix such that x2'Fx1 = 0.
%inliers - An array of indices of the elements of x1, x2 that were the inliers for the best model.

if nargin == 3
    feedback = 1; 
end    


[X1,T1] = normalize(x1);
[X2,T2] = normalize(x2);

s = 8; 

%create handles to functions
modelfunct = @fundmatrix;
distfunct = @funddist;
degenfunct = @isdegenerate;

[F, inliers] = ransac([X1;X2],modelfunct, distfunct, degenfunct,s,t,feedback)



if isempty(F) 
    warning('No F computed by the ransac code')
    return;
end 
F = fundmatrix([X1(:,inliers);X2(:,inliers)]);

%Denormalize
f = T2'*F*T1;

display('F')
F
abs(X2'*F*X1)
det(X2'*F*X1)
display('f')
f
abs(x2'*f*x1)
det(x2'*f*x1)


if nargout == 4  	% Solve for epipoles
	[U,D,V] = svd(f,0);
	e1 = hnormalise(V(:,3));
	e2 = hnormalise(U(:,3));
    
    end
    
    
    
 function nx = hnormalise(x)
    [rows,npts] = size(x);
    nx = x;

    % Find the indices of the points that are not at infinity
    finiteind = find(abs(x(rows,:)) > eps);

    if length(finiteind) ~= npts
        warning('Some points are at infinity');
 end
 for r = 1:rows-1
	nx(r,finiteind) = x(r,finiteind)./x(rows,finiteind);
    end
    nx(rows,finiteind) = 1;
    




