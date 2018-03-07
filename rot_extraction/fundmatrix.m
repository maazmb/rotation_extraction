function [f,e1,e2] = fundmatrix(x)
%computes the fundamental matrix from 8 or more points in a stereo pair of images. 
%The normalized 8 point algorithm has accurate results with 12 or more points.

if size(x,1) ~= 6
     error('x must be 6 by N');
else
    x1 = x(1:3,:);
    x2 = x(4:6,:);
end

npts = size(x,2)

if npts<8
    error('Atleast 8 points are needed to compute the fundamental matrix');
end

%[X1, T1] = normalize(x1);
%[X2, T2] = normalize(x2);
X1 = x1;
X2 = x2;
%build the constraint matrix of 9 columns

A = [X2(1,:)'.*X1(1,:)'   X2(1,:)'.*X1(2,:)'  X2(1,:)' ...
     X2(2,:)'.*X1(1,:)'   X2(2,:)'.*X1(2,:)'  X2(2,:)' ...
     X1(1,:)'             X1(2,:)'            ones(npts,1) ];       

	[U,D,V] = svd(A,0); 
    

   
    F = reshape(V(:,9),3,3)';
    
   
    [U,D,V] = svd(F,0);
    F = U*diag([D(1,1) D(2,2) 0])*V';
   
   
    f=F;
    
    %with normalized coordinates..det is zero.
    if nargout == 3  	% Solve for epipoles
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
    
    
    
    
    
