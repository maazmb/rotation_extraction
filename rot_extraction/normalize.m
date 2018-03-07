function [Xn,T] = normalize(x)
%function translates 2d homogenous points (3*n) and normalizes them. the process typically improves the 
%conditioning of any equation used to solve H,F,etc.
%points that have already been normalized will have no effect and the T
%will just be a diagonal of ones.
%Xn are new points
%T is a transformation matrix such that newpoints = T*x(old points)

%the scaling parameter is normalized to 1 unless the point is at infinty.

%find the indices of points not at infinity

finiteind = find(abs(x(3,:) > eps));


if length(finiteind)~= size(x,2)
    warning('some points are at infinity');
end

%for the finite points ensure homogenous coordinates have scale of 1.
x(1,finiteind) = x(1,finiteind)./x(3,finiteind);
x(2,finiteind) = x(2,finiteind)./x(3,finiteind);
x(3,finiteind) = 1;

centroid = mean(x(1:2,finiteind)')';
%shift origin to centroid
Xn(1,finiteind) = x(1,finiteind)-centroid(1);
Xn(2,finiteind) = x(2,finiteind)-centroid(2);

dist = sqrt(Xn(1,finiteind).^2 + Xn(2,finiteind).^2);
meandist = mean(dist(:)); %ensure dist is a column vector

scale = sqrt(2)/meandist;

T = [scale 0 -scale*centroid(1); 0 scale -scale*centroid(2); 0 0 1];



Xn = T*x;
y=1;


