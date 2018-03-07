function [Inliers] = funddist(F,x,t);
%function to evaluate the first order approximation of the geometric error
%(sampson distance) of the fit of a fundamental matrix w.r.t a set of
%matched points as needed by RANSAC.

x1 = x(1:3,:);
x2 = x(4:6,:);

x2tFx1 = zeros(1,length(x1));
%transpose
for n = 1:length(x1)
    x2tFx1(n) = x2(:,n)'*F*x1(:,n);
end

Fx1 = F*x1;
Ftx2 = F'*x2;

d = (x2tFx1.^2)./(Fx1(1,:).^2 + Fx1(2,:).^2 + Ftx2(1,:).^2 + Ftx2(2,:).^2);

Inliers = find(abs(d)<t);





