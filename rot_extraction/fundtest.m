%variable already from main
x=[x1h(:,inliers);x2h(:,inliers)];

[f,e1,e2] = fundmatrix(x);

display(x2h(:,inliers)'*f*x1h(:,inliers))
display('a')
det(x2h(:,inliers)'*f*x1h(:,inliers))


Fm = estimateFundamentalMatrix(x1(:,inliers)',x2(:,inliers)')
display(x2h(:,inliers)'*Fm*x1h(:,inliers))

