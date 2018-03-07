function Xhat = worldcords (x1,x2,P1,P2,homogeneous)
%this function finds world coordinates using triagulation as described
%Reference : Andrew Zisserman's Multiple View Geometry pg312 
%Xhat is homogenous too

Xhat = ones(4,size(x1,2));
for i = 1 : size(x1,2)
    if homogeneous
        A =[x1(1,i)*P1(3,:) - P1(1,:);
            x1(2,i)*P1(3,:) - P1(2,:);
            x2(1,i)*P2(3,:) - P2(1,:);
            x2(2,i)*P2(3,:) - P2(2,:)];    
        [u,s,v] = svd(A);
        Xhat(:,i) = v(:,end);
    else
        A =[x1(1,i)*P1(3,1:3) - P1(1,1:3);
            x1(2,i)*P1(3,1:3) - P1(2,1:3);
            x2(1,i)*P2(3,1:3) - P2(1,1:3);
            x2(2,i)*P2(3,1:3) - P2(2,1:3)];
        b = -[x1(1,i)*P1(3,4) - P1(1,4);
            x1(2,i)*P1(3,4) - P1(2,4);
            x2(1,i)*P2(3,4) - P2(1,4);
            x2(2,i)*P2(3,4) - P2(2,4)];    
        [u,s,v] = svd(A);
        bprime = u'*b;
        y = bprime(1:3) ./ (diag(s));
        Xhat(1:3,i) = v*y;
    end
end
Xhat = Xhat./repmat(Xhat(4,:),4,1);
end