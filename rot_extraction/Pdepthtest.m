   function [P2final] = Pdepthtest(P2E, K, x1,x2)
    xtestpoint1 = x1(:,1); 
    xtestpoint2 = x2(:,1);
    P1 = [eye(3),zeros(3,1)];
    P1K = K*P1;
    xhat1 = inv(K)*xtestpoint1;
    X = zeros(4,4);
    Depth = zeros(4,2);
    for i=1:4
        xhat2 = inv(K)*xtestpoint2;
 
        A = [P1(3,:).*xhat1(1,1)-P1(1,:);
             P1(3,:).*xhat1(2,1)-P1(2,:);
             P2E(3,:,i).*xhat2(1,1)-P2E(1,:,i);
             P2E(3,:,i).*xhat2(2,1)-P2E(2,:,i)];
        % Normalize matrix A
        A1n = sqrt(sum(A(1,:).*A(1,:)));
        A2n = sqrt(sum(A(2,:).*A(2,:)));
        A3n = sqrt(sum(A(1,:).*A(1,:)));
        A4n = sqrt(sum(A(1,:).*A(1,:))); 
        Anorm = [A(1,:)/A1n;
                 A(2,:)/A2n;
                 A(3,:)/A3n;
                 A(4,:)/A4n];
             
        [U,S,V] = svd(Anorm);
        X(:,i) = V(:,end);
        % Check depth on second camera
        xi = P2E(:,:,i)*X(:,i);
        w = xi(3);
        T = X(end,i);
        mn = sqrt(sum(P2E(3,1:3,i).*P2E(3,1:3,i)));
        Depth(i,1) = (sign(det(P2E(:,1:3,i)))*w)/(T*mn);
        
        % Check depth on first camera
        xi = P1(:,:)*X(:,i);
        w = xi(3);
        T = X(end,i);
        mn = sqrt(sum(P1(3,1:3).*P1(3,1:3)));
        Depth(i,2) = (sign(det(P1(:,1:3)))*w)/(T*mn);
    end;
 
    % Check which solution is the right one and return
    if(Depth(1,1)>0 && Depth(1,2)>0)
        P2final = P2E(:,:,1);
    elseif(Depth(2,1)>0 && Depth(2,2)>0)
        P2final = P2E(:,:,2);    
    elseif(Depth(3,1)>0 && Depth(3,2)>0)
        P2final = P2E(:,:,3);
    elseif(Depth(4,1)>0 && Depth(4,2)>0)
        P2final = P2E(:,:,4);
    end;
