%CODE TO CALCULATE THE ROTATION ANGLE DIFFERENCE BETWEEN TWO IMAGES (WHILE THE SCALE IS KEPT CONSTANT)
%written by Maaz But

% Load images
 I1=imread('twodeg3.jpg');
 I2=imread('twodeg4.jpg');
 
 %This code can be altered easily to calculate angle differences between
 %frames in a video stream using a loop.
 
 imshow(I1);


 ig1 = rgb2gray(I1);
 ig2 = rgb2gray(I2);
 
 image1 = im2double(ig1);
 image2 = im2double(ig2);

 
 % Get the Key Points
  Options.upright=true;
  Options.tresh=0.0001; 
  Ipts1=OpenSurf(I1,Options);
  Ipts2=OpenSurf(I2,Options);
  
 % points1 = detectSURFFeatures(image1);
 % points2 = detectSURFFeatures(image2);
 %If you want to use these basic matlab function..alter the code accordingly
  
  %[feature1,validpoints1] = extractFeatures(image1,points1);
  %[feature2,validpoints2] = extractFeatures(image2,points2);
 
% inputs,
%   I : The 2D input image color or greyscale
%   (optional)
%   Options : A struct with options (see below)
%
% outputs,
%   Ipts : A structure with the information about all detected Landmark points
%     Ipts.x , ipts.y : The landmark position
%     Ipts.scale : The scale of the detected landmark
%     Ipts.laplacian : The laplacian of the landmark neighborhood
%     Ipts.orientation : Orientation in radians
%     Ipts.descriptor : The descriptor for corresponding point matching
%
% options,
%   Options.verbose : If set to true then useful information is 
%                     displayed (default false)
%   Options.upright : Boolean which determines if we want a non-rotation
%                       invariant result (default false)
%   Options.extended : Add extra landmark point information to the
%                   descriptor (default false)
%   Options.tresh : Hessian response treshold (default 0.0002)
%   Options.octaves : Number of octaves to analyse(default 5)
%   Options.init_sample : Initial sampling step in the image (default 2)


% Put the landmark descriptors in a matrix
  D1 = reshape([Ipts1.descriptor],64,[]); 
  D2 = reshape([Ipts2.descriptor],64,[]);
  
  err=zeros(1,length(Ipts1));
  cor1=1:length(Ipts1); 
  cor2=zeros(1,length(Ipts1));
  
  for i=1:length(Ipts1),
      distance=sum((D2-repmat(D1(:,i),[1 length(Ipts2)])).^2,1);
      [err(i),cor2(i)]=min(distance);
  end
  
% Sort matches on vector distance in ascending order with indexes
  [err, ind]=sort(err); 
  cor1=cor1(ind); 
  cor2=cor2(ind);
  
% Show both images
  I = zeros([size(I1,1) size(I1,2)*2 size(I1,3)]);
  I(:,1:size(I1,2),:)=I1; I(:,size(I1,2)+1:size(I1,2)+size(I2,2),:)=I2;
 
  

  
 figure (2)
 imshow(I/255);
 hold on;
  
  Nofmatches = 30;%number of best matches required
  x1=zeros(2,Nofmatches);
  x2=zeros(2,Nofmatches);

 % Show the best matches
  
for i=1:Nofmatches,
      c=rand(1,3);
      x1(:,i) = [Ipts1(cor1(i)).x; Ipts1(cor1(i)).y];
      x2(:,i) = [Ipts2(cor2(i)).x; Ipts2(cor2(i)).y];
      plot([Ipts1(cor1(i)).x Ipts2(cor2(i)).x+size(I1,2)],[Ipts1(cor1(i)).y Ipts2(cor2(i)).y],'-','Color',c)
      plot([Ipts1(cor1(i)).x Ipts2(cor2(i)).x+size(I1,2)],[Ipts1(cor1(i)).y Ipts2(cor2(i)).y],'o','Color',c)
  
  end
  
  
  %making homogenous coordinates for fundamental matrix
  x1h = [x1;ones(1,length(x1))];
  x2h = [x2;ones(1,length(x1))];
  
  
  
  
  t = .001; 
%distance threshold for deciding outliers in ransac Note that point coordinates are normalised to that their
%mean distance from the origin is sqrt(2). The value of t should be set relative to this, say in the range 0.001 - 0.01.
  
[F,inliers,e1,e2] = ransacfund(x1h,x2h,t);
  
%displaying epipolar lines
l2 = F*x1h(:,inliers);
l1 = F'*x2h(:,inliers); 

%make sure homogenous points lie in z=1 plane
l1 = l1./l1(3);
l2 = l2./l2(3);




nofelines = 4

for i = 1:nofelines
    figure(i+2)
     subplot(1,2,1),imshow(I1)
     hold on

plot(x1h(1,inliers(i)),x1h(2,inliers(i)),'+','Color','r');
if abs(l1(1,i))> abs(l1(2,i)) %line is more vertical
  ylim = get(gca,'Ylim');
  p1 = cross(l1(:,i),[0 -1 ylim(1)]');
  p1 = p1./p1(3);
  p2 = cross(l1(:,i),[0 -1 ylim(2)]');
  p2 = p2./p2(3);
else
    xlim =get(gca,'Xlim');
    p1 = cross(l1(:,i),[-1 0 xlim(1)]');
    p1 = p1./p1(3);
    p2 = cross(l1(:,i),[-1 0 xlim(2)]');
    p2 = p2./p2(3);
end
%+size(I1,2)
 line([p1(1) p2(1)], [p1(2) p2(2)],'color', 'blue');

 subplot(1,2,2),imshow(I2)
 hold on
 plot(x2h(1,inliers(i)),x2h(2,inliers(i)),'+','Color','r');
 if abs(l2(1,i))> abs(l2(2,i)) %line is more vertical
  ylim = get(gca,'Ylim');
  p1 = cross(l2(:,i),[0 -1 ylim(1)]');
  p1 = p1./p1(3);
  p2 = cross(l2(:,i),[0 -1 ylim(2)]');
  p2 = p2./p2(3);
else
    xlim =get(gca,'Xlim');
    p1 = cross(l2(:,i),[-1 0 xlim(1)]');
    p1 = p1./p1(3);
    p2 = cross(l2(:,i),[-1 0 xlim(2)]');
    p2 = p2./p2(3);
 end
 %+size(I1,2)
line([p1(1) p2(1)], [p1(2) p2(2)],'color', 'blue');

end  

%find rotation and angle between frames

%for finding E..sample camera intrinsic/calibration matrix
K = [2321.4068 0 654.6514; 0 2325.5211 459.3121; 0 0 1];
E = K'*F*K

abs(x2h(:,inliers)'*E*x1h(:,inliers))  
det(x2h(:,inliers)'*E*x1h(:,inliers))  
  

[Uinit,Dinit,Vinit] = svd(E,0); 

d=[1,1,0];
newE = Uinit*diag(d)*Vinit';
[U D V] = svd(newE);
W = [0 -1 0; 1 0 0;0 0 1];
Z = [0 1 0;-1 0 0; 0 0 0];
S = U*Z*U'
T1 = U(:,3)
T2 = -U(:,3)
R1 = U*W*V'
R2 = U*W'*V'
P1 = [R1 T1];
P2 = [R1 T2];
P3 = [R2 T1];
P4 = [R2 T2];

P2E(: ,:, 1) = P1;
P2E(: ,:, 2) = P2;
P2E(: ,:, 3) = P3;
P2E(: ,:, 4) = P4;

P2final = Pdepthtest(P2E, K, x1h(:,inliers), x2h(:,inliers));

Angle1(1)=atan2d(R1(3,2),R1(3,3));
Angle1(2)=atan2d(-1*R1(3,1)  , sqrt(R1(3,2)^2 + R1(3,3)^2));
Angle1(3)=atan2d(R1(2,1),R1(1,1));

Angle2(1)=atan2d(R2(3,2),R2(3,3));
Angle2(2)=atan2d(-1*R2(3,1)  , sqrt(R2(3,2)^2 + R2(3,3)^2));
Angle2(3)=atan2d(R2(2,1),R2(1,1));

Anglef(1)=atan2d(P2final(3,2),P2final(3,3));
Anglef(2)=atan2d(-1*P2final(3,1)  , sqrt(P2final(3,2)^2 + P2final(3,3)^2));
Anglef(3)=atan2d(P2final(2,1),P2final(1,1));

%PA and PB camera projection matrix of two cameras

PA = eye(3,4); %[I 0]
% By result 9.14, pg. 256 p from F.
% P2 = [ [e2]*cross*F | e2 ]
e2cross = get_cross(e2); 
e2crossF = e2cross*F;
PB = [e2crossF e2]

%triangulate to get world coordinates...based on principle (x)*cross*(PX)=0.   

Xhat = worldcords(x1h,x2h,PA,PB,1);


sprintf('The Rotation difference between the two images is = %f', abs(Anglef(2)))
  

  