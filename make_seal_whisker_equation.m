for ii = 1:60
    x(ii) = file(1,ii*3 -2)
    y(ii) = file(1,ii*3 -1)
    z(ii) = file(1,ii*3 -0)
end

hold on 
plot(x)
plot(y)
plot(z)

%%
clf
clf
t = linspace(0,100);
r = sin(-linspace(0,2*pi)) +1 ;
th = linspace(0,2*pi);
h = cos(linspace(0,pi)) ;
j = linspace(0,1);
[R,TH] = meshgrid(r,th);
[H,J] = meshgrid(h,j); 
X = R.*cos(TH) * sin(H) ;
Y = R.*sin(TH) * sin(H);
Z = R ;
figure(1)
surf(X,Y,Z)
xlabel('X')
ylabel('Y')
zlabel('Z')

figure(2)
hold on
plot(r)
plot(th)
plot(r.*cos(th).*(1-h))

%%
D = 1 ; % 1mm
d = 0.010 ; % 10um 
s = 30 ; %arclength
r_osc = (D-d)/2; 
amp = 0.095/2; 
n = 100;
n_cir = 8;
period = 10*2;
t = linspace(0,period*pi, n);
taper = linspace(1,0, length(t));
r = taper.*(r_osc + amp*sin(t)); r = r(:);
m = length(r);
theta = (0:n_cir)/n_cir*2*pi;
sintheta = sin(theta); sintheta(n_cir+1) = 0;

y0 = zeros(1,n);
z0 = linspace(0,s,n);
x0 = 0.005*z0.^2; 
curve = [x0;y0;z0]; 
c_s=repmat(cos(theta),[3,1]);
s_s=repmat(sintheta,[3,1]);


%deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);
  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;
  curve_2=repmat([0],[3,n_cir+1,n+2]);

  %r=1;
  for k = 1:n
      convec=cross(nvec,dv(:,k));
      convec=convec./norm(convec);
      nvec=cross(dv(:,k),convec);
      nvec=nvec./norm(nvec);
      curve_2(:,:,k+1)=repmat(curve(:,k),[1,n_cir+1])+...
          c_s.*repmat(r(k)*nvec,[1,n_cir+1])+...
          s_s.*repmat(r(k)*convec,[1,n_cir+1]);
  end
 %finally, cap the ends:
  curve_2(:,:,1)=repmat(curve(:,1),[1,n_cir+1]);
  curve_2(:,:,end)=repmat(curve(:,end),[1,n_cir+1]);
  
  %,extract results:
  X=squeeze(curve_2(1,:,:));
  Y=squeeze(curve_2(2,:,:));
  Z=squeeze(curve_2(3,:,:));
  
% X = r.*cos(theta).*nvec + curve(1,k);
% Y = r.*sintheta.*convec + curve(2,k);
% Z = (0:m-1)'/(m-1) * ones(1,n_cir+1) + curve(3,k);


%Z = Z*s; 
figure
surf(X,Y,Z)
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal

figure
hold on
plot3(squeeze(curve_2(1,1,:)), squeeze(curve_2(2,1,:)),squeeze(curve_2(3,1,:)),'-')
plot3(squeeze(curve_2(1,2,:)), squeeze(curve_2(2,2,:)),squeeze(curve_2(3,2,:)),'-')
plot3(squeeze(curve_2(1,3,:)), squeeze(curve_2(2,3,:)),squeeze(curve_2(3,3,:)),'-')
axis square
%%

n = 100; 
y0 = zeros(1,n/2);
z0 = linspace(0,50,n/2);
x0 = 0.005*z0.^2; 
curve = [x0;y0;z0];
r = 1;
n =8; 
%ct=0.5*r;

clf
% Usage: [x,y,z]=tubeplot(curve,r,n,ct)
% 
% Tubeplot constructs a tube, or warped cylinder, along
% any 3D curve, much like the build in cylinder function.
% If no output are requested, the tube is plotted.
% Otherwise, you can plot by using surf(x,y,z);
%
% Example of use:
% t=linspace(0,2*pi,50);
% tubeplot([cos(t);sin(t);0.2*(t-pi).^2],0.1);
% daspect([1,1,1]); camlight;
%
% Arguments:
% curve: [3,N] vector of curve data
% r      the radius of the tube
% n      number of points to use on circumference. Defaults to 8
% ct     threshold for collapsing points. Defaults to r/2 
%
% The algorithms fails if you have bends beyond 90 degrees.
% Janus H. Wesenberg, july 2004
  
  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct;
      npoints=npoints+1;
      curve(:,npoints)=curve(:,k);
    end
  end
  %Always include endpoint
  if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
  end
  
  t = linspace(0,period*pi, npoints);
  r_2= sin(t);
  
  
  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);
  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;
  xyz=repmat([0],[3,n+1,npoints+2]);
  
  %precalculate cos and sing factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
  
  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r_2(k)*nvec,[1,n+1])...
        +sfact.*repmat(r_2(k)*convec,[1,n+1]);
  end;
  
  %finally, cap the ends:
  xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
  xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
  
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  
  %... and plot:
  surf(x,y,z);
 xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal

