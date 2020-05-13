%%


%% Method 1 (Pizza): Sin/Cos Waves around Outside of Whisker
clear
close all

% basic parameters
D = 0.875 ; % 1mm
d = 0.025 ; % 10um 
s = 30 ; %arclength
A = 0.005; %curvature

% build whisker
n = 100; %along the whisker length
[z0,y0,x0] = xyz_whisker(s,A,n);
curve = [x0;y0;z0]; 

n_cir = 8; %around outside of whisker

theta = (0:n_cir)/n_cir*pi*2;
sintheta = sin(theta); sintheta(n_cir+1) = 0;
c_s=repmat(cos(theta),[3,1]);
s_s=repmat(sintheta,[3,1]);

r_osc = (D-d)/2; %base radius for oscillation
amp = 0.10; %percentage undulation
period = 0.91*2; % period of undulation 
t = 2*pi*z0/period;
taper = linspace(1,d/D, length(t));
r_maj = taper.*r_osc.*(1 + amp*sin(t)); r_maj = r_maj(:);
r_min = taper.*r_osc.*(1 + amp*cos(t))/2; r_min = r_min(:);

%deltavecs: average for internal points. first strecth for endpoitns.
dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);
nvec=zeros(3,1); %make nvec not parallel to dv(:,1)
[buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;
curve_2=repmat([0],[3,n_cir+1,n+1]);

for k = 1:n
  convec=cross(nvec,dv(:,k));
  convec=convec./norm(convec);
  nvec=cross(dv(:,k),convec);
  nvec=nvec./norm(nvec);
  curve_2(:,:,k+1)=repmat(curve(:,k),[1,n_cir+1])+...
      c_s.*repmat(r_min(k)*nvec,[1,n_cir+1])+...
      s_s.*repmat(r_maj(k)*convec,[1,n_cir+1]);
end
curve_2(:,:,1)=repmat(curve(:,1),[1,n_cir+1]);
curve_2(:,:,end-1)=repmat(curve(:,end),[1,n_cir+1]);
curve_2 = curve_2(:,:,1:end-1);

X=squeeze(curve_2(1,:,:));
Y=squeeze(curve_2(2,:,:));
Z=squeeze(curve_2(3,:,:));


alpha = 15.27*pi/180; %peak offset
beta = 17.60*pi/180; %trough offset
pi_slices = linspace(0,sin(2*alpha), length(X(:,1)));
theta_shift = repmat(sin(((beta-alpha)/2)*cos(t) +((beta-alpha)/2))-sin(alpha) ,[length(X(:,1)),1]) + repmat(pi_slices, [length(X(1,:)),1])';

temp = arrayfun(@(k) circshift(theta_shift(k,:),k),1:length(X(:,1)),'uni',0);
theta_shift_circ = [];
for ii = 1:length(X(:,1))
    theta_shift_circ(ii,:) = temp{1,ii}';
end

Z_shift = Z + theta_shift_circ;
a = find(Z_shift<0 | Z_shift>s); 
X_f = X; Y_f = Y; Z_f = Z_shift; 

X_f(a) = NaN;
Y_f(a) = NaN;
Z_f(a) = NaN;

figure(3)
subplot(2,1,1)
surf(Z_f,Y_f,X_f)
xlabel('Z (mm)')
ylabel('Y (mm)')
zlabel('X (mm)')
axis equal

% figure
% hold on
% plot(Y(1,:)/max(Y(1,:)),'-')
% plot((X(3,:)-X(7,:))/max((X(3,:)-X(7,:))),'-')
% axis square

%% Method 2 (Ellipses): Stack Ellipses along Z-axis

% basic parameters
D = 0.875 ; % 1mm
d = 0.025 ; % 10um 
s = 30 ; %arclength
A = 0.005; %curvature

r_osc = (D-d)/2; %base radius for oscillation
amp = 0.1; %percentage undulation
n = 100; %along the whisker length
period = 0.91*2; % period of undulation 
taper = linspace(1,d/D, n);

% build whisker
[z0,y0,x0] = xyz_whisker(s,A,n);
curve = [x0;y0;z0];

alpha = 15.27; %peak offset
beta = 17.60; %trough offset
shift = ((beta-alpha)/2)*cos(2*pi*curve(3,:)./period) +alpha +((beta-alpha)/2);

majAxWide = r_osc*(1 + amp*sin(2*pi*curve(3,:)./period));
majAxThin = r_osc*(1/2)*(1 + amp*cos(2*pi*curve(3,:)./period));
undWh=[];
h = figure;
hold on

for i=1:length(curve(1,:))
   undWh=[undWh; drawEllipse3d(curve(:,i)',taper(i)*majAxThin(i),taper(i)*majAxWide(i),-shift(i),0)];
end

close(h)

figure(3)
subplot(2,1,2)
plot3(undWh(:,3), undWh(:,2), undWh(:,1), '.')
xlabel('Z (mm)')
ylabel('Y (mm)')
zlabel('X (mm)')
axis equal


