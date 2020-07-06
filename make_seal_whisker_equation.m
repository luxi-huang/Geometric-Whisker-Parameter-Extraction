%%


%% Method 1 (Pizza): Sin/Cos Waves around Outside of Whisker
clear
close all

% basic parameters
D = 0.875 ; % 1mm
d = 0.025 ; % 10um 
s = 30 ; %arclength
A = 0.0025; %curvature

% build whisker
n = 1000; %along the whisker length
[z0,y0,x0] = xyz_whisker(s,A,n);
curve = [x0;y0;z0]; 
%curve = [zeros(1,n); zeros(1,n); linspace(0,s,n)];

n_cir = 8; %around outside of whisker

theta = 2*pi*((0:n_cir)/(n_cir));
sinhalftheta = cos(theta*2); sinhalftheta(n_cir+1) =1; 
sintheta = sin(theta); sintheta(n_cir+1) = 0;
c_s=repmat(cos(theta),[3,1]);
s_s=repmat(sintheta,[3,1]);

r_osc = (D-d)/2; %base radius for oscillation
amp = 0.20; %percentage undulation
period = 0.91*2; % period of undulation 
t = 2*pi*z0/period;
taper = linspace(1,d/D, length(t));
r_maj = taper.*r_osc.*(1 + amp*sin(t)); r_maj = r_maj(:);
r_min = 0.75*taper.*r_osc.*(1 - amp*sin(t)); r_min = r_min(:);

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
  
  r_step = (r_maj(k)+r_min(k))/2;
  amp_shift = [repmat([r_maj(k), r_step, r_min(k),r_step],[1,2]),r_maj(k)];
%   shift = circshift(amp_shift,k*period/(2*pi))
  
  curve_2(:,:,k+1)=repmat(curve(:,k),[1,n_cir+1])+...
        c_s.*amp_shift.*nvec+... 
        s_s.*amp_shift.*convec;
%         c_s.*repmat(r_maj(k)*nvec,[1,n_cir+1])+...
%         s_s.*repmat(r_min(k)*convec,[1,n_cir+1]);

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

figure(1)
%subplot(2,1,1)
surf(Z_f,Y_f,X_f); %, 'EdgeColor', 'none')
xlabel('Z (mm)')
ylabel('Y (mm)')
zlabel('X (mm)')
axis equal

figure(2)
%subplot(2,1,1)
surf(Z,Y,X)
xlabel('Z (mm)')
ylabel('Y (mm)')
zlabel('X (mm)')
axis equal

%%
Y_ff = []; Z_ff = []; X_ff = [];
for ii = 1:length(Y_f(:,1))
    h = fitlm(1:n, Y_f(ii,:)); 
    g = fitlm(1:n, Z_f(ii,:), 'quadratic');
    f = fitlm(1:n, X_f(ii,:), 'quadratic');
    yy = h.Coefficients.Estimate(2)*[1:n] + h.Coefficients.Estimate(1); 
    zz = g.Coefficients.Estimate(3)*[1:n].^2 + g.Coefficients.Estimate(2)*[1:n] + g.Coefficients.Estimate(1); 
    xx = f.Coefficients.Estimate(3)*[1:n].^2 + f.Coefficients.Estimate(2)*[1:n] + f.Coefficients.Estimate(1); 
    Y_ff(ii,:) = Y_f(ii,:) - yy; 
    Z_ff(ii,:) = Z_f(ii,:);% - zz(1,1000); 
    X_ff(ii,:) = X_f(ii,:) - xx; 
end
close all
figure(5)
hold on 
ax = gca;
ax.ColorOrder = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;];
%    0.4660    0.6740    0.1880;];
ax.LineStyleOrder = {'-','--'};
%plot(X_ff([1,2,4,5],:)'); plot(X_f([1,2,4,5],:)'); 
%ylabel('Z (mm)')
plot(Y_ff([1,2,4,5],:)'); plot(Y_f([1,2,4,5],:)');
ylabel('Y (mm)')
set(gcf,'color','w')
xlabel('n')
%xlim([50,850])
%%
xlim([50,150])
%%
figure(3)
%subplot(2,1,1)
surf(Z_f,Y_f,X_ff)
xlabel('Z (mm)')
ylabel('Y (mm)')
zlabel('X (mm)')
axis equal
%%
close all
figure(8)
hold on
ax = gca;
ax.ColorOrder = [
         0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    %0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    %0.6350    0.0780    0.1840;
    0           0           0;
    %      0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    %0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0           0           0;];
plot((Y_ff([1,2,4,5,6,8],:)')); 
plot((X_ff([2,3,4,6,7,8],:)' + X_f([2,3,4,6,7,8],55)')); 
legend(string([1,2,4,5,6,8,2,3,4,6,7,8]))

% ax.ColorOrder = [
%          0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560;
%     0.4660    0.6740    0.1880;
%     0.3010    0.7450    0.9330;
%     0.6350    0.0780    0.1840;
%     0           0           0;
%           0    0.4470    0.7410;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560;
%     0.4660    0.6740    0.1880;
%     0.3010    0.7450    0.9330;
%     0.6350    0.0780    0.1840;
%     0           0           0;];
% plot((Y_f([1:8],:)')); 
% plot(X_ff([1:8],:)' + X_ff([1:8],55)'); 
%legend(string([1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8]))
set(gcf,'color','w')
xlabel('n')
ylabel('amp (mm)')
%ylim([-0.1,0.1])
%xlim([55,150])


%%
close all
figure(8)
hold on
ax = gca;
ax.ColorOrder = [
         0    0.4470    0.7410;
    0.4660    0.6740    0.1880;
    0.9290    0.6940    0.1250;
    0.6350    0.0780    0.1840;  ];
plot(Z_ff([1,5],:)',(Y_ff([1,5],:)')); 
plot(Z_ff([1,5],:)',((X_ff([3,7],:)'))); % + X_f([3,7],55)'))); 
legend(string([1,5,3,7]))
set(gcf,'color','w')
xlabel('Z (mm)')
ylabel('')
%ylim([0.005,0.07])
xlim([8,10])

%%
figure(9)
hold on
plot(X_f(1:n_cir/4,:)')
xlim([0,150])
legend(string([1:n_cir/4]))

% figure(6)
% hold on
% for ii = 1:length(Y_ff(:,1))
%     plot3(Z_ff(ii,:),Y_ff(ii,:),X_ff(ii,:))
% end

figure(7)
hold on
% for ii = [1,2,3,4,5,6,7,8,9] %,1:length(Y_f(:,1))
%     plot3(abs(Z_f(ii,:)),abs(Y_f(ii,:)),abs(X_f(ii,:)))
% end
% for ii = [1,4] %n_cir/4] %,1:length(Y_f(:,1))
%     %plot3(abs(Z(ii,:)),abs(Y(ii,:)),abs(X(ii,:)))
%     plot(abs(Z(ii,:)),abs(Y(ii,:)))
%     plot(abs(Z(ii,:)),abs(X(ii,:)))
% end
% for ii = [3,4,7,8] %,1:length(Y_f(:,1))
%     plot3(abs(Z_f(ii,:)),abs(Y_f(ii,:)),abs(X_f(ii,:)))
% end
plot((Z_ff(1,:)),(Y_ff(1,:)))
plot((Z_ff(4,:)),(X_ff(4,:)))
%xlim([0,4.5])
xlabel('Z (mm)')
ylabel('(mm)')
set(gcf,'color','w')
%legend('3','4','7','8')
legend(string([1:n_cir/2]))

%% Method 2 (Ellipses): Stack Ellipses along Z-axis

% basic parameters
D = 0.875 ; % 1mm
d = 0.025 ; % 10um 
s = 30 ; %arclength
A = 0.0025; %curvature

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


