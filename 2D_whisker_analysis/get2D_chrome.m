% Extract whisker from 2D scan image
%
% Follow the instruction:
% 1. select 11 bars from the ruler, which is the span of 1cm. The selected
% area shouldn't contain tip of the bar.
% 2. select whisker area, which should cover the whole whisker.
%
%
% By Yifu
% 2018/07/16
%
close all
clear
warning off
nBar = 11;
animal = 'GBL03';
whisker = 'LB06';
filename = 'Harbor7_008';
dir = '/home/luxi/Spring_2020/Whiskers/matlab/';

%% Image processing
% Load image
im0 = imread([dir,filename,'.tif']);
% isa: determine the data type,imu2unit8 convert image to integer (Note by LH)
if isa(im0,'uint16'), im0 = im2uint8(im0); end
% im0 = im0(round(end/3):round(end*3/4),round(end/4):round(end*3/4));
imshow(im0);
zoom off

% % %% code by luxi, find the mean value of G channel
% % Ig = im0(:,:,2);
% % % Extract the background (black) region
% % % Igray = rgb2gray(im0);
% % % idx = Igray == 0;
% % Gave = uint8(mean(Ig));
%  % fprintf('!!!!! %d,\n',Gave);
%% Crop image
fprintf('%s: \n',whisker)
fprintf('(1) Crop the metric now. Please include %d lines. (1cm)\n\n', nBar)
im_metric = imcrop(im0);
fprintf('(2) Crop the whisker now.\n\n')
im = imcrop(im0);
%% code by luxi, find the mean value of G channel
% % fprintf('(3) Crop the background now.\n\n')
% % im3 = imcrop(im0);
% 
% Ir = im0(:,:,1);
% Ig = im0(:,:,2);
% Ib = im0(:,:,3);
% % Extract the background (black) region
% % Igray = rgb2gray(im);
% % idx = Igray == 0;
% % Gave = uint8(mean(Ig(~idx)));
% Rave = uint8(mean(Ir));
% Gave = uint8(mean(Ig));
% Bave = uint8(mean(Ib));
% fprintf('rr!!!!! %d,\n',Rave)
% % fprintf('gg!!!!! %d,\n',Gave)
% % fprintf('bb!!!!! %d,\n',Bave)
%%
close

% Changed cropped image to vertical
if ~issorted(size(im), 'descend'), im = permute(im,[2,1,3]); end
if ~issorted(size(im_metric), 'descend'), im_metric = permute(im_metric,[2,1,3]); end

%% Fix whisker
f1 = figure('Position',[100 100 900 800]);
subplot(1,3,2); imshow(im); 
% % code by luxi, find the mean value of G channel
% fprintf('(3) Crop the background now.\n\n')
% im3 = imcrop(im); 
% Ig = im3(:,:,2);
% Extract the background (black) region
% Igray = rgb2gray(im);
% 
% %idx = Igray == 0;
% % Gave = uint8(mean(Ig(~idx)));
% Gave = uint8(mean(Ig));
% fprintf('!!!!! %d,\n',Gave);
%%

hold on
fprintf('Now please manually select the base.\n')
fprintf('You can zoom in/out before you''re ready to select the base, then press enter.\n')
fprintf('Click once, and press enter again.\n')
waitfor(f1, 'CurrentCharacter', char(13))
[xb, yb] = ginput();
zoom OUT
xb = round(xb(end)); yb = round(yb(end));
plot(xb,yb,'ro');
set(f1,'CurrentCharacter', '1')
fprintf('When you''re ready to select the tip, press enter.\n')
fprintf('Click once, and press enter again.\n')
waitfor(f1, 'CurrentCharacter', char(13))
[xt, yt] = ginput();
zoom OUT
xt = round(xt(end)); yt = round(yt(end));
plot(xt,yt,'ro');

% perform chrome key 130 110
%[im1,~] = chromeKey(im,'green',[153 115]);
[im1,~] = chromeKey(im,'green',[140 100]);
subplot(1,3,1); hold on
imshow(im1);
plot([xb,xt],[yb,yt],'ro')
% % Cut in the middle
% set(f1,'CurrentCharacter', '1')
% waitfor(f1, 'CurrentCharacter', char(13))
% [xm, ym] = ginput();


%% new line add by luxi
% change white color area to black, and change black area to white 
im1 = bwareaopen(im1, 8000);
im1 = imcomplement(bwareaopen(imcomplement(im1),8000));
%% Generate whisker
radius = 40;
x = zeros(abs(yb-yt)+1,1);
y = zeros(abs(yb-yt)+1,1);
center = xb;
jj = 1;
for j = yb:(-1)^(yb>yt):yt
%     if j <= ym(1) || j >= ym(2)
%         if j == ceil(ym(2))
%             center = xm(2); 
%         end
        y(jj) = j;
        x(jj) = center - radius + expIndex(im1(j,center-radius:center+radius));
        center = x(jj);
        jj = jj + 1;
%     end  
end
subplot(1,3,3);
imshow(im1); 
hold on;
%% 
% plot(x,y,'g-',x(1),y(1),'co');
% suptitle('Processed <- origin -> Traced')
%% newline by luxi (find x and y area along the perpendicular line)
% smooth x and y value 
x1 = x;
y1 = y;
x = smooth(x,2000,'loess');
y = smooth(y,2000,'loess');
% find the perpendicular line
% fprintf ('%d,\n', y(6000));
radius1 = 60;
% array1 = [];
array_value = []; % white color area along the perpendicular line.
lower_line_x = [];
lower_line_y = [];
upper_line_x = [];
upper_line_y = [];
new_center_line_x = [];
count = 0;
for j = 1:2:length(y)-10
    count = count +1;
    k = (y(j)-y(j+1))/(x(j)-x(j+1));
    b =  y(j)-x(j)*k;
    k1 = (x(j)-x(j+1))/(y(j+1)-y(j));
    b1 = y(j)-x(j)*k1;
    sum_value =0; 
    check_value1 = 0;
    check_value2 = 0; 
    % find the lowest white pixel point on perpendicualr line 
    for m = x(j)-radius1:x(j)+radius1;
% %     for m = x(j)-radius1:x(j);
        m = round(m);
        y1 = m*k1+b1; % function of perpendicular line. 
        y1 = round(y1);
%         fprintf('%d,\n',m);
        sum_value = im1(y1,m) + sum_value;  
        % find the lower line position:
        if (check_value1 == 0) && (im1(y1,m) == 1)
            lower_line_x (count) = m;
            lower_line_y (count) = y1;
            check_value1 = 1;
            
        end
            
    end
    % find the highest white pixel point on perpendicualr line 
    for m = x(j)+radius1:-1:x(j)-radius1;
        m = round(m);
        y1 = m*k1+b1; % function of perpendicular line. 
        y1 = round(y1);
        
        y2 = (m +1)*k1+b1;
        y2 = round(y2);
        m2 = m +1;
        if (check_value2 == 0) && (im1(y1,m) == 1) && (im1(y2,m2) == 0)
            check_value2 = 1;
            upper_line_x (count) = m;
            upper_line_y (count) = y1;            
        end
        
    end    

    new_center_line_x(count) = (upper_line_x(count) + lower_line_x(count))/2.0;
    new_center_line_y(count) = (upper_line_y(count) + lower_line_y(count))/2.0;
    
%     fprintf('%d,\n',sum_value);
    array_value(count) = sum_value;
end 
%% smooth data by luxi
% new_center_line_x_smooth = smooth(new_center_line_x,2000,'rloess');
% new_center_line_y_smooth = smooth(new_center_line_y,2000,'rloess'); 
%% white color area on the new centerline by luxi
radius2 = 50;
% array1 = [];
new_array_value_upper_line = [];% white color area along the perpendicular line.
new_array_value_lower_line = [];
new_array_value = [];
count = 0;
new_y = new_center_line_y;
new_x = new_center_line_x;

for j = 1:length(new_y)-11
    count = count +1;
    k = (new_y(j)-new_y(j+5))/(new_x(j)-new_x(j+5));
    b =  new_y(j)-new_x(j)*k;
    k1 = -1/k;
    b1 = new_y(j+1)-new_x(j+1)*k1;
    fprintf('k1: %f,\n',k1);
    fprintf('b1: %f,\n',b1);
    sum_value1 =0; 
    sum_value2 =0;
%     for m = new_x(j)-radius1:new_x(j)+radius1;
    for m = new_x(j+3)-radius1:new_x(j+3);
        y1 = m*k1+b1; % function of perpendicular line. 
%         fprintf('m: %d,\n',m);
%         fprintf('y1!!!: %f,\n',y1);
        y1 = round(y1);
        m = round(m);
%         fprintf('m: %d,\n',m);
%         fprintf('y1~~~: %f,\n',y1);
%         fprintf('im1(y1,m): %d,\n',im1(y1,m));
        sum_value1 = im1(y1,m) + sum_value1;
            
    end
    
%     for m = new_x(j)+radius1:-1:new_x(j)-radius1;

    for m = new_x(j+3):new_x(j+3)+radius1;
        %for m = x(j)-radius1:x(j);
        m = round(m);
        y1 = m*k1+b1; % function of perpendicular line. 
        y1 = round(y1);
%         fprintf('%d,\n',m);
        sum_value2 = im1(y1,m) + sum_value2;
    end    
    new_array_value_lower_line(count) = sum_value1;
    new_array_value_upper_line(count) = sum_value2;
    
    new_array_value(count) = sum_value1 + sum_value2; 
end 
%% modified by luxi, draw center lines 
% plot(x,y,'g-',x(1),y(1),'co',new_center_line_x, new_center_line_y,'r', new_center_line_x_smooth, new_center_line_y_smooth,'b');
plot(x(1),y(1),'co',new_center_line_x, new_center_line_y,'r');

% plot(x,y,'g-',x(1),y(1),'co');
suptitle('Processed <- origin -> Traced')
%% plot by luxi
% x_value = 1:10000;
fig2 = figure('Position',[1050,100,600,400]); hold on
% plot(x_value(1:numel(y_value)));
plot(new_array_value);
%array_value_smooth = smooth(array_value,150,'rloess');
% get smooth line convex array
array_value_smooth = smooth(new_array_value,200,'rloess');
convex_xx = []; 
convex_yy = [];
count = 0;
for i = 2: 1: length(array_value_smooth)-1
%     fprintf('inside loop %d\n', i);
    if (array_value_smooth(i)< array_value_smooth(i+1)) && (array_value_smooth(i)<array_value_smooth(i-1))
        count = count +1;
        convex_xx(count) = i;
        convex_yy(count) = array_value_smooth(i); 
%         fprintf('x value is %f\n',convex_xx(count));
    end     
end
% 
plot(array_value_smooth);
plot(convex_xx , convex_yy, 'go');
grid on

%% Code by luxi
% compare convex point
x_value = 0; 
peaks_xx = [];
peaks_yy = [];

count = 1;
peaks_xx(count) = convex_xx(1);
peaks_yy(count) = convex_yy(1)
for i = 2: 1: length(convex_xx)
    if convex_yy(i) < convex_yy(i-1)
        count = count + 1;
        peaks_xx(count) = convex_xx(i);
%         fprintf('real!!! %d\n',peaks_xx(count));
        peaks_yy(count) = convex_yy(i)
    else
        peaks_xx = [];
        peaks_yy = [];
        count = 1;
        peaks_xx(count) = convex_xx(i);
%         fprintf('update!!! %d\n',peaks_xx(count));
        peaks_yy(count) = convex_yy(i)        
    end 
end
plot(peaks_xx , peaks_yy, 'r*');

%% lower_line_plot by luxi
% x_value = 1:10000;
fig5 = figure('Position',[1050,100,600,400]); hold on
plot(new_array_value_lower_line);
%array_value_smooth = smooth(array_value,150,'rloess');
% get smooth line convex array
new_array_value_lower_line_smooth = smooth(new_array_value_lower_line,200,'rloess');
convex_xx = []; 
convex_yy = [];
count = 0;
for i = 2: 1: length(new_array_value_lower_line_smooth)-1
%     fprintf('inside loop %d\n', i);
    if (new_array_value_lower_line_smooth(i)< new_array_value_lower_line_smooth(i+1)) && (new_array_value_lower_line_smooth(i)<new_array_value_lower_line_smooth(i-1))
        count = count +1;
        convex_xx(count) = i;
        convex_yy(count) = new_array_value_lower_line_smooth(i); 
%         fprintf('x value is %f\n',convex_xx(count));
    end     
end
% 
plot(new_array_value_lower_line_smooth);
plot(convex_xx , convex_yy, 'go');
grid on

%% Code by luxi
% compare convex point
x_value = 0; 
peaks_xx = [];
peaks_yy = [];

count = 1;
peaks_xx(count) = convex_xx(1);
peaks_yy(count) = convex_yy(1)
for i = 2: 1: length(convex_xx)
    if convex_yy(i) < convex_yy(i-1)
        count = count + 1;
        peaks_xx(count) = convex_xx(i);
%         fprintf('real!!! %d\n',peaks_xx(count));
        peaks_yy(count) = convex_yy(i)
    else
        peaks_xx = [];
        peaks_yy = [];
        count = 1;
        peaks_xx(count) = convex_xx(i);
%         fprintf('update!!! %d\n',peaks_xx(count));
        peaks_yy(count) = convex_yy(i)        
    end 
end
plot(peaks_xx , peaks_yy, 'r*');


%% upper line plot by luxi
% x_value = 1:10000;
fig6 = figure('Position',[1050,100,600,400]); hold on
plot(new_array_value_upper_line);
%array_value_smooth = smooth(array_value,150,'rloess');
% get smooth line convex array
new_array_value_upper_line_smooth = smooth(new_array_value_upper_line,200,'rloess');
convex_xx = []; 
convex_yy = [];
count = 0;
for i = 2: 1: length(new_array_value_upper_line_smooth)-1
%     fprintf('inside loop %d\n', i);
    if (new_array_value_upper_line_smooth(i)< new_array_value_upper_line_smooth(i+1)) && (new_array_value_upper_line_smooth(i)<new_array_value_upper_line_smooth(i-1))
        count = count +1;
        convex_xx(count) = i;
        convex_yy(count) = new_array_value_upper_line_smooth(i); 
%         fprintf('x value is %f\n',convex_xx(count));
    end     
end
% 
plot(new_array_value_upper_line_smooth);
plot(convex_xx , convex_yy, 'go');
grid on

%% Code by luxi
% compare convex point
x_value = 0; 
peaks_xx = [];
peaks_yy = [];

count = 1;
peaks_xx(count) = convex_xx(1);
peaks_yy(count) = convex_yy(1)
for i = 2: 1: length(convex_xx)
    if convex_yy(i) < convex_yy(i-1)
        count = count + 1;
        peaks_xx(count) = convex_xx(i);
%         fprintf('real!!! %d\n',peaks_xx(count));
        peaks_yy(count) = convex_yy(i)
    else
        peaks_xx = [];
        peaks_yy = [];
        count = 1;
        peaks_xx(count) = convex_xx(i);
%         fprintf('update!!! %d\n',peaks_xx(count));
        peaks_yy(count) = convex_yy(i)        
    end 
end
plot(peaks_xx , peaks_yy, 'r*');
%% Generate metric
bw_metric = imbinarize(im_metric);
[y_metric, x_metric] = find(bw_metric == 0);
% Initial location for the centroids
C_start = [ones(nBar,1)*size(bw_metric,1), linspace(min(y_metric), max(y_metric), nBar)'];
[~, C] = kmeans([x_metric, y_metric], nBar, 'Start', C_start);
fig3 = figure('Position',[1050,500,600,400]); hold on
plot(x_metric, y_metric,'k.'); 
plot(C(:,1),C(:,2),'w*')

% The pixel length of 1 mm
[~, iC] = sort(C(:,2));
C = C(iC,:);

% distance
mm1 = mean(sqrt(sum((C(2:end,:) - C(1:end-1,:)).^2,2)));
% vertical
mm2 = mean(C(2:end,2) - C(1:end-1,2));


%%
x = x(x~=0)/mm2;
y = y(y~=0)/mm2;
Px = smooth(x,2000,'loess');
Py = smooth(y,2000,'loess');

fig4 = figure('Position',[1050,100,600,400]); hold on
plot(x,y,'k-',Px,Py,'r.',Px(1),Py(1),'co');
[Pxstd, Pystd, th] = stdOrientation2D(Px', Py', 'up');
plot(Pxstd, Pystd, 'r-', Pxstd(1), Pystd(1), 'co')
text(Pxstd(end),Pystd(end),'Std')
axis equal
grid on
fprintf('Whisker length: %.2fmm\n', whiskerLength([Px, Py]'));
eval([whisker,'=[Pxstd;Pystd]'';']);
% save([dir,animal,' Extractions/',animal,'_',whisker,'.mat'],whisker)







