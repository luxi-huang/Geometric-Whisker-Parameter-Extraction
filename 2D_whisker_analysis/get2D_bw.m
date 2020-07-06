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

close all
clear
warning off
nBar = 11;
animal = 'SLN01';
whisker = 'LC03';
filename = 'img005';
dir = '..\_Raw_Data\SLN\2D\SLN01\';


%% Image processing
% Load image
im0 = imread([dir,filename,'.tif']);
if isa(im0,'uint16'), im0 = im2uint8(im0); end
if size(im0,3) == 3, im0 = rgb2gray(im0); end
% im0 = im0(round(end/3):round(end*3/4),round(end/4):round(end*3/4));
imshow(im0);
zoom off
% Crop image
fprintf('(1) Crop the metric now. Please include %d lines. (1cm)\n\n', nBar)
im_metric = imcrop(im0);
fprintf('(2) Crop the whisker now.\n\n')
im = imcrop(im0);
close

% Changed cropped image to vertical
if ~issorted(size(im), 'descend'), im = im'; end
if ~issorted(size(im_metric), 'descend'), im_metric = im_metric'; end


%% Fix whisker
fig1 = figure('Position',[100 100 900 800]);
subplot(1,3,2); imshow(im); hold on

fprintf('%s: \n',whisker)
fprintf('Now please manually select the base.\n')
fprintf('You can zoom in/out before you''re ready to select the base, then press enter.\n')
fprintf('Click once, and press enter again.\n')
waitfor(fig1, 'CurrentCharacter', char(13))
[xb, yb] = ginput();
zoom OUT
xb = round(xb(end)); yb = round(yb(end));
plot(xb,yb,'ro');
set(fig1,'CurrentCharacter', '1')
fprintf('When you''re ready to select the tip, press enter.\n')
fprintf('Click once, and press enter again.\n')
waitfor(fig1, 'CurrentCharacter', char(13))
[xt, yt] = ginput();
zoom OUT
xt = round(xt(end)); yt = round(yt(end));
plot(xt,yt,'ro');

% Change image to white background
iswhite = mean(mean(im)) > 255/3;
if ~iswhite, im = 255 - im; end

rowmin = min(im,[],2);
rowmax = max(im,[],2);
%rescale eceryrow from [min, max] to [0, 255], and clip at 127
im1 = int8(255*rescale(double(im),'InputMin',rowmin,'InputMax',rowmax));
%rescale [0, 127] to [0,255], and convert to uint8
im1 = uint8(rescale(im1,0,255));
im1 = 255 - im1;

% im1 = 255 - im;
subplot(1,3,1); imshow(im1); hold on
plot([xb,xt],[yb,yt],'ro')


% % Cut in the middle
cut = 0;
if cut
    set(fig1,'CurrentCharacter', '1')
    waitfor(fig1, 'CurrentCharacter', char(13))
    [xm, ym] = ginput();
end



%% Generate whisker
% radius = uint16(size(im,2)/40);
radius = 100;
% if whisker is horizontal-ish

x = zeros(abs(yb-yt)+1,1);
y = zeros(abs(yb-yt)+1,1);
center = xb;
jj = 1;
for j = yb:(-1)^(yb>yt):yt
    if cut, if j <= ym(1) || j >= ym(2)
        if j == floor(ym(1))
            center = xm(1); 
        end
    end; end
    y(jj) = j;
    x(jj) = center - radius + expIndex(im1(j,center-radius:center+radius));
    center = x(jj);
    jj = jj + 1;
end

subplot(1,3,3);
imshow(im1); hold on;
plot(x,y,'g-',x(1),y(1),'co')
suptitle('Processed  <-  origin  ->  Traced')

%% Generate metric
bw_metric = imbinarize(im_metric);
[y_metric, x_metric] = find(bw_metric == 0);
% Initial location for the centroids
if ~issorted(size(bw_metric),'ascend')
    C_start = [ones(nBar,1)*size(bw_metric,2)/2, linspace(min(y_metric), max(y_metric), nBar)'];
else
    C_start = [linspace(min(x_metric), max(x_metric), nBar)', ones(nBar,1)*size(bw_metric,1)/2];
end
% k-means
[~, C] = kmeans([x_metric, y_metric], nBar, 'Start', C_start);
fig2 = figure('Position',[1050,500,600,400]); hold on
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
% [Px, Py] = equidistlen(x,y,0.5);
Px = smooth(x,2000,'loess');
Py = smooth(y,2000,'loess');

fig3 = figure('Position',[1050,100,600,400]); hold on
plot(x,y,'k-',Px,Py,'r.',Px(1),Py(1),'co');
[Pxstd, Pystd, th] = stdOrientation2D(Px', Py', 'up');
plot(Pxstd, Pystd, 'r-', Pxstd(1), Pystd(1), 'co')
text(Pxstd(end),Pystd(end),'Std')
axis equal
grid on
fprintf('Whisker length: %.2fmm\n', arclength(Px,Py));

% g = fittype(@(a,b,x) a*x.^3+b*x.^2,...
%         'coefficients', {'a','b'});
% b = fit(Pxstd',Pystd',g,'StartPoint',[0.01 0.01]);
% Pxstdfit = linspace(0,max(Pxstd),1000);
% Pystdfit = b.a*Pxstdfit.^3 + b.b*Pxstdfit.^2;
% plot(Pxstdfit,Pystdfit,'b-')
% eval([whisker,'=[Pxstdfit;Pystdfit]'';']);
% save([pwd,'/',animal,'/',whisker,'.mat'],whisker)


%%
eval([whisker,'=[Pxstd;Pystd]'';']);
% save([dir,animal,' Extractions/',animal,'_',whisker,'.mat'],whisker)
