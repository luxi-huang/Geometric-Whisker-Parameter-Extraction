function [foreground, mask] = chromeKey(im,varargin)
%chromeKey takes a unit8 rgb image as input, and filters out the
%background specified as R/G(default)/B. The function outputs the foreground
%and mask(background) as a logical matrix.
%
%   [foreground, mask] = chromeKey(im) % case 1 
%   [foreground, mask] = chromeKey(im, CHANNEL) % case 2
%   [foreground, mask] = chromeKey(im, CHANNEL, bounds)  takes a unit8 rgb
%   % case three 
%   image and returns the foreground and background as a logical matrix.
%   The CHANNEL should be specified as a string containing R/G(default)/B 
%   or r/g(default)/b. Mask is supposed to have pixels which are stronger
%   than bounds(1) in specified channel, and weaker than bounds(2) in other
%   channels.
%
%   Example: [foreground, mask] = chromeKey(im,'green',[130 110])
%
% By Yifu
% 2018/04/18

switch nargin
    case 1
        channel = 'g';
        bounds = [130 110];
    case 2
        %  "ischar" determine if input is character array; 
        if ~ischar(varargin{1}), error('Please specify a channel for filtering.'); end
        if contains(varargin{1},'g') || contains(varargin{1},'G')
            channel = 'g';
        elseif contains(varargin{1},'r') || contains(varargin{1},'R')
            channel = 'r';
        elseif contains(varargin{1},'b') || contains(varargin{1},'B')
            channel = 'b';
        end
        bounds = [130 110];
    case 3
        if ~ischar(varargin{1}), error('Please correctly specify a channel for filtering.'); end
        if contains(varargin{1},'g') || contains(varargin{1},'G')
            channel = 'g';
        elseif contains(varargin{1},'r') || contains(varargin{1},'R')
            channel = 'r';
        elseif contains(varargin{1},'b') || contains(varargin{1},'B')
            channel = 'b';
        end
        % isvector determine it is vector
        if ~isvector(varargin{2}), error('Please correctly specify the filtering bounds'); end
        bounds = varargin{2};
    otherwise
        error('Too many input arguments.')
end


%% apply chrome key trick to fidderent color channels
imr = im(:,:,1);  img = im(:,:,2);  imb = im(:,:,3); 
if strcmp(channel,'g')
    mask = imr < bounds(2) & img > bounds(1) & imb < bounds(2);
elseif strcmp(channel,'r')
    mask = imr > bounds(1) & img < bounds(2) & imb < bounds(2);
elseif strcmp(channel,'b')
    mask = imr < bounds(2) & img < bounds(2) & imb > bounds(1);
end
foreground = ~mask;

end

