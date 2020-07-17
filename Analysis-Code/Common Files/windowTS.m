function [x_mat] = windowTS(x,winSize,winStep,pad,winType)
%Function windows a time series and returns a mxn matrix where each row (m)
%is a window of length winSize. 
%winStep = the advance between consecutive windows
%pad = NaN pad x(start) by winSize/2 and x(end by up to winSize/2 ('none' the default)
%winType = determines the type of window to apply. Options include:
%   'boxcar' (the default)
%   'hanning'
%   'hamming'

%Check input structure and fill in default values
if nargin < 3
    error('Not enough input arguments.')
    return
elseif nargin == 3
    pad = 'none';
    winType = 'boxcar';
elseif nargin == 4
    winType = 'boxcar';
elseif nargin > 5
    error('Too many input arguments.')
    return
end

%Check size of array and force to row vector
[m n] = size(x);
if (m~=1 && n~=1) || (m==1 && n==1)
    error('x must be a vector')
    return
elseif m==1 && n~=1
    x = x;
elseif m~=1 && n==1
    x = x';
end

%Create window functions as specified
if strcmp(winType,'boxcar')
    winVect = rectwin(winSize)';
elseif strcmp(winType,'hanning')
    winVect = hann(winSize)';
elseif strcmp(winType,'hamming')
    winVect = hamming(winSize)';
else
    error('Invalid value for input ''winType''.')
    return
end

%Add padding to start of x, if selected
if strcmp(pad,'pad')
%     %Padding should be no more than 50% of winSize, but winStep affactes
%     %how these accumulate
%     s = floor(winSize/(2*winStep))*winStep;

    %Centers the first window on the x(start)
    s = floor(winSize/2);
    fpad = NaN(1,s);
    x_prep = [fpad, x];
elseif strcmp(pad,'none')
    x_prep = x;
else
    error('Invalid value for input ''pad''. Options are ''none'' and ''pad''.')
    return
end

%Set up indices for the start and stops of the windows
winStarts = 1:winStep:(length(x_prep)-winSize+1);
winStops = winSize:winStep:length(x_prep);

%Add padding to the end of x, if selected
if strcmp(pad,'pad')
    %Append the trailing NaN padding such that all 
    s = floor(winSize/(2*winStep))*winStep;
    extra = length(x_prep)-winStops(end); %all the leftover points
    t = s-extra;
    while t+winStep <= (winSize/2)
        t = t+winStep;
    end
    epad = NaN(1,t);
    x_prep = [x_prep, epad];
    
    %Recalculate indices for the start and stops of the windows
    winStarts = 1:winStep:(length(x_prep)-winSize+1);
    winStops = winSize:winStep:length(x_prep);
end

%Populate matrix (can this be vectorized?)
x_mat = [];
for i = 1:length(winStarts)
    %Multiply indexed segment by the windowing function
    x_mat(i,:) = x_prep(winStarts(i):winStops(i)).*winVect;
end












