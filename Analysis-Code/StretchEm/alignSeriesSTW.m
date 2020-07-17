function [alignV] = alignSeriesSTW(tSeries,p)
% This script aligns two series given a DTW path by averaging over the
% repeated DTW paths. Given is the tSeries, and DTW path with first column template
% index and the second column template is the series index. Returned is a
% row vector of size of template with entries that map to that point.
%
%Updated with comments 7/8/2015 TMO

%Predefine the output matrix
alignV = [];

%Check for the size of the output matrix (not necessary for the script to run, just diagnostic)
tSize = max(p(:,1));       %size in the warping dimension (i.e., time)
[a, b] = size(tSeries);     %size in the orthogonal dimension (i.e., frequency for spectrograms)

%Determine the number of steps in the alignement
[pX, ~] = size(p);

%Cycle through steps in the path (save the last), linearly warping between defined anchors
for i = 1:pX-1
    %Warp between adjacent anchor points
    dum = linearWarp(tSeries(:,p(i,2):p(i+1,2)),p(i+1,1)-p(i,1)+1);
    
    %Build up the warped series from the individual warped snips
    alignV = [alignV dum(:,1:end-1)];
end

%The last point in the aligned series just is the last point of the input series
alignV = [alignV tSeries(:,p(pX,2))];

end