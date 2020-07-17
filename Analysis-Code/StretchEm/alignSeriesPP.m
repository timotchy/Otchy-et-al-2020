function [alignPP] = alignSeriesPP(ppSeries,p)
%This script aligns two series given a DTW path by averaging over the repeated DTW paths. Main inputs are
%a DTW path, p, with first column the template index and the second column is the series index, and ppSeries,
%a point process that spans the time 0:max(p(:,2)).  Returned is a row vector of point process times (of size
%ppSeries) that span the time 0:max(p(:,1)).
%
%Updated with comments 7/9/2015 TMO

%Predefine the output matrix
alignPP = [];

%Identify the parameters for warping
ppSize = length(ppSeries);      %size of the point process
tSpan = max(p(:,1));               %max span of the output series
ppSpan = max(p(:,2));            %max span of the output series

%Step through each event in the point process and remap to the new time
for i = 1:ppSize
    if ppSeries(i) >= p(end,2) % if sp time is at or beyond last rendition sample
        alignPP(i) = (ppSeries(i) - p(end,2)) + p(end,1);
        
    elseif ppSeries(i) < p(1,2) % if sp time is before first rendition sample
        alignPP(i) = (ppSeries(i) - p(1,2)) + p(1,1);
        
    else
        %Time warp based on linear interpolation between neighboring points 
        int = find(p(:,2) <= ppSeries(i), 1, 'last');
        alignPP(i) = p(int,1) + (ppSeries(i) - p(int,2)) * diff(p(int:int+1,1)) / diff(p(int:int+1,2));
        
    end
end