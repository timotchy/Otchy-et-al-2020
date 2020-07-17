% This script aligns two series given a DTW path by averaging over the
% repeated DTW paths. Given is the tSeries, and DTW path with first column template
% index and the second column template is the series index. Returned is a
% row vector of size of template with entries that map to that point.

function [path alignV] = alignSeriesPath(tSeries,p)

tSize = max(p(:,1));
tStartInd = p(1,2);
[a b] = size(tSeries);
alignV = zeros(a,tSize);

path =[];

[pX pY] = size(p);
alignV = zeros(a,tSize);

DTWcount = 1;

for i = 1:tSize
    
   
   count = 0;
   dum = 0;
   pVector = [];
   
   while p(DTWcount,1) == i
       
       pVector = [pVector (p(DTWcount,2)-tStartInd+1)];
       dum = dum + tSeries(:,p(DTWcount,2));
       count = count+1;
       DTWcount = DTWcount+1;
   
       if DTWcount > pX
           break;
       end
   end
    
   path(i) = mean(pVector);
   alignV(:,i) = dum/count;
    
end






end