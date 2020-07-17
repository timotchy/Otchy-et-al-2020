% This script aligns two series given a DTW path by averaging over the
% repeated DTW paths. Given is the tSeries, and DTW path with first column template
% index and the second column template is the series index. Returned is a
% row vector of size of template with entries that map to that point.

function [alignV] = alignSeries(tSeries,p)

tSize = max(p(:,1));
[a b] = size(tSeries);
alignV = zeros(a,tSize);


[pX pY] = size(p);
alignV = zeros(a,tSize);

DTWcount = 1;

for i = 1:tSize
    
   
   count = 0;
   dum = 0;
   
   
   while p(DTWcount,1) == i
       
       dum = dum + tSeries(:,p(DTWcount,2));
       count = count+1;
       DTWcount = DTWcount+1;
   
       if DTWcount > pX
           break;
       end
   end
    
   alignV(:,i) = dum/count;
    
end






end