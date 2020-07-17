function [leadingEdgeNdx, fallingEdgeNdx] = detectThresholdCrossings(sig, fThres, bAbove)
%Returns indices of values above or below threshold.  LeadingEdge indices
%are those that first cross the thres.  FallingEdge contains the indices of
%the last value to be above the thres.

leadingEdgeNdx = [];
fallingEdgeNdx = [];
if(bAbove)
    exceedsThres = find(sig>fThres);
else
    exceedsThres = find(sig<fThres);
end
    
if(length(exceedsThres)>0)
    ndx = find(diff(exceedsThres)>1);
    
    leadingEdgeNdx = exceedsThres(ndx + 1);
    leadingEdgeNdx = [exceedsThres(1); leadingEdgeNdx];
    fallingEdgeNdx = exceedsThres(ndx);  
    fallingEdgeNdx = [fallingEdgeNdx; exceedsThres(end)];
end
    
    
   