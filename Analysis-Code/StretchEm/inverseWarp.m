function [futureStr,basepoints,futurepoints] = inverseWarp(paths,future_mNeuro)
%This script will take in the mean path for alignments from future to base
%and apply that path to the future trace.  The result is that the temporal
%structure of the original future renditions can be recovered.

paths.pathMean;
templatesyllBreaks = paths.templatesyllBreaks;
future_mNeuro;

halfpath = round((paths.pathMean + (1:length(paths.pathMean))));
fullpath = [(1:length(paths.pathMean))',halfpath'];

anchorpoints = [];
for ind = 1:size(templatesyllBreaks,1)     
     anchorpoints = [anchorpoints, templatesyllBreaks(ind,:)];
end
basepoints = [1, anchorpoints, length(paths.pathMean)];%Include the buffer regions
%basepoints = [anchorpoints];%Disregard buffers

[futurepoints] = getRemapAnchorsShort(fullpath,basepoints);

futureStr = [];
futureTime = [];
for i = 2:length(basepoints)
    segStops = ((basepoints(i-1)*44)-43):((basepoints(i)*44));
    futStops = ((futurepoints(i-1)*44)-43):((futurepoints(i)*44));
    seg = future_mNeuro(segStops);
    %remSeg = interp1(1:length(segStops),seg,1:length(futStops));
    remSeg = interp1(1:length(segStops),seg,linspace(1,length(segStops),length(futStops)));
    if i~=length(basepoints)
        futureStr = [futureStr, remSeg(1:end-44)];
        futureTime = [futureTime,futStops(1:end-44)];
    else
        futureStr = [futureStr, remSeg];
        futureTime = [futureTime,futStops];
    end
end
