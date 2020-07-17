function [alignV] = alignRawTS(tSeries,bins,p)
% Updated on 5/26/14 at 11:25pm by TMO
% This script aligns two series given a DTW path. Given is the tSeries, and DTW path with first column template
% index and the second column template is the series index. Returned is a
% row vector of size of template with entries that map to that point.

%Determine the extents of the series to warp
tSize = max(p(:,1));
[a b] = size(tSeries);
[pX pY] = size(p);

alignV = []; %Main Output
anchors = round(linspace(1,b,bins+1)); %Set the boundaries for the signal segments

DTWcount = 1;
for i = 1:tSize
    %Reset values for this signal segment
    count = 0;
    dum = [];
    prev = 0; ahead = 0;
    
    %Check for p(:,2) series previous
    if DTWcount~=1
        if p(DTWcount-1,2)==p(DTWcount,2) 
            prev = 1;
            r = 2;
            if DTWcount~=2
                while p(DTWcount-r,2)==p(DTWcount,2) && (DTWcount-r>0) %count up how many matching counts there are previously
                    prev = prev + 1;
                    r = r+1;
                    if DTWcount-r<=0
                        break;
                    end
                end
            end
        end
    end
    
    %Check for p(:,2) series ahead
    if DTWcount~=pX
        if p(DTWcount,2)==p(DTWcount+1,2) 
            ahead = 1;
            r = 2;
            if DTWcount ~= pX-1
                while p(DTWcount,2)==p(DTWcount+r,2) && (DTWcount+r<pX) %count up how many matching counts there are ahead
                    ahead = ahead + 1;
                    r = r+1;
                    if DTWcount+r>length(p)
                        break;
                    end
                end
            end
        end
    end
    
    %Calculate the start of this segment and how long it is
    seriesPos = prev+1;
    seriesSize = prev+1+ahead;
    
    %Build up the series of points to stretch in this segment
    while p(DTWcount,1) == i
        if seriesSize == 1
            segs = anchors(p(DTWcount,2)):anchors(p(DTWcount,2)+1);
        else
            steps = round(linspace(anchors(p(DTWcount,2)),anchors(p(DTWcount,2)+1),seriesSize+1));
            segs = steps(seriesPos):steps(seriesPos+1)-1;
            
            %The following 2 conditionals were added to catch errors that come from rounding 'steps' two lines above so that
            %the result can be used as indices.(as 'segs')
            if isempty(segs)
                segs = steps(seriesPos+1)-1;
            end
            if segs == 0
                segs = 1;
            end
            
        end
        dum = [dum tSeries(:,segs)];
        
        %Update advance/exit conditions
        count = count+1;
        DTWcount = DTWcount+1;
        
        if DTWcount > pX
            break;
        end
    end
    
    %Use linear interpolation to rescale the segment
    if length(dum)>1
        temp{i} = interp1(1:length(dum),dum,linspace(1,length(dum),44),'linear');
    else
        temp{i} = ones(1,44).*dum;
    end
end

%Build up the completely warped time series from the individual snippets stored in temp.
for j = 1:tSize
    alignV = [alignV,temp{j}];
end

end