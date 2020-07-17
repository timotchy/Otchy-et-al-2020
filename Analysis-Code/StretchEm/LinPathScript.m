 

tempAnchors = basepoints*44;
%tempAnchors = round(basepoints);
tempAnchors(1) = 1;
rendAnchors = round(longHVCfuturepoints);

 %Create linear path between the chosen anchor points
 LinPath = [];
 for j = 2:length(tempAnchors)
     path_t = linspace(rendAnchors(j-1),rendAnchors(j),tempAnchors(j)-tempAnchors(j-1)+1);
     if j~=length(tempAnchors)
        LinPath = [LinPath,path_t(1:end-1)];
     else
        LinPath = [LinPath,path_t];
     end
 end
 