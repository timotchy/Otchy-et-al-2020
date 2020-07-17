function [warpedOut] = getWarpedStarts(path,anchors)
[m, n] = size(anchors);
warpedOut = zeros(m,n);

for i = 1:m
    for j = 1:n
        ind = find(path(:,1)==anchors(i,j));
        warpedOut(i,j) = round(mean(path(ind,2)));
    end
end

end