function [warpedOut] = getRemapAnchorsShort(path,anchors)
[m, n] = size(anchors);
warpedOut = zeros(m,n);

for i = 1:m
    for j = 1:n
        ind = find(path(:,1)==anchors(i,j));
        %ind = find(path(:,1)>=(anchors(i,j)-43) & path(:,1)<=anchors(i,j));
        warpedOut(i,j) = round(mean(path(ind,2)));
        %warpedOut(i,j) = round(min(path(ind,2)));
    end
end