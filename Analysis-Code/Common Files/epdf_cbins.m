function [X,Ys] = epdf_cbins(data,step,minXin,maxXin)
%Take in data vector and create probability distribution function given the
%bin size and min/max values

xin = reshape(data, numel(data), 1 );
if ~isreal( xin )
    xin = abs( xin );
end

% if floor( nbins ) ~= nbins
%     error( 'Number of bins should be integer value' );
% end
% if nbins < 2
%     error( 'Number of bins should be positive integer greater than 1 ' );
% end

% step = (maxXin - minXin) / (nbins-1);
binc = minXin : step : maxXin;     
[N, X] = hist(xin, binc);
Ys = N/sum(N);