function [tSx, tSy]  = rasterLine(times, ymin, ymax)
%This function takes in the spike times for a rendition (times) and converts it to a line represenation for fast plotting

%Constants for ease of testing (delete later)
% ymin = 0.3;
% ymax = 1;
% times = [0.5, 1, 1.75, 2, 3, 4, 5, 5.5]';

%Initialize the output
tSx = [];
tSy = [];

%Total number of spikes to plot
numSpikes = numel(times);

if numSpikes > 0
    %Reformat times
    rTimes = [times, times, times];
    tSx = reshape(rTimes', 3*numSpikes,1);
    
    %Build height lines
    base = [ymin, ymax, NaN]';
    tSy = repmat(base, numSpikes, 1);
end









