% function sapDistDistr
%This is a self-contained function for generating the short- and long-distance distributions that are needed to calculate
%similarity p-values. The proceedure is described in Ofer's Animal Behavior (2000) paper on automated similarity measures.
%The basic gist is that I'm going to collect the bin-wise distances between example song from unrelated birds. The sample
%dataset is loaded from data TMO recorded. The output will be summary cumulative distributions of the Ds and Dl matrices
%which can then be used to convert distances to p-values.

%Load the test data from file
% load('sapTestSet.mat', 'TestMotifs')

%Reset Vars
Dshort = [];
Dlong = [];

%Collect Distances
numSongs = numel(TestMotifs);
idx = 1;
for i = 1:numSongs
    for j = (i+1):numSongs
     [Dshort{idx}, Dlong{idx}] = sapSimilarity(TestMotifs(i).wav, TestMotifs(j).wav, 44150);
     idx = idx+1;
    end
end

%Concatenate cells to pool distance cross all comparisons
Ds_all = cell2mat(Dshort);
Dl_all = cell2mat(Dlong);

%Calculate cumulative distribution functions for each distance set
[pDs,xDs] = ecdf(Ds_all);
[pDl,xDl] = ecdf(Dl_all);
xDs(1) = xDs(1) - 0.000001;
xDl(1) = xDl(1) - 0.000001;

%Down sample the distribution
xDs_small = linspace(xDs(1), xDs(end), 10000);
pDs_small = interp1(xDs, pDs, xDs_small);

xDl_small = linspace(xDl(1), xDl(end), 10000);
pDl_small = interp1(xDl, pDl, xDl_small);

save('sapCDFs.mat', 'xDs_small', 'pDs_small', 'xDl_small', 'pDl_small')
