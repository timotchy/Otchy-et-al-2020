function [accScores, simScores, indx] = batchSapSimilarity(sounds, dur)
%The function batch processes a cell array of syllables/motifs to return SAP-style similarity scores. Similarity is
%calculated by a separate script (sapSimilarity.m), but this higher-level script takes care of file handling for doing
%pairwise similarity for all sound samples. If sounds is an array of length m, this function will make m(m-1) similarity
%comparisons with outputs to match.
% INPUTS:
%   sounds = a one-dimensional cell array with each cell containing the sound sample to be analyzed
%
% OUTPUTS:
%   accScores = an array of length m(m-1) containing the "accuracy" (i.e., local similarity)
%   simScores = an array of length m(m-1) containing the "% similarity" (i.e., global similarity)
%   indx = an array of size (m(m-1),2) showing the indices for sounds used to make those comparisons.
%
% Last updated 7/27/2015 by TMO

%Find number of sound samples to compare
numSounds = numel(sounds);

%Predefine output vars
outputSize = (numSounds^2 - numSounds)/2;
accScores = zeros(outputSize,1);
simScores = zeros(outputSize,1);
indx = zeros(outputSize,2);

%Collect Distances
idx = 1;
for i = 1:numSounds
    s1 = sounds{i};
    for j = (i+1):numSounds
        try
        %Similarity metrics
%         [accScores(idx), simScores(idx)] = sapSimilarity(s1, sounds{j}, 44150, false);
        [accScores(idx), simScores(idx)] = SapSimilarity2(s1, sounds{j}, 44150, dur, false);

        catch
        %Place holder
        accScores(idx) = NaN;
        simScores(idx) = NaN;  
        
        end
        
        %Capture the sounds used
        indx(idx,:) = [i,j];
        
        %Increment the pointer
        idx = idx+1;
    end
end







