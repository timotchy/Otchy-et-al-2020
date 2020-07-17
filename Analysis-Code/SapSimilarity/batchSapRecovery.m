function [accScores, simScores, indx] = batchSapRecovery(sounds, dur, REFsounds)
%The function batch processes a cell array of syllables/motifs to return SAP-style similarity scores that are relative
%to a specified baseline. Similarity is calculated by a separate script (sapSimilarity.m), but this higher-level
%script takes care of file handling for doing pairwise similarity for all sound samples. If REFsounds is an array of
%length m and sounds is of length n, this function will make m*n similarity comparisons with outputs to match.
% INPUTS:
%   sounds = a one-dimensional cell array with each cell containing the sound sample to be analyzed
%   dur = a one-dimensional integer array with the length of the syllable during each rendition (in ms) 
%   REFsounds = a one-dimensional cell array with each cell containing the sound of the reference syllable
%
% OUTPUTS:
%   accScores = an array of length m*n containing the "accuracy" (i.e., local similarity)
%   simScores = an array of length m*n containing the "% similarity" (i.e., global similarity)
%   indx = an array of size (m*n,2) showing the indices for sounds used to make those comparisons.
%
% Last updated 12/15/2016 by TMO

%Find number of sound samples to compare
numSounds = numel(sounds);
numREF = numel(REFsounds);

%Predefine output vars
outputSize = (numSounds * numREF);
accScores = NaN(outputSize,1);
simScores = NaN(outputSize,1);
indx = zeros(outputSize,2);

%Collect Distances
idx = 1;
for i = 1:numREF
    s1 = REFsounds{i};
    for j = 1:numSounds
        
        %Check to avoid same-syllable comparisons
        zeroCHK = s1(1) - sounds{j}(1);
        if zeroCHK ~= 0
            
            try
                %Similarity metrics
                [accScores(idx), simScores(idx)] = SapSimilarity3(s1, sounds{j}, 44150, dur);
                
            catch
                %Place holder
                accScores(idx) = NaN;
                simScores(idx) = NaN;
                display('Something looks wrong here... betta check it out!')
            end
            
            %Capture the sounds used
            indx(idx,:) = [i,j];
            
            %Increment the pointer
            idx = idx+1;
        end
    end
end







