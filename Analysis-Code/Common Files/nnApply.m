function [nn_classification] = nnApply(nn_classifier,new_data,classes,classification_thresh)
%Function to apply a given neural network classifier to new data 

% [nn_classification] = nnApply(new_data,classes) returns a vector of numeric
% labels based on the nn_classifier object, new_data matrix and classes label

% Input, output and dimensions are as follows:

%    nn_classifier (required)
%       Matlab object containing an existing trained network. Use output
%       from nnTrain
%    new_data  - N x S matrix (required)
%       Use output from nnFormatInputSpec. Each entry is power in log space. 
%       Each row is a different sample. Each column corresponds to the 
%       time-frequency bins from the vector output of nnFormatInputSpec.
%    classes - 1 x C vector (required)
%      Use output from nnFormatTargetClass. Each columns represents the 
%      sorted numeric labels (lowest to highest).
%    classification_thresh - a single scalar (optional)
%      A single scalar between 0 and 1 to use as a threshold for
%      classification. If none given, defaults to 0.8.
%    nn_classification - N x 1 vector
%      Each row represents a sample. Each entry corresponds to classification 
%      produced by the neural network in terms of the numeric labels in classes
%      Any sample that is not classifiable based on the
%      classification_thres is given NaN.
% 
%  Where:
%    N  = number of training samples
%    C  = number of classes to classify into based on the number of discrete
%         unique numeric labels given in classes

if nargin < 4
classification_thresh = 0.8;
end

%Apply neural network to new data
[output]=nn_classifier(new_data');
output=output';
% %Output classification
empty_count=0;classify_count=0;
for i=1:size(output,1)
    idx=find(output(i,:)>classification_thresh);
    if isempty(idx) %no output > classification_thresh
    empty_count = empty_count + 1;
    nn_classification(i,1) = NaN;    
    elseif length(idx)>1 %more than one output > classification_thresh, to assign to unknown class
    empty_count = empty_count + 1;
%     nn_classification(i,1) = NaN;    %Assign NaN
    idx = find(output(i,:) == max(output(i,:))); %Assign the highest class
    nn_classification(i,1) = classes(idx);
    else %only one output > classification_thresh
    classify_count = classify_count+1;
    nn_classification(i,1) = classes(idx);        
    end
    idx=[];
end

% classification_rate = classify_count./(empty_count+classify_count);


end