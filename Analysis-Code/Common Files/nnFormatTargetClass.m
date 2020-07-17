function [classes target_classification] = nnFormatTargetClass(input_classes)
% Function to format an input vector of classes to a matrix for use in 
% neural network training

% [classes target_classification] = nnFormatTargetClass(input_classes) returns a
% matrix target_classification and vector classes based on the input vector 
% input_classes

% Input, outputs and dimensions are as follows:
%    input_classes  - N  x 1 vector (required)
%       Each row is a sample with a numeric label given to indicate its class.
%       Numeric labels do not need to start from any specific number or be 
%       consecutive.
%    classes - 1 x C vector
%       Each columns represents the sorted numeric labels.
%       (lowest to highest)
%    target_classification  - N  x C matrix
%       Each row of the matrix is a sample and has a column with a one to 
%       indicate the class the sample is from, with zeros in the other columns.  
%       Columns correspond to the numeric labels given in the columns of output 
%       classes. The ordering of the rows is maintained based on the input. 
%       See example below

%  Where:
%    N  = number of training samples based on the input
%    C  = number of classes to classify into based on the number of discrete
%         unique numeric labels given in the input

% Example
% input_classes = [1;2;2;5;5;6;6;6]
% target_classification = nnFormatTargetClass(input_classes)
% classes = [ 1 2 5 6 ];
% target_classification = [ 1 0 0 0
%                           0 1 0 0
%                           0 1 0 0
%                           0 0 1 0
%                           0 0 1 0
%                           0 0 0 1
%                           0 0 0 1
%                           0 0 0 1 ];

    

%Error checking
if isrow(input_classes)
        error('Input needs to be N X 1 single column vector');    
end


%find number of discrete classes based on unique numbers
classes = unique(input_classes);
num_input_classes = length(classes); 

%initialize output
target_classification = zeros(length(input_classes),num_input_classes);

%create output matrix by iterating through each row of input
for i=1:length(input_classes)
    column = find(input_classes(i,1)== classes); %find the column position
    target_classification(i,column)=1; %and set it to one       
end

classes=classes';

end