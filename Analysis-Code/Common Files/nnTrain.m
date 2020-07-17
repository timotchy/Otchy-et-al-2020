function [nn_scores,nn_success_rate,nn_classifier] = nnTrain(data,targets)
%Function to create neural network classifier with the following parameters and
%defaults otherwise

% [nn_scores,nn_success_rate,nn_classifier] = nnTrain(data,targets) returns the matrix
% nn_classification, scalar nn_success_rate and Matlab object nn_classifier 

% Input, output and dimensions are as follows:

%    data  - N x S matrix (required)
%       Use output from nnFormatInputSpec. Each entry is power in log space. 
%       Each row is a different sample.Each column corresponds to the 
%       time-frequency bins from the vector output of nnFormatInputSpec.
%    targets  - N x C matrix (required)
%       Use output from nnFormatTargetClass. Each row is a sample. Each column 
%       of the matrix has a column with a one to indicate the class 
%       the sample is from, with zeros in the other columns.
%    nn_scores  - N x C matrix (required)
%       Output scores from the classifier. Numbers ranging from 1
%       to zero. The column with the number closest to 1 represents the 
%       class to be assigned to for each sample.
%    nn_success_rate  - single scalar
%       % of successful classification.
%    nn_classifier
%       Matlab object containing all the information about the trained network.
%       Type nn_classifier for detailed information.


% Set up neural network and parameters

% Parameter for number of neurons in hidden layer. This seems to work well 
% thus far. It is also within the common advice of setting the number of neurons 
% in the hidden layer to be between the input size (5000 neurons given the
% current spectrogram resampling parameters) and output size (~2-10 output
% classes). This is a rule of thumb and can potentially be further optimized
% in more systematic ways including finding optimal number of hidden layers
layer_size = 100; 

%Create neural network for pattern recognition
nn_classifier = patternnet(layer_size);

%Parameters for proportion of samples for training and validation
%since the data being fed is exclusively used for training, no testing is
%done in this function
nn_classifier.divideParam.trainRatio = 70/100;
nn_classifier.divideParam.valRatio = 30/100;

%Set some other parameters of the training
nn_classifier.trainParam.min_grad = 1e-8;
nn_classifier.trainParam.max_fail = 5;
nn_classifier.trainParam.showWindow = 0; %suppress training GUI

%Train and test neural network using parameters above
[nn_classifier,training_record] = train(nn_classifier,data',targets');

%Test neural network
nn_scores = nn_classifier(data');
misclassification_rate = perform(nn_classifier,targets,nn_scores');
nn_success_rate = (1-misclassification_rate).*100;

end