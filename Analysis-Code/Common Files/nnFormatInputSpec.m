function [log_spec_formatted] = nnFormatInputSpec(log_spec)
% Function to format an input spectrogram matrix for use in 
% neural network training

% [log_spec_formatted] = nnFormatInputSpec(log_spec) returns a
% vector spec_formatted based on the input matrix spec 

% Input, output and dimensions are as follows:

%    log_spec  - F x T logarithmic transformed spectrogram matrix (required)
%       Each entry is power in log space; Each row is a frequency bin 
%       and each column is a time bin.
%    log_spec_formatted - 1 x (F' x T') vector
%       Each entry is power in log space.

%  Where:
%    F  = number of frequency bins in the spectrogram
%    T  = number of time bins in the spectrogram
%    F' = number of resampled frequency bins based on parameters below
%    T' = number of resampled tim bins based on parameters below  



% Parameters for resampling spectrogram. Spectrogram needs to be resampled
% to ensure that all data are of identical sizes across classes and samples 
% as the neural network is currently not set up to handle unequal data
% sizes.
% These parameters have worked well so far based on the balance between size
% of input data to neural network and amount of information lost but can 
% potentially be further optimized. Currently produces data of 5000
% elements.
resample_len_freq = 50; %number of resampled frequency bins
resample_len_time = 100; %number of resampled tim bins based on parameters below

%resample spectrogram using above params
log_spec_resampled = imresize(log_spec,[resample_len_freq, resample_len_time]);

%Convert spectrogram to a vector. 
log_spec_formatted = reshape(log_spec_resampled,[],1)';

% In essence, it is concatenating the different frequencies to make full use of 
% the spectral information. This has worked well so far, but can be potentially be
% further optimized, say by contenating in time. However, because a neural 
% network input need not be ordered (meaning whether Var1 comes before Var2 or after),
% concatenating it in time should produce similar results.
% Use code below to concatenate in time
% log_spec_formatted = reshape(log_spec_resampled',1,[])';


end