function [PP] = createPreProcess
% This function creates a PreProcess file that specifies a default
% configuration for TweetVision and (possibly) future programs.  This can
% be used to create a variable for system testing or to define particular
% configurations

% Process steps:
%   1 = None
%   2 = Local Detrend
%   3 = Denoising
%   4 = High Pass
%   5 = Low Pass
%   6 = Amplify
%   7 = Smooth
%   8 = Down Sample
PP.Steps = [1 1 1 1 1 1 1];

% General Params
PP.params.Fs = 44150;       % Sampling rate in Hz
PP.params.tapers = [5 9];   % Tapers for spectral analysis
PP.params.pad = 5;          % Zero padding (next 2^x)
PP.params.fpass = [0 15000];% Frequencies of interest
PP.params.err = [0 0.05];   % Error handling
PP.params.trialave = 0;     % Trial averaging

% Local Detrending
PP.Detrend.window = 0.1;    %window length in seconds
PP.Detrend.winstep = 0.05;  %window overlap in seconds

% Denoising
PP.Denoise.LP = 0;          %lower bound of region to denoise
PP.Denoise.HP = 300;        %upper bound of region to denoise
PP.Denoise.p = 0.05;        %p-value for finding line noise and harmonics
PP.Denoise.plt = 'n';       %plotting on/off
PP.Denoise.f0 = 60;         %array of specific frequencies to remove

% High Pass Filtering
PP.HP.type = 'butter';      %filter type
PP.HP.cutoff = 300;         %cutoff frequency
PP.HP.order = 4;            %filter order
PP.HP.rpass_ripple = 0.5;   %response ripple int he passband
PP.HP.rstop_ripple = 0.5;   %response ripple int he stopband

% Low Pass Filtering
PP.LP.type = 'butter';
PP.LP.cutoff = 15000;
PP.LP.order = 4;
PP.LP.rpass_ripple = 0.5;
PP.LP.rstop_ripple = 0.5;

% Amplification
PP.Amp.gain = 10;           %gain

% Smoothing
PP.Smooth.type = 'sgolay';  %smoothing algorithm
PP.Smooth.window = 0.1;     %smoothing window in seconds
PP.Smooth.degree = 4;       %polynomial degree

% Downsampling
PP.DSample.target = 500;    %desired sampling rate

