function [channels,fs] = getChannels(file)
% Function reads in raw multi-channel recording from LABVIEW and outputs
% the channels stripped to rows the rows of a matrix.  Assumes the
% structure of the encoding, but not the content.

%Preallocate variables
chunks = [];
channels = [];

%Read the file
fid = fopen(file);
%dataset = fread(fid, inf,'single',0,'b');
%fclose(fid);

%Parse header information
fs = fread(fid, 1,'single',0,'b');
fs = 44150; %Hardcode this for now as many of the headers are incorrect.
crap = fread(fid, 1,'single',0,'b');
numChan = fread(fid, 1,'single',0,'b');
%numChan = 6;

%The first one second of data from each channel is read sequentially
for i=1:numChan
   chunks(i,:) =  fread(fid, fs,'single',0,'b');
end

%The the data is interleaved to the end of the file
MultiPlexed = fread(fid, inf,'single',0,'b');
for i=1:numChan
   channels(i,:) = [chunks(i,:), MultiPlexed(i:numChan:end)'];
%    channels(i,:) = [chunks(i,:), MultiPlexed(i:7:end)'];  %Correction for 6-chan mistake
end
fid = fclose('all');
end