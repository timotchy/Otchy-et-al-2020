function [accScore, simScore] = SapSimilarity3(audio1, audio2, fs, dur)
%An attempt to recapituate the similarity measures of SAP, as described in the Animal Behavior (2000) article. In keeping
%with changes mad to SAP in later versions (i.e., 2011), Amplitude Modulation and Pitch Goodness are used instead of Spectral
%Continuity. This version differs from method described in that it
%implements linear temporal rescaling of syllables to match mean syllable
%duration.
%
%INPUTS:
%   audio1 = raw audio of first sound (typically of only a single syllable)
%   audio2 = raw audio of second sound (typically of only a single syllable)
%   fs = sampling rate (Hz)
%   dur = alignment duration (ms)
%
% Last updated 01/16/2017 by TMO

%Fill in default values
if nargin < 3
    fs = 44150; %default: 44150Hz sampling rate
end
if nargin < 4
    dur = []; %default: use the mean length of audio1 and audio2 (calculated below)
end

%Set constants for the program (these are based on default SAP values)
numBins = 30; %number of milliseconds in large windows

%Define spectral feature constants (from SAP: http://soundanalysispro.com/manual-1/)
m_logPitch = 6.534; mad_logPitch = 0.687;
m_FM = 44.10; mad_FM = 22.60;
m_AM = 0.00; mad_AM = 0.127;
m_entropy = -1.810; mad_entropy = 0.937;
m_logGood = 5.386; mad_logGood = 0.540;

%Format incoming data
audio1 = audio1(:);
audio2 = audio2(:);

%Lightly filter and mean subtract to remove noise
HP_fNorm = 300/(fs/2); % High pass @ 300 Hz
LP_fNorm = 6500/(fs/2); % %Low pass @ 8500 Hz
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]); %2-pole butterworth filter
audio1 = filtfilt(BP_b,BP_a,audio1); audio1 = audio1-mean(audio1); %zero-phase filtering
audio2 = filtfilt(BP_b,BP_a,audio2); audio2 = audio2-mean(audio2);

%Calculate the spectral features for each sound
feats1 = koenigSpectral(audio1, fs);
feats2 = koenigSpectral(audio2, fs);

%Normalize desired features to mad and format as a matrix for easier coding (these are basically z-scores)
fMat1(1,:) = (log(feats1.Pitch_chose)-m_logPitch)./mad_logPitch;
fMat1(2,:) = (feats1.FM-m_FM)./mad_FM;
fMat1(3,:) = (feats1.AM-m_AM)./mad_AM;
fMat1(4,:) = (feats1.Entropy-m_entropy)./mad_entropy;
fMat1(5,:) = (log(feats1.PitchGoodness)-m_logGood)./mad_logGood;

fMat2(1,:) = (log(feats2.Pitch_chose)-m_logPitch)./mad_logPitch;
fMat2(2,:) = (feats2.FM-m_FM)./mad_FM;
fMat2(3,:) = (feats2.AM-m_AM)./mad_AM;
fMat2(4,:) = (feats2.Entropy-m_entropy)./mad_entropy;
fMat2(5,:) = (log(feats2.PitchGoodness)-m_logGood)./mad_logGood;

%At this step, these two matrices need to be made the same length
%If duration wasn't specified, use the mean of the two passed files
if isempty(dur)
    dur = round(mean([size(fMat1,2), size(fMat2,2)]));
end

for i = 1:5
    fMat1_sL(i,:) = interp1(1:size(fMat1,2), fMat1(i,:), linspace(1, size(fMat1,2), dur));
    fMat2_sL(i,:) = interp1(1:size(fMat2,2), fMat2(i,:), linspace(1, size(fMat2,2), dur));
end

%Construct short and long-range Euclidean distance matrices time windows in the the two sounds; convert to p-values
%[Ds, Dl, Ps, Pl] = createDistMatrices(fMat1_sL ,fMat2_sL, numBins);
[~, ~, Ps, Pl] = createDistMatrices(fMat1_sL ,fMat2_sL, numBins);

%Create similarity matrices
localSimMatrix = 1 - Ps;
globalSimMatrix = 1 - Pl;

%Get the size of the square distance matrix
[m, ~] = size(Ps);

% Calculate accuracy score based on the short/local distances (i.e., the
% mean of the diagonal elements)
accScore = trace(localSimMatrix)/m;

% Calculate similarity score based on the long distances
simScore = trace(globalSimMatrix)/m;

%Plot a check-in figure if selected
% if bCheckPlot
%     load('sapColormap.mat', 'cmap');
%     figure(666); clf
%     subplot(5,5,2:5) %SpecX
%     displaySpecgramQuick(audio2, 44150, [0, 8000], [], 0)
%     set(gca, 'TickDir', 'out', 'XTick', [], 'YTick', [], 'Xlabel', [], 'Ylabel', [])
%     
%     subplot(5,5,[6,11,16,21]) %SpecY
%     displaySpecgramQuick(audio1, 44150, [0, 8000], [], 0)
%     set(gca, 'TickDir', 'out', 'XTick', [], 'YTick', [], 'Xlabel', [], 'Ylabel', [])
%     view(90, 90)
%     
%     subplot(5,5,[7:10, 12:15, 17:20, 22:25]) %SpecX
%     imagesc(localSimMatrix); colormap(cmap); hold on
%     for i = 1:numel(warpPath)
%         plot(accWarpPath{i}(:,2), accWarpPath{i}(:,1), 'xw')
%         plot(warpPath{i}(:,2), warpPath{i}(:,1), 'g')
%     end
%     set(gca, 'TickDir', 'out', 'YTick', [], 'Xlabel', [], 'Ylabel', [])
% end

function [Ds, Dl, Ps, Pl] = createDistMatrices(fMat1,fMat2, numBins)
%Load the conversion table for distance-to-pval conversion
load('sapCDFs.mat');

%Construct Euclidean distance matrix between all 1ms timebins in the the two sounds
Ds = pdist2(fMat1', fMat2', 'euclidean'); %MxN matrix, where M is the duration of audio1; N length of audio2

%Preallocate variables for the other distance/p-val matrices
[m, n] = size(Ds);
Dl = zeros(m, n);
%d = nan(1, n);

%Create the indices that will window the short-scale distances
ind1 = windowTS(1:size(fMat1,2), numBins, 1, 'pad', 'boxcar');
ind2 = windowTS(1:size(fMat2,2), numBins, 1, 'pad', 'boxcar');

%Cycle through all bins of the matrixs (super slow... improve if possible)
for i = 1:m
    idx1 = ind1(i,:);
    use1 = ~isnan(idx1);
    d = [];
    
    idx2 = ind2(i,:);
    use2 = ~isnan(idx2);
    
    %Use only the bins that match
    use = use1 & use2;
    k1 = idx1(use);
    k2 = idx2(use);
    
    %Long-range distances are simply the mean of short-range distances
    %Compile all elements of the row into an array
    Dl(i,i) = sum(diag(Ds(k1, k2)).^2)/numel(k1);
end

%Finish up Dl
Dl = Dl.^0.5;

%Do the conversion to p-values while you're here; the conversion function was determined empirically as described in
%the SAP document
Ps = interp1(xDs_small, pDs_small, Ds);
Ps(Ds < xDs_small(1)) = 0; %deal with out-of-bounds cases
Ps(Ds > xDs_small(end)) = 1;

Pl = interp1(xDl_small, pDl_small, Dl);
Pl(Dl < xDl_small(1)) = 0; %deal with out-of-bounds cases
Pl(Dl > xDl_small(end)) = 1;







