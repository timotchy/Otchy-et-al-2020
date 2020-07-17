function [accScore, simScore] = SapSimilarity(audio1, audio2, fs, bCheckPlot)
%An attempt to recapituate the similarity measures of SAP, as described in the Animal Behavior (2000) article. In keeping
%with changes mad to SAP in later versions (i.e., 2011), Amplitude Modulation and Pitch Goodness are used instead of Spectral Continuity.
%INPUTS:
%   audio1 = raw audio of first sound (typically of only a single syllable)
%   audio2 = raw audio of second sound (typically of only a single syllable)
%   fs = sampling rate
%   bCheckPlot = logical flag for plotting the summary figure at end
%
% Last updated 7/24/2015 by TMO

%Fill in default values
if nargin < 3
    fs = 44150; %default: 44150Hz sampling rate
end
if nargin < 4
    bCheckPlot = false; %default: no summary figure
end

%Set constants for the program (these are based on default SAP values)
numBins = 30; %number of milliseconds in large windows
globalThresh = 0.25; %The p-val threshold for interval setting 

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
HP_fNorm = 300/(fs/2);
LP_fNorm = 6500/(fs/2);
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);
audio1 = filtfilt(BP_b,BP_a,audio1); audio1 = audio1-mean(audio1);
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

%Construct short and long-range Euclidean distance matrices time windows in the the two sounds; convert to p-values
[Ds, Dl, Ps, Pl] = createDistMatrices(fMat1,fMat2, numBins);

%Create similarity matrices
localSimMatrix = zeros(size(Ds));
indx = Pl <= globalThresh;
localSimMatrix(indx) = 1 - Ps(indx);
globalSimMatrix = 1 - Pl;

%Format for calibration output
% Dshort = Ds(:);
% Dlong = Dl(:);

% Calculate accuracy score based on the short/local distances
[accScore, accPartialSim, sT, accWarpPath] = getScore(localSimMatrix);

% Calculate similarity score based on the long distances
 for i = 1:size(sT,1)
    [partialSim(i), num(i), warpPath{i}] = getPartialSim(globalSimMatrix, sT(i,:), 0.03, 0);
 end
 
 %Scale partial similiarity scores and compute final score
simScore = sum(partialSim.*(num./min(size(globalSimMatrix))));

%Plot a check-in figure if selected
if bCheckPlot
    load('sapColormap.mat', 'cmap');
    figure(666); clf
    subplot(5,5,2:5) %SpecX
    displaySpecgramQuick(audio2, 44150, [0, 8000], [], 0)
    set(gca, 'TickDir', 'out', 'XTick', [], 'YTick', [], 'Xlabel', [], 'Ylabel', [])
    
    subplot(5,5,[6,11,16,21]) %SpecY
    displaySpecgramQuick(audio1, 44150, [0, 8000], [], 0)
    set(gca, 'TickDir', 'out', 'XTick', [], 'YTick', [], 'Xlabel', [], 'Ylabel', [])
    view(90, 90)
    
    subplot(5,5,[7:10, 12:15, 17:20, 22:25]) %SpecX
    imagesc(localSimMatrix); colormap(cmap); hold on
    for i = 1:numel(warpPath)
        plot(accWarpPath{i}(:,2), accWarpPath{i}(:,1), 'xw')
        plot(warpPath{i}(:,2), warpPath{i}(:,1), 'g')
    end
    set(gca, 'TickDir', 'out', 'YTick', [], 'Xlabel', [], 'Ylabel', [])
end

function [partialSim, num, p] = getPartialSim(simMatrix, rect, width, bScaleIt)
%Given a similarity matrix and an enclosing rectangle, theis function returns a partial similarity score, an interval
%length, and a scaling factor. By default (i.e., width == 0), the similarity path is constrained to the diagonal of the
%enclosing rectangle; however, if width > 0, the path is allowed to deviate from the diagonal by up to +/-width # of
%bins.

%Limit test input
width = max([0, width]); %width should be input as a % stretch; no negatives; width > 1 is possible, but path is always constrained to the enclosing rectangle
% width = 0.04;

%Define the rectanle diagonal
maxSize =  min(rect(3:4));
diagPnts = [(rect(2):(rect(2)+maxSize-1))', (rect(1):(rect(1)+maxSize-1))'];

%Preallocate the warped path output
p = nan(size(diagPnts));

%Adjust points off the diagonal as specified
simScores = [];
if width ~=0
    %Pre calculate the search areas
    searchDist = ceil([1:size(diagPnts,1)].*width);
    searchEnd = diagPnts+ [searchDist',searchDist'];
    
    %Limit the warping to the extents of the enclosing rectangle
    searchEnd(searchEnd(:,1) > (rect(2)+rect(4)-1),1) = (rect(2)+rect(4)-1);
    searchEnd(searchEnd(:,2) > (rect(1)+rect(3)-1),2) = (rect(1)+rect(3)-1);

    %Find the optimal/max similarity path
    for i = 1:size(diagPnts,1)
        %Row search indices
        r = (diagPnts(i,2):searchEnd(i,2))';
        rIdx = [diagPnts(i,1).*ones(size(r)), r];
        
        %Column search indices
        c = (diagPnts(i,1)+1:searchEnd(i,1))';
        cIdx = [c, diagPnts(i,2).*ones(size(c))];
        
        %Combine search indices
        rcIdx = [rIdx; cIdx];
        
        %Retrieve candidate similarities
        sims = diag(simMatrix(rcIdx(:,1),rcIdx(:,2)));
        
        %Select the point of maximum similarity and record the bin coordinates
        [simScores(i), idx] = max(sims);
        p(i,:) = rcIdx(idx,:);
    end
    
else
    p = diagPnts;
    simScores = diag(simMatrix(p(:,1), p(:,2)));
end

%Partial similarity is simply the mean of scores along the path
num = numel(simScores);
meanSim = mean(simScores);

%Length scaling factor
[m, n] = size(simMatrix);
scale = log(sqrt(rect(3)^2 + rect(4)^2)) / log(sqrt(m^2 + n^2));

%Final partial similarity score for the interval
if bScaleIt
    partialSim = meanSim * scale;
else
    partialSim = meanSim;
end

function [simScore,fPartialSim, sT, fWarpPath] = getScore(simMat)
%Cache original similarity matrix
simMatrix = simMat;

%Find all enclosing regions in the current similarity matrix
s = regionprops(simMat>0, 'BoundingBox');
t = ceil(reshape(struct2array(s), 4, [])'); %[y1, x1, dy, dx]

%Iteratively find the most significant regions of similarity and eliminate redundancy until there are no unallocated regions
sT = []; sSim = []; sNum = []; sPath = [];
while ~isempty(t)
    %Cycle through regions and calculate the partial similarity scores for each
    partialSim = []; num = []; warpPath = [];
    for i = 1:size(t,1)
        [partialSim(i), num(i), warpPath{i}] = getPartialSim(simMat, t(i,:), 0.03, 1);
    end
    
    %Sort all data by partial similarity
    [partialSim, idx] = sort(partialSim, 'descend');
    t = t(idx,:);
    num = num(idx);
    warpPath = warpPath(idx);
    
    %Retain the data for the highest partial sim window
    sT(end+1, :) = t(1,:);
    sSim(end+1) = partialSim(1);
    sNum(end+1) = num(1);
    sPath{end+1} = warpPath{1};
    
    %Update the similarity matrix to exclude all bins that overlap with the chosen window
    rIdx = warpPath{1}(1,1):warpPath{1}(end,1);
    cIdx = warpPath{1}(1,2):warpPath{1}(end,2);
    simMat(rIdx, :) = 0;
    simMat(:, cIdx) = 0;
    
    %Find all enclosing regions in the current similarity matrix
    s = regionprops(simMat>0, 'BoundingBox');
    t = ceil(reshape(struct2array(s), 4, [])'); %[y1, x1, dy, dx]
end

%Calculate final similarity values for the remaining windows using the cached sim matrix
fPartialSim = []; fNum = []; fWarpPath = [];
for i = 1:size(sT,1)
    [fPartialSim(i), fNum(i), fWarpPath{i}] = getPartialSim(simMatrix, sT(i,:), 0.02, 0);
end

%Scale and compute final simailrity scores
simScore = sum(fPartialSim.*(fNum./min(size(simMatrix))));

function [Ds, Dl, Ps, Pl] = createDistMatrices(fMat1,fMat2, numBins)
%Load the conversion table for distance-to-pval conversion
load('sapCDFs.mat');

%Construct Euclidean distance matrix between all 1ms timebins in the the two sounds
Ds = pdist2(fMat1', fMat2', 'euclidean'); %MxN matrix, where M is the duration of audio1; N length of audio2

%Preallocate variables for the other distance/p-val matrices
[m, n] = size(Ds);
Dl = zeros(m, n);
d = nan(1, n);

%Create the indices that will window the short-scale distances
ind1 = windowTS(1:size(fMat1,2), numBins, 1, 'pad', 'boxcar');
ind2 = windowTS(1:size(fMat2,2), numBins, 1, 'pad', 'boxcar');

%Cycle through all bins of the matrixs (super slow... improve if possible)
for i = 1:m
    idx1 = ind1(i,:);
    use1 = ~isnan(idx1);
    %Serial version
    for j = 1:n
        idx2 = ind2(j,:);
        use2 = ~isnan(idx2);
        
        %Use only the bins that match
        use = use1 & use2;
        k1 = idx1(use); 
        k2 = idx2(use);
        
        %Long-range distances are simply the mean of short-range distances
        %Compile all elements of the row into an array
        d(j) = sum(diag(Ds(k1, k2)).^2)/numel(k1);
    end
    
    %Parallel version (turned out to be very slow; too much overhead?)
%     parfor j = 1:n
%         idx2{j} = ind2(j,:);
%         use2(j,:) = ~isnan(idx2{j});
%         
%         %Use only the bins that match
%         use(j,:) = use1 & use2(j,:);
%         k1{j} = idx1(use(j,:)); 
%         k2{j} = idx2{j}(use(j,:));
%         
%         %Long-range distances are simply the mean of short-range distances
%         %Compile all elements of the row into an array
%         d(j) = sum(diag(Ds(k1{j}, k2{j})).^2)/numel(k1{j});
%     end


    %Copy the row array to the main data structure (saves a few % time)
    Dl(i,:) = d;
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







