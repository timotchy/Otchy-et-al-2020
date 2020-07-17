function [songSnips, scores] =  batchProcessFolder
%Function runs SAP similarity analysis on all wav-files in a folder. The batch file handling and actual similarity scoring is
%done by external scripts. Here I just handles loading from disk, formating, and calling subfunctions

%Select folders containing song snips to process
folders = uipickfiles('FilterSpec', 'c:\Users\Tim\Desktop\', 'Prompt', 'Select folders to process');
numFolders = numel(folders);

%Predefine vars
songSnips = [];

%Cycle through folders and extract out songs from each
for i = 1:numFolders
    %retrieve names of all wav-files in the folder
    files = dir([folders{i}, filesep, '*.wav']);
    
    %Sequentially load each fileand parse away the data 
    for j = 1:numel(files)
        [data, ~] = audioread([folders{i}, filesep, files(j).name]);
        
        %Sort into data structure
        songSnips(i).folder{j} = folders{i};
        songSnips(i).name{j} = files(j).name;
        songSnips(i).wavs{j} = data;
    end
    
    %Batch process all songs in the folder with MxN comparisons
    [scores(i).acc, scores(i).sim, scores(i).indx] = batchSapSimilarity(songSnips(i).wavs);
    
end


%Plot something simple
figure;
for i = 1:numFolders
    mSim(i) = nanmean(scores(i).acc);
    sSim(i) = nanstd(scores(i).acc);
    sp = regexp(folders{i}, '\', 'split');
    day(i) = str2double(sp{end}(4:end));
    
    %Scatter individual points + x-jitter
    jitter = -0.4 + (0.8).*rand(numel(scores(i).acc),1);
    scatter(day(i)+jitter, scores(i).acc, '.'); hold on

end
errorbar(day, mSim, sSim, 'ok')
ylabel('Similarity (%)')
xlabel('Age (dph)')
xlim([40, 120]); ylim([0.5, 1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 40:20:120, 'YTick', [0.5:0.25:1])

%Plot something simple
figure;
for i = 1:numFolders
    mSim(i) = nanmean((1-scores(i).acc)/0.5);
    sSim(i) = nanstd((1-scores(i).acc)/0.5);
    
    %Scatter individual points + x-jitter
    jitter = -.4 + (0.8).*rand(numel(scores(i).acc),1);
    scatter(day(i)+jitter+jitter, (1-scores(i).acc)/0.5, '.'); hold on
    
end
errorbar(day, mSim, sSim, 'ok')
ylabel('Variability Score')
xlabel('Age (dph)')
xlim([40, 120]); ylim([0, 0.5])
set(gca, 'Box', 'off', 'TickDir', 'out', 'XTick', 40:20:120, 'YTick', [0.0:0.25:1])
















