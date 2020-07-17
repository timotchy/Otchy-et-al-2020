function analyzeCV
% Does full cyclic voltammery anaylsis for selected files

% Pad dimensions for CIC calculations of Smania arrays
smth = 5; %smoothing elements

pW = 0.007; %in cm
pH = 0.0035; %in cm
pArea = pH*pW; %in cm^2

SweepRate = 50; %in mV/s

%Ask user which file to load; check for suffix
[filename, pathname, filterindex] = uigetfile('X:\CV*.DTA', 'Pick a DTA file', 'MultiSelect', 'on');
if iscell(filename)
    numFiles = numel(filename);
else
    numFiles = 1;
end


%Read in the data from each of the Gamry files
if numFiles > 1
    %Sequentially read in the nultiple files
    for i = 1:numFiles
        fullName = [pathname, filesep, filename{i}];
        dataBlock(i) = DTAreader(fullName);
    end
else
    %Read in a single file
    fullName = [pathname, filesep, filename];
    dataBlock = DTAreader(fullName);
end

%Designate the figure for output
figure(100); clf
lgnd = []; h = []; hs = []; t = [];
for i = 1:numFiles
    %Get how many curves in file
    numCurves = numel(dataBlock(i).cvcurve);
    
    V = []; I = [];
    for j = 2:numCurves
        %Retrieve color
        cT = colorTint(i, j/(numCurves-1));
        
        %Plot this curve
        h(end+1) = plot(dataBlock(i).cvcurve(j).Vf, smooth(dataBlock(i).cvcurve(j).Im, smth)*1000./pArea, 'Color', cT); hold on
        
        
        %Copy out, minus first and second
        if j~=numCurves
            V(end+1,:) = dataBlock(i).cvcurve(j).Vf;
            I(end+1,:) = dataBlock(i).cvcurve(j).Im;
        end
        
        %In the legend, only show the first trace for each file
        if j ~= 2
            set(get(get(h(end),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
    
    %Calculate the mean CV sweeps for this file
    mV(i, :) = mean(V,1);
    mI(i, :) = mean(I,1);
    
    %Calculate the charge injection capacity
    cic(i) = cicCalc(mV(i, :), mI(i, :), pArea, SweepRate); %in mA/cm^2 (for an average sweep)
%     cic(i) = cicCalc(V(3,:), I(3,:), pArea, SweepRate); %in mA/cm^2 (just for one sweep)
    
    %Eliminate the positive half
    current = mI(i, :);
    current(current>0) = 0;

    %Plot the residual area
    p(i) = patch(mV(i, :), current*1000/pArea, cT, 'FaceAlpha', 0.5);
    set(get(get(p(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    %Assemble the CIC label
    t{end+1} = ['CIC' num2str(i) ' = ' num2str(cic(i),3) ' mC/cm^2'];
end

%Format the figure
fls = getFieldVectorCell(dataBlock, 'filename');
legend(fls, 'interpreter', 'none', 'Location','southeast', 'FontSize', 8)
xlim([-1, 1]); %ylim([-.2e-6, .2e-6]);
set(gca,'Box', 'off', 'TickDir', 'out', 'FontSize', 12);
ylabel('I/Area (mA/cm^2)');
xlabel('V vs Ag|AgCl');

%Plot CIC
ys = ylim; y = mean(ys);
text(0, y, t, 'interpreter', 'none', 'FontSize', 8)

function cic = cicCalc(V, I, pArea, SweepRate)
%Using the mean CV trace, calculate the charge injection capacity as the
%area of the negative portion of the CV sweep

%Voltage series
volt = V;

%Current series
smth = 5;
current = smooth(I, smth)';

%Eliminate the positive half
current(current>0) = 0;

%Calculate the residual area
cic = polyarea(volt,current*1000/pArea)*SweepRate/2;

function cT = colorTint(fileNum, rendFrac)
%Based on the file number and rendition fraction, return a tint/shade of a
%pre-specified color

%Define colors (hg2 standards)
colorBase = [0    0.4470    0.7410;...
            0.8500    0.3250    0.0980;...
            0.9290    0.6940    0.1250;...
            0.4940    0.1840    0.5560;...
            0.4660    0.6740    0.1880;...
            0.3010    0.7450    0.9330;...
            0.6350    0.0780    0.1840];
        
%We don't want to get too light (i.e., white), so scale with max = 0.5
fraction = rendFrac * 0.5;

%Generate tint by scaling color vector
cT = (1-colorBase(fileNum, :)) * fraction + colorBase(fileNum, :);






















