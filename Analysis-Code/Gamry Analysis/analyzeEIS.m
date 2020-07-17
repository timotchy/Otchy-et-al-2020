function analyzeEIS
% Does full electrode impedance spectroscopy anaylsis for selected files

% Pad dimensions for CIC calculations of Smania arrays
smth = 5; %smoothing elements

pW = 0.007; %in cm
pH = 0.0035; %in cm
pArea = pH*pW; %in cm^2

%Ask user which file to load; check for suffix
[filename, pathname, filterindex] = uigetfile('X:\EIS*.DTA', 'Pick a DTA file', 'MultiSelect', 'on');
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
figure(101); %clf
F = []; PH = []; Z = [];
h = []; t = [];
for i = 1:numFiles
    %Retrieve color
    cT = colorTint(i, 0);
    
    %Plot phase vs frequency
    subplot(2,1,1)
    semilogx(dataBlock(i).eis.freq, dataBlock(i).eis.Zph, 'Color', cT); hold on
    
    %Log data
    F(end+1,:) = dataBlock(i).eis.freq;
    PH(end+1,:) = dataBlock(i).eis.Zph;
    
    %Plot impedance vs frequency
    subplot(2,1,2)
    loglog(dataBlock(i).eis.freq, dataBlock(i).eis.Zmod, 'Color', cT); hold on
    
    %Log data
    Z(end+1,:) = dataBlock(i).eis.Zmod;
    
end

%Calculate the mean sweeps across files
mF = mean(F,1);
mPH = smooth(mean(PH,1), 5); sPH = smooth(std(PH, 1, 1), 5);
mZ = smooth(mean(Z,1), 5); sZ = smooth(std(Z, 1, 1), 5);

% Calculate nominal phase and impedance @ 1kHz
phase = interp1(mF, mPH, 1000);
phaseStr = num2str(phase, 4);

imp = interp1(mF, mZ, 1000);
impStr = num2str(round(imp)/1000, 5);

%Plot meanF vs PH
subplot(2,1,1)
h1 = shadedErrorBar(mF, mPH, sPH, '-k', 1);
h1.patch.FaceColor = colorTint(1, 0);
h1.mainLine.Color = colorTint(1, 0);

%Format
xlim([10^2, 10^5]);
set(gca,'Box', 'off', 'TickDir', 'out', 'XTickLabels', [], 'FontSize', 10);
ylabel('Phase (\Theta)');

%Plot Phase label
ys = ylim; y = mean(ys);
t = ['Phase @ 1kHz = ' phaseStr ' deg'];
text(10^2.5, y, t, 'FontSize', 10)

%Plot meanF vs PH
subplot(2,1,2)
Q = shadedErrorBar(mF, mZ, sZ, '-k', 1);
Q.patch.FaceColor = colorTint(1, 0);
Q.mainLine.Color = colorTint(1, 0);

% %Format
% set(get(get(Q.mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(Q.patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(Q.edge(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% set(get(get(Q.edge(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% fls = getFieldVectorCell(dataBlock, 'filename');
% legend(fls, 'interpreter', 'none', 'Location','southwest', 'FontSize', 7)

xlim([10^2, 10^5]);
set(gca,'Box', 'off', 'TickDir', 'out', 'FontSize', 10);
xlabel('Frequency (Hz)')
ylabel('Impedance (\Omega)');

%Plot Imp label
ys = ylim; y = mean(ys);
t = ['Z @ 1kHz = ' impStr 'k\Omega'];
text(10^3.5, 10^6, t, 'FontSize', 10)

set(gcf, 'Units', 'Inches', 'Position', [0.2, 0.5, 4.5, 7])






function cT = colorTint(fileNum, rendFrac)
%Based on the file snumber and rendition fraction, return a tint/shade of a
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






















