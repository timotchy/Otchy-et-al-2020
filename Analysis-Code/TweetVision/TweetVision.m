
function varargout = TweetVision(varargin)
% TWEETVISION M-file for TweetVision.fig
%      TWEETVISION, by itself, creates a new TWEETVISION or raises the existing
%      singleton*.
%
%      H = TWEETVISION returns the handle to a new TWEETVISION or the handle to
%      the existing singleton*.
%
%      TWEETVISION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWEETVISION.M with the given input arguments.
%
%      TWEETVISION('Property','Value',...) creates a new TWEETVISION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TweetVision_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TweetVision_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TweetVision

% Last Modified by GUIDE v2.5 23-Jun-2018 12:45:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TweetVision_OpeningFcn, ...
                   'gui_OutputFcn',  @TweetVision_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before TweetVision is made visible.
function TweetVision_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TweetVision (see VARARGIN)

% Choose default command line output for TweetVision
handles.output = hObject;

load('Default PreProcess.mat');
handles.PreProcess = PP;

handles.filenum = -1;
handles.selectedSyll = -1;
handles.startNdx = -1;
handles.endNdx = -1;
handles.txtHandles = [];
handles.rectHandles = [];
handles.navRect = [];
handles.bWaitingForAddClick = false;
handles.bWaitingForSeparationClick = false;
handles.annotation = mhashtable;
handles.annotationLoaded = false;

% set numSylls default to avoid error in network training
handles.numSylls = 1;

%%% This part is for the .DAT functionality of TweetVision
handles.cell_no = [0, 0, 0, 0];
%Create the baseline cell structure
cell.birdname = [];
cell.filename = [];
cell.channel = [];
cell.cell_no = [];
cell.filenum = [];
cell.spikes = [];
cell.noise = [];
cell.spont = [];

%Copy to the permenant structure
% creates four empty cells
handles.cell{1} = cell;
handles.cell{2} = cell;
handles.cell{3} = cell;
handles.cell{4} = cell;

contents = {'No Drug'; 'TTX'; 'Lidocaine'; 'Muscimol'; 'Ibotanic'};
set(handles.popup_drugstatus,'String', contents);
set(handles.popup_drugstatus,'Value', 1);
handles.drugstatus = 'No Drug';

contents = {'Undirected';'Directed'};
set(handles.popup_directstatus,'String', contents);
set(handles.popup_directstatus,'Value', 1);
handles.directstatus = 'Undirected';

% This is also not in TVL; probably used for the E-phys graphs
set(handles.axes_1,'XTick',[]);
set(handles.axes_2,'XTick',[]);
set(handles.axes_3,'XTick',[]);
set(handles.axes_4,'XTick',[]);
set(handles.axes_raster1,'XTick',[],'YTick',[]);
set(handles.axes_raster2,'XTick',[],'YTick',[]);
set(handles.axes_raster3,'XTick',[],'YTick',[]);
set(handles.axes_raster4,'XTick',[],'YTick',[]);

handles.Cedit_edgeThresh.value = '1.0';
handles.Cedit_contThresh.value = '0.8';
handles.Cedit_boutThresh.value = '75';
handles.Cedit_minSyl.value = '16';
handles.Cedit_maxSyl.value = '500';
handles.Cedit_minGap.value = '7';
handles.Cedit_boutOverhang.value = '150';
handles.Cedit_minBoutSep.value = '300';
handles.Cedit_minBout.value = '100';
handles.Cedit_silenceStart.value = '1000';
%handles.Cedit_axisLims.value = '0.5';

%Set slider control values
handles.Cslider_edgeThresh.max = 5;
handles.Cslider_edgeThresh.min = 0;
handles.Cslider_edgeThresh.sliderStep = [0.1, 0.1];
handles.Cslider_edgeThresh.value = str2double(handles.Cedit_edgeThresh.value);

handles.Cslider_contThresh.max = 5;
handles.Cslider_contThresh.min = 0;
handles.Cslider_contThresh.sliderStep = [0.1, 0.1];
handles.Cslider_contThresh.value = str2double(handles.Cedit_contThresh.value);

handles.Cslider_boutThresh.max = 250;
handles.Cslider_boutThresh.min = 0;
handles.Cslider_boutThresh.sliderStep = [0.25, 0.25];
handles.Cslider_boutThresh.value = str2double(handles.Cedit_boutThresh.value);

handles.Cslider_minSyl.max = 200;
handles.Cslider_minSyl.min = 0;
handles.Cslider_minSyl.sliderStep = [1, 1];
handles.Cslider_minSyl.value = str2double(handles.Cedit_minSyl.value);

handles.Cslider_maxSyl.max = 1500;
handles.Cslider_maxSyl.min = 0;
handles.Cslider_maxSyl.sliderStep = [50, 50];
handles.Cslider_maxSyl.value = str2double(handles.Cedit_maxSyl.value);

handles.Cslider_minGap.max = 50;
handles.Cslider_minGap.min = 0;
handles.Cslider_minGap.sliderStep = [1, 1];
handles.Cslider_minGap.value = str2double(handles.Cedit_minGap.value);

handles.Cslider_boutOverhang.max = 1000;
handles.Cslider_boutOverhang.min = 0;
handles.Cslider_boutOverhang.sliderStep = [25, 25];
handles.Cslider_boutOverhang.value = str2double(handles.Cedit_boutOverhang.value);

handles.Cslider_minBoutSep.max = 1500;
handles.Cslider_minBoutSep.min = 0;
handles.Cslider_minBoutSep.sliderStep = [50, 50];
handles.Cslider_minBoutSep.value = str2double(handles.Cedit_minBoutSep.value);

handles.Cslider_minBout.max = 1500;
handles.Cslider_minBout.min = 0;
handles.Cslider_minBout.sliderStep = [50, 50];
handles.Cslider_minBout.value = str2double(handles.Cedit_minBout.value);

handles.Cslider_silenceStart.max = 2500;
handles.Cslider_silenceStart.min = 0;
handles.Cslider_silenceStart.sliderStep = [100, 100];
handles.Cslider_silenceStart.value = str2double(handles.Cedit_silenceStart.value);

%handles.Cslider_axisLims.max = 20;
%handles.Cslider_axisLims.min = 0.01;
%handles.Cslider_axisLims.sliderStep = [.25, .25];
%handles.Cslider_axisLims.value = str2double(handles.Cedit_axisLims.value);

%populate slider values
handles.sliderObj = {'slider_edgeThresh';'slider_contThresh';'slider_boutThresh';'slider_minSyl';'slider_maxSyl';'slider_minGap';'slider_boutOverhang';'slider_minBoutSep';'slider_minBout';'slider_silenceStart'};
for i = 1:length(handles.sliderObj)
    handles = populateSlider(handles.sliderObj{i}, handles);
end

%populate edit values
handles.editObj = {'edit_edgeThresh';'edit_contThresh';'edit_boutThresh';'edit_minSyl';'edit_maxSyl';'edit_minGap';'edit_boutOverhang';'edit_minBoutSep';'edit_minBout';'edit_silenceStart'};
for i = 1:length(handles.editObj)
    handles = populateEdit(handles.editObj{i}, handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TweetVision wait for user response (see UIRESUME)
% uiwait(handles.TweetVision);
setappdata(0  , 'hTweetVision'    , gcf);
setappdata(gcf,   'PreProcess'    , handles.PreProcess);
setappdata(gcf, 'fhUpdateDataAxes', @updateDataAxes);


function handles = populateSlider(sliderName, handles)
%Set the specified values
eval(['set(handles.' sliderName ',''Max'', handles.C' sliderName '.max)'])
eval(['set(handles.' sliderName ',''Min'', handles.C' sliderName '.min)'])
eval(['set(handles.' sliderName ',''SliderStep'', handles.C' sliderName '.sliderStep./(handles.C' sliderName '.max - handles.C' sliderName '.min))'])
eval(['set(handles.' sliderName ',''Value'', handles.C' sliderName '.value)'])

function handles = populateEdit(editName, handles)
%Set the specified values
eval(['set(handles.' editName ',''String'', handles.C' editName '.value)'])

function handles = edit2slider(Obj, handles)
%Get the value from the edit box
value = str2num(get(Obj,'String'));

%Pump the value to the slider
editName = get(Obj,'Tag');
sliderName = ['slider' editName(5:end)];

eval(['set(handles.' sliderName ',''Value'', value)'])

function handles = slider2edit(Obj, handles)
value = num2str(get(Obj,'Value'));
    
%Update the value in the edit box by user selection
sliderName = get(Obj,'Tag');
editName = ['edit' sliderName(7:end)];

eval(['set(handles.' editName ',''String'', value)'])

function [currAnnot, handles] = multiSegment(handles, audio, boutThresh)
%If there has been segmentation before delete all boxes and numbers
if(size(handles.rectHandles,1)>1)
    handles=clearAll(handles);
end
handles.thresSyll = []; 
handles.thresEdge = []; 
handles.thresBout = []; 
handles.edgeBout = []; 

%Create low-passed sound amplitude
winSize = 4400; %100ms
winStep = 44; %1ms
[pow, ~] = aSAP_getLogPower(audio, handles.fs);
smPow = mean(windowTS_spec(pow,winSize,winStep))';

%The last section of recording should be "silence", so use that to measure the background noise
boutRatio = str2double(get(handles.edit_boutThresh,'String'));
silenceOffset = str2double(get(handles.edit_silenceStart,'String'));

%This is the switch for the old and new methods.

% mSilence = mean(smPow(end-silenceOffset:end));
% stdSilence = std(smPow(end-silenceOffset:end));
% ampThresh_new = (boutRatio * stdSilence) + mSilence;

Silence = smPow(end-silenceOffset:end);
% winSize = round(length(Silence)/5); % 20% of the total silence period
winSize = round(length(Silence)/20); % 10% of the total silence period
winStep = winSize; % no overlap
mSilence = mean(windowTS_spec(Silence,winSize,winStep))';
stdSilence = std(windowTS_spec(Silence,winSize,winStep))';
%cvSilence = abs(stdSilence./mSilence);
[~, indx] = min(mSilence);
ampThresh_new = (boutRatio * stdSilence(indx)) + mSilence(indx);

if isempty(boutThresh)
    ampThresh = ampThresh_new;  
else
    ampThresh = boutThresh;
end

%Find crossing points
X = smPow>ampThresh;
boutStart = find(diff(X)==1); %start of each period of interest
boutEnd = find(diff(X)==-1); %end of each period of interest
if isempty(boutStart)
    boutStart = 1;
end
if boutEnd(1) < boutStart(1) && (length(boutEnd)-length(boutStart) == 1)
    if boutEnd(1) < 400
        boutEnd(1) = [];
    else
        boutStart = [1; boutStart];
    end
end
if boutEnd(end) < boutStart(end) && (length(boutStart) - length(boutEnd) == 1)
    boutStart(end) = [];
end
    
%Refine points to excluse noise/waste
minBoutSep = str2double(get(handles.edit_minBoutSep,'String')); %in ms
minBoutLength = str2double(get(handles.edit_minBout,'String')); %in ms
boutGap = boutStart(2:end) - boutEnd(1:(end-1));
tooClose = find(boutGap < minBoutSep);
if ~isempty(tooClose)
    boutEnd(tooClose) = [];
    boutStart(tooClose+1) = [];
end

boutLength = boutEnd - boutStart;
tooShort = find(boutLength < minBoutLength);
if ~isempty(tooShort)
    boutEnd(tooShort) = [];
    boutStart(tooShort) = [];
end

%Retrieve syllable segmentation parameters
boutBuffer = str2double(get(handles.edit_boutOverhang,'String')); %in ms
edges = str2double(get(handles.edit_edgeThresh,'String'));
threshold = str2double(get(handles.edit_contThresh,'String'));
minSyl = str2double(get(handles.edit_minSyl,'String'))/1000; %in sec
maxSyl = str2double(get(handles.edit_maxSyl,'String'))/1000; %in sec
minGap = str2double(get(handles.edit_minGap,'String'))/1000; %in sec

%Loop through the bouts and segment each into syllables
syllStartTimes = [];
syllEndTimes = [];
thresEdge = [];
thresSyll = [];
numBouts = length(boutStart);
for i = 1:numBouts
    %Calculate time points
    startTime = max([boutStart(i)-boutBuffer, 1])/1000; %in sec
    startNdx = floor(startTime * handles.fs);
    endTime = min([boutEnd(i)+boutBuffer, length(audio)/44.15])/1000; %in sec
    endNdx = ceil(endTime * handles.fs);
    audio_select = audio(startNdx:endNdx);
    startB(i) = startTime*1000;
    endB(i) = endTime*1000;
    %Segment the bout

    [boutSyllStarts, boutSyllEnds, ~, ~, ~, thresEdge(i),thresSyll(i), ~] = aSAP_segSyllablesFromRawAudioBence(audio_select', handles.fs, edges, threshold, minSyl, maxSyl, minGap);
    boutSyllStarts = boutSyllStarts + startTime;
    boutSyllEnds = boutSyllEnds + startTime;

    if isempty(syllStartTimes)
        %If this is the first set of points added to the count...
        syllStartTimes = boutSyllStarts;
        syllEndTimes = boutSyllEnds;
    else
        %...otherwise, find the location to insert the syllable so everything remains in order
        insertStart = find(syllStartTimes < startTime,1,'last');
        insertEnd = find(syllStartTimes > endTime,1,'first');
%         insertStart = find(syllStartTimes < boutSyllStarts(1),1,'last');
%         insertEnd = find(syllStartTimes > boutSyllStarts(end),1,'first');
        if isempty(insertStart) && ~isempty(insertEnd)
            syllStartTimes = [boutSyllStarts, syllStartTimes];
            syllEndTimes = [boutSyllEnds, syllEndTimes];
        elseif ~isempty(insertStart) && isempty(insertEnd)
            syllStartTimes = [syllStartTimes(1:insertStart), boutSyllStarts];
            syllEndTimes = [syllEndTimes(1:insertStart), boutSyllEnds];
        else
            syllStartTimes = [syllStartTimes(1:insertStart), boutSyllStarts, syllStartTimes(insertEnd:end)];
            syllEndTimes = [syllEndTimes(1:insertStart), boutSyllEnds, syllEndTimes(insertEnd:end)];
        end
    end
end

%Remove any remaining overlaps
overLaps = find((syllStartTimes(2:end) - syllEndTimes(1:end-1)) < 0);
if ~isempty(overLaps)
    syllStartTimes(overLaps+1) = [];
    syllEndTimes(overLaps) = [];
end

%Recover guidelines for plotting
handles.thresSyll = thresSyll'; %in dB
handles.thresEdge = thresEdge'; %in dB
handles.thresBout = ampThresh; %in dB
handles.edgeBout = [startB', endB']; %in ms from start

%Generate Absolute times for the syllables
time = getFileTime(handles.filelist(handles.curfile).name);
segType=repmat(-1,length(syllStartTimes), 1);

%Create the annotation record
currAnnot.exper = handles.exper;
currAnnot.filenum = handles.filenum;
currAnnot.segAbsStartTimes = time + (syllStartTimes/(24*60*60));
currAnnot.segFileStartTimes = syllStartTimes;
currAnnot.segFileEndTimes = syllEndTimes;
currAnnot.segType = segType;
currAnnot.fs = handles.fs;
currAnnot.drugstatus=handles.drugstatus;
currAnnot.directstatus=handles.directstatus;

handles.lastMode = 'multi';

function [currAnnot, handles] = frameSegment(handles, audio)
%if there has been segmentation before delete all boxes and numbers
if(size(handles.rectHandles,1)>1)
    handles=clearAll(handles);
end
handles.thresSyll = []; 
handles.thresEdge = []; 
handles.thresBout = []; 
handles.edgeBout = []; 

fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/handles.fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;
audio_select = audio(startNdx:endNdx);

%Segment
%[syllStartTimes, syllEndTimes, noiseEst, noiseStd, soundEst, thresEdge,thresSyll, soundStd] = aSAP_segSyllablesFromRawAudioBence(audio_select, fs, handles.edges, handles.threshold);
edges = str2double(get(handles.edit_edgeThresh,'String'));
threshold = str2double(get(handles.edit_contThresh,'String'));
minSyl = str2double(get(handles.edit_minSyl,'String'))/1000; %in sec
maxSyl = str2double(get(handles.edit_maxSyl,'String'))/1000; %in sec
minGap = str2double(get(handles.edit_minGap,'String'))/1000; %in sec

[syllStartTimes, syllEndTimes, noiseEst, noiseStd, soundEst, thresEdge,thresSyll, soundStd] = aSAP_segSyllablesFromRawAudioBence(audio_select', handles.fs, edges, threshold, minSyl, maxSyl, minGap);

syllStartTimes=syllStartTimes+startTime;
syllEndTimes=syllEndTimes+startTime;
handles.thresSyll=thresSyll;
handles.thresEdge=thresEdge;
%set(handles.editSyllAbs, 'String',num2str(thresSyll,4));
%set(handles.editEdgeAbs, 'String',num2str(thresEdge,4));

time = getFileTime(handles.filelist(handles.curfile).name);
segType=repmat(-1,length(syllStartTimes), 1);


currAnnot.exper = handles.exper;
currAnnot.filenum = handles.filenum;
currAnnot.segAbsStartTimes = time + (syllStartTimes/(24*60*60));
currAnnot.segFileStartTimes = syllStartTimes;
currAnnot.segFileEndTimes = syllEndTimes;
currAnnot.segType = segType;
currAnnot.fs = fs;
currAnnot.drugstatus=handles.drugstatus;
currAnnot.directstatus=handles.directstatus;

% --- Outputs from this function are returned to the command line.
function varargout = TweetVision_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axis interaction functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function cb_keypress(hObject, evnt)
% hObject    handle to buttonAutoSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Keyboard shortcuts!! :)
handles = guidata(hObject);

%get modifier booleans
bShift = false;
bControl = false;
if(length(evnt.Modifier) >= 1)
    for(nMod = 1:length(evnt.Modifier))
        bShift = bShift || strcmp(evnt.Modifier{nMod}, 'shift');
        bControl = bControl || strcmp(evnt.Modifier{nMod}, 'control');
    end
end

[label, ok] = str2num(evnt.Key);
if(ok)
    if(bShift)
        label = label + 200;
    elseif(bControl)
        label = label + 300;
    end
    handles = labelSyllable(handles, label);
else
    if(strcmp(evnt.Key,'u')) 
        handles = labelSyllable(handles, 103);
    elseif(strcmp(evnt.Key,'s')) 
        handles = labelSyllable(handles, 100);
    elseif(strcmp(evnt.Key,'n')) 
        handles = labelSyllable(handles, 101);
    elseif(strcmp(evnt.Key,'c')) 
        handles = labelSyllable(handles, 102);
    elseif(strcmp(evnt.Key,'period')) 
        width = handles.endNdx - handles.startNdx;
        if(handles.endNdx + width/2 > length(handles.audio))
            handles.startNdx = length(handles.audio) - width;
            handles.endNdx = length(handles.audio);
        else
            handles.startNdx = round(handles.startNdx + width/1.5);
            handles.endNdx = round(handles.endNdx + width/1.5);
        end
        handles = updateZoomOrPan(handles);
    elseif(strcmp(evnt.Key,'comma'))
        width = handles.endNdx - handles.startNdx;
        if(handles.startNdx - width/2 < 1)
            handles.startNdx = 1;
            handles.endNdx = width;
        else
            handles.startNdx = round(handles.startNdx - width/2);
            handles.endNdx = round(handles.endNdx - width/2);
        end
        handles = updateZoomOrPan(handles);
    elseif(strcmp(evnt.Key,'space')) 
        if(bShift)
            handles = setSelectedSyllable(handles, handles.selectedSyll - 1);   
        else
            handles = setSelectedSyllable(handles, handles.selectedSyll + 1);
        end
    elseif(strcmp(evnt.Key,'p')) 
        soundsc(handles.audio(handles.startNdx:handles.endNdx), handles.exper.desiredInSampRate);
    elseif(strcmp(evnt.Key,'delete'))
        guidata(hObject, handles);
        buttonDeleteSyll_Callback(hObject, [], handles)
        handles = guidata(hObject);
    elseif(strcmp(evnt.Key,'m'))
        guidata(hObject, handles);
        buttonDeletePause_Callback(hObject, [], handles)
        handles = guidata(hObject);
    end
end
guidata(hObject, handles);

function cb_navbar_click(hObject, evnt)
% graphical input from mouse/cursor??
handles = guidata(hObject);
fs = handles.fs;
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axes_navBar, 'CurrentPoint');
axes(handles.axes_navBar);
if(strcmp(mouseMode, 'open'))
    % DOES NOT WORK; GETS LOOPED WITH ZOOM FUNCTION
    %if double click then zoom out
    %disp('open');
    %handles.startNdx = 1;
    %handles.endNdx = length(handles.audio);
    %handles = updateAudAxes(handles);
    %handles = updateDataAxes(handles);
elseif(strcmp(mouseMode, 'alt'))
    %pan length of drag.
elseif(strcmp(mouseMode, 'normal'))
    %zoom 
    rect = rbbox;
    endPoint = get(gca,'CurrentPoint'); 
    point1 = clickLocation(1,1:2);              % extract x and y
    point2 = endPoint(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    l = xlim;
    if((offset(1) / l(2))< .005)
        %if didn't drag a rectangle...
        %recenter on the click location.
        width = handles.endNdx - handles.startNdx;
        mid = (p1(1)*fs) + 1;
        handles.startNdx = round(mid - width/2);
        handles.endNdx = round(mid + width/2);
        if(handles.startNdx<1)
            handles.startNdx = 1;
            handles.endNdx = width + handles.startNdx;
        end
        if(handles.endNdx > length(handles.audio))
            handles.startNdx = length(handles.audio) - width;
            handles.endNdx = length(handles.audio);
        end        
        handles = updateAudAxes(handles);
        handles = updateZoomOrPan(handles);
        handles = updateDataAxes(handles);
    else
        handles.startNdx = max(1,round((p1(1)*fs) + 1));
        handles.endNdx = min(length(handles.audio), round(((p1(1) + offset(1))*fs) + 1));
        handles = updateAudAxes(handles);
        handles = updateZoomOrPan(handles);
        handles = updateDataAxes(handles);
    end
end

guidata(hObject, handles);

function cb_specgram_click(hObject, evnt)
handles = guidata(hObject);
fs = handles.exper.desiredInSampRate;
filenum = handles.filenum;
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axes_spec, 'CurrentPoint');
axes(handles.axes_spec);

if(handles.bWaitingForAddClick)
    %handle add syllable click...
    handles.bWaitingForAddClick = false;
    set(handles.push_addsyllable, 'BackgroundColor', 'white');
    
    %Check to see if there is an annotation for the current file
    if ~handles.annotation.containsKey(handles.filelist(handles.curfile).name)
        if ~isfield(handles,'exper')
            handles.exper = getExperInfo;
        end
        currAnnot.exper = handles.exper;
        currAnnot.filenum = filenum;
        currAnnot.segAbsStartTimes = [];
        currAnnot.segFileStartTimes = [];
        currAnnot.segFileEndTimes = [];
        currAnnot.segType = [];
        currAnnot.fs = fs;
        currAnnot.drugstatus=handles.drugstatus;
        currAnnot.directstatus=handles.directstatus;
        %currAnnot.syllthresh=handles.threshold;
        %currAnnot.edgethresh=handles.edges;
        %currAnnot.abssyllthresh=handles.thresSyll;
        %currAnnot.absedgethresh=handles.thresEdge;
        %currAnnot.drugindex=handles.drugindex;
        %currAnnot.drugstatus=handles.drugstatus;
        
    else
        currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
    end

    %get region for addition
    rect = rbbox;
    endPoint = get(gca,'CurrentPoint'); 
    point1 = clickLocation(1,1:2);       % extract x and y
    point2 = endPoint(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    
    if(get(handles.check_edgefind, 'Value'))
        %adjust p1 and offset accordingly.
        startNdx = round((p1(1)*fs) + 1);
        endNdx = round(((p1(1)+ offset(1))*fs) + 1);
        edgeThres = '0.6';
        [edgeThres,ok] = str2num(edgeThres);
          [pow, filtAud] = aSAP_getLogPower(handles.audio, fs);
        if(ok)
            ndx = find(pow(startNdx:endNdx) > edgeThres);
            if(length(ndx)>=2)
                p1(1) = ((startNdx + ndx(1) - 1) - 1) / fs;
                offset(1) = (((startNdx + ndx(end) - 1) - 1) / fs) - p1(1);
            end
        end
    end
    
    %be sure the region doesn't overlap another syllable.
    syllStartTimes = currAnnot.segFileStartTimes;
    syllEndTimes = currAnnot.segFileEndTimes;
    nSyllAfter = find(p1(1) + offset(1) < [syllStartTimes,Inf]); nSyllAfter = nSyllAfter(1);
    nSyllBefore = find(p1(1) > [-Inf, syllEndTimes]); nSyllBefore = nSyllBefore(end)-1;
    if(nSyllBefore - nSyllAfter == -1)
        time = getFileTime(handles.filelist(handles.curfile).name);
        currAnnot.segAbsStartTimes = [currAnnot.segAbsStartTimes([1:nSyllBefore]),time + (p1(1)/(24*60*60)),currAnnot.segAbsStartTimes([nSyllAfter:end])];
        currAnnot.segFileStartTimes = [currAnnot.segFileStartTimes([1:nSyllBefore]),p1(1),currAnnot.segFileStartTimes([nSyllAfter:end])];
        currAnnot.segFileEndTimes = [currAnnot.segFileEndTimes([1:nSyllBefore]),p1(1) + offset(1),currAnnot.segFileEndTimes([nSyllAfter:end])];
        currAnnot.segType = [currAnnot.segType([1:nSyllBefore]);-1;currAnnot.segType([nSyllAfter:end])];
        handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
        handles = drawSyllableRect(handles, nSyllBefore + 1, true);

    else
        warndlg('Specified segment overlaps another segment.');
        uiwait;
    end    
elseif(handles.bWaitingForSeparationClick)
    handles.bWaitingForSeparationClick = false;
elseif(handles.annotation.containsKey(handles.filelist(handles.curfile).name))
    %Otherwise check for syllable to select.
    currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
    nSelect = find(clickLocation(1,1) >= currAnnot.segFileStartTimes & ...
                   clickLocation(1,1) <= currAnnot.segFileEndTimes);
    if(length(nSelect) == 1)
        if(nSelect ~= handles.selectedSyll)
            handles = setSelectedSyllable(handles, nSelect);
        end
        if(strcmp(mouseMode, 'open'))
            if(~isfield(handles, 'templates'))
                handles.templates = [];
                handles.templates.wavs = [];
                %handles.templateFeatures=[];
            else
                    if (size(handles.templates.wavs,2)~=0)
                        segTypes=getFieldVector(handles.templates.wavs,'segType');
                        if (~isempty(find(segTypes==currAnnot.segType(handles.selectedSyll))))
                            warndlg('Delete the existing syllable with the same Annotation first');
                        return;
                        end
                    end
            end


            startTime = currAnnot.segFileStartTimes(handles.selectedSyll);
            endTime = currAnnot.segFileEndTimes(handles.selectedSyll);
            if (isempty(handles.templates.wavs))
                handles.templates.wavs(1).filename = handles.filelist(handles.curfile).name;
                handles.templates.wavs(1).startTime = startTime;
                handles.templates.wavs(1).endTime = endTime;
                handles.templates.wavs(1).fs = fs;
                handles.templates.wavs(1).wav = handles.audio(round((startTime*fs)+1):round((endTime*fs)+1));
                handles.templates.wavs(1).segType = currAnnot.segType(handles.selectedSyll); 
            else
                segType=getFieldVector(handles.templates.wavs,'segType');
                index=find(segType==currAnnot.segType(handles.selectedSyll));
                if (isempty(index))
                    index=find(segType<currAnnot.segType(handles.selectedSyll),1,'last')+1;
                end
                if (isempty(index))
                    index=1;
                end
                for i=length(handles.templates.wavs)+1:-1:index+1
                    handles.templates.wavs(i).filename = handles.templates.wavs(i-1).filename;
                    handles.templates.wavs(i).startTime = handles.templates.wavs(i-1).startTime;
                    handles.templates.wavs(i).endTime = handles.templates.wavs(i-1).endTime;
                    handles.templates.wavs(i).fs = handles.templates.wavs(i-1).fs;
                    handles.templates.wavs(i).wav = handles.templates.wavs(i-1).wav; 
                    handles.templates.wavs(i).segType = handles.templates.wavs(i-1).segType;
                end
                handles.templates.wavs(index).filename = handles.filelist(handles.curfile).name;
                handles.templates.wavs(index).startTime = startTime;
                handles.templates.wavs(index).endTime = endTime;
                handles.templates.wavs(index).fs = fs;
                handles.templates.wavs(index).wav = handles.audio(round((startTime*fs)+1):round((endTime*fs)+1));
                handles.templates.wavs(index).segType = currAnnot.segType(handles.selectedSyll); 
            end
            
            handles = updateTemplates(handles);
            %handles=updateTemplateFeatures(handles,size(handles.templates.wavs,2),'add');
        
            
        end
    end
end
guidata(hObject, handles);

function cb_templates_click(hObject, evnt)
hTweetVision = getappdata(0, 'hTweetVision');
handles = guidata(hTweetVision);
fs = handles.exper.desiredInSampRate;
filenum = handles.filenum;
mouseMode = get(get(hObject,'Parent'), 'SelectionType');
clickLocation = get(handles.axesTemplates, 'CurrentPoint');
axes(handles.axesTemplates);

if(strcmp(mouseMode, 'open'))
    if(isfield(handles, 'templates'))
        stdfs = fs;
        edges = [0];
        allwavs = [];
        for(nWav = 1:length(handles.templates.wavs))
            fs = handles.templates.wavs(nWav).fs;
            wav = handles.templates.wavs(nWav).wav;
            label = handles.templates.wavs(nWav).segType;
            wav = resample(wav, stdfs, fs);
            allwavs = [allwavs; wav'];
            edges = [edges, edges(end) + (length(wav) - 1)/stdfs];
        end
        nSelect = find(clickLocation(1,1) < edges);
        if(length(nSelect) > 0)
            nSelect = nSelect(1) - 1;
            %handles=updateTemplateFeatures(handles,nSelect,'delete');
            handles.templates.wavs(nSelect) = [];
            handles = updateTemplates(handles);
            
        end
    end
end
guidata(hTweetVision, handles);

% specific to TV
function cb_rasters_click(hObject, evnt)
%Capture the click point and axes name
handles = guidata(hObject);
click = get(gca,'CurrentPoint');
tag = get(gca,'Tag');
ind = str2num(tag(end));

%Extract the raster index
if isfield(handles,'spikemarker')
    if (~isempty(handles.spikemarker))
        delete(handles.spikemarker);
       handles.spikemarker = [];
    end
end

%Record and pllot the marker position
handles.spiketime = click(1,1);
hold on
handles.spikemarker = plot(handles.spiketime, -0.1,'g:*'); 
hold off

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axis update functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = drawNavBarSyllableRects(handles)
filenum = handles.filenum;
fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/fs;
endTime = (handles.endNdx-1)/fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

if(handles.annotation.containsKey(handles.filelist(handles.curfile).name))
    currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
    syllStartTimes = currAnnot.segFileStartTimes;
    syllEndTimes = currAnnot.segFileEndTimes;
    syllType = currAnnot.segType;
    x = [syllStartTimes; syllStartTimes; syllEndTimes; syllEndTimes; syllStartTimes];   
    color = zeros(3,length(syllStartTimes));
    color(3,syllType~=-1) = 1;
    color(1,syllType==-1) = 1;
    if(handles.selectedSyll~=-1)
        color(:,handles.selectedSyll) = [0,1,0];
    end
    
    %draw nav bar rectangles.
    axes(handles.axes_navBar)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
    for(nSyll = 1:length(syllStartTimes))
        handles.rectHandles(nSyll, 1) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', 'none', 'FaceAlpha', .2);
    end
    set(handles.rectHandles(:,1), 'HitTest', 'off');
end

function handles = drawSyllableRect(handles, nSyll, bInsert)
fs = handles.fs;
windowStartTime = (handles.startNdx-1)/fs;
windowEndTime = (handles.endNdx-1)/fs;
currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
startTime = currAnnot.segFileStartTimes(nSyll);
endTime = currAnnot.segFileEndTimes(nSyll);
label = currAnnot.segType(nSyll);
if(label ~= -1)
    color = [0 0 1];
else
    color = [1 0 0];
end

if(handles.selectedSyll == nSyll)
    color = [0 1 0];
end

if(~exist('bInsert') || ~bInsert)
    handles = deleteSyllableRectAndTxt(handles, nSyll);
end

x = [startTime, startTime, endTime, endTime, startTime];
axes(handles.axes_navBar)
lims = ylim;
y = [lims(1), lims(2), lims(2), lims(1), lims(1)];

if(exist('bInsert') && bInsert)
    handles.rectHandles = [handles.rectHandles([1:nSyll-1],:);repmat(-1,1,size(handles.rectHandles,2));handles.rectHandles([nSyll:end],:)];
    handles.txtHandles = [handles.txtHandles([1:nSyll-1],:);repmat(-1,1,size(handles.txtHandles,2));handles.txtHandles([nSyll:end],:)];
end

handles.rectHandles(nSyll, 1) = patch(x,y, color, 'EdgeColor', 'none', 'FaceAlpha', .2);
set(handles.rectHandles(nSyll,1), 'HitTest', 'off'); 

if(startTime<=windowEndTime && endTime>=windowStartTime)    
    axes(handles.axes_spec)
    lims = ylim;
    y = [lims(1), lims(2), lims(2), lims(1), lims(1)];
    handles.rectHandles(nSyll, 2) = patch(x,y, color, 'EdgeColor', color, 'FaceAlpha', 0);

%     axes(handles.axesPower)
%     lims = ylim;
%     y = [lims(1), lims(2), lims(2), lims(1), lims(1)];
%     handles.rectHandles(nSyll, 3) = patch(x,y, color, 'EdgeColor', 'none', 'FaceAlpha', .2);
% 
%     axes(handles.axesSignal)
%     lims = ylim;
%     y = [lims(1), lims(2), lims(2), lims(1), lims(1)];
%     handles.rectHandles(nSyll, 4) = patch(x,y, color, 'EdgeColor', 'none', 'FaceAlpha', .2);
    
    set(handles.rectHandles(nSyll,2), 'HitTest', 'off');
    %set(handles.rectHandles(nSyll,2:4), 'HitTest', 'off');

    if(label ~= -1) 
        midSyll = mean([startTime,endTime]);
        %create a new text handles
        axes(handles.axes_spec);
        handles.txtHandles(nSyll, 1) = text(midSyll,mean(ylim),num2str(label));
        
%         axes(handles.axesPower);
%         handles.txtHandles(nSyll, 2) = text(midSyll,mean(ylim),num2str(label));
%         axes(handles.axesSignal);
%         handles.txtHandles(nSyll, 3) = text(midSyll,mean(ylim),num2str(label));
        set(handles.txtHandles(nSyll,:), 'Color', 'black');
        set(handles.txtHandles(nSyll,:), 'FontSize', 14);
        set(handles.txtHandles(nSyll,:), 'FontWeight', 'bold');
        set(handles.txtHandles(nSyll,:), 'HitTest', 'off');   
        set(handles.txtHandles(nSyll,:), 'HorizontalAlignment', 'center');         
    end
end

function handles = drawZoomSyllableRects(handles)
filenum = handles.filenum;
fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/fs;
endTime = (handles.endNdx-1)/fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

if(handles.annotation.containsKey(handles.filelist(handles.curfile).name))
    currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
    syllStartTimes = currAnnot.segFileStartTimes;
    syllEndTimes = currAnnot.segFileEndTimes;
    midSyll = (syllStartTimes + syllEndTimes) / 2;
    syllType = currAnnot.segType;
    x = [syllStartTimes; syllStartTimes; syllEndTimes; syllEndTimes; syllStartTimes];   
    color = zeros(3,length(syllStartTimes));
    color(3,syllType~=-1) = 1;
    color(1,syllType==-1) = 1;
    if(handles.selectedSyll~=-1)
        color(:,handles.selectedSyll) = [0,1,0];
    end
 
    %Only the visible rectangles need be drawn in the other axes.
    ndx = find(syllStartTimes < endTime & syllEndTimes> startTime);
    
    % Added from TVL; may not be necessary
    if x(1,ndx(1)) < startTime
         x([1,2,5],ndx(1)) = startTime+0.005;
    end
    if x(3,ndx(end)) > endTime
         x([3,4],ndx(end)) = endTime-0.005;
    end
    
    axes(handles.axes_spec)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
    for(nSyll = ndx)
        handles.rectHandles(nSyll, 2) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', color(:,nSyll)', 'FaceAlpha', 0);
        if(syllType(nSyll) ~= -1)
            handles.txtHandles(nSyll, 1) = text(midSyll(nSyll),mean(ylim),num2str(syllType(nSyll)));
            set(handles.txtHandles(nSyll,1), 'Color', 'black');
            set(handles.txtHandles(nSyll,1), 'FontSize', 14);
            set(handles.txtHandles(nSyll,1), 'FontWeight', 'bold');
            set(handles.txtHandles(nSyll,1), 'HitTest', 'off'); 
            set(handles.txtHandles(nSyll,1), 'HorizontalAlignment', 'center'); 
        else
            handles.txtHandles(nSyll,1) = -1;
        end
    end

%     axes(handles.axesPower)
%     lims = ylim;
%     y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
%     for(nSyll = ndx)
%         handles.rectHandles(nSyll, 3) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', 'none', 'FaceAlpha', .2);
%         if(syllType(nSyll) ~= -1)
%             handles.txtHandles(nSyll, 2) = text(midSyll(nSyll),mean(ylim),num2str(syllType(nSyll)));
%             set(handles.txtHandles(nSyll,2), 'Color', 'black');
%             set(handles.txtHandles(nSyll,2), 'FontSize', 14);
%             set(handles.txtHandles(nSyll,2), 'FontWeight', 'bold');
%             set(handles.txtHandles(nSyll,2), 'HitTest', 'off'); 
%             set(handles.txtHandles(nSyll,2), 'HorizontalAlignment', 'center'); 
%         else
%             handles.txtHandles(nSyll,2) = -1;
%         end
%     end
%     
%     axes(handles.axesSignal)
%     lims = ylim;
%     y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
%     for(nSyll = ndx)
%         handles.rectHandles(nSyll, 4) = patch(x(:,nSyll),y, color(:,nSyll)', 'EdgeColor', 'none', 'FaceAlpha', .2);
%         if(syllType(nSyll) ~= -1)
%             handles.txtHandles(nSyll, 3) = text(midSyll(nSyll),mean(ylim),num2str(syllType(nSyll)));
%             set(handles.txtHandles(nSyll,3), 'Color', 'black');
%             set(handles.txtHandles(nSyll,3), 'FontSize', 14);
%             set(handles.txtHandles(nSyll,3), 'FontWeight', 'bold');
%             set(handles.txtHandles(nSyll,3), 'HitTest', 'off'); 
%             set(handles.txtHandles(nSyll,3), 'HorizontalAlignment', 'center'); 
%         else
%             handles.txtHandles(nSyll,3) = -1;
%         end
%     end
    
    set(handles.rectHandles(ndx,2), 'HitTest', 'off'); 
    handles.rectHandles(handles.rectHandles==0) = -1;
    handles.rectHandles(handles.txtHandles==0) = -1;
    
    
    
    %%%%% THIS NEXT PART MIGHT BE WRONG
        if get(handles.check_audPower,'Value') 
        %Add here all of the code for drawing the segmentation guides on the amplitude
        handles = deleteSegmentGuides(handles);
        colorBout = [0 0 0; 1 0 1; 0 1 1; 0 1 0]; 
        colorBout = [colorBout; colorBout; colorBout; colorBout; colorBout; colorBout];

        %Check for the existence of the prereq values (b/c files not all lcalls of this function will already have right info)
        if get(handles.check_showGuides,'Value') && isfield(handles,'thresBout') && isfield(handles,'edgeBout') && ~isempty(handles.thresBout) && ~isempty(handles.edgeBout) %For the bout guides
            
            %Plot patches for the bout start/stops
            boutStartTimes = handles.edgeBout(:,1)'/1000; %in sec
            boutEndTimes = handles.edgeBout(:,2)'/1000; %in sec
            x = [boutStartTimes; boutStartTimes; boutEndTimes; boutEndTimes; boutStartTimes];
            y = ylim;
            y = [y(1), y(2), y(2), y(1), y(1)];
            
            %Only the visible objects need be drawn in the axes.
            ndx = find(boutStartTimes < endTime & boutEndTimes> startTime);
            if x(1,ndx(1)) < startTime
                x([1,2,5],ndx(1)) = startTime;
            end
            if x(3,ndx(end)) > endTime
                x([3,4],ndx(end)) = endTime;
            end
            
            for(i = ndx)
                %Patches for bout edges
                handles.boutHandles(i, 2) = patch(x(:,i),y, color(:,i)', 'EdgeColor', colorBout(i,:), 'FaceAlpha', 0, 'LineWidth',2);
            end
            
            %Line for bout threshold
            handles.boutThreshHandles(1,2) = line([startTime endTime], [handles.thresBout handles.thresBout],'Color','k','LineStyle','-');
        end

        if get(handles.check_showGuides,'Value') && isfield(handles,'thresSyll') && isfield(handles,'thresEdge') && ~isempty(handles.thresSyll) && ~isempty(handles.thresEdge) %For the syllable guides
            %Set start and end points for the syllThresh Line
            if isfield(handles,'edgeBout') && ~isempty(handles.edgeBout)
                
            else
                boutStartTimes = min(syllStartTimes); %in sec
                boutEndTimes = max(syllEndTimes); %in sec
            end
            
            %Only the visible objects need be drawn in the axes.
            ndx = find(boutStartTimes < endTime & boutEndTimes> startTime);
            
            if boutStartTimes(ndx(1)) < startTime
                boutStartTimes(ndx(1)) = startTime;
            end
            if boutEndTimes(ndx(end)) > endTime
                boutEndTimes(ndx(end)) = endTime;
            end
            
            %Plot the lines for the intervals at the given thresholds
            for i = 1:length(boutStartTimes)
                handles.syllEdgeHandles(i,2) = line([boutStartTimes(i) boutEndTimes(i)], [handles.thresEdge(i) handles.thresEdge(i)],'Color',colorBout(i,:),'LineStyle','--');
                
                handles.syllContHandles(i,2) = line([boutStartTimes(i) boutEndTimes(i)], [handles.thresSyll(i) handles.thresSyll(i)],'Color',colorBout(i,:),'LineStyle',':');
            end
        end
    end
    %%%%%%%%%%%%%%%%
    
end

function handles = deleteSyllableRectAndTxt(handles, nSyll)
if(size(handles.txtHandles,1) >= nSyll)
    for(nTxt = 1:size(handles.txtHandles,2))
        if(ishandle(handles.txtHandles(nSyll,nTxt))) && (handles.txtHandles(nSyll,nTxt) ~= 0)
            delete(handles.txtHandles(nSyll,nTxt));
            handles.txtHandles(nSyll,nTxt) = -1;
        end
    end
end
if(size(handles.rectHandles,1) >= nSyll)
    for(nRect = 1:size(handles.rectHandles,2))
        if(ishandle(handles.rectHandles(nSyll,nRect))) && (handles.rectHandles(nSyll,nRect) ~= 0)
            delete(handles.rectHandles(nSyll,nRect));
            handles.rectHandles(nSyll,nRect) = -1;
        end
    end
end

%%% Pulled from TVL; might be uneccessary
function handles = deleteSegmentGuides(handles)
if isfield(handles, 'boutHandles')
    for(i = 1:size(handles.boutHandles,1))
        if(ishandle(handles.boutHandles(i,2))) && (handles.boutHandles(i,2) ~= 0)
            delete(handles.boutHandles(i,2));
        end
    end
end

if isfield(handles, 'boutThreshHandles')
    for i = 1:size(handles.boutThreshHandles,1)
        if(ishandle(handles.boutThreshHandles(i,2))) && (handles.boutThreshHandles(i,2) ~= 0)
            delete(handles.boutThreshHandles(i,2));
        end
    end
end

if isfield(handles, 'syllEdgeHandles')
    for(i = 1:size(handles.syllEdgeHandles,1))
        if(ishandle(handles.syllEdgeHandles(i,2))) && (handles.syllEdgeHandles(i,2) ~= 0)
            delete(handles.syllEdgeHandles(i,2));
        end
    end
end

if isfield(handles, 'syllContHandles')
    for(i = 1:size(handles.syllContHandles,1))
        if(ishandle(handles.syllContHandles(i,2))) && (handles.syllContHandles(i,2) ~= 0)
            delete(handles.syllContHandles(i,2));
        end
    end
end

function handles = labelSyllable(handles, label)
if(handles.filenum > 0 && handles.selectedSyll > 0)
    filenum = handles.filenum;
    nSyll = handles.selectedSyll;
    fs = handles.exper.desiredInSampRate;
    
    %set the label value
    currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
    currAnnot.segType(nSyll) = label;
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %delete any old text handles.
    if((size(handles.txtHandles,1)) >= nSyll)
        for(nTxt = 1:size(handles.txtHandles,2))
            if ishandle(handles.txtHandles(nSyll,nTxt))
                delete(handles.txtHandles(nSyll,nTxt));
                handles.txtHandles(nSyll,nTxt) = -1;
            end
        end
    end
    
    %Compute middle of syllable
    midSyll = (currAnnot.segFileEndTimes(nSyll) + currAnnot.segFileStartTimes(nSyll)) / 2;
    
    %create a new text
    axes(handles.axes_spec);
    handles.txtHandles(nSyll, 1) = text(midSyll,mean(ylim),num2str(label));
    
%     axes(handles.axesPower);
%     handles.txtHandles(nSyll, 2) = text(midSyll,mean(ylim),num2str(label));
%     axes(handles.axesSignal);
%     handles.txtHandles(nSyll, 3) = text(midSyll,mean(ylim),num2str(label));
    set(handles.txtHandles(nSyll,1), 'Color', 'black');
    set(handles.txtHandles(nSyll,1), 'FontSize', 14);
    set(handles.txtHandles(nSyll,1), 'FontWeight', 'bold');
    set(handles.txtHandles(nSyll,1), 'HitTest', 'off');
    set(handles.txtHandles(nSyll,1), 'HorizontalAlignment', 'center'); 
 
    %change rect colors to blue
    if(ishandle(handles.rectHandles(nSyll,:)))
        set(handles.rectHandles(nSyll,:),'FaceColor','blue','EdgeColor','none');
        set(handles.rectHandles(nSyll,2),'EdgeColor','blue');
    end
    
    %select next syllable and shift specgram axis.
    handles = setSelectedSyllable(handles, handles.selectedSyll + 1);
    if(currAnnot.segFileEndTimes(nSyll) > (handles.endNdx-1)/fs)
        width = handles.endNdx - handles.startNdx;
        handles.startNdx = max(0,round(round((midSyll*fs)+1) - width/2));
        handles.endNdx = min(length(handles.audio), round(round((midSyll*fs)+1) + width/2));
        handles = updateZoomOrPan(handles);
    end
end

function handles = updateAudAxes(handles)
% Update the pressure time series and spectrogram of the audio signal
axes(handles.axes_spec);
startTime = (handles.startNdx-1)/handles.fs;
if ~get(handles.check_audPower,'Value')
    displaySpecgramQuick(handles.audio(handles.startNdx:handles.endNdx), handles.fs,[0,10000],[],startTime);
    xlabel('');
    ylabel('');
else
    [pow, filtAud] = aSAP_getLogPower(handles.audio, handles.fs);
    plotTimeSeriesQuick(handles.time(handles.startNdx:handles.endNdx), pow(handles.startNdx:handles.endNdx), 'Color', [0, 0.4470, 0.7410]);
end

function handles  = updateDataAxes(handles)
% Update the data plots
% Get workspace name
hTweetVision = getappdata(0, 'hTweetVision');
handles.PreProcess = getappdata(hTweetVision, 'PreProcess');
handles.plotdata = [];
minval = 0;
maxval = -10;

for i = 1:4
    % Select a different axes on each pass through the loop
    axes(eval(['handles.axes_' num2str(i)]))
    
    lFreq = 0; hFreq = 10000;
    
    %Grab CommonMode (if desired)
    if get(handles.check_CMS, 'Value')
        ref = commonMode(handles.data, i);
        tmp = handles.data(i,:);
        data = tmp-ref;
    else
        data = handles.data(i,:);
    end
    
    % PreProcess data, if selected
    if get(handles.check_PreProcess,'Value')
%         PPdata = PreProcess(handles.data(i,:), handles.PreProcess(i));
        PPdata = PreProcess(data, handles.PreProcess(i));
        %The below is require to temporarily reshape the data for display
        %after downsampling; these changes are reversed at the end of this
        %routine
        
        if get(handles.check_ScaleFreq1,'Value')
            if ~isempty(find(handles.PreProcess(i).Steps==4, 1))
                lFreq = handles.PreProcess(i).HP.cutoff;
            end
            if ~isempty(find(handles.PreProcess(i).Steps==5, 1))
                hFreq = handles.PreProcess(i).LP.cutoff;
            end
        end
        
        if ~isempty(find(handles.PreProcess(i).Steps==8, 1))
            factor = floor(handles.PreProcess(i).params.Fs/handles.PreProcess(i).DSample.target);
            tempSNdx = handles.startNdx; tempENdx = handles.endNdx; tempTime = handles.time;
            handles.startNdx = max(floor(tempSNdx/factor),1);
            handles.endNdx = floor(tempENdx/factor);
            handles.time = downsample(handles.time,factor);
        end
    else
%         PPdata = handles.data(i,:);
        PPdata = data;
    end
    
    % Determine and plot according to the selected display option
    if get(handles.radio_power1,'Value')
        if ~get(handles.check_powDiff1,'Value')
            handles.plotdata(i,:) = smooth(PPdata.^2,441,'sgolay',4);
            %handles.plotdata(i,:) = smooth(PPdata.^2,4410,'moving');
        else
            plottemp = diff(smooth(PPdata.^2,441,'sgolay',4));
           % plottemp = diff(smooth(PPdata.^2,4410,'moving'));
            handles.plotdata(i,:) = [plottemp',0];
        end
        plotTimeSeriesQuick(handles.time(handles.startNdx:handles.endNdx), handles.plotdata(i,(handles.startNdx:handles.endNdx)), 'Color', [0, 0.4470, 0.7410]);
    elseif get(handles.radio_spectrogram1,'Value')
        handles.plotdata(i,:) = PPdata;
        startTime = (handles.startNdx-1)/handles.fs;
        displaySpecgramQuick(handles.plotdata(i,(handles.startNdx:handles.endNdx)), handles.fs,[lFreq,hFreq],[],startTime);
    elseif get(handles.radio_powSpec1,'Value')
        %handles.plotdata(i,:) = [diff(PPdata),0];
        %plotTimeSeriesQuick(handles.time(handles.startNdx:handles.endNdx), handles.plotdata(i,(handles.startNdx:handles.endNdx)));
        handles.plotdata(i,:) = PPdata;
        params = handles.PreProcess(i).params;
        params.fpass = [0 str2double(get(handles.edit_maxFreq1,'String'))];
        [S,f]=mtspectrumc(handles.plotdata(i,(handles.startNdx:handles.endNdx)),params);
        plot_vector(S,f,[])
    else
        handles.plotdata(i,:) = PPdata;
        plotTimeSeriesQuick(handles.time(handles.startNdx:handles.endNdx), handles.plotdata(i,(handles.startNdx:handles.endNdx)), 'Color', [0, 0.4470, 0.7410]);
    end
    minval = min(minval,min(handles.plotdata(i,(handles.startNdx:handles.endNdx))));
    maxval = max(maxval,max(handles.plotdata(i,(handles.startNdx:handles.endNdx))));
    
    %Return the values to match the orginal data
    if get(handles.check_PreProcess,'Value') && ~isempty(find(handles.PreProcess(i).Steps==8, 1))
        handles.startNdx = tempSNdx;
        handles.endNdx = tempENdx;
        handles.time = tempTime;
    end
end

set(handles.axes_1,'ButtonDownFcn','');set(handles.axes_1,'XTick',[]);
set(handles.axes_2,'ButtonDownFcn','');set(handles.axes_2,'XTick',[]);
set(handles.axes_3,'ButtonDownFcn','');set(handles.axes_3,'XTick',[]);
set(handles.axes_4,'ButtonDownFcn','');set(handles.axes_4,'XTick',[]);
set(handles.axes_raster1,'XTick',[],'YTick',[]);
set(handles.axes_raster2,'XTick',[],'YTick',[]);
set(handles.axes_raster3,'XTick',[],'YTick',[]);
set(handles.axes_raster4,'XTick',[],'YTick',[]);


if get(handles.check_syncScale,'value')
    range = [str2double(get(handles.edit_syncScale_min,'String')) str2double(get(handles.edit_syncScale_max,'String'))];
    axes(handles.axes_1); ylim([range]);
    axes(handles.axes_2); ylim([range]);
    axes(handles.axes_3); ylim([range]);
    axes(handles.axes_4); ylim([range]);
else
    set(handles.edit_syncScale_min,'String',num2str(minval));
    set(handles.edit_syncScale_max,'String',num2str(maxval));
end

function handles = updateZoomOrPan(handles)

fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/fs;
endTime = (handles.endNdx-1)/fs;
startNdx = handles.startNdx;
endNdx = handles.endNdx;

%fix navigator rect
axes(handles.axes_navBar);
if(ishandle(handles.navRect))
    delete(handles.navRect);
end
x = [startTime, startTime, endTime, endTime, startTime];
y = ylim;
y = [y(1), y(2), y(2), y(1), y(1)];
handles.navRect = line(x,y);
set(handles.navRect,'Color', 'red');
set(handles.navRect,'LineWidth', 2);
set(handles.navRect,'HitTest','off');

handles = updateAudAxes(handles);
xlabel('');
ylabel('');
%set(handles.axes_spec,'XTick',[]);

% Update data plots
handles = updateDataAxes(handles);

%Link all the x-axes together
linkaxes([handles.axes_1, handles.axes_raster1],'x');
linkaxes([handles.axes_2, handles.axes_raster2],'x');
linkaxes([handles.axes_3, handles.axes_raster3],'x');
linkaxes([handles.axes_4, handles.axes_raster4],'x');

set(get(handles.axes_navBar,'Parent'), 'KeyPressFcn', @cb_keypress);
set(handles.axes_navBar, 'ButtonDownFcn', @cb_navbar_click);
set(handles.axes_spec, 'ButtonDownFcn', @cb_specgram_click);
if get(handles.check_showSpikes,'Value')
    set(handles.axes_raster1, 'ButtonDownFcn', @cb_rasters_click);
    set(handles.axes_raster2, 'ButtonDownFcn', @cb_rasters_click);
    set(handles.axes_raster3, 'ButtonDownFcn', @cb_rasters_click);
    set(handles.axes_raster4, 'ButtonDownFcn', @cb_rasters_click);
end

handles.rectHandles(:,2) = -1;
handles.txtHandles(:,1) = -1;
handles = drawZoomSyllableRects(handles);

function handles = updateTemplates(handles)
h = figure(6969);
position = get(h,'Position');
handles.axesTemplates = axes;
cla(handles.axesTemplates);

if(isfield(handles, 'templates'))
    stdfs = handles.exper.desiredInSampRate;
    allwavs = [];
    edges = [0];
    for(nWav = 1:length(handles.templates.wavs))
        fs = handles.templates.wavs(nWav).fs;
        wav = handles.templates.wavs(nWav).wav;
        label = handles.templates.wavs(nWav).segType;
        wav = resample(wav, stdfs, fs);
        allwavs = [allwavs; wav'];
        edges = [edges, edges(end) + (length(wav) - 1)/stdfs];
    end
    
    %zero pad to 1 sec in length.
    if(edges(end) < 1)
        allwavs(end+1:stdfs) = 0;
    end
    
    %display the templates
    axes(handles.axesTemplates);
    displaySpecgramQuick(allwavs, stdfs, [0,10000]);
    
    %seperate with lines and label.
    for(nWav = 1:length(edges)-1)
        l = line([edges(nWav),edges(nWav)], ylim);
        set(l,'Color','k');
        set(l,'LineWidth', 2);
        mid = mean(edges([nWav,nWav+1]));
        if(handles.templates.wavs(nWav).segType ~= -1)
            t = text(mid, mean(ylim), num2str(handles.templates.wavs(nWav).segType), 'FontWeight', 'bold', 'FontSize', 20);
        end
    end
    
    set(handles.axesTemplates, 'ButtonDownFcn', @cb_templates_click);
    position(1,3) = 1024; position(1,4) = 180;
    set(h,'Position',position);
end

function handles = clearAll(handles)

for(nSyll = 1:max(size(handles.txtHandles,1),size(handles.rectHandles,1)))
    handles = deleteSyllableRectAndTxt(handles, nSyll);
end
handles.txtHandles = [];
handles.rectHandles = [];
handles.annotation.remove(handles.filelist(handles.curfile).name);
handles.selectedSyll = -1;

function handles = setSelectedSyllable(handles, newSelectedSyll)
%make sure selection is in bounds
currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
if(newSelectedSyll == -1)
    handles.selectedSyll = -1;
    return;
elseif(newSelectedSyll < 1 || newSelectedSyll > length(currAnnot.segAbsStartTimes))
    return
end

if(handles.selectedSyll ~= -1 && handles.selectedSyll <= length(currAnnot.segAbsStartTimes))
    %unhighlight the old selected syll
    if(currAnnot.segType(handles.selectedSyll) == -1)
        color = 'red';
    else
        color = 'blue';
    end
    if(size(handles.rectHandles,1) >= handles.selectedSyll)
        for(nRect = 1:size(handles.rectHandles,2))
            if(ishandle(handles.rectHandles(handles.selectedSyll,nRect)))
                set(handles.rectHandles(handles.selectedSyll,nRect),'FaceColor',color,'EdgeColor','none');
                if (nRect==2)
                    set(handles.rectHandles(handles.selectedSyll,nRect),'EdgeColor',color);
                end
            end
        end
    end
end

%set the selected syllable to return.
selectedSyll = newSelectedSyll;
handles.selectedSyll = selectedSyll;

%highlight the new selected syllable
if(size(handles.rectHandles,1) >= selectedSyll)
    for(nRect = 1:size(handles.rectHandles,2))
        if(ishandle(handles.rectHandles(selectedSyll,nRect)))
            set(handles.rectHandles(selectedSyll,nRect),'FaceColor','green','EdgeColor','none');
            if (nRect==2)
                set(handles.rectHandles(selectedSyll,nRect),'EdgeColor','green');
            end
        end
    end
end

%recenter if necessary
filenum = handles.filenum;
fs = handles.exper.desiredInSampRate;
startTime = (handles.startNdx-1)/fs;
endTime = (handles.endNdx-1)/fs;
syllEndTime = currAnnot.segFileEndTimes(selectedSyll);
syllStartTime = currAnnot.segFileStartTimes(selectedSyll);
if(syllEndTime > endTime || syllStartTime < startTime)
    width = handles.endNdx - handles.startNdx;
    mid = (syllStartTime*fs) + 1;
    handles.startNdx = round(mid - width/2);
    handles.endNdx = round(mid + width/2);
    if(handles.startNdx<1)
        handles.startNdx = 1;
        handles.endNdx = width + handles.startNdx;
    end
    if(handles.endNdx > length(handles.audio))
        handles.startNdx = length(handles.audio) - width;
        handles.endNdx = length(handles.audio);
    end        
    handles = updateZoomOrPan(handles);
end
    
function handles = loadfile(handles)
% Loads the file passed to it and returns the data in the handles structure

% Save the current annotation so no data is lost
handles.elements = saveAnnotation(handles);

% Save cell files, if they're loaded
for i = 1:4
    handles = saveCell(handles, i);
    handles.spikes{i}.sorted_spikes = [];
    handles.spikes{i}.unsorted_spikes = [];
    handles.spikes{i}.spont=[];
end

% Strip the channels out of the raw recording file
[rawdata, handles.fs] = getChannels(handles.filelist(handles.curfile).name);
%rawdata = rawdata(1,:)'; think i want all the channels

fprintf("Data frequency is %dHz\n", handles.fs);

if handles.fs ~= handles.exper.desiredInSampRate;
    warndlg('Experiment parameters and data file do not have the same sampling rate.');
    uiwait;
end
set(handles.text_curFilename,'String',handles.filelist(handles.curfile).name);
strparts = regexp(handles.filelist(handles.curfile).name,'_', 'split');
handles.filenum = str2double(strparts{2});

%Constants for Bandpass Audio (300-8000Hz)
HP_fNorm = 300/(44150/2);
LP_fNorm = 8000/(44150/2); %Changed 7/26/12 to match Farhan and Cengiz
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);

% Parse the data into the storage structures
handles.audio = filtfilt(BP_b,BP_a,rawdata(1,:));
handles.data = rawdata(2:end,:);
handles.time = 0:1/handles.fs:(length(handles.audio)-1)/handles.fs;

%Reset pltting indices for the new file
handles.startNdx = 1;
handles.endNdx = length(handles.audio);

% Update audio plots
handles.selectedSyll = -1;
handles.rectHandles = [];
handles.txtHandles = [];
axes(handles.axes_navBar);
plotTimeSeriesQuick(handles.time, handles.audio, 'Color', [0.8500, 0.3250, 0.0980]);

% handles = updateAudAxes(handles);
% 
% % Update data plots
% handles = updateDataAxes(handles);

%Draw rectangles on the nav bar...
if(handles.annotation.containsKey(handles.filelist(handles.curfile).name))
    currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
    handles.drugstatus=currAnnot.drugstatus;
    handles.directstatus=currAnnot.directstatus;
%     handles.threshold=currAnnot.syllthresh;
%     handles.edges=currAnnot.edgethresh;
%     handles.thresSyll=currAnnot.abssyllthresh;
%     handles.thresEdge=currAnnot.absedgethresh;
%     handles.drugindex=currAnnot.drugindex;
%     set(handles.editEdgeThreshold,'String',num2str(handles.edges));
%     set(handles.editPowerThreshold,'String',num2str(handles.threshold));
%     set(handles.editSyllAbs, 'String',num2str(handles.thresSyll));
%     set(handles.editEdgeAbs, 'String',num2str(handles.thresEdge));
    contents = get(handles.popup_drugstatus,'String');
    set(handles.popup_drugstatus,'Value',strmatch(handles.drugstatus,contents));
    
    contents = get(handles.popup_directstatus,'String');
    set(handles.popup_directstatus,'Value',strmatch(handles.directstatus,contents));
%      contents = get(handles.popupDrug,'String');% returns popupDrug contents as cell array
%      if (~strcmp(currAnnot.drugstatus,contents{handles.drugindex}))
%          errordlg('drug definitions have changed')
%      end
 
    syllStartTimes = currAnnot.segFileStartTimes;
    syllEndTimes = currAnnot.segFileEndTimes;
    x = [syllStartTimes; syllStartTimes; syllEndTimes; syllEndTimes; syllStartTimes];
    
    axes(handles.axes_navBar)
    lims = ylim;
    y = [lims(1); lims(2); lims(2); lims(1); lims(1)];
    for nSyll = 1:length(syllStartTimes)
        if(currAnnot.segType(nSyll) == -1)
            color = 'red';
        else
            color = 'blue';
        end
        handles.rectHandles(nSyll, 1) = patch(x(:,nSyll),y, color, 'EdgeColor', 'none', 'FaceAlpha', .2);
    end
    %if ~isempty(handles.rectHandles)
        set(handles.rectHandles(:,1), 'HitTest', 'off');
    %end
end
handles = updateZoomOrPan(handles);

if get(handles.check_showSpikes,'Value')
    for i = 1:4
        handles = ShowSpikes(handles, i);
        handles.cell{i}.filenum = handles.filenum;
        handles.cell{i}.filename = handles.filelist(handles.curfile).name;
    end
end

function elements = saveAnnotation(handles)
% Save the currently Annotation to file
if(isfield(handles,'annotFilename'))
    aaSaveHashtable(handles.annotFilename, handles.annotation);
    
    %Copy out the two main structures of the annotation file
    annotation = aaLoadHashtable(handles.annotFilename);

    if isempty(annotation)
        elements = [];
    else
        elements = annotation.elements;
    end
else
    elements = [];
    warndlg('You have not yet created an annotation.');
    uiwait;
end

function filenums = getStructFilenums(presCell)
    
 numRends =  size(presCell,2);
 
 if numRends ==1 && isempty(getfield(presCell(1), 'filenum'))
     filenums = [];
 else
     for i = 1:size(presCell,2)
        if isempty(getfield(presCell(i), 'filenum'))
            filenums(i) = 0;
        else
            filenums(i) = getfield(presCell(i), 'filenum');
        end
    end       
 end
 
 function filenums = getStructFields(presCell, field)
    
 numRends =  size(presCell,2);
 
 if numRends ==1 && isempty(getfield(presCell(1), field))
     filenums = [];
 else
     for i = 1:size(presCell,2)
        if isempty(getfield(presCell(i), field))
            filenums(i) = 0;
        else
            filenums(i) = getfield(presCell(i), field);
        end
    end       
end    
        
function handles = saveCell(handles,chan)
%Determine the cell file pntr, if any that selected for this channel
eval(['contents =  get(handles.popup_cellFile' num2str(chan) ',''String'');']);
eval(['select =  get(handles.popup_cellFile' num2str(chan) ',''Value'');']);
pntr = str2num(contents(select));

if(isfield(handles,'cell_allrecs')) && pntr > 0
    eval(['contents =  get(handles.popup_cellNum' num2str(chan) ',''String'');']);
    eval(['select =  get(handles.popup_cellNum' num2str(chan) ',''Value'');']);
    num = str2num(contents(select));
    presentCell=handles.cell_allrecs{pntr}; 
     
    %Check for empty records in the cell file and remove them
    Idx = [];
    for i = 1:numel(presentCell)
        Idx(i) = ~isempty(presentCell(i).filenum);
    end
    presentCell = presentCell(logical(Idx));
    
    % if there are no spikes for this record don't save anything
    if ~isfield(handles.spikes{chan},'sorted_spikes') || (isempty(handles.spikes{chan}.sorted_spikes) && isempty(handles.spikes{chan}.unsorted_spikes))
       
    else%if there are spikes for that record, save to file
        %insert the present record and then sort the structure based on filenumber before saving   
        handles.cell{chan}.spikes = handles.spikes{chan}.sorted_spikes;
        handles.cell{chan}.noise = handles.spikes{chan}.unsorted_spikes;
        handles.cell{chan}.channel = num2str(chan);
        handles.cell{chan}.cell_no = num;
        
        %Look up if current file is saved in the cell file
        if ~isempty(presentCell)
            filenumVector = getStructFilenums(presentCell);
            index = find(filenumVector == handles.filenum);
        else
            filenumVector = [];
            index = [];
        end
        
        if numel(index) > 1 %Is this really necessary? not in the original...
            index = index(end);
        end
                
        %if there is no record for the present file number then insert it appropriately
        if (isempty(index)) 
            %Find the insert location
            index = find(filenumVector>handles.filenum,1,'first');
            if (isempty(index))%present file is the latest for the cell; simply append
                presentCell(size(presentCell,2)+1) = handles.cell{chan};
            elseif (index == 1) % present file if the first one for the cell
               newpresentCell(1,1) = handles.cell{chan};
               newpresentCell(1,2:size(presentCell,2)+1) = presentCell;
               presentCell = newpresentCell;
            else
               newpresentCell(1,1:index-1) = presentCell(1:index-1);
               newpresentCell(1,index) = handles.cell{chan};
               newpresentCell(1,index+1:size(presentCell,2)+1) = presentCell(index:end);
               presentCell = newpresentCell;
            end
        else %override the old cell with the new one
            presentCell(1,index) = handles.cell{chan};        
        end   
            
    save([handles.cellPath{pntr},handles.cellName{pntr}], 'presentCell', '-v7');
    end
    handles.cell_allrecs{pntr} = presentCell; 
end

function [info] = getExperInfo
% Gather user input for the parameters
prompt = {'Enter bird name:','Experiment Description:','Sampling Rate:', 'Audio Channel:', 'Signal Channels:', 'Date Created:', 'Researcher:'};
dlg_title = 'Input experiment parameters';
num_lines = 1;
def = {'','Recording from Zebra Finch','44150', '0', '[1,2,3,4]', date, 'TMO'};
options.Resize='on';
answer = inputdlg(prompt,dlg_title,num_lines,def,options);
info.birdname = answer{1};
info.expername = answer{2};
info.desiredInSampRate = str2double(answer{3});
info.audioCh = str2double(answer{4});
info.sigCh = str2num(answer{5});
info.datecreated = answer{6};
info.researcher = answer{7};

function [filtAudio] = Prep(handles,audio)
if ~isempty(audio)

    %Get signal processing values from the GUI
    %these fields dont yet exist in the GUI
    gain = 1;%str2double(get(handles.edit_sigGain, 'String'));
    hp = 300;%str2double(get(handles.edit_sigHP, 'String'));
    lp = 8500;%str2double(get(handles.edit_sigLP, 'String'));

    %Constants for Bandpass Audio (300-8000Hz)
    HP_fNorm = hp/(handles.fs/2);
    LP_fNorm = lp/(handles.fs/2);
    [BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);

    renditions = size(audio,2);
    for i = 1:renditions
        if ~isempty(audio)
            filtAudio{i} = filtfilt(BP_b,BP_a,audio{i}.*gain);
        end
    end
end

function [syll, motif, audio] = getSubSet(handles,buffer,~)
%Given all of the previous filtering and selecting, this function will
%parse the data from the dataset, perform the requested alignments, and
%display it on both the audio and neuro axes.

%Parse data from the dataset into two matrices (audio and neuro)
%Filter the annotation and dataset
filt_elements = handles.elements; %(handles.filtInd==1);
filt_dataset = handles.rawDataset; %(handles.filtInd==1);

for i=1:size(handles.chosenStartSeq,1)
    %Grab the start and end times of each syllable in the motif 
    syll(1:2:handles.numSylls*2-1,i)=filt_elements{handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numSylls-1);
    syll(2:2:handles.numSylls*2,i)=filt_elements{handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2):handles.chosenStartSeq(i,2)+handles.numSylls-1);
    
    %Grab the start and end times of each motif
    motif(1,i)= filt_elements{handles.chosenStartSeq(i,1)}.segFileStartTimes(handles.chosenStartSeq(i,2));
    motif(2,i)= filt_elements{handles.chosenStartSeq(i,1)}.segFileEndTimes(handles.chosenStartSeq(i,2)+handles.numSylls-1);

    startT = max(floor(((motif(1,i)-buffer)*handles.fs)),1);
    endT = min(ceil(((motif(2,i)+buffer)*handles.fs)), length(filt_dataset{handles.chosenStartSeq(i,1)}(:)));
    
    
    %Grab the segment of data that corresponds to each selected motif
    audio{i} = filt_dataset{handles.chosenStartSeq(i,1)}(startT:endT);

end

function [SpecCube] = calcPower(data,type)
% Grab the PreProcessing data in the cloud
%hStretchEm = getappdata(0, 'hStretchEm');
%PP = getappdata(hStretchEm, 'PreProcess');

%Run the Chronux spec to grab power functions
movingwin=[0.005 0.001]; % FFT window size/ overlap

%Parameters for the audio
audPP.Fs = 44150;           % Sampling rate in Hz
audPP.tapers = [3 5];       % Tapers for spectral analysis
audPP.pad = 2;              % Zero padding (next 2^x)
audPP.fpass = [300 8000];   % Frequencies of interest
audPP.err = [0 0.05];       % Error handling
audPP.trialave = 0;         % Trial averaging

renditions = size(data,2);
if strcmp(type,'audio')
    for i=1:renditions
        [S,F,T,P] = spectrogram((data{i}/(sqrt(mean(data{i}.^2)))),220,220-44,512,44150);
        strt = 4; stp = 94; %These bins correspond to 300-8000Hz
        SpecCube{i} = abs(P(strt:stp,:));
    end
elseif strcmp(type,'neuro')
    for i=1:renditions
        % Calculate STFT features for both sounds (90% window overlap)
        [S,~,f]=mtspecgramc(data{i}',movingwin,PP(1).params);
        SpecCube{i} = S';
    end
end

function CM = commonMode(data,leaveOut)
%This function calculates the common mode of a set of simultaneously
%recorded ephys signals. Becuase there are expected to be few channels in
%the dataset (i.e., <10), I specify one channel to be left out of the
%calculations

%Dataset channels
chans = 2:size(data,1);

%Select the common channels
commonChans = ~(chans == leaveOut);

%Calculate common mode
% CM = mean(data(commonChans,:),1);
CM = median(data(commonChans,:),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% GUI Controls...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Panel controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_quarentine_Callback(hObject, eventdata, handles)
% hObject    handle to push_quarentine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'quarfiles')
    handles.quarfiles = handles.filelist_str(get(handles.list_filenames,'Value'));
    handles.filelist_str(get(handles.list_filenames,'Value')) = cellstr(['**' char(handles.filelist_str(get(handles.list_filenames,'Value'))) '**']);
    set(handles.text_message,'String',['File ' char(handles.filelist_str(get(handles.list_filenames,'Value'))) ' moved to quarentine.'])
else
    if sum(strcmp(handles.filelist_str(get(handles.list_filenames,'Value')),handles.quarfiles)) == 0
        handles.quarfiles = [handles.quarfiles; handles.filelist_str(get(handles.list_filenames,'Value'))];
        handles.filelist_str(get(handles.list_filenames,'Value')) = cellstr(['**' char(handles.filelist_str(get(handles.list_filenames,'Value'))) '**']);
        set(handles.text_message,'String',['File ' char(handles.filelist_str(get(handles.list_filenames,'Value'))) ' moved to quarentine.'])
    else
        set(handles.text_message,'String',['File ' char(handles.filelist_str(get(handles.list_filenames,'Value'))) ' is already in the quarentine.'])
    end
    
end

set(handles.list_filenames,'string',handles.filelist_str);

guidata(hObject, handles);

function push_unquarentine_Callback(hObject, eventdata, handles)
% hObject    handle to push_unquarentine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'quarfiles') || isempty(handles.quarfiles)
    set(handles.text_message,'String','There are no files to clear in quarentine.')
else
    handles.quarfiles = {};
    cd(handles.dirname)
    handles.filelist = dir('*.dat');
    
    % Populate listbox of filenames
    handles.filelist_str = cell(1,length(handles.filelist));
    for c = 1:length(handles.filelist)
        handles.filelist_str{c} = handles.filelist(c).name(1:end-4);
    end
    set(handles.list_filenames,'string',handles.filelist_str);
    set(handles.list_filenames,'value',handles.curfile);
    set(handles.text_message,'String','Quarentine cleared of files.  Nothing Deleted.')
end

guidata(hObject, handles);

function push_delete_Callback(hObject, eventdata, handles)
% hObject    handle to push_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles, 'quarfiles')
    set(handles.text_message,'String','There are no quarentined files to delete')
    return
else
    if length(handles.quarfiles)>0
        ansButton = questdlg('This will delete all quarentined files.  Are you sure?','Confirm Delete','Yep!','Shit! No!','Shit! No!');
    else
        set(handles.text_message,'String','There are no quarentined files to delete')
        return
    end
end

%Turn delete files to the recycle bin for possible recovery (if needed)
recycle('on');
switch ansButton
case {'Shit! No!'}
	return
case 'Yep!'
    for i=1:length(handles.quarfiles)
        delete([char(handles.quarfiles(i)) '.dat']);
    end
end

handles.quarfiles = [];

cd(handles.dirname)
handles.filelist = dir('*.dat');
handles.curfile = 1;

% Populate listbox of filenames
handles.filelist_str = cell(1,length(handles.filelist));
for c = 1:length(handles.filelist)
    handles.filelist_str{c} = handles.filelist(c).name(1:end-4);
end
set(handles.list_filenames,'string',handles.filelist_str);
set(handles.list_filenames,'value',handles.curfile);
set(handles.text_message,'String','Quarentined files deleted. You can still recover them from the recycle bin.')
[handles] = loadfile(handles);

guidata(hObject, handles);

function push_load_Callback(hObject, eventdata, handles)

directory_name = uigetdir;
%directory_name = 'C:\Users\Tim\Desktop\MATLAB Scripts\Tim''s Scripts\TweetVision';
if ~isstr(directory_name)
    return
else
    handles.dirname = directory_name;
    cd(handles.dirname)
    
    listDatNames = [];
    listStmNames = [];
    
    if get(handles.check_loadDATs, 'Value')
        listDat =dir('*.dat');
        for i = 1:length(listDat)
            % makes a cell array
            listDatNames{i} = listDat(i).name;
        end
    else
        listDat = [];
    end
    
    if get(handles.check_loadSTMs, 'Value')
        listStm = dir('*.stm');
        for i = 1:length(listStm)
            listStmNames{i} = listStm(i).name;
        end
    else
        listStm = [];
    end
    
    handles.filelist = [listDat; listStm];
    [~,indx] = sort([listDatNames,listStmNames]);
    handles.filelist = handles.filelist(indx);
    
    %handles.filelist = dir('*.dat');
    handles.curfile = 1;
    
    % Grab needed information from the filenames
    filelist_all = cell(1,length(handles.filelist));
    filelist_str = cell(1,length(handles.filelist));
    filelist_time = zeros(1,length(handles.filelist));
    filelist = handles.filelist;
    for c = 1:length(handles.filelist)
        filelist_all{c} = filelist(c).name;
        filelist_str{c} = filelist(c).name(1:end-4);
        sp = regexp(filelist(c).name(1:end-4),'_', 'split');
        filelist_time(c) = str2double(sp(6)) + str2double(sp(7))/60;
    end
    handles.filelist_all = filelist_all;
    handles.filelist_str = filelist_str;
    handles.filelist_time = filelist_time;
    
    %Place the conditional followed by java/HTML formatting call here
    richFileList = formatFileList(handles);
    
    %Push the values to the listbox
    set(handles.list_filenames,'string',richFileList);
    set(handles.list_filenames,'value',handles.curfile);
    
%     % Populate listbox of filenames
%     handles.filelist_str = cell(1,length(handles.filelist));
%     for c = 1:length(handles.filelist)
%         handles.filelist_str{c} = handles.filelist(c).name(1:end-4);
%     end
%     set(handles.list_filenames,'string',handles.filelist_str);
%     set(handles.list_filenames,'value',handles.curfile);
    
    handles.exper = getExperInfo;
    [handles] = loadfile(handles);
    set(handles.text_message,'String',['Files loaded from ' directory_name])
end
guidata(hObject, handles);

function richList = formatFileList(handles)
%This function will provide rich formatting for the filelist contents based on recent actions and contents of the annotation file

%Set the constants for all lines
numFiles = length(handles.filelist_str);
pre = '<html> ';
post = '</html>';

%Generate masks for marking up filelist
isAnnotated = [];
isErrored = [];
isBadSegmented = [];
isST = [];
stList = [];

if isfield(handles, 'annotation') && handles.annotationLoaded && ~isempty(handles.elements)
    isAnnotated = ismember(handles.filelist_all, handles.annotation.keys);
end

if isfield(handles, 'errorFiles')
    isErrored = ismember(handles.filelist_all, handles.errorFiles);
end

if isfield(handles, 'badSegFiles')
    isBadSegmented = ismember(handles.filelist_all, handles.badSegFiles);
end

%Append filelist with spikeTimes annotation
filelist_annot = annotateSpikeList(handles);

%Append filelist with completed/sorted data in the cell files (if those files are present)
if isfield(handles, 'cellName') && ~isempty(handles.cellName)
    filelist_annot = annotateCellList(handles, filelist_annot);
end

%Loop through the filelist and format each line by the markup mask
% richList = handles.filelist_str; %Original code; no additional annotation
richList = filelist_annot;
for i = 1:numFiles
    if ~isempty(isBadSegmented) && isBadSegmented(i) %Make segmentation errors italic
            richList{i} = ['<Font color="blue">'  richList{i} '</Font>'];
    end
    
    if ~isempty(isAnnotated) && isAnnotated(i) %Make annotated files bold and green highlighted
         richList{i} = ['<Body bgcolor="green"><b>' richList{i} '</Body>'];
    end
    
    if ~isempty(isErrored) && isErrored(i) %Make errored files red
         richList{i} = ['<Font color="red">'  richList{i} '</Font>'];
    end

    richList{i} = [pre  richList{i} post];
    
end

function filelist_annot = annotateSpikeList(handles)
%Copyover output
filelist_annot = handles.filelist_str;

%Retrieve the list of spikeTime files
spikelist = dir([handles.dirname filesep 'times_*.mat']);
if ~isempty(spikelist)
spikelist_base = cell(1,length(spikelist)); %base filename
spikelist_str = cell(1,length(spikelist)); %chan list to add
%For each, parse the base name and the sorted channel
idx = 1;
for c = 1:length(spikelist)
    if ~strcmp(spikelist_base{idx},spikelist(c).name(7:end-9)) && c>1
        %If this is a new base name, increment counter and copynew data
        idx = idx + 1;
        spikelist_base{idx} = spikelist(c).name(7:end-9);
        spikelist_str{idx} = spikelist(c).name(end-4);
    elseif c==1
        spikelist_base{idx} = spikelist(c).name(7:end-9);
        spikelist_str{idx} = spikelist(c).name(end-4);
    else
        %If this is not new, add the new channel to the current list
        spikelist_str{idx} = [spikelist_str{idx} '+'  spikelist(c).name(end-4)];
    end
end

%Eliminate extraneous, empty arrays
spikelist_base = spikelist_base(1:idx);
spikelist_str = spikelist_str(1:idx);

%Look up each base string the filelist and append the sorted annotation
for c = 1:length(spikelist_base)
    pntr = ~cellfun(@isempty,strfind(handles.filelist_str,spikelist_base{c}));
    target = find(pntr==1,1,'first');
    if ~isempty(target)
        filelist_annot{target} = [handles.filelist_str{target} '     ( ' spikelist_str{c} ' )'];
    end
end
end

function filelist_annot = annotateCellList(handles, filelist_annot)
fnames = []; chans = [];

%Retrieve the list of completed files from the cell files
numCells = numel(handles.cell_allrecs);
for i = 1:numCells
    if ~isempty(handles.cell_allrecs{i})
        %Get vectors from this cell file
        filename = getFieldVectorCell(handles.cell_allrecs{i},'filename');
        chan = getFieldVectorCell(handles.cell_allrecs{i},'channel');

        %Add to the running lists
        fnames = [fnames, filename];
        chans = [chans, chan];
    end
end
[fnames, idx] = sort(fnames);
chans = chans(idx);

if isempty(fnames)
    filelist_annot = filelist_annot;
else
    cellList_base = cell(1,numel(fnames)); %base filename
    chanList_str = cell(1,numel(fnames)); %chan list to add

    %For each, parse the base name and the sorted channel
    idx = 1;
    for c = 1:length(fnames)
        if ~strcmp(cellList_base{idx}, fnames{c}(1:end-4)) && c>1
            %If this is a new base name, increment counter and copynew data
            idx = idx + 1;
            cellList_base{idx} = fnames{c}(1:end-4);
            chanList_str{idx} = chans{c};
        elseif c==1
            cellList_base{idx} = fnames{c}(1:end-4);
            chanList_str{idx} = chans{c};
        else
            %If this is not new, add the new channel to the current list
            chanList_str{idx} = [chanList_str{idx} '+' chans{c}];
        end
    end

    %Eliminate extraneous, empty arrays
    cellList_base = cellList_base(1:idx);
    chanList_str = chanList_str(1:idx);

    %Look up each base string the filelist and append the sorted annotation
    for c = 1:length(cellList_base)
        pntr = ~cellfun(@isempty,strfind(handles.filelist_str, cellList_base{c}));
        target = find(pntr==1,1,'first');
        if ~isempty(target)
            filelist_annot{target} = [filelist_annot{target} '     = ' chanList_str{c} ' ='];
        end
    end
end

function push_previous_Callback(hObject, eventdata, handles)
% hObject    handle to push_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'filelist')
    return
elseif isempty(handles.filelist)
    return
else
    handles.curfile = handles.curfile - 1;
    if handles.curfile < 1
        handles.curfile = 1;
    end
    [handles] = loadfile(handles);
    
    %Place the conditional followed by java/HTML formatting call here
    richFileList = formatFileList(handles);
    
    %Push the values to the listbox
    set(handles.list_filenames,'string',richFileList);
    set(handles.list_filenames,'value',handles.curfile);
end
guidata(hObject, handles);

function push_next_Callback(hObject, eventdata, handles)
% hObject    handle to push_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'filelist')
    return
elseif isempty(handles.filelist)
    return
else
    handles.curfile = handles.curfile + 1;
    if handles.curfile > length(handles.filelist)
        handles.curfile = length(handles.filelist);
    end
     [handles] = loadfile(handles);
    
    %Place the conditional followed by java/HTML formatting call here
    richFileList = formatFileList(handles);
    
    %Push the values to the listbox
    set(handles.list_filenames,'string',richFileList);
    set(handles.list_filenames,'value',handles.curfile);
end
guidata(hObject, handles);

function list_filenames_Callback(hObject, eventdata, handles)

if strcmp(get(gcf,'selectiontype'),'open')
    handles.curfile = get(handles.list_filenames,'value');
else
    return
end

handles = loadfile(handles);

%Place the conditional followed by java/HTML formatting call here
richFileList = formatFileList(handles);

%Push the values to the listbox
set(handles.list_filenames,'string',richFileList);
set(handles.list_filenames,'value',handles.curfile);

guidata(hObject, handles);

function list_filenames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_filenames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annotation Panel controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_deletesyll_Callback(hObject, eventdata, handles)
% hObject    handle to push_deletesyll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.selectedSyll ~= -1)
    nSyll = handles.selectedSyll;
    handles = deleteSyllableRectAndTxt(handles, nSyll);
    if(size(handles.txtHandles,1) >= nSyll)
        handles.txtHandles(nSyll,:) = [];
    end
    if(size(handles.rectHandles,1) >= nSyll)
        handles.rectHandles(nSyll,:) = [];
    end
    currAnnotation = handles.annotation.get(handles.filelist(handles.curfile).name);
    currAnnotation.segAbsStartTimes(nSyll) = [];
    currAnnotation.segFileStartTimes(nSyll) = [];
    currAnnotation.segFileEndTimes(nSyll) = [];
    currAnnotation.segType(nSyll) = [];
    selectedSyll = handles.selectedSyll;
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnotation);
    if(length(currAnnotation.segAbsStartTimes)== 0)
        selectedSyll = -1;
    elseif(selectedSyll > length(currAnnotation.segAbsStartTimes));
        selectedSyll = selectedSyll - 1;
    end    
    
    %Remove the file from the errorlist, if it's in there.
    if isfield(handles, 'errorFiles')
        handles.errorFiles = removeFile(handles.filelist(handles.curfile).name, handles.errorFiles);
    end
    
    %Remove the file from the bad segmentation list, if it's in there.
    if isfield(handles, 'badSegFiles') && ~ismember(104,currAnnotation.segType)
        handles.badSegFiles = removeFile(handles.filelist(handles.curfile).name, handles.badSegFiles);
    end
    
    handles = setSelectedSyllable(handles, selectedSyll);    
    guidata(hObject, handles);
end

%function push_clearall_Callback(hObject, eventdata, handles)
% hObject    handle to push_clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles=clearAll(handles);
%guidata(hObject, handles);

function push_deletepause_Callback(hObject, eventdata, handles)
% hObject    handle to push_deletepause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.selectedSyll ~= -1)
    filenum = handles.filenum;
    nSyll = handles.selectedSyll;
    nSyllNext = nSyll + 1;
    currAnnotation = handles.annotation.get(handles.filelist(handles.curfile).name);
    numSyll = length(currAnnotation.segAbsStartTimes);
    if(nSyll<numSyll)
        %if there is a syllable beyond the selected one, then delete
        %the pause in between.
        %delete the rectangles and text for both syllables.
        handles = deleteSyllableRectAndTxt(handles,nSyll);
        handles = deleteSyllableRectAndTxt(handles,nSyllNext);
        if((size(handles.txtHandles,1)) >= nSyllNext)
            handles.txtHandles(nSyllNext,:) = [];
        end
        if((size(handles.rectHandles,1)) >= nSyllNext)
            handles.rectHandles(nSyllNext,:) = [];
        end       
        
        currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
        currAnnot.segAbsStartTimes(nSyllNext) = [];
        currAnnot.segFileStartTimes(nSyllNext) = [];
        currAnnot.segFileEndTimes(nSyll) = [];
        currAnnot.segType(nSyll) = -1;
        currAnnot.segType(nSyllNext) = [];
        handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
        handles = drawSyllableRect(handles, nSyll);
        if(ishandle(handles.rectHandles(nSyll,:)))
            set(handles.rectHandles(nSyll,:),'FaceColor','green','EdgeColor','green');
        end
        guidata(hObject, handles);
    end
end

function push_newannot_Callback(hObject, eventdata, handles)
% hObject    handle to push_newannot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%If annotation is already open check whether they want to append or clear.
if isfield(handles, 'annotFilename')
    ansButton = questdlg('Would you like to save the current annotations and begin a new one, or append to the existing annotation?','Annotation File Already Loaded','Save & Clear','Append','Append');
end

% Get the new location for the annotation file
[file,path] = uiputfile('*.mat', 'Create a .mat file for the annotation:');

if ~(isequal(file,0) || isequal(path,0))
    if isfield(handles,'annotFilename') && strcmp(ansButton, 'Save & Clear')
        %If there is already have an annotation object loaded, then save & 
        %clear current annotation.
        saveAnnotation(handles);
        guidata(hObject, handles);
        push_clearall_Callback(hObject, eventdata, handles);
        handles = guidata(hObject);
        
        %Delete the currently loaded hashtable
        handles.annotation.delete;        
        
    elseif(~isfield(handles,'annotFilename'))
        %Create new annotation hashtable
        handles.annotation = mhashtable;
    end
    handles.annotFilename = [path,file];
    saveAnnotation(handles);
end
guidata(hObject, handles);

function push_loadannot_Callback(hObject, eventdata, handles)

[file,path] = uigetfile('*.mat', 'Create a .mat file for the audioAnnotation:');

if ~(isequal(file,0) || isequal(path,0))
    if isfield(handles,'annotFilename') 
        ansButton = questdlg('There is already an annotation loaded. You want to save this and load another or continue with the current one?','Annotation File Already Loaded','Save & Load','Use Current','Use Current');
        if strcmp(ansButton, 'Save & Load')
            %Save the current copy of annotation to file
            handles.elements = saveAnnotation(handles);
            guidata(hObject, handles);
            push_clearall_Callback(hObject, eventdata, handles);
            handles = guidata(hObject);

            %Delete old and create new hashtables
            handles.annotation.delete;
            handles.annotation = mhashtable;
            handles.annotationLoaded = true;
            
            handles = loadfile(handles);
            richFileList = formatFileList(handles);
            
            %Push the values to the listbox
            set(handles.list_filenames,'string',richFileList);
            set(handles.list_filenames,'value',handles.curfile);
            
        elseif strcmp(ansButton, 'Use Current')
            %Do nothing (for now)
        end
    else
        ansButton ='';
    end
        %Update path to annotation file and load in to memory
        handles.annotFilename = [path file];
        annotation = aaLoadHashtable(handles.annotFilename);
        handles.annotation = annotation;
        handles.annotationLoaded = true;
        
%         %If the currently loaded file has a record in the annotation, load
%         %the record
%         if(handles.annotation.containsKey(handles.filelist(handles.curfile).name))
%             currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
%             %set(handles.popupDrug,'Value',currAnnot.drugindex);
%         end

        %Update the GUI plots to show the syllable boundaries on the Nav
        %and spec bars
        handles = loadfile(handles);

        %Create a list of files with segmentation problems (likely only present if the AnnotDisplay program was used to
        %create the annotation file)
        keys = handles.annotation.keys;
        handles.badSegFiles = {};
        for i = 1:length(handles.annotation.keys)
            if ismember(104, handles.elements{i}.segType)
                handles.badSegFiles{end+1} = keys{i};
            end
        end

        richFileList = formatFileList(handles);
        
        %Push the values to the listbox
        set(handles.list_filenames,'string',richFileList);
        set(handles.list_filenames,'value',handles.curfile);
        
        handles = drawNavBarSyllableRects(handles);
        handles = drawZoomSyllableRects(handles);
        
    if strcmp(ansButton, 'Use Current')
        set(handles.text_message,'String',['Continuing to use existing annotation from ' handles.annotFilename]);
    else
        set(handles.text_message,'String',['Annotation loaded from ' handles.annotFilename]);
    end
end
guidata(hObject, handles);

function push_savetemplate_Callback(hObject, eventdata, handles)
% hObject    handle to push_savetemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(isfield(handles,'templates'))
    templates=handles.templates;
    %templateFeatures=handles.templateFeatures;
    [FileName,PathName]=uiputfile('*.mat','Save Template As:','template');
    save([PathName FileName],'templates');
else
    warndlg('No templates to save.');
    uiwait;
end

function push_viewtemplate_Callback(hObject, eventdata, handles)
% hObject    handle to push_viewtemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'templates')
    ansButton = questdlg('A template file is already loaded. Keep this one or load another?','Template Already Loaded','Keep & View','Load Another','Keep & View');
    if strcmp(ansButton,'Load Another')
        uiload;
        handles.templates = templates;
    end
else
    uiload;
    handles.templates = templates;
end
handles = updateTemplates(handles);
guidata(hObject, handles);

function push_addsyllable_Callback(hObject, eventdata, handles)
% hObject    handle to push_addsyllable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'annotFilename')
    handles.bWaitingForAddClick = true;
    handles.bWaitingForSeparationClick = false;
    set(handles.push_addsyllable, 'BackgroundColor', 'red');
else
    set(handles.text_message,'String','There is no annotation file! Create or load annotation file to add syllable.')
end
guidata(hObject, handles);

function check_edgefind_Callback(hObject, eventdata, handles)

function check_audPower_Callback(hObject, eventdata, handles)

handles = updateZoomOrPan(handles)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Panel controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  PRETTY SURE A LOT OF THESE AREN'T USED
function check_syncScale_Callback(hObject, eventdata, handles)
% hObject    handle to check_syncScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_syncScale

handles = updateDataAxes(handles);
guidata(hObject, handles);

function edit_syncScale_max_Callback(hObject, eventdata, handles)

function edit_syncScale_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_syncScale_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_syncScale_min_Callback(hObject, eventdata, handles)

function edit_syncScale_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_syncScale_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_Update1_Callback(hObject, eventdata, handles)
% hObject    handle to push_Update1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get workspace name
hTweetVision = getappdata(0, 'hTweetVision');
handles.PreProcess = getappdata(hTweetVision, 'PreProcess');

handles = updateDataAxes(handles);

guidata(hObject, handles);

function push_ProcWind1_Callback(hObject, eventdata, handles)
% hObject    handle to push_ProcWind1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get workspace name
hTweetVision = getappdata(0, 'hTweetVision');
setappdata(hTweetVision, 'PreProcess', handles.PreProcess);

TweetVisProcOptions(handles.PreProcess,1,'hTweetVision');

guidata(hObject, handles);

function push_saveData_Callback(hObject, eventdata, handles)
% hObject    handle to push_saveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Setup the data to be saved
if isfield(handles,'plotdata')
    saveAnnotation(handles);
    dataset.filename = handles.filelist(handles.curfile).name;
    dataset.audio = handles.audio;
    dataset.neural = handles.plotdata;
    
    if isfield(handles,'PreProcess') && get(handles.check_PreProcess,'Value')
        dataset.PP = handles.PreProcess;
    else
        dataset.PP = [];
    end
    
    if handles.annotation.containsKey(handles.filelist(handles.curfile).name)
        dataset.annotation = handles.annotation.get(handles.filelist(handles.curfile).name);
    else
        dataset.annotation = [];
    end
    
    [file,path] = uiputfile('*.mat','Save data as...');
    save([path file],'dataset')
    set(handles.text_message,'String',['Data saved to ' path file])
else
    set(handles.text_message,'String','No data available to save. Load files and retry.')
end

function push_corr_Callback(hObject, eventdata, handles)
% hObject    handle to push_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Gather user input for the parameters
prompt = {'Signal 1:','Signal 2:'};
dlg_title = 'Correlation of which channels?';
num_lines = 1;
def = {'audio','1'};
options.Resize='on';
answer = inputdlg(prompt,dlg_title,num_lines,def,options);

if strcmp(lower(answer{1}),'audio')
    [sig1, filtAud] = aSAP_getLogPower(handles.audio, handles.fs);
else
    pntr = str2double(answer{1});
    sig1 = handles.plotdata(pntr,:);
end
if strcmp(lower(answer{2}),'audio')
    [sig2, filtAud] = aSAP_getLogPower(handles.audio, handles.fs);
else
    pntr = str2double(answer{2});
    sig2 = handles.plotdata(pntr,:);
end
figure;
if strcmp(answer{1},answer{2})
    [c,lags] = xcorr(sig1);
else
    [c,lags] = xcorr(sig1,sig2);
end
plot((lags/handles.fs),(c/max(abs(c))))
title(['Correlation of channel' answer{1} ' and channel ' answer{2} '.'])
xlabel('Time lag (in sec)')

function check_PreProcess_Callback(hObject, eventdata, handles)
% This isnt used??
handles = updateDataAxes(handles);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axis  Panel controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E-Phys Data Display
function radio_voltage1_Callback(hObject, eventdata, handles)
% hObject    handle to radio_voltage1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_voltage1

if get(handles.radio_voltage1,'Value')
    set(handles.radio_power1,'Value',0);
    set(handles.radio_spectrogram1,'Value',0);
    set(handles.radio_powSpec1,'Value',0);
    
    handles = updateDataAxes(handles);
end

guidata(hObject, handles);

function radio_power1_Callback(hObject, eventdata, handles)
% hObject    handle to radio_power1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_power1
if get(handles.radio_power1,'Value')
    set(handles.radio_voltage1,'Value',0);
    set(handles.radio_spectrogram1,'Value',0);
    set(handles.radio_powSpec1,'Value',0);
    
    handles = updateDataAxes(handles);
end

guidata(hObject, handles);

function radio_spectrogram1_Callback(hObject, eventdata, handles)
% hObject    handle to radio_spectrogram1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_spectrogram1
if get(handles.radio_spectrogram1,'Value')
    set(handles.radio_power1,'Value',0);
    set(handles.radio_voltage1,'Value',0);
    set(handles.radio_powSpec1,'Value',0);
    
    handles = updateDataAxes(handles);
end

guidata(hObject, handles);

function radio_powSpec1_Callback(hObject, eventdata, handles)
% hObject    handle to radio_powSpec1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_powSpec1
if get(handles.radio_powSpec1,'Value')
    set(handles.radio_power1,'Value',0);
    set(handles.radio_voltage1,'Value',0);
    set(handles.radio_spectrogram1,'Value',0);
    
    handles = updateDataAxes(handles);
end

guidata(hObject, handles);

function check_ScaleFreq1_Callback(hObject, eventdata, handles)
% hObject    handle to check_ScaleFreq1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_ScaleFreq1
handles = updateDataAxes(handles);

guidata(hObject, handles);

function check_powDiff1_Callback(hObject, eventdata, handles)

handles = updateDataAxes(handles);

guidata(hObject, handles);

function edit_maxFreq1_Callback(hObject, eventdata, handles)
    
if str2double(get(handles.edit_maxFreq1,'String'))>(handles.fs/2)
    warndlg(['Upper limit cannot be higher than Nyquist Frequency: ' int2str(floor(handles.fs/2)) '.']);
    uiwait;
    set(handles.edit_maxFreq1,'String',int2str(floor(handles.fs/2)))
end
handles = updateDataAxes(handles);
guidata(hObject, handles);

function edit_maxFreq1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial Params Panel controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popup_drugstatus_Callback(hObject, eventdata, handles)

contents = get(handles.popup_drugstatus,'String');
handles.drugstatus = contents{get(handles.popup_drugstatus,'Value')};

%Check to see if there is an annotation for the current file
if handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
    currAnnot.drugstatus=handles.drugstatus;
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
end
guidata(hObject, handles);

function popup_drugstatus_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_directstatus_Callback(hObject, eventdata, handles)

contents = get(handles.popup_directstatus,'String');
handles.directstatus = contents{get(handles.popup_directstatus,'Value')};

%Check to see if there is an annotation for the current file
if handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
    currAnnot.directstatus=handles.directstatus;
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
end
guidata(hObject, handles);

function popup_directstatus_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike Sorting controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_cellNames_Callback(hObject, eventdata, handles)

function listbox_cellNames_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    
function check_showSpikes_Callback(hObject, eventdata, handles)

if get(handles.check_showSpikes,'Value')
    for i = 1:4
        handles = ShowSpikes(handles, i);
    end
end
guidata(hObject, handles);

function push_newCell_Callback(hObject, eventdata, handles)
%Function creates a new cell structure and saves to user specified location

%Determine which cell structure index is the next to be filled
 if ~isfield(handles, 'cellName')
     handles.cellName = [];
 end
 if isempty(handles.cellName)
     pntr = 1;
 else
        target = cellfun(@isempty,handles.cellName);
        pntr = find(target==1,1,'first');
        if isempty(pntr)
            pntr = length(handles.cellName)+1;
        end
 end
 
 %Copy to the permenant structure
handles.cell{1}.birdname = handles.exper.birdname;
handles.cell{2}.birdname = handles.exper.birdname;
handles.cell{3}.birdname = handles.exper.birdname;
handles.cell{4}.birdname = handles.exper.birdname;
 
 %Create the baseline cell structure
cell.birdname = handles.exper.birdname;
cell.filename = [];
cell.channel = [];
cell.cell_no = [];
cell.filenum = [];
cell.spikes = [];
cell.noise = [];
cell.spont = [];

%Request the location for the cell file to be saved
 [handles.cellName{pntr},handles.cellPath{pntr}]=uiputfile('cell.mat','Save cell as:');
 
 %Find the currently loaded cell pntrs
handles.cellSlots = find(~cellfun(@isempty,handles.cellName));

%Copy the base structure to the "completed" spike array
presentCell(1) = cell; 
handles.cell_allrecs{pntr} = presentCell;

%Save the file to the specified location
save ([handles.cellPath{pntr},handles.cellName{pntr}], 'presentCell');

%Display the name of the new cell file in the list box
contents = genCellList(handles.cellName);
set(handles.listbox_cellNames,'String',contents);

%Push new dropdown options to the file popups without changing the selection
for i = 1:4
    %Retrieve the old values
    cont = eval(['get(handles.popup_cellFile' num2str(i) ',''String'')']);
    select = eval(['get(handles.popup_cellFile' num2str(i) ',''Value'')']);
    if isempty(select)
        pick = 0;
    else
        pick = str2num(cont(select));
    end

    %Find the new values and pointers
    newCont = ([0; handles.cellSlots'])';
    newVal = find(newCont==pick,1,'first');

    %Pump them to the GUI
    eval(['set(handles.popup_cellFile' num2str(i) ',''String'', newCont)'])
    eval(['set(handles.popup_cellFile' num2str(i) ',''Value'', newVal)'])
end

guidata(hObject, handles);

function contents = genCellList(cellNames)
%From the stored list of cell names, generate a cell array of names for the listbox
contents = [];
if ~iscell(cellNames)
    contents = ['1)   ' char(cellNames)];
else
    ind = length(cellNames);
    j = 1;
    for i = 1:ind
        if ~isempty(cellNames{i})
            contents{j} = [num2str(i) ')   ' char(cellNames{i})];
            j = j+1;
        end
    end
end

function spikesName = getSpikeTimesFname(curFilename,chan)
%Generate the spike times filename and check if it exists
spikesName = ['times_' curFilename(1:end-4) 'chan' num2str(chan) '.mat'];

if exist(spikesName,'file') ~= 2
    spikesName = [];
end

function handles=loadSpikes(handles,chan)
%This function handles loading and formatting spikes data from spike times files

%This figures out what the matching times file name would be for the current file and returns 0 if it doesn't exist
%handles.spike_filename = getspikeDatafile(handles.exper,handles.filenum,handles.chan);
handles.spikes{chan}.filename = getSpikeTimesFname(handles.filelist(handles.curfile).name,chan);

%Clear the insertion point
handles.spikes{chan}.marker=[];

%If it doesn't exist, bail
 if isempty(handles.spikes{chan}.filename)
     %Clear CellNum popup box
%     eval(['set(handles.popup_cellNum' num2str(chan) ',''String'' ,''-'');'])
%     eval(['set(handles.popup_cellNum' num2str(chan) ',''Value'' ,1);'])

 else
    %Determine which of the loaded cell files has been user selected for the channel
    %pntr = handles.cellSlots(get(handles.popup_cellFile4, 'Value'));
    
    %Otherwise, load the file into memory and log filename
    load(handles.spikes{chan}.filename);
    handles.cell{chan}.filename=handles.spikes{chan}.filename;
    if isempty(cluster_class)
        set(handles.text_message, 'String',['The Spike Times file ' handles.spikes{chan}.filename ' is empty... you should delete it.']);
        beep
        return
    end
    
    %Covert spike times into seconds
    cluster_class(:,2)=cluster_class(:,2)/1000;
    
    %Determine how many cell cluster classes are in the times files and push to popup box
    popup_string=(1:1:(max(cluster_class(:,1))));
    eval(['set(handles.popup_cellNum' num2str(chan) ',''String'' ,popup_string);'])
    if ~ismember(handles.cell_no(chan), popup_string)
        handles.cell_no(chan) = 1;
    end
    eval(['set(handles.popup_cellNum' num2str(chan) ',''Value'', handles.cell_no(chan));'])
    %For the user selected cell (from this times file), copy over the good spikes to the cell file...
    cellindex=find(cluster_class(:,1)==handles.cell_no(chan));
    handles.spikes{chan}.sorted_spikes=cluster_class(cellindex,2);
    handles.cell{chan}.spikes=handles.spikes{chan}.sorted_spikes;
    
    %...and the potential spikes...
    cellindex=find(cluster_class(:,1)~=handles.cell_no(chan)); %fix this
    cellindex=sort(cellindex);
    handles.spikes{chan}.unsorted_spikes=cluster_class(cellindex,2);
    handles.cell{chan}.noise=handles.spikes{chan}.unsorted_spikes;
    
    %Zero out the spontaneous spikes... don't know if I'll use this.
    handles.cell{chan}.spont=[];
    
    %Plot the "good spikes" on the raster
    eval(['axes(handles.axes_raster' num2str(chan) ')'])
%     for i=1:length(handles.spikes{chan}.sorted_spikes)
%         line([handles.spikes{chan}.sorted_spikes(i) handles.spikes{chan}.sorted_spikes(i)],[0.3 1],'color','r');
%     end
    %New single-line version of plot
    [tSx, tSy]  = rasterLine(handles.spikes{chan}.sorted_spikes, 0.3, 1);
    line(tSx, tSy, 'color', 'r');   
    hold on;
    
    %Plot the "potential spikes" on the raster
%     y = zeros(length(handles.spikes{chan}.unsorted_spikes));
%     plot(handles.spikes{chan}.unsorted_spikes,y,'b+');
    y = zeros(length(handles.spikes{chan}.unsorted_spikes),1);
    props = {'LineStyle','none','Marker','+','MarkerEdge','b','MarkerSize',6};
    line([handles.spikes{chan}.unsorted_spikes, handles.spikes{chan}.unsorted_spikes],[y,y],props{:});
    ylim([-0.3 1]);
    hold off;
end

function handles=ShowSpikes(handles,chan)
%This is the primary function for controlling spike plotting... it calls everything else.

%Select and clear the rasters axis
eval(['axes(handles.axes_raster' num2str(chan) ')'])
cla;

%If checkbox is on...
if get(handles.check_showSpikes,'Value')
    %Determine which of the loaded cell files has been user selected for the channel
    %eval(['pntr = handles.cellSlots(get(handles.popup_cellFile' num2str(chan) ', ''Value''));'])
    
    eval(['pntr = get(handles.popup_cellFile' num2str(chan) ', ''Value'');'])
    pntr = pntr - 1; %Adjust for the null file
    if pntr > 0
        pntr = handles.cellSlots(pntr);
    end
    
    if (isfield(handles,'cell_allrecs')) && pntr > 0
        
        if ~isempty(handles.cell_allrecs{pntr})
            %Get the list of filenumbers that have already been sorted, and look up index for currently loaded DAT file
            if size(handles.cell_allrecs{pntr},2)>=1
                %fieldVector = getFieldVector(handles.cell_allrecs{pntr},'filenum');
                presentCell = handles.cell_allrecs{pntr};
                fieldVector = getStructFilenums(handles.cell_allrecs{pntr});
                index = find(fieldVector == handles.filenum');
                if length(index) > 1
                    index = index(end);
                end
                
                %If the cell record for this file, pop the spike file over
                if ~isempty(index)
                    handles.cell{chan}=presentCell(index);
                end
                
            end

            %If spikes have been sorted, check the cell number to make it match what's in the file; if not, update the popup
            if (~isempty(index) && (str2double(handles.cell{chan}.channel) == chan) && ~isempty(handles.cell_allrecs{pntr}(index).spikes)) && ~(handles.cell{chan}.cell_no == handles.cell_no(chan))
                %Update the cell pointer
                handles.cell_no(chan) = handles.cell{chan}.cell_no;
                
                %Update the popup menu
                eval(['set(handles.popup_cellNum' num2str(chan) ', ''String'', 1:handles.cell_no(chan));'])
                eval(['set(handles.popup_cellNum' num2str(chan) ', ''Value'', handles.cell_no(chan));'])
%                 a = 1
            end
            
            %If the spikes in this file have already been sorted
            if (~isempty(index) && (str2double(handles.cell{chan}.channel) == chan) && (handles.cell{chan}.cell_no == handles.cell_no(chan)) && ~isempty(handles.cell_allrecs{pntr}(index).spikes))
                
                %Load the cell file spike/noise times
                handles.spikes{chan}.sorted_spikes = handles.cell_allrecs{pntr}(index).spikes;
                handles.spikes{chan}.unsorted_spikes = handles.cell_allrecs{pntr}(index).noise;
                handles.cell{chan}.spikes = handles.spikes{chan}.sorted_spikes;
                handles.cell{chan}.noise = handles.spikes{chan}.unsorted_spikes;

                %If there's no spikes in the file, throw a warning and exit
                 if (isempty(handles.cell{chan}.spikes))
                    warndlg('no spikes!');
                    return;
                 end

                %Plot the "good spikes" on the raster
                %axes(handles.axesRasters);
                eval(['axes(handles.axes_raster' num2str(chan) ')'])
                
                %Old version of raster plot
%                 for i=1:length(handles.spikes{chan}.sorted_spikes)
%                     line([handles.spikes{chan}.sorted_spikes(i) handles.spikes{chan}.sorted_spikes(i)],[0.3 1],'color','k');
%                 end
                
                %New single-line version of plot
                [tSx, tSy]  = rasterLine(handles.spikes{chan}.sorted_spikes, 0.3, 1);
                line(tSx, tSy, 'color', 'k');
                
                hold on;

                %Plot the "potential spikes" on the raster
%                 y = zeros(length(handles.spikes{chan}.unsorted_spikes));
%                 plot(handles.spikes{chan}.unsorted_spikes,y,'b+');
                y = zeros(length(handles.spikes{chan}.unsorted_spikes),1);
                props = {'LineStyle','none','Marker','+','MarkerEdge','b','MarkerSize',6};
                line([handles.spikes{chan}.unsorted_spikes, handles.spikes{chan}.unsorted_spikes],[y,y],props{:});
                ylim([-0.3 1]);
                hold off;
            else 
                %If the file has not yet been sorted, procede to loading spikes from the spike times file (if present)
                handles = loadSpikes(handles,chan);
            end
        else 
            %If the file has not yet been sorted, procede to loading spikes from the spike times file (if present)
            handles = loadSpikes(handles,chan);
        end
    else
        %If there is nothing in the cell file or the file is not selected, then procede to loading spikes from spike times file
        %handles = loadSpikes(handles,chan); 
     end
else
    % Checkbox is not checked-take appropriate action
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell File Handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_loadCell_Callback(hObject, eventdata, handles)
%Function loads an already created cell structure from file

 %Determine which cell structure index is the next to be filled
 if ~isfield(handles, 'cellName')
     handles.cellName = [];
 end
 if isempty(handles.cellName)
     pntr = 1;
 else
        target = cellfun(@isempty,handles.cellName);
        pntr = find(target==1,1,'first');
        if isempty(pntr)
            pntr = length(handles.cellName)+1;
        end
 end

%   %See if the basic structure exists yet
%  if ~isfield(handles, 'cell')
%     %Create the baseline cell structure
%     cell.birdname = handles.exper.birdname;
%     cell.filename = [];
%     cell.channel = [];
%     cell.cell_no = [];
%     cell.filenum = [];
%     cell.spikes = [];
%     cell.noise = [];
%     cell.spont = [];
% 
%     %Copy to the permenant structure
%     handles.cell{1} = cell;
%     handles.cell{2} = cell;
%     handles.cell{3} = cell;
%     handles.cell{4} = cell;
%  end
 
 %Copy to the permenant structure
handles.cell{1}.birdname = handles.exper.birdname;
handles.cell{2}.birdname = handles.exper.birdname;
handles.cell{3}.birdname = handles.exper.birdname;
handles.cell{4}.birdname = handles.exper.birdname;
 
%Load the file into memory
[handles.cellName{pntr},handles.cellPath{pntr}]=uigetfile('*.mat');
load([handles.cellPath{pntr}, handles.cellName{pntr}]);

if ~isvarname('presentCell')
        %If this isn't a cell file, boot the process
        set(handles.text_message,'String','Selected cell file is not of the expected structure... check the file and try again.')
        beep
else
    %Check to see if this is a cell file for the bird currently loaded
    if ~strcmp(presentCell(1).birdname,handles.exper.birdname)
        %If they aren't the same, request permission from the user to continue loading file
        ansButton = questdlg('The bird names specified in the experiment file and celll file don''t match. Load anyway?','Birdname mismatch','Load it anyway','Nevermind','Nevermind');
        switch ansButton
            case 'Nevermind'
                %Jump out
                return
            case 'Load it anyway'
                %Then procede with the below code
        end
    end
    
    
%     if presentCell(1).channel==handles.chan || presentCell(1).cell_no==handles.cell_no
         %make the old cell files up to date with the spont field
%          if ~isfield(presentCell,'spont')
%             for i=1:size(presentCell,2) 
%                 presentCell(i).spont=[];
%             end
%          end
         
         %Standardize the cell file order and fields (old cell files had fucked up structure and need to be reordered)
         q = [];
         for i=1:numel(presentCell)
             q(i).birdname = presentCell(i).birdname;
             q(i).filename = presentCell(i).filename;
             q(i).channel = presentCell(i).channel;
             q(i).cell_no = presentCell(i).cell_no;
             q(i).filenum = presentCell(i).filenum;
             q(i).spikes = presentCell(i).spikes;
             q(i).noise = presentCell(i).noise;
             if ~isfield(presentCell,'spont')
                 q(i).spont = [];
             else
                 q(i).spont = presentCell(i).spont;
             end
         end
         presentCell = q;
        
        %Check for empty records in the cell file and remove them
         Idx = [];
         for i = 1:numel(presentCell)
             Idx(i) = ~isempty(presentCell(i).filenum);
         end
         presentCell = presentCell(logical(Idx));
         
        %Copy the loaded spike data to the "completed" structure
        handles.cell_allrecs{pntr} = presentCell;
 
         %Determine whether there is already a spike set for the current DAT file, if so, load it to
         filenum = getFieldVector(presentCell,'filenum');
         index=find(filenum==handles.filenum);
         if length(index) > 1
            index = index(end);
         end
         if ~isempty(index)
             handles.cell{pntr}=presentCell(index);
             for i = 1:4
                handles=ShowSpikes(handles,i);
             end
         end
    
    %Copy the loaded spike data to the "completed" structure
    handles.cell_allrecs{pntr} = presentCell;
    
    %Find the currently loaded cell pntrs
    handles.cellSlots = find(~cellfun(@isempty,handles.cellName));
    
    %Display the name of the new cell file in the list box
    contents = genCellList(handles.cellName);
    set(handles.listbox_cellNames,'String',contents);

    %Push new dropdown options to the file popups without changing the selection
    for i = 1:4
        %Retrieve the old values
        cont = eval(['get(handles.popup_cellFile' num2str(i) ',''String'')']);
        select = eval(['get(handles.popup_cellFile' num2str(i) ',''Value'')']);
        if isempty(select) || select == 0
            pick = 0;
        else
            pick = str2num(cont(select));
        end

        %Find the new values and pointers
        newCont = ([0; handles.cellSlots'])';
        newVal = find(newCont==pick,1,'first');
        if isempty(newVal)
            newVal = 0;
        end

        %Pump them to the GUI
        eval(['set(handles.popup_cellFile' num2str(i) ',''String'', newCont)'])
        eval(['set(handles.popup_cellFile' num2str(i) ',''Value'', newVal)'])
    end
end


if get(handles.check_showSpikes,'Value')
    for i = 1:4
        handles = ShowSpikes(handles, i);
    end
end

guidata(hObject, handles);

function push_unloadCell_Callback(hObject, eventdata, handles)
%Function unloads the cell structure highlighted in the listbox

%Get the highlighted index and filename
usrSelect = get(handles.listbox_cellNames,'Value');
aux = get(handles.listbox_cellNames,'String');
usrFilename = char(aux(usrSelect));
pntr = str2num(usrFilename(1));

%Check to see to make sure that this really oought to be unloaded
ansButton = questdlg(['Sure you want to unload ' usrFilename(5:end)],'Confirm unload','Hell yes!','Nevermind','Nevermind');
switch ansButton
    case 'Nevermind'
        %Jump out
        return
    case 'Hell yes!'
        %Then procede with the below code
end

%Save the cell info to file
%presentCell(end)=handles.cell{pntr}; 
%handles.cell_allrecs=presentCell;

%Save the file to the specified location
presentCell = handles.cell_allrecs{pntr};
save ([handles.cellPath{pntr},handles.cellName{pntr}], 'presentCell');

%Remove all traces of the cell in memory
%handles.cell{pntr} = [];
handles.cellPath{pntr} = [];
handles.cellName{pntr} = [];
handles.cell_allrecs{pntr} = [];

%Find the currently loaded cell pntrs
handles.cellSlots = find(~cellfun(@isempty,handles.cellName));

%Display the name of the new cell file in the list box
contents = genCellList(handles.cellName);
set(handles.listbox_cellNames,'String',contents,'Value',1);

%Push new dropdown options to the file popups
for i = 1:4
    %Retrieve the old values
    cont = eval(['get(handles.popup_cellFile' num2str(i) ',''String'')']);
    select = eval(['get(handles.popup_cellFile' num2str(i) ',''Value'')']);
    if isempty(select)  || select == 0
        pick = 0;
    else
        pick = str2num(cont(select));
    end

    %Find the new values and pointers
    newCont = ([0; handles.cellSlots'])';
    newVal = find(newCont==pick,1,'first');
    if isempty(newVal)
        newVal = 0;
    end

    %Pump them to the GUI
    eval(['set(handles.popup_cellFile' num2str(i) ',''String'', newCont)'])
    eval(['set(handles.popup_cellFile' num2str(i) ',''Value'', newVal)'])
end

guidata(hObject, handles);

function handles = cellFileSelect(handles, chan)
%Only execute if the spikes are displayed
if get(handles.check_showSpikes,'Value')
    
    %Get the file selection fromt he popup box
    eval(['contents = get(handles.popup_cellFile' num2str(chan) ', ''String'');'])
    eval(['pntr = str2num(contents(get(handles.popup_cellFile' num2str(chan) ', ''Value'')));'])
    if pntr == 0
        %Do somthing to clear the channel
        eval(['set(handles.popup_cellNum' num2str(chan) ',''String'',''0'',''Value'',1);'])
        handles.cell_no(chan) = 1;

    else
        %Check to see if the cell file exists,
        if isfield(handles, 'cell')
             %Copy the loaded spike data to temp location
            presentCell = handles.cell_allrecs{pntr};
            
            %Determine whether there is already a spike set for the current DAT file, if so, load it to
            filenum = getFieldVector(presentCell,'filenum');
            index=find(filenum==handles.filenum);
            if length(index) > 1
                index = index(end);
            end
            if ~isempty(index)
                handles.cell{chan}=presentCell(index);
            end
            
           %Push the set of unique cell numbers in the file to the popup
            cell_no = unique(getFieldVector(presentCell,'cell_no'));
            if isempty(cell_no)
                eval(['set(handles.popup_cellNum' num2str(chan) ',''String'',''0'',''Value'',1);'])
            else
                eval(['set(handles.popup_cellNum' num2str(chan) ',''String'', [' num2str(cell_no) '] ,''Value'',1);'])
            end
            handles.cell_no(chan) = 1;
            
        end
    end
    
    %Update spikes
    handles = ShowSpikes(handles, chan);
    
end

function handles = addSpike(handles, chan)
%If there's a spike marker and unsorted spikes to add, find the unsorteds that you're closeted to
if (~isempty(handles.spiketime) && ~isempty(handles.spikes{chan}.unsorted_spikes))
    added_spike1 = find (handles.spikes{chan}.unsorted_spikes<handles.spiketime, 1,'last');
    added_spike2 = find (handles.spikes{chan}.unsorted_spikes>handles.spiketime, 1,'first');
    if isempty(added_spike2)
        added_spike=added_spike1;
    elseif (abs(handles.spikes{chan}.unsorted_spikes(added_spike1)-handles.spiketime)<abs(handles.spikes{chan}.unsorted_spikes(added_spike2)-handles.spiketime))
        added_spike=added_spike1;
    else
        added_spike=added_spike2;
    end

    %Find the insertion point in the sorted spikes array
    added_spike_index=find(handles.spikes{chan}.sorted_spikes<handles.spikes{chan}.unsorted_spikes(added_spike), 1,'last');

    %Insert and remove where required
    handles.spikes{chan}.sorted_spikes=[handles.spikes{chan}.sorted_spikes(1:added_spike_index); handles.spikes{chan}.unsorted_spikes(added_spike); handles.spikes{chan}.sorted_spikes((added_spike_index+1):end)];
    handles.spikes{chan}.unsorted_spikes=[handles.spikes{chan}.unsorted_spikes(1:added_spike-1); handles.spikes{chan}.unsorted_spikes(added_spike+1:end)];
    handles.spiketime=[];

    %Plot the results
    eval(['subplot(handles.axes_raster' num2str(chan) ')'])
    cla;
%     for i=1:length(handles.spikes{chan}.sorted_spikes)
%         line([handles.spikes{chan}.sorted_spikes(i) handles.spikes{chan}.sorted_spikes(i)],[0.3 1],'color','r');
%     end
    %New single-line version of plot
    [tSx, tSy]  = rasterLine(handles.spikes{chan}.sorted_spikes, 0.3, 1);
    line(tSx, tSy, 'color', 'r');
    hold on;
    
%     y = zeros(length(handles.spikes{chan}.unsorted_spikes));
%     plot(handles.spikes{chan}.unsorted_spikes,y,'b+');
    y = zeros(length(handles.spikes{chan}.unsorted_spikes),1);
    props = {'LineStyle','none','Marker','+','MarkerEdge','b','MarkerSize',6};
    line([handles.spikes{chan}.unsorted_spikes, handles.spikes{chan}.unsorted_spikes],[y,y],props{:});
    hold off;

    %Update the record cell files
    handles.spikemarker=[];
    handles.cell{chan}.spikes=handles.spikes{chan}.sorted_spikes;
    handles.cell{chan}.noise=handles.spikes{chan}.unsorted_spikes;
end

function handles = deleteSpike(handles, chan)
%Only run if spikesmarker time is loaded
if (~isempty(handles.spiketime))
    %Find the closest spike before and aft
    added_spike1 = find (handles.spikes{chan}.sorted_spikes<handles.spiketime, 1,'last');
    added_spike2 = find (handles.spikes{chan}.sorted_spikes>handles.spiketime, 1,'first');

    %Choose the closer one as the one to be deleted
    if (abs(handles.spikes{chan}.sorted_spikes(added_spike1)-handles.spiketime)<abs(handles.spikes{chan}.sorted_spikes(added_spike2)-handles.spiketime));
        added_spike=added_spike1;
    else
        added_spike=added_spike2;
    end

    %Find where in the unsorted list this newly deleted spike will be inserted in to
    added_spike_index=find(handles.spikes{chan}.unsorted_spikes<handles.spikes{chan}.sorted_spikes(added_spike), 1,'last');
    if isempty(added_spike_index)
        added_spike_index=1;
    end

    %Add to the unsorted list and subtract from the sorted
    if isempty(handles.spikes{chan}.unsorted_spikes)
        handles.spikes{chan}.unsorted_spikes=handles.spikes{chan}.sorted_spikes(added_spike);
    else
        handles.spikes{chan}.unsorted_spikes=[handles.spikes{chan}.unsorted_spikes(1:added_spike_index); handles.spikes{chan}.sorted_spikes(added_spike); handles.spikes{chan}.unsorted_spikes((added_spike_index+1):end)];
    end
    handles.spikes{chan}.sorted_spikes=[handles.spikes{chan}.sorted_spikes(1:added_spike-1); handles.spikes{chan}.sorted_spikes(added_spike+1:end)];

    %Clear out the marker position
    handles.spiketime=[];

    %Plot the results
    eval(['subplot(handles.axes_raster' num2str(chan) ')'])
    cla;
%     for i=1:length(handles.spikes{chan}.sorted_spikes)
%         line([handles.spikes{chan}.sorted_spikes(i) handles.spikes{chan}.sorted_spikes(i)],[0.3 1],'color','r');
%     end
    %New single-line version of plot
    [tSx, tSy]  = rasterLine(handles.spikes{chan}.sorted_spikes, 0.3, 1);
    line(tSx, tSy, 'color', 'r');
    hold on;
    
%     y = zeros(length(handles.spikes{chan}.unsorted_spikes));
%     plot(handles.spikes{chan}.unsorted_spikes,y,'b+');
    y = zeros(length(handles.spikes{chan}.unsorted_spikes),1);
    props = {'LineStyle','none','Marker','+','MarkerEdge','b','MarkerSize',6};
    line([handles.spikes{chan}.unsorted_spikes, handles.spikes{chan}.unsorted_spikes],[y,y],props{:});
    hold off;
end

%Update the current record cell structures
handles.cell{chan}.spikes=handles.spikes{chan}.sorted_spikes;
handles.cell{chan}.noise=handles.spikes{chan}.unsorted_spikes;
handles.spikemarker=[];

function handles = forceSpike(handles, chan)
%If there's a spike marker and unsorted spikes to add, find the unsorteds that you're closeted to
if ~isempty(handles.spiketime) %&& ~isempty(handles.spikes{chan}.unsorted_spikes))
    
    %Find the insertion point in the sorted spikes array
    added_spike_index=find(handles.spikes{chan}.sorted_spikes<handles.spiketime, 1,'last');
    
    %Insert and remove where required
    if ~isempty(added_spike_index)
        handles.spikes{chan}.sorted_spikes = [handles.spikes{chan}.sorted_spikes(1:added_spike_index); handles.spiketime; handles.spikes{chan}.sorted_spikes((added_spike_index+1):end)];
    else
        %If it's the first spike, just append to the beginning
        handles.spikes{chan}.sorted_spikes = [handles.spiketime; handles.spikes{chan}.sorted_spikes];
    end
    
    handles.spiketime=[];
    
    %Plot the results
    eval(['subplot(handles.axes_raster' num2str(chan) ')'])
    cla;
%      for i=1:length(handles.spikes{chan}.sorted_spikes)
%          line([handles.spikes{chan}.sorted_spikes(i) handles.spikes{chan}.sorted_spikes(i)],[0.3 1],'color','r');
%      end
    %New single-line version of plot
    [tSx, tSy]  = rasterLine(handles.spikes{chan}.sorted_spikes, 0.3, 1);
    line(tSx, tSy, 'color', 'r');
                
    hold on;
%     y = zeros(length(handles.spikes{chan}.unsorted_spikes));
%     plot(handles.spikes{chan}.unsorted_spikes,y,'b+');
    y = zeros(length(handles.spikes{chan}.unsorted_spikes),1);
    props = {'LineStyle','none','Marker','+','MarkerEdge','b','MarkerSize',6};
    line([handles.spikes{chan}.unsorted_spikes, handles.spikes{chan}.unsorted_spikes],[y,y],props{:});
    hold off;
    
    %Update the record cell files
    handles.spikemarker=[];
    handles.cell{chan}.spikes=handles.spikes{chan}.sorted_spikes;
    handles.cell{chan}.noise=handles.spikes{chan}.unsorted_spikes;
end

function handles = reloadSpikes(handles, chan)
%Reload the spiketimes from the default spiketimes file location
handles=loadSpikes(handles,chan);

%Plot the results
eval(['subplot(handles.axes_raster' num2str(chan) ')'])
cla;
% for i=1:length(handles.spikes{chan}.sorted_spikes)
%     line([handles.spikes{chan}.sorted_spikes(i) handles.spikes{chan}.sorted_spikes(i)],[0.3 1],'color','r');
% end
%New single-line version of plot
[tSx, tSy]  = rasterLine(handles.spikes{chan}.sorted_spikes, 0.3, 1);
line(tSx, tSy, 'color', 'r');
                
hold on;
% y = zeros(length(handles.spikes{chan}.unsorted_spikes));
% plot(handles.spikes{chan}.unsorted_spikes,y,'b+');
y = zeros(length(handles.spikes{chan}.unsorted_spikes),1);
props = {'LineStyle','none','Marker','+','MarkerEdge','b','MarkerSize',6};
line([handles.spikes{chan}.unsorted_spikes, handles.spikes{chan}.unsorted_spikes],[y,y],props{:});
hold off;

%Update the record cell files
handles.spikemarker=[];
handles.cell{chan}.spikes=handles.spikes{chan}.sorted_spikes;
handles.cell{chan}.noise=handles.spikes{chan}.unsorted_spikes;

function handles = clearSpike(handles, chan)

%Add all the "sorted spikes" to the unsorted list and subtract from the sorted
if isempty(handles.spikes{chan}.unsorted_spikes)
    handles.spikes{chan}.unsorted_spikes=handles.spikes{chan}.sorted_spikes;
else
    handles.spikes{chan}.unsorted_spikes=sort([handles.spikes{chan}.unsorted_spikes; handles.spikes{chan}.sorted_spikes]);
end
handles.spikes{chan}.sorted_spikes=NaN;

%Clear out the marker position
handles.spiketime=[];

%Plot the results
eval(['subplot(handles.axes_raster' num2str(chan) ')'])
cla;

%New single-line version of plot
[tSx, tSy]  = rasterLine(handles.spikes{chan}.sorted_spikes, 0.3, 1);
line(tSx, tSy, 'color', 'r');
hold on;

y = zeros(length(handles.spikes{chan}.unsorted_spikes),1);
props = {'LineStyle','none','Marker','+','MarkerEdge','b','MarkerSize',6};
line([handles.spikes{chan}.unsorted_spikes, handles.spikes{chan}.unsorted_spikes],[y,y],props{:});
hold off;

%Update the current record cell structures
handles.cell{chan}.spikes=handles.spikes{chan}.sorted_spikes;
handles.cell{chan}.noise=handles.spikes{chan}.unsorted_spikes;
handles.spikemarker=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell File Pushbutton calls (real control code above)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_addSpike4_Callback(hObject, eventdata, handles)
%Channel constant
chan = 4;

%Common control file
handles = addSpike(handles, chan);

guidata(hObject, handles);
    
function push_deleteSpike4_Callback(hObject, eventdata, handles)
%Channel constant
chan = 4;

%Common control file
handles = deleteSpike(handles, chan);

guidata(hObject, handles);

function popup_cellNum4_Callback(hObject, eventdata, handles)
%Set cell number (with indirect address)
t = get(hObject,'Tag');
handles.cell_no(str2num(t(end))) = get(hObject,'Value');

handles=ShowSpikes(handles,str2num(t(end)));

guidata(hObject, handles);

function popup_cellNum4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_cellFile4_Callback(hObject, eventdata, handles)
    %ID chan number (with indirect address)
t = get(hObject,'Tag'); 
chan = str2num(t(end));

%Common control code
handles = cellFileSelect(handles, chan);

guidata(hObject, handles);

function popup_cellFile4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_addSpike3_Callback(hObject, eventdata, handles)
%Channel constant
chan = 3;

%Common control file
handles = addSpike(handles, chan);

guidata(hObject, handles);
    
function push_deleteSpike3_Callback(hObject, eventdata, handles)
%Channel constant
chan = 3;

%Common control file
handles = deleteSpike(handles, chan);

guidata(hObject, handles);

function popup_cellNum3_Callback(hObject, eventdata, handles)
%Set cell number (with indirect address)
t = get(hObject,'Tag');
handles.cell_no(str2num(t(end))) = get(hObject,'Value');

handles=ShowSpikes(handles,str2num(t(end)));

guidata(hObject, handles);

function popup_cellNum3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_cellFile3_Callback(hObject, eventdata, handles)
%ID chan number (with indirect address)
t = get(hObject,'Tag'); 
chan = str2num(t(end));

%Common control code
handles = cellFileSelect(handles, chan);

guidata(hObject, handles);

function popup_cellFile3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_addSpike2_Callback(hObject, eventdata, handles)
%Channel constant
chan = 2;

%Common control file
handles = addSpike(handles, chan);

guidata(hObject, handles);
    
function push_deleteSpike2_Callback(hObject, eventdata, handles)
%Channel constant
chan = 2;

%Common control file
handles = deleteSpike(handles, chan);

guidata(hObject, handles);
    
function popup_cellNum2_Callback(hObject, eventdata, handles)
%Set cell number (with indirect address)
t = get(hObject,'Tag');
handles.cell_no(str2num(t(end))) = get(hObject,'Value');

handles=ShowSpikes(handles,str2num(t(end)));

guidata(hObject, handles);

function popup_cellNum2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_cellFile2_Callback(hObject, eventdata, handles)
%ID chan number (with indirect address)
t = get(hObject,'Tag'); 
chan = str2num(t(end));

%Common control code
handles = cellFileSelect(handles, chan);

guidata(hObject, handles);

function popup_cellFile2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_addSpike1_Callback(hObject, eventdata, handles)
%Channel constant
chan = 1;

%Common control file
handles = addSpike(handles, chan);

guidata(hObject, handles);
    
function push_deleteSpike1_Callback(hObject, eventdata, handles)
%Channel constant
chan = 1;

%Common control file
handles = deleteSpike(handles, chan);

guidata(hObject, handles);
    
function popup_cellNum1_Callback(hObject, eventdata, handles)
%Set cell number (with indirect address)
t = get(hObject,'Tag');
handles.cell_no(str2num(t(end))) = get(hObject,'Value');

handles=ShowSpikes(handles,str2num(t(end)));

guidata(hObject, handles);

function popup_cellNum1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_cellFile1_Callback(hObject, eventdata, handles)
%ID chan number (with indirect address)
t = get(hObject,'Tag'); 
chan = str2num(t(end));

%Common control code
handles = cellFileSelect(handles, chan);

guidata(hObject, handles);
    
function popup_cellFile1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_reload1_Callback(hObject, eventdata, handles)
%Reload the raw spiking information from the spikeTimes file (generated by wav_clus)

%Channel constant
chan = 1;

%Common control code
handles = reloadSpikes(handles, chan);
    
guidata(hObject, handles);

function push_force1_Callback(hObject, eventdata, handles)
%Force the addition of a spike exactly where the marker is, instead of searching for the nearest "suspected" spike

%Channel constant
chan = 1;

%Common control code
handles = forceSpike(handles, chan);

guidata(hObject, handles);

function push_reload2_Callback(hObject, eventdata, handles)
%Reload the raw spiking information from the spikeTimes file (generated by wav_clus)

%Channel constant
chan = 2;

%Common control code
handles = reloadSpikes(handles, chan);
    
guidata(hObject, handles);

function push_force2_Callback(hObject, eventdata, handles)
%Force the addition of a spike exactly where the marker is, instead of searching for the nearest "suspected" spike

%Channel constant
chan = 2;

%Common control code
handles = forceSpike(handles, chan);

guidata(hObject, handles);

function push_reload3_Callback(hObject, eventdata, handles)
%Reload the raw spiking information from the spikeTimes file (generated by wav_clus)

%Channel constant
chan = 3;

%Common control code
handles = reloadSpikes(handles, chan);
    
guidata(hObject, handles);

function push_force3_Callback(hObject, eventdata, handles)
%Force the addition of a spike exactly where the marker is, instead of searching for the nearest "suspected" spike

%Channel constant
chan = 3;

%Common control code
handles = forceSpike(handles, chan);

guidata(hObject, handles);

function push_reload4_Callback(hObject, eventdata, handles)
%Channel constant
chan = 4;

%Common control code
handles = reloadSpikes(handles, chan);

guidata(hObject, handles);

function push_force4_Callback(hObject, eventdata, handles)
%Force the addition of a spike exactly where the marker is, instead of searching for the nearest "suspected" spike

%Channel constant
chan = 4;

%Common control code
handles = forceSpike(handles, chan);

guidata(hObject, handles);

function push_clear1_Callback(hObject, eventdata, handles)
%Clear all the spikes intot he 'noise' bin

%Channel constant
chan = 1;

%Common control code
handles = clearSpike(handles, chan);

guidata(hObject, handles);

function push_clear2_Callback(hObject, eventdata, handles)
%Clear all the spikes intot he 'noise' bin

%Channel constant
chan = 2;

%Common control code
handles = clearSpike(handles, chan);

guidata(hObject, handles);

function push_clear3_Callback(hObject, eventdata, handles)
%Clear all the spikes intot he 'noise' bin

%Channel constant
chan = 3;

%Common control code
handles = clearSpike(handles, chan);

guidata(hObject, handles);

function push_clear4_Callback(hObject, eventdata, handles)
%Clear all the spikes intot he 'noise' bin

%Channel constant
chan = 4;

%Common control code
handles = clearSpike(handles, chan);

guidata(hObject, handles);


% --- Executes on button press in check_CMS.
function check_CMS_Callback(hObject, eventdata, handles)
% hObject    handle to check_CMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_CMS

% --- Executes on button press in push_multiSegment.
%function push_multiSegment_Callback(hObject, eventdata, handles)

%fprintf("multiseg pushed\n");
%Call the multisegmenting routine
%[currAnnot, handles] = multiSegment(handles, handles.audio, []);

%fprintf("multiseg routine successfully run");

%Update the annotation folder
%handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);

%frintf("annotation folder updated"); 

%Remove the file from the errorlist, if it's in there.
%if isfield(handles, 'errorFiles')
%    handles.errorFiles = removeFile(handles.filelist(handles.curfile).name, handles.errorFiles);
%end

%Remove the file from the bad segmentation list, if it's in there.
%if isfield(handles, 'badSegFiles') && ~ismember(104,currAnnot.segType)
%    handles.badSegFiles = removeFile(handles.filelist(handles.curfile).name, handles.badSegFiles);
%end

%draw syllable information
%handles = drawNavBarSyllableRects(handles);
%handles = drawZoomSyllableRects(handles);

%flag
%handles.lastMode = 'multi';

%guidata(hObject, handles);

function listOut = removeFile(file, listIn)
listOut = listIn;
t = ismember(listIn, file);
if ~isempty(t)
    indx = find(t==1);
    listOut(indx) = [];
end

function push_clearall_Callback(hObject, eventdata, handles)
% hObject    handle to push_clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=clearAll(handles);
guidata(hObject, handles);

% --- Executes on button press in push_resetParams.
function push_resetParams_Callback(hObject, eventdata, handles)
%populate edit values
for i = 1:length(handles.editObj)
    handles = populateEdit(handles.editObj{i}, handles);
end
    
    %populate slider values
for i = 1:length(handles.sliderObj)
    handles = populateSlider(handles.sliderObj{i}, handles);
end

% Update handles structure
guidata(hObject, handles);
% This should open a second GUI in which the segmenting parameters can be
% edited


% --- Executes on button press in push_frameSegment.
function push_frameSegment_Callback(hObject, eventdata, handles)
%Call the frame segmenting routine
[currAnnot, handles] = frameSegment(handles, handles.audio);

%Update the annotation folder
handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);

%Remove the file from the errorlist, if it's in there.
if isfield(handles, 'errorFiles')
    handles.errorFiles = removeFile(handles.filelist(handles.curfile).name, handles.errorFiles);
end

%Remove the file from the bad segmentation list, if it's in there.
if isfield(handles, 'badSegFiles') && ~ismember(104,currAnnot.segType)
    handles.badSegFiles = removeFile(handles.filelist(handles.curfile).name, handles.badSegFiles);
end


%draw syllable information
handles.rectHandles = [];
handles.txtHandles = [];
handles = drawNavBarSyllableRects(handles);
handles = drawZoomSyllableRects(handles);

%flag
handles.lastMode = 'frame';

guidata(hObject, handles);


% --- Executes on button press in check_powerEnv.
function check_powerEnv_Callback(hObject, eventdata, handles)
handles = updateZoomOrPan(handles);
guidata(hObject, handles);


% --- Executes on button press in check_dynamic.
function check_dynamic_Callback(hObject, eventdata, handles)
% hObject    handle to check_dynamic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_dynamic


% --- Executes on button press in check_segGuide.
function check_segGuide_Callback(hObject, eventdata, handles)
handles = updateZoomOrPan(handles);
guidata(hObject, handles);


% --- Executes on slider movement.
function slider_silenceStart_Callback(hObject, eventdata, handles)
    %Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);  


% --- Executes during object creation, after setting all properties.
function slider_silenceStart_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_boutOverhang_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);  


% --- Executes during object creation, after setting all properties.
function slider_boutOverhang_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_boutOverhang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_minBoutSep_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider_minBoutSep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_minBoutSep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_minBout_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);  


% --- Executes during object creation, after setting all properties.
function slider_minBout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_minBout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_boutThresh_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_boutThresh_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_silenceStart_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);    

% --- Executes during object creation, after setting all properties.
function edit_silenceStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_silenceStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_boutOverhang_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);    

% --- Executes during object creation, after setting all properties.
function edit_boutOverhang_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_boutOverhang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minBoutSep_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function edit_minBoutSep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minBoutSep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minBout_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function edit_minBout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minBout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_boutThresh_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_boutThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_boutThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_edgeThresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_edgeThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of edit_edgeThresh as text
%        str2double(get(hObject,'String')) returns contents of edit_edgeThresh as a double


% --- Executes during object creation, after setting all properties.
function edit_edgeThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_contThresh_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something if dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end
% Hints: get(hObject,'String') returns contents of edit_contThresh as text
%        str2double(get(hObject,'String')) returns contents of edit_contThresh as a double


% --- Executes during object creation, after setting all properties.
function edit_contThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minSyl_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);   

% --- Executes during object creation, after setting all properties.
function edit_minSyl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxSyl_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);    
    
% --- Executes during object creation, after setting all properties.
function edit_maxSyl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minGap_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Only do something is dynamic mode is selected
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);   

% --- Executes during object creation, after setting all properties.
function edit_minGap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minGap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_edgeThresh_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);    

% --- Executes during object creation, after setting all properties.
function slider_edgeThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_edgeThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_contThresh_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function slider_contThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_contThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_minSyl_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider_minSyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_minSyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_maxSyl_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);  


% --- Executes during object creation, after setting all properties.
function slider_maxSyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_maxSyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_minGap_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Take additional action only if you're in dynamic mode    
if get(handles.check_dynamic,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
    
    if strcmp(handles.lastMode, 'multi')
        %Call the multisegmenting routine
        [currAnnot, handles] = multiSegment(handles, handles.audio, []);
    else
        [currAnnot, handles] = frameSegment(handles, handles.audio);
    end

    %Update the annotation folder
    handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
    
    %Draw updated syllable information
    handles = drawNavBarSyllableRects(handles);
    handles = drawZoomSyllableRects(handles);
end

guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function slider_minGap_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_multisegment.
function push_multiSegment_Callback(hObject, eventdata, handles)

%fprintf("multiseg pushed\n");
%Call the multisegmenting routine
[currAnnot, handles] = multiSegment(handles, handles.audio, []);

%fprintf("multiseg routine successfully run");

%Update the annotation folder
handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);

%fprintf("annotation folder updated"); 

%Remove the file from the errorlist, if it's in there.
if isfield(handles, 'errorFiles')
    handles.errorFiles = removeFile(handles.filelist(handles.curfile).name, handles.errorFiles);
end

%Remove the file from the bad segmentation list, if it's in there.
if isfield(handles, 'badSegFiles') && ~ismember(104,currAnnot.segType)
    handles.badSegFiles = removeFile(handles.filelist(handles.curfile).name, handles.badSegFiles);
end

%draw syllable information
handles = drawNavBarSyllableRects(handles);
handles = drawZoomSyllableRects(handles);

%flag
handles.lastMode = 'multi';

guidata(hObject, handles);


% --- Executes on button press in push_applyNetwork.
function push_applyNetwork_Callback(hObject, eventdata, handles)
%Function applies a previously trained neural network to un-annotated song file. The main prerequisites for the 
%function are (1) a trained network and (2) a segmented song file.

%Check to see if the network is trained; if not, exit out
if ~isfield(handles, 'nn_classifier')
    set(handles.text_message,'String','Your must have a neural net trained and in memory before using Automated ID functions.')
    return
end

%Check to see that song file is segmented
currAnnot = handles.annotation.get(handles.filelist(handles.curfile).name);
if isempty(currAnnot) | isempty(currAnnot.segFileStartTimes)
    set(handles.text_message,'String','Song file must be segmented prior to applying Automated ID functions. Segment and try again.')
    return
end

%Because you have the prerequisites, procede to format the data for classifying.

%Get the Syl duration padding range from user box
weightPercent = str2num(get(handles.edit_durWeight, 'String'))/100;
durPadVectLength =  floor((5000 * weightPercent)/(1 - weightPercent));
durPadVect = ones(1,durPadVectLength);

%Constants for Bandpass Audio (300-10000kHz)
HP_fNorm = 300/(44150/2);
LP_fNorm = 6500/(44150/2);
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);
PPaudio = filtfilt(BP_b,BP_a,handles.audio);
PPaudio = (PPaudio/(sqrt(mean(PPaudio.^2))));

%Parse the spectrogram using the segmentation boundaries
buffer = str2num(get(handles.edit_overhang, 'String')); %buffer in ms (10ms by default)
startBins = max(floor((currAnnot.segFileStartTimes - (buffer/1000)) * handles.fs), 1);
endBins = min(ceil((currAnnot.segFileEndTimes + (buffer/1000)) * handles.fs), length(PPaudio));
strt = 4; stp = 94; %These bins correspond to 300-8000Hz
numSyls = length(startBins);
for i = 1:numSyls
    [~,~,~,P{i}] = spectrogram(PPaudio(startBins(i):endBins(i)),220,220-44,512,44150);
    logSpecs{i,1} = -1*log10(abs(P{i}(strt:stp,:)));
    
    durPad(i,:) = (size(logSpecs{i,1},2)/10) * durPadVect;
    vectorSpec(i,:) = nnFormatInputSpec(logSpecs{i,1});
    log_spec_formatted(i,:) = [vectorSpec(i,:) , durPad(i,:)];

%     log_spec_formatted(i,:) = nnFormatInputSpec(logSpecs{i});
%     figure(1)
%     imagesc(-logSpecs{i}); axis xy
%     pause(1)
end

% Apply network to the current dataset
[nn_classification] = nnApply(handles.nn_classifier,log_spec_formatted,handles.nn_classes,str2double(get(handles.edit_decisionThresh,'String')));

%Transform network output to the annotation segTypes 
ind = isnan(nn_classification);
nn_classification(ind) = 103; %label as unknown
currAnnot.segType = nn_classification;

%Finish up
handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);
handles = updateZoomOrPan(handles);
set(handles.text_message,'String','Neural Net based syllable identification complete.')
    
guidata(hObject, handles);

% --- Executes on button press in push_trainNetwork.
function push_trainNetwork_Callback(hObject, eventdata, handles)

fprintf("button was pushed\n");
%Get the list of syllables to train on
if isfield(handles, 'templates')
    sylList = getFieldVector(handles.templates.wavs,'segType');
    sylList = [sylList, 101, 102, 103]; %Add the defaults for calls, unknown, etc.
    %sylList = [sylList, 101, 102]; %Add the defaults for calls, and noise; NOT the unknown...
else
    set(handles.text_message, 'String', 'Template must be loaded before training network. Load and try again.')
    beep
    return
end
fprintf("got syll list\n");

%Create dataset of audio files to parse
load(handles.annotFilename);
handles.rawDataset = [];
for i=1:length(keys)
    fprintf("%d\n",i);
    n = char(keys{i});
    % Strip the channels out of the raw recording file
    if strcmp(n(end-3:end), '.wav')
        %[r, fs] = wavread(keys{i});
        [r, fs] = audioread(keys{i});
        if size(r,2)>1
            handles.rawDataset{i} = r(:,1); %Take only the first channel from stereo audio recordings
        else
            handles.rawDataset{i} = r;
        end
    elseif strcmp(n(end-3:end), '.dat') || strcmp(n(end-3:end), '.stm')
        [r, fs] = getChannels(keys{i});
        handles.rawDataset{i} = r(1,:); %Take only the first channel from multiplexed A&N recordings
    end
    %[handles.rawDataset{i}, fs] = wavread(keys{i});
end
set(handles.text_message,'String','Dataset creation completed')

%Preallocate the main output structures
logSpecs = cell(length(sylList),1);
sylTypeIndx = [];

%Chose to include calls and other 'trash' in the training set.
if get(handles.check_trainSylsOnly,'Value')
    [handles.allSeq, handles.startSeq] = sortSequences(handles, 1);
else
    [handles.allSeq, handles.startSeq] = sortALLSequences(handles, 1);
end
fprintf("gettin there\n");

%Cycle through each of the syllable types, finding syllable indices constructing a lookup table for each annotated file.
buffer = str2num(get(handles.edit_overhang, 'String'))/1000; %buffer in ms (10ms by default)
for i = 1:length(sylList)
    fprintf("Syll number %d\n",i);
    %Generate indices for each syllable
    handles.chosenStartSeq = parseSequences(handles,sylList(i));
    
    if ~isempty(handles.chosenStartSeq)
        %Given the index, parse the loaded dataset
        fprintf("I dont get it\n");
        [~, ~, audio] = getSubSet(handles,buffer,[]);
    
        %Generate logSpecs for each audio snip
        PPaudio = Prep(handles,audio);
        audioSpecs = calcPower(PPaudio,'audio');
        renditions = size(audio,2);
        for j=1:renditions
            logSpecs{i,1}{j,1} = -1*log10(audioSpecs{j});%./mean(mean(log10(handles.data.audioSpecs{i})));
        end
    
        %Create the sylTypeIndx for the block
        blockTypes = sylList(i)*ones(renditions,1);
    
        %Then add it to the running summary list
        sylTypeIndx = [sylTypeIndx; blockTypes];
    else
        display(['Syllable type ' num2str(sylList(i)) ' not present in annotation'])
    end
end

%Get the Syl duration padding range from user box
weightPercent = str2num(get(handles.edit_durWeight, 'String'))/100;
durPadVectLength =  floor((5000 * weightPercent)/(1 - weightPercent));
durPadVect = ones(1,durPadVectLength);
fprintf("got syll duration padding range\n");

%Loop through the entire logSpec structure to format the spectrograms for the training routine
log_spec_formatted={};
for i = 1:length(logSpecs)
     for k = 1:length(logSpecs{i,1})
         durPad(k,:) = ((size(logSpecs{i,1}{k,1},2))/10) * durPadVect ;
         vectorSpec(k,:) = nnFormatInputSpec(logSpecs{i,1}{k,1});
         [log_spec_formatted{i,1}(k,:)] = [vectorSpec(k,:) , durPad(k,:)];
         
%          [log_spec_formatted{i,1}(k,:)] = nnFormatInputSpec(logSpecs{i,1}{k,1}); %Original code
         
     end
end
log_spec_formatted=cell2mat(log_spec_formatted);

%Format the TypeIndx
[handles.nn_classes, target_classification] = nnFormatTargetClass(sylTypeIndx);

%Once the audio is extracted, run the training routine and copy the classifier to the handles structure
[nn_scores, nn_success_rate, handles.nn_classifier] = nnTrain(log_spec_formatted, target_classification);
handles.nn_buffer = str2num(get(handles.edit_overhang, 'String'));
handles.nn_durWeighting = str2num(get(handles.edit_durWeight, 'String'));

fprintf("Did we get there?\n");
%Signal training complete
set(handles.text_message, 'String', 'Syllable classifying neural net trained and ready for use.');
set(handles.push_trainNetwork,'BackgroundColor',[0,1,0])
    
guidata(hObject, handles);
    
function chosenStartSeq = parseSequences(handles,chosenSeq)
%Returns the numerical sequence that the user has selected
%chosenSeq = handles.seq(get(handles.popup_seqSylls, 'Value')); %chosen sequence as number

%Cycle through all of the syllables listed in the template file by pumping the syllable label (as a double)  into chosenSeq
%handles.Seq = chosenSeq;

%Finds the index of the chosen sequence as they appear in the
%handles.all_seq catalog
startIndex = find(handles.allSeq==chosenSeq);

%Gives the filtered record/annotation number and the syllable number within
%the annotation
chosenStartSeq = [];
for i=1:length(startIndex)
    chosenStartSeq(i,:) = handles.startSeq(startIndex(i),:); 
end

% --- Executes on button press in push_saveNetwork.
function push_saveNetwork_Callback(hObject, eventdata, handles)
if ~isfield(handles,'nn_classifier')
    warndlg('No trained network available to save. Train one and try again.');
else
    nn_classes = handles.nn_classes;
    nn_classifier = handles.nn_classifier;
    nn_buffer = handles.nn_buffer;
    nn_durWeighting = handles.nn_durWeighting;
    [FileName,PathName]=uiputfile('*.mat','Save Network As:','network');
    save([PathName FileName],'nn_classes','nn_classifier','nn_buffer','nn_durWeighting');
    set(handles.text_message, 'String', ['Network parameters saved to ' FileName]);
end

guidata(hObject, handles);

% --- Executes on button press in push_loadNetwork.
function push_loadNetwork_Callback(hObject, eventdata, handles)
if isfield(handles,'nn_classifier')
    ansButton = questdlg('A trained classifier network is already loaded. Keep this one or load another?','Network Already Loaded','Keep this One','Load Another','Keep this One');
    if strcmp(ansButton,'Load Another')
        uiload;
        handles.nn_classes = nn_classes;
        handles.nn_classifier = nn_classifier;
        handles.nn_durWeighting = nn_durWeighting;
        handles.nn_buffer = nn_buffer;
        set(handles.edit_overhang,'String',num2str(nn_buffer));
        set(handles.edit_durWeight,'String',num2str(nn_durWeighting));
        set(handles.text_message, 'String', 'Network parameters loaded from file.');
        set(handles.push_trainNetwork,'BackgroundColor',[0,1,0])
    end
else
    uiload;
    handles.nn_classes = nn_classes;
    handles.nn_classifier = nn_classifier;
    handles.nn_durWeighting = nn_durWeighting;
    handles.nn_buffer = nn_buffer;
    set(handles.edit_overhang,'String',num2str(nn_buffer));
    set(handles.edit_durWeight,'String',num2str(nn_durWeighting));
    set(handles.text_message, 'String', 'Network parameters loaded from file.');
    set(handles.push_trainNetwork,'BackgroundColor',[0,1,0])
end
    
guidata(hObject, handles);

% --- Executes on button press in check_trainSylsOnly.
function check_trainSylsOnly_Callback(hObject, eventdata, handles)

% --- Executes on button press in check_skipAnnotated.
function check_skipAnnotated_Callback(hObject, eventdata, handles)

% --- Executes on button press in push_batchProcess.
function push_batchProcess_Callback(hObject, eventdata, handles)
%Test function for batch processing the segmentations for the day
if ~isfield(handles,'filelist')
    return
elseif isempty(handles.filelist)
    return
else
    tStart = tic;
    h = waitbar(0,'Processing all files in the folder. Please wait...');
    handles.errorFiles = {};
    numFiles = length(handles.filelist);
    
    %Set the time control bounds
    sp = regexp(get(handles.edit_batchStartTime,'String'), ':','split');
    startTime = str2double(sp(1)) + str2double(sp(2))/60;
    sp = regexp(get(handles.edit_batchEndTime,'String'), ':','split');
    endTime = str2double(sp(1)) + str2double(sp(2))/60;
    
    %Get signal processing values from the GUI
    gain = 1; %str2double(get(handles.edit_sigGain, 'String'));
    hp = 300; %str2double(get(handles.edit_sigHP, 'String'));
    lp = 8500; %str2double(get(handles.edit_sigLP, 'String'));
    
    %Constants for Bandpass Audio
    HP_fNorm = hp/(handles.fs/2);
    LP_fNorm = lp/(handles.fs/2);
    [BP_b,BP_a] = butter(4,[HP_fNorm LP_fNorm]);
    
    %Debugging switch
    endPoint = numFiles;
    %     endPoint = 5;
    for i = 1:endPoint
        handles.curfile = i;
        
        if  get(handles.check_skipAnnotated,'Value') && handles.annotation.containsKey(handles.filelist(handles.curfile).name)
            display(['File  ' handles.filelist(handles.curfile).name ' already annotated.'])
        
        elseif (handles.filelist_time(handles.curfile) < startTime) || (handles.filelist_time(handles.curfile) > endTime)
            %Do nothing, I suppose...
            display(['File  ' handles.filelist(handles.curfile).name ' out of range.'])
        else
            % Save the current annotation so no data is lost
            handles.elements = saveAnnotation(handles);

            % Strip the channels out of the raw recording file
            if strcmp(handles.filelist(handles.curfile).name(end-3:end), '.wav')
%                 [rawdata] = wavread(handles.filelist(handles.curfile).name);
                    [rawdata] = audioread(handles.filelist(handles.curfile).name);
                if size(rawdata,2)>1
                    rawdata = rawdata(:,1); %Take only the first channel from stereo audio recordings
                end
            elseif strcmp(handles.filelist(handles.curfile).name(end-3:end), '.dat') || strcmp(handles.filelist(handles.curfile).name(end-3:end), '.stm')
                [rawdata, ~] = getChannels(handles.filelist(handles.curfile).name);
                rawdata = rawdata(1,:)'; %Take only the first channel from multiplexed A&N recordings
            end
            
            if handles.fs ~= handles.exper.desiredInSampRate;
                warndlg('Experiment parameters and data file do not have the same sampling rate.');
                uiwait;
            end

            strparts = regexp(handles.filelist(handles.curfile).name,'_', 'split');
            handles.filenum = str2double(strparts{2});
            
            %             %Constants for Bandpass Audio (300-6500Hz)
            %             HP_fNorm = 300/(44150/2);
            %             LP_fNorm = 6500/(44150/2);
            %             [BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);

            % Do the filtering
            PPaudio = filtfilt(BP_b,BP_a,rawdata'.*gain);
            
            try
                thresh = [];
                %Enable this snippet for multiSegment
                 [currAnnot, handles] = multiSegment(handles, PPaudio, thresh);

                %Enable this snippet for batch Frame Segmenting
%                 handles.startNdx = 1;
%                 handles.endNdx = length(PPaudio);
%                 [currAnnot, handles] = frameSegment(handles, PPaudio);

                %If a neural network classifier exists, attempt to classify syllable types
                if isfield(handles, 'nn_classifier')
                    PPaudio = (PPaudio/(sqrt(mean(PPaudio.^2))));
                    
                    %Get the Syl duration padding range from user box
                    weightPercent = str2num(get(handles.edit_durWeight, 'String'))/100;
                    durPadVectLength =  floor((5000 * weightPercent)/(1 - weightPercent));
                    durPadVect = ones(1,durPadVectLength);
                    
                    %Parse the spectrogram using the segmentation boundaries
                    buffer = str2num(get(handles.edit_overhang, 'String')); %buffer in ms (10ms by default)
                    startBins = floor((currAnnot.segFileStartTimes - (buffer/1000)) * handles.fs);
                    startBins(startBins <= 0) = 1;
                    
                    endBins = ceil((currAnnot.segFileEndTimes + (buffer/1000)) * handles.fs);
                    endBins(endBins > length(PPaudio)) = length(PPaudio);

                    strt = 4; stp = 94; %These bins correspond to 300-8000Hz
                    numSyls = length(startBins);
                    P = cell(numSyls,1);
                    logSpecs = cell(numSyls,1);
                    log_spec_formatted = zeros(numSyls,(5000+durPadVectLength));
                    durPad = zeros(numSyls,durPadVectLength);
                    vectorSpec = zeros(numSyls,5000);
                    parfor i = 1:numSyls
                        [~,~,~,P{i}] = spectrogram(PPaudio(startBins(i):endBins(i)),220,220-44,512,44150);
                        logSpecs{i} = -1*log10(abs(P{i}(strt:stp,:)));
                        
                        durPad(i,:) = (size(logSpecs{i},2)/10) * durPadVect;
                        vectorSpec(i,:) = nnFormatInputSpec(logSpecs{i});
                        log_spec_formatted(i,:) = [vectorSpec(i,:) , durPad(i,:)];
                        
                        %log_spec_formatted(i,:) = nnFormatInputSpec(logSpecs{i});
                    end

                    % Apply network to the current dataset
                    [nn_classification] = nnApply(handles.nn_classifier,log_spec_formatted,handles.nn_classes,str2double(get(handles.edit_decisionThresh,'String')));

                    %Transform network output to the annotation segTypes 
                    ind = isnan(nn_classification);
                    nn_classification(ind) = 103; %label as unknown
                    currAnnot.segType = nn_classification;
                end

                %Update the annotation file
                handles.annotation.put(handles.filelist(handles.curfile).name, currAnnot);

            catch
                display(['Error generated for file ' handles.filelist(handles.curfile).name])
                handles.errorFiles{end+1} = handles.filelist(handles.curfile).name;
            end
        end
        
        %Update the user on progress at discrete 
        waitbar(handles.curfile/endPoint)
    end

    %Update the 
    elapsedMin = round((toc(tStart)/60)*10)/10;
    set(handles.text_message, 'String',['Completed processing ' num2str(handles.curfile) ' files in ' num2str(elapsedMin) ' minutes.']);
    close(h) 
end

guidata(hObject, handles);


function edit_overhang_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_overhang_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_batchStartTime_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_batchStartTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_batchEndTime_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_batchEndTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_decisionThresh_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_decisionThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_durWeight_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_durWeight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_numSylls.
function popup_numSylls_Callback(hObject, eventdata, handles)
handles.numSylls = get(handles.popup_numSylls, 'Value');
handles = extractSequences(handles);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popup_numSylls_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_numSylls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popup_seqSylls.
function popup_seqSylls_Callback(hObject, eventdata, handles)
%empty in TV.... 

% --- Executes during object creation, after setting all properties.
function popup_seqSylls_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_seqStats.
function push_seqStats_Callback(hObject, eventdata, handles)
%go through the annotations and whenever there is a stretch of contiguous
%sequence assign the transition to the matrix. display the figure.
numberSyll=2;
[allSeq startSeq counter]=getSequences(handles,numberSyll);
numberSyllTypes=max(max(allSeq));
seqMatrix=zeros(numberSyllTypes,numberSyllTypes);   
for i=1:counter-1
    from=allSeq(i,1);
    to=allSeq(i,2);
    seqMatrix(from,to)=seqMatrix(from,to)+1;
end
figure;
h=imagesc([1 numberSyllTypes],[1 numberSyllTypes],seqMatrix);

function handles = extractSequences(handles)
%Run through the now refined set of records to extract the sequences
if ~isfield(handles,'numSylls')
    handles.numSylls = 1;
end
[handles.allSeq handles.startSeq]=sortSequences(handles, handles.numSylls); %sort the sequences of handles.numSyll in terms of how often they appear

%Sort and display sequences by prevalance
histx=(1:10^handles.numSylls);
hist_seq=hist(handles.allSeq,histx);

%Choose the 10 most abundant sequences to display in the sequence pop-up window
[sortSeq, sortIndex]=sort(hist_seq,'descend');
numSeq=length(find(sortSeq>0));
if isempty(numSeq)
    warndlg('No such sequences found.');
    uiwait;
    return;
elseif (numSeq<15)
    popupNum=numSeq;
else
    popupNum=15;
end
popupSeqString(1,:)=[num2str(sortIndex(1)),' (n=', ' ',num2str(sortSeq(1)),')'];
histSeqString(1,:) = [num2str(sortIndex(1))];

handles.seq=sortIndex(1);

%Clear axes for new data
axes(handles.axes_annotationHist)
cla

if (popupNum>1)
    for i=2:popupNum %display the frequency of the various sequences)
        padLength=length(num2str(sortSeq(1)))-length(num2str(sortSeq(i)));
        if (padLength==1)
            freq(i,:)=['  ',num2str(sortSeq(i))];
        elseif (padLength==2)
            freq(i,:)=['   ',num2str(sortSeq(i))];
        elseif (padLength==3)
            freq(i,:)=['    ',num2str(sortSeq(i))];
        elseif (padLength==4)
            freq(i,:)=['    ',num2str(sortSeq(i))];
        else
            freq(i,:)=[' ',num2str(sortSeq(i))];
        end
        handles.seq(i)=sortIndex(i);
        popupSeqString(i,:)=[num2str(sortIndex(i)),' (n=', freq(i,:),')'];
        histSeqString = [histSeqString; num2str(sortIndex(i))];
    end
    
    %Plot the histogram to the axes
    xs = 1:popupNum;
    xL = [0, 1+popupNum];
    hold on
    bar(xs,sortSeq(xs));
    line([0,popupNum+1],[100,100],'LineStyle',':','Color',[1 0 0])
    axis tight
    set(gca,'Box','off','TickDir','out')
    set(gca,'XTick',xs,'XLim',xL,'XTickLabel','')
    ax = axis;
    axis(axis)
    yL = ax(3:4);
    
    % Place the text labels
    t = text(xs,yL(1)*ones(1,length(xs)),histSeqString);
    set(t,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);
    
    % Get the Extent of each text object.
    for i = 1:length(t)
      ext(i,:) = get(t(i),'Extent');
    end
    
    % Determine the lowest point.  The X-label will be placed so that the top is aligned with this point.
    LowYPoint = min(ext(:,2));
    
    % Place the axis label at this point
    XMidPoint = xL(1)+abs(diff(xL))/2;
    tl = text(XMidPoint,LowYPoint,'Sequence','VerticalAlignment','top','HorizontalAlignment','center');
    xlabel('');
    ylabel('');
    hold off
end

%Set the values in the 
set(handles.popup_seqSylls,'String',popupSeqString);
set(handles.popup_seqSylls,'Value',1)

function [allSeq startSeq counter]=getSequences(handles,numberSyll)
%Gets all the sequences containing numSyll.  Note that this is entirely
%dependent on the current filter settings.

%Reset placeholders
counter=1;
song_interruptions=[];

%Restrict sequence search to those files filtered by the user
%filt_elements = handles.elements(handles.filtInd == 1);
filt_elements = handles.elements;

%Find all pauses greater than 100ms to determine boundaries of song motifs
for i = 1:length(filt_elements)
    for j = 1:length(filt_elements{i}.segFileStartTimes)-1
        pauses(j) = filt_elements{i}.segFileStartTimes(j+1)-filt_elements{i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.3);
    end

    if length(filt_elements{i}.segFileStartTimes)>=numberSyll
        for j = 1:length(filt_elements{i}.segFileStartTimes)-numberSyll+1
            current_seq(1:numberSyll)=filt_elements{i}.segType(j:j+numberSyll-1);
            if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
                allSeq(counter,1:numberSyll)=current_seq;
                startSeq(counter,:)=[i,j];%index (index in annotation file,number of starting syllable of sequences)
                counter=counter+1;
            end
        end
    end
end

if ~exist('allSeq')
    allSeq=0;
    startSeq=0;
    counter=0;
end

function [allSeqNumbers startSeq]=sortALLSequences(handles,numberSyll)
%make a cell array with all the song sequences

[allSeq startSeq counter]=getALLSequences(handles,numberSyll);
if (allSeq==0)
    warndlg('There are no sequences in the specified range.');
    uiwait;
    return;
end
%go through this cell array and build up an array of sequences
%numberSyll long
for i=1:numberSyll
    multiplier(i)=10^(numberSyll-i);
end
allSeqNumbers=[];
for i=1:counter-1
    allSeqNumbers(i)=dot(allSeq(i,:),multiplier);
end

function [allSeq startSeq counter]=getALLSequences(handles,numberSyll)
%Gets all the sequences containing numSyll.  Note that this is entirely
%dependent on the current filter settings.

%Reset placeholders
counter=1;
song_interruptions=[];

%Restrict sequence search to those files filtered by the user
%filt_elements = handles.elements(handles.filtInd == 1);
filt_elements = handles.elements;

%Find all pauses greater than 100ms to determine boundaries of song motifs
for i = 1:length(filt_elements)
    for j = 1:length(filt_elements{i}.segFileStartTimes)-1
        pauses(j) = filt_elements{i}.segFileStartTimes(j+1)-filt_elements{i}.segFileEndTimes(j);
        song_interruptions=find(pauses>0.1);
    end

    if length(filt_elements{i}.segFileStartTimes)>=numberSyll
        for j = 1:length(filt_elements{i}.segFileStartTimes)-numberSyll+1
            current_seq(1:numberSyll)=filt_elements{i}.segType(j:j+numberSyll-1);
            %if (isempty(find(current_seq<1 | current_seq>9)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
            if (isempty(find(current_seq<1)) && isempty(find(song_interruptions>j-1 & song_interruptions<j+numberSyll-1)))
                allSeq(counter,1:numberSyll)=current_seq;
                startSeq(counter,:)=[i,j];%index (index in annotation file,number of starting syllable of sequences)
                counter=counter+1;
            end
        end
    end
end

if ~exist('allSeq')
    allSeq=0;
    startSeq=0;
    counter=0;
end

function [allSeqNumbers startSeq]=sortSequences(handles,numberSyll)
%make a cell array with all the song sequences

[allSeq startSeq counter]=getSequences(handles,numberSyll);
if (allSeq==0)
    warndlg('There are no sequences in the specified range.');
    uiwait;
    return;
end
%go through this cell array and build up an array of sequences
%numberSyll long
for i=1:numberSyll
    multiplier(i)=10^(numberSyll-i);
end
allSeqNumbers=[];
for i=1:counter-1
    allSeqNumbers(i)=dot(allSeq(i,:),multiplier);
end


% --- Executes on button press in check_loadDATs.
function check_loadDATs_Callback(hObject, eventdata, handles)
% hObject    handle to check_loadDATs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_loadDATs


% --- Executes on button press in check_loadSTMs.
function check_loadSTMs_Callback(hObject, eventdata, handles)
fprintf("Ya checked the box!\n");
% hObject    handle to check_loadSTMs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_loadSTMs
