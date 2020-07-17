function varargout = TemplateBender(varargin)
% TEMPLATEBENDER MATLAB code for TemplateBender.fig
%      Created by TMO 3/22/2013; Last update 3/31/13 1:15am
%
%      TEMPLATEBENDER, by itself, creates a new TEMPLATEBENDER or raises the existing
%      singleton*.
%
%      H = TEMPLATEBENDER returns the handle to a new TEMPLATEBENDER or the handle to
%      the existing singleton*.
%
%      TEMPLATEBENDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEMPLATEBENDER.M with the given input arguments.
%
%      TEMPLATEBENDER('Property','Value',...) creates a new TEMPLATEBENDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TemplateBender_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TemplateBender_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TemplateBender

% Last Modified by GUIDE v2.5 27-Mar-2014 16:22:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TemplateBender_OpeningFcn, ...
                   'gui_OutputFcn',  @TemplateBender_OutputFcn, ...
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

function TemplateBender_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TemplateBender (see VARARGIN)

handles.addFlag = false;
handles.deleteFlag = false;

%Setup the slider controls
set(handles.slider_edgeThresh, 'Max', 5, 'Min', -5, 'SliderStep', [0.025 0.025]./10)

% Choose default command line output for TemplateBender
handles.output = [];
handles.template = [];
handles.templatesyllBreaks = [];
handles.threshold = [];
handles.renditionTraces = [];
        
%If data is passed in, strip it out in the handles structure and start processing it
if ~isempty(varargin)
    %Copy the passed data from the input cell structure to the local handles structure
    handles.template = varargin{1}.template;
    %handles.templatesyllBreaks = varargin{1}.templatesyllBreaks;
    if isfield(varargin{1},'renditionTraces') %Not always passed?
        handles.renditionTraces = varargin{1}.renditionTraces;
    end
    
    set(handles.edit_edgeThresh,'String',num2str(varargin{1}.threshold))
    handles = edit2slider(handles.edit_edgeThresh, handles);
    
    %Copy untouched data to the perm structure so that original values can always be recalled
    handles.templatePerm = handles.template;
    
    handles.thresholdPerm = varargin{1}.threshold;
    
%     %Plot the basic data to the GUI access
%     handles = drawAudioAxes(handles);
    
    %Locate edges in the power envelop 
    [handles.powerM, handles.powerStd, handles.specPowTS] = getPowerDistr(handles.template);
    [handles.templatesyllBreaks, handles.powThresh] = getEdges(handles.template, handles.powerM, handles.powerStd, get(handles.slider_edgeThresh,'Value'), handles.specPowTS);
    handles.templatesyllBreaksPerm = handles.templatesyllBreaks;
    
        %Plot the basic data to the GUI access
    handles = drawAudioAxes(handles);
    
%     set(handles.text_message,'String',num2str(handles.templatesyllBreaks)); %for debug...
    
    %Draw lines on axes for syllable boundaries
    handles = drawSylLines(handles);
    
    %Draw edit boxes at the boundaries
    handles = drawSylBoxes(handles);

end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SylClus wait for user response (see UIRESUME)
uiwait(handles.figure1);

function varargout = TemplateBender_OutputFcn(hObject, eventdata, handles) 

%Assign all outputs to the 'output' variable as a structure.
output.template = handles.template;
output.templatesyllBreaks = handles.templatesyllBreaks;
output.threshold = get(handles.slider_edgeThresh, 'Value');

%This is the main output to the calling function
varargout{1} = output;

%Terminates the figure after pushing data to the calling function
delete(handles.figure1);

function figure1_CloseRequestFcn(hObject, eventdata, handles)
% Hint: delete(hObject) closes the figure

%Prevents the output function from running until the user is finished with the sub-GUI
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent GUI functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = drawAudioAxes(handles)

%Plot the spectrogram
axes(handles.axes_spec)
cla;
imagesc(-1*handles.template)
axis xy; axis tight;
set(gca,'Box','off','TickDir', 'out','XTick',[], 'YTick',[])
set(handles.axes_spec, 'ButtonDownFcn', @cb_axes_click);
hold on

%Plot the power envelope
axes(handles.axes_powEnv)
cla;
plot(handles.specPowTS,'k','LineWidth',2)
axis tight
set(gca,'Box','off','TickDir', 'out','YTick',[0, 1])
set(handles.axes_powEnv, 'ButtonDownFcn', @cb_axes_click);
xlabel('Template Time (ms)')
ylabel('Sound Amplitude (AU)')
hold on

function handles = drawSylLines(handles)
%Draw lines on axes for syllable boundaries
axes(handles.axes_spec)
y1 = ylim;
numSyls = size(handles.templatesyllBreaks,1);
for i = 1:numSyls
    line([handles.templatesyllBreaks(i,1), handles.templatesyllBreaks(i,1)], [y1(1), y1(2)],'Color', 'k', 'LineWidth',1.5);
    line([handles.templatesyllBreaks(i,2), handles.templatesyllBreaks(i,2)], [y1(1), y1(2)],'Color', 'k', 'LineWidth',1.5);
end

axes(handles.axes_powEnv)
x1 = xlim;
y1 = ylim;
line(x1, [handles.powThresh, handles.powThresh],'Color', 'r', 'LineStyle', ':', 'LineWidth',0.5);
for i = 1:numSyls
    line([handles.templatesyllBreaks(i,1), handles.templatesyllBreaks(i,1)], [y1(1), y1(2)],'Color', 'g', 'LineWidth',1.5);
    line([handles.templatesyllBreaks(i,2), handles.templatesyllBreaks(i,2)], [y1(1), y1(2)],'Color', 'b', 'LineWidth',1.5);
end

function handles = drawSylBoxes(handles)
%Get the parameters of the current axes
axes(handles.axes_powEnv)
x1 = xlim;
location = get(handles.axes_spec, 'Position');

%Set the box position constants
boxWidth = 30;
boxHeight = 20;
scale = location(3)/(x1(2)-x1(1));
start = location(1)-(boxWidth/2);

%Draw the boxes on the figure
numSyls = size(handles.templatesyllBreaks,1);
for i = 1:numSyls
    %On Box
    eval(['handles.edit_OnSyl' num2str(i) '  = uicontrol(handles.figure1,''Tag'', ''edit_OnSyl' num2str(i) ''', ''Style'', ''edit'', ''String'', ' num2str(handles.templatesyllBreaks(i,1)) ', ''Position'', [' num2str(start+scale*handles.templatesyllBreaks(i,1)) ', ' num2str(location(2)-40) ', ' num2str(boxWidth) ', ' num2str(boxHeight) '],  ''BackgroundColor'', [0, 1, 0], ''Callback'', @syllBoxCallback);'])

    %Off Box
    eval(['handles.edit_OffSyl' num2str(i) '  = uicontrol(handles.figure1,''Tag'', ''edit_OffSyl' num2str(i) ''', ''Style'', ''edit'', ''String'', ' num2str(handles.templatesyllBreaks(i,2)) ', ''Position'', [' num2str(start+scale*handles.templatesyllBreaks(i,2)) ', ' num2str(location(2)-70) ', ' num2str(boxWidth) ', ' num2str(boxHeight) '],  ''BackgroundColor'', [0, 0, 1], ''Callback'', @syllBoxCallback);'])
end

function handles = deleteSylBoxes(handles)
%Draw the boxes on the figure
numSyls = size(handles.templatesyllBreaks,1);
for i = 1:numSyls
    %On Box
    eval(['delete(handles.edit_OnSyl' num2str(i) ');'])

    %Off Box
    eval(['delete(handles.edit_OffSyl' num2str(i) ');'])
    
end

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

function syllBoxCallback(hObject, eventdata)
%This function controls the behavior of all the syllBoxes
handles = guidata(hObject);

%Get the object name and value from the edit box
editName = get(hObject, 'Tag');
value = str2num(get(hObject, 'String'));

%Refresh the display
handles = drawAudioAxes(handles);
handles = deleteSylBoxes(handles);

%Based on the name, determine the Syllable Number and On/Off
sylNum = str2num(editName(end));
if strcmp(editName(7),'n') %It's an 'On' box
    handles.templatesyllBreaks(sylNum,1) = value;
else %It's an 'Off' box
    handles.templatesyllBreaks(sylNum,2) = value;
end

%Update the axes
handles = drawSylLines(handles);

%Update the syllable boxes
handles = drawSylBoxes(handles);

% guidata(gco, handles);

%Because we deleted the original calling hObject, we have to store the handles structure in the GUI using
%the newly created hObject that just so happens to have exactly the same name. Easiest to do it with string evaluation:
eval(['guidata(handles.' editName ', handles)'])

function cb_axes_click(hObject, evnt)
handles = guidata(hObject);
clickLocation = get(hObject, 'CurrentPoint');

axes(hObject)
if(handles.addFlag)
    %Reset the button and flag
    handles.addFlag = false;
    set(handles.push_add, 'BackgroundColor', [0.941, 0.941, 0.941]);
    
    rect = rbbox;
    endPoint = get(gca,'CurrentPoint'); 
    point1 = clickLocation(1,1);       % x coordinate (time)
    point2 = endPoint(1,1);
    onset = round(min(point1, point2));             % calculate locations
    offset = round(max(point1, point2));
    
    %Delete the current boxes
    handles = deleteSylBoxes(handles);
    
    %Contained within existing syllable
    ind_within = find(handles.templatesyllBreaks(:,1) < onset & handles.templatesyllBreaks(:,2) > offset);
    
    %Overlaps one or more syllable onsets
    ind_overlap = find((handles.templatesyllBreaks(:,1) > onset & handles.templatesyllBreaks(:,1) < offset) | (handles.templatesyllBreaks(:,2) > onset & handles.templatesyllBreaks(:,2) < offset));
    
    %Update the breaks structure
    if ~isempty(ind_within) %simply replace
        handles.templatesyllBreaks(ind_within, :) = [onset, offset];
        
    elseif ~isempty(ind_overlap) %simply replace
        if length(ind_overlap) == 1
            handles.templatesyllBreaks(ind_overlap, :) = [onset, offset];
        else
            handles.templatesyllBreaks(ind_overlap(1), :) = [onset, offset];
            handles.templatesyllBreaks(ind_overlap(2:end), :) = [];
        end
    else %Insert the new syllable in the breaks structure 
        ind_on = find(handles.templatesyllBreaks(:,1) < onset, 1, 'last');
        if isempty(ind_on)
            %Add to the beginning
            handles.templatesyllBreaks = [onset, offset; handles.templatesyllBreaks];
        elseif ind_on == size(handles.templatesyllBreaks,1)
            %Add to the end
            handles.templatesyllBreaks = [handles.templatesyllBreaks; onset, offset];
        else
            %Insert in the middle
            handles.templatesyllBreaks = [handles.templatesyllBreaks(1:ind_on,:); onset, offset; handles.templatesyllBreaks((ind_on+1):end,:)];
        end
    end
    
elseif(handles.deleteFlag)
    %Reset the button and flag
    handles.deleteFlag = false;
    set(handles.push_delete, 'BackgroundColor', [0.941, 0.941, 0.941]);

    point = clickLocation(1,1);       % x coordinate (time)
    
    %Find the syllable boundary that's closest to the click point
    [~, ind_onoff] = min(min(abs(handles.templatesyllBreaks-point)));
    [~, ind_syl] = (min(abs(handles.templatesyllBreaks-point)));
    ind_syl = ind_syl(ind_onoff);
    
    %Delete the current boxes
    handles = deleteSylBoxes(handles);
    
    if ind_onoff == 1 %it's an onset
        % so only delete it and the offset that follows
        handles.templatesyllBreaks(ind_syl,:) = [];
        
    elseif ind_onoff == 2 && ind_syl < size(handles.templatesyllBreaks,1)
        %so only delete it and the onset that follows
        handles.templatesyllBreaks(ind_syl,2) = handles.templatesyllBreaks(ind_syl+1,2);
        handles.templatesyllBreaks(ind_syl+1,:) = [];
    end
end

%Plot the basic data to the GUI access
handles = drawAudioAxes(handles);

%Draw lines on axes for syllable boundaries
handles = drawSylLines(handles);

%Draw edit boxes at the boundaries
handles = drawSylBoxes(handles);

guidata(hObject, handles);

function [powerM, powerStd, specPowTS] = getPowerDistr(data)
%Sum across frequency bins to create power envelop
specPowTS2 = -1*sum(data);
specPowTS = mat2gray(specPowTS2); %Added 9/1/2015

%Rerun EM until a solution is found
powerM = [];
powerStd = [];
iterations = 1;
while (isempty(powerM) || isempty(powerStd))
    gMix = gmdistribution.fit(specPowTS',2);
    [powerM,pnt] = sort(gMix.mu);
    temp = sqrt(squeeze(gMix.Sigma));
    powerStd = temp(pnt);
    powerM = powerM(1);
    powerStd = powerStd(1);
    
    if iterations > 5
        print(['gmdistribution.fit has been run ' num2str(iterations) ' without converging on a solution.'])
        return
    end
    iterations = iterations + 1;
end

function [syllBreaks, thresh] = getEdges(data, powerM, powerStd, edgeThresh, specPowTS)
Crossings.Up = [];
Crossings.Down = [];

%Sum across frequency bins to create power envelop
% specPowTS = mat2gray(-1*sum(data));

%Calculate Threshold
thresh = powerM+(edgeThresh*powerStd); % num of SDs above the silent Gaussian mean

%Find the crossing points to estimate the onset and offsets of
%motifs
ind = 2:length(specPowTS);
Crossings.Up = find(specPowTS(ind)>thresh & specPowTS(ind-1)<thresh);
Crossings.Down = find(specPowTS(ind)<thresh & specPowTS(ind-1)>thresh);

if numel(Crossings.Up) ~= numel(Crossings.Down)
    if Crossings.Up(end) > Crossings.Down(end)
        Crossings.Up(end) = [];
    end
    if Crossings.Up(1) > Crossings.Down(1)
        Crossings.Down(1) = [];
    end
end

%Probably need to do something fancier here, but ok for now.
syllBreaks = [];
if ~isempty(Crossings.Up) && ~isempty(Crossings.Down)
    syllBreaks(:,1) = Crossings.Up';
    syllBreaks(:,2) = Crossings.Down';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main GUI controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_load_Callback(hObject, eventdata, handles)
%Retrieve the template file
[fname, pathname] = uigetfile('*.mat','Select *.mat file containing template information');
load([pathname fname]);

%Copy the loaded variables to the handles structure
handles.template = template;
handles.templatesyllBreaks = templatesyllBreaks;
handles.templatemotifBreaks = templatemotifBreaks;
if exist('threshold', 'var')
    handles.threshold = threshold;
else
    handles.threshold = 2;
end

%Copy untouched data to the perm structure so that original values can always be recalled
handles.templatePerm = handles.template;
handles.templatesyllBreaksPerm = handles.templatesyllBreaks;
handles.thresholdPerm = handles.threshold;

%Setup up axes for new information


%Plot the basic data to the GUI access
handles = drawAudioAxes(handles);

%Update the axes
handles = drawSylLines(handles);

%Update the syllable boxes
handles = drawSylBoxes(handles);

%Set messages for state change
set(handles.text_message,'String',['Template loaded from ' pathname fname '.']);

guidata(hObject, handles);

function push_save_Callback(hObject, eventdata, handles)
if ~isfield(handles,'template')

else
    %Prep variables to save
    template = handles.template;
    templatesyllBreaks = handles.templatesyllBreaks;
    templatemotifBreaks = [handles.templatemotifBreaks(1,1), handles.templatemotifBreaks(end,end)];
    threshold = get(handles.slider_threshold,'Value');
    
    %Get location and save file
    [fname, pathname] = uiputfile('*.mat','Select location and file name to save template.');
    save([pathname fname],'template','templatesyllBreaks','templatemotifBreaks','threshold')
    
    %Update display
    set(handles.text_message,'String',['Current template saved to file: ' pathname fname]);
end

guidata(hObject, handles);

function push_reset_Callback(hObject, eventdata, handles)
%Update plots accordingly
%Refresh the display
handles = drawAudioAxes(handles);
handles = deleteSylBoxes(handles);

%Back copy data fromt he Perm Structure
handles.templatesyllBreaks = handles.templatesyllBreaksPerm;
set(handles.edit_edgeThresh,'String',num2str(handles.thresholdPerm))
handles = edit2slider(handles.edit_edgeThresh, handles);

%Calculate new syllable breaks
[handles.templatesyllBreaks, handles.powThresh] = getEdges(handles.template, handles.powerM, handles.powerStd, get(handles.slider_edgeThresh,'Value'),handles.specPowTS);

% set(handles.text_message,'String',num2str(handles.templatesyllBreaks)); %for debug...

%Update the axes
handles = drawSylLines(handles);

%Update the syllable boxes
handles = drawSylBoxes(handles);

guidata(hObject, handles);

function push_close_Callback(hObject, eventdata, handles)

%Prevents the output function from running until the user is finished with the sub-GUI
if isequal(get(handles.figure1,'waitstatus'),'waiting')
    uiresume(handles.figure1);
else
    delete(handles.figure1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edge definition controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_delete_Callback(hObject, eventdata, handles)
%Given the state of the flag, change button color and change out flag
if handles.deleteFlag
    set(hObject,'BackgroundColor',[0.941, 0.941, 0.941]); %Set it to gray
    handles.deleteFlag = false;
elseif ~handles.deleteFlag
    set(hObject,'BackgroundColor',[1, 0, 0]); %Set it to red
    handles.deleteFlag = true;
    
    %Make sure the opposing button is off
    set(handles.push_add,'BackgroundColor',[0.941, 0.941, 0.941]); %Set it to gray
    handles.addFlag = false;
end

guidata(hObject, handles);    

function push_add_Callback(hObject, eventdata, handles)
%Given the state of the flag, change button color and change out flag
if handles.addFlag
    set(hObject,'BackgroundColor',[0.941, 0.941, 0.941]); %Set it to gray
    handles.addFlag = false;
elseif ~handles.addFlag
    set(hObject,'BackgroundColor',[0, 1, 0]); %Set it to red
    handles.addFlag = true;
    
    %Make sure the opposing button is off
    set(handles.push_delete,'BackgroundColor',[0.941, 0.941, 0.941]); %Set it to gray
    handles.deleteFlag = false;
end

guidata(hObject, handles);     

function slider_edgeThresh_Callback(hObject, eventdata, handles)
%Pump the slider value to the edit box
handles = slider2edit(hObject, handles);
    
%Refresh the display
handles = drawAudioAxes(handles);
handles = deleteSylBoxes(handles);

%Calculate new syllable breaks
[handles.templatesyllBreaks, handles.powThresh] = getEdges(handles.template, handles.powerM, handles.powerStd, get(handles.slider_edgeThresh,'Value'),handles.specPowTS);

% set(handles.text_message,'String',num2str(handles.templatesyllBreaks)); %for debug...

%Update the axes
handles = drawSylLines(handles);

%Update the syllable boxes
handles = drawSylBoxes(handles);

guidata(hObject, handles);    

function slider_edgeThresh_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_edgeThresh_Callback(hObject, eventdata, handles)
%Pump the new value to the slider
handles = edit2slider(hObject, handles);
    
%Refresh the display
handles = drawAudioAxes(handles);
handles = deleteSylBoxes(handles);

%Calculate new syllable breaks
[handles.templatesyllBreaks, handles.powThresh] = getEdges(handles.template, handles.powerM, handles.powerStd, get(handles.slider_edgeThresh,'Value'), handles.specPowTS);

% set(handles.text_message,'String',num2str(handles.templatesyllBreaks)); %for debug...

%Update the axes
handles = drawSylLines(handles);

%Update the syllable boxes
handles = drawSylBoxes(handles);

guidata(hObject, handles);    

function edit_edgeThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
