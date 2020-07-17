function varargout = FinchShocker(varargin)
% FINCHSHOCKER MATLAB code for FinchShocker.fig
%      FINCHSHOCKER, by itself, creates a new FINCHSHOCKER or raises the existing
%      singleton*.
%
%      H = FINCHSHOCKER returns the handle to a new FINCHSHOCKER or the handle to
%      the existing singleton*.
%
%      FINCHSHOCKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FINCHSHOCKER.M with the given input arguments.
%
%      FINCHSHOCKER('Property','Value',...) creates a new FINCHSHOCKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FinchShocker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FinchShocker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FinchShocker

% Last Modified by GUIDE v2.5 20-Dec-2017 18:02:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FinchShocker_OpeningFcn, ...
                   'gui_OutputFcn',  @FinchShocker_OutputFcn, ...
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

% --- Executes just before FinchShocker is made visible.
function FinchShocker_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for FinchShocker
handles.output = hObject;

global hardware
hardware = [];

%Plexon
hardware.plexon.id = 1;   % Assume (hardcode) 1 Plexon box
hardware.plexon.open = false;
hardware.plexon.chans = [3, 5, 7, 9, 11, 13]; %This is for the 6-ch Nanoclip
% hardware.stim_trigger = 'ni';

%NI DAQ
hardware.ni.session = [];

%TCP/IP Comms
hardware.tcp.host = 'localhost';
hardware.tcp.port = 2057; %This code use for communication "port" and "port+1"
hardware.tcp.sr = 150; %in ms -- I think we're cool with setting this much lower than tha orginal 1000ms
hardware.tcp.numRetries = 1; %Set to -1 for infinite

global in_stim_loop
in_stim_loop = false;

global loop_stopped
loop_stopped = false;

%Set handles constants
handles.fs = 44150; %STM file sampling rate (in Hz)
handles.chans = 1:6;%This is for the 6-ch Nanoclip
handles.stimFileLoaded = false;

%Load default values into the GUI
def.I1 = 100;
def.T1 = 500;
def.TI = 5;
def.I2 = -100;
def.T2 = 500;
for i = handles.chans
    handles = updateGUIChan(handles, def, i);
end

set(handles.popup_monChan, 'String', num2cell(handles.chans), 'Value', 1)

%Initialize hardware
init_plexon(hObject, handles);
updatePlexonMain(handles);
init_ni(hObject, handles);

% Update handles structure
guidata(hObject, handles);

function varargout = FinchShocker_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = updateGUIChan(handles, vals, chan)
%This function updates the GUI controls for stim params for the channel
%specified in 'chan'. This is GUI only and has no effect on the hardware.
cats = {'I1'; 'T1'; 'TI'; 'I2'; 'T2'};
for i = 1:numel(cats)
   c =  ['handles.edit_ch' num2str(chan) '_' cats{i} ', ''String'', num2str(vals.' cats{i} ')'];
   eval(['set(' c ');'])
end

function params = scrapeChan(handles, chan)
%This function scrapes the GUI controls for stim params for the channel
%specified in 'chan' and puts them in a structure to be sent to the Plexon
cats = {'I1'; 'I2'; 'T1'; 'T2'; 'TI'}; %This order is necessary for the Plexon
params = [];
for i = 1:numel(cats)
   %Collect the values from the GUI
   c =  ['handles.edit_ch' num2str(chan) '_' cats{i} ', ''String'''];
   eval(['q = get(' c ');'])
   
   %Copy to output structure
   params(i) = str2double(q);
end

function init_plexon(hObject, handles)
%Initializes the Plexon stimulator for use with the GUI

%Retrieve hardware globals
global hardware;
global plexon_newly_initialised;

% Close/clear the Plexon
PS_CloseAllStim;
if hardware.plexon.open
    hardware.plexon.open = false;
end

%Initialize Plexon and evaluate status message
err = PS_InitAllStim;
switch err
    case 1
        msgbox({'Error: Could not open the Plexon box.', ' ', 'POSSIBLE CAUSES:', '* Device is not attached', '* Device is not turned on', '* Another program is using the device', ...
            '* Device needs rebooting', '', 'TO REBOOT:', '1. DISCONNECT THE BIRD!!!', '2. Power cycle', '3. Reconnect bird.'});
        error('plexon:init', 'Plexon initialisation error: %s', PS_GetExtendedErrorInfo(err));
    case 2
        msgbox({'Error: Could not open the Plexon box.', ' ', 'POSSIBLE CAUSES:', '* Device is not attached', '* Device is not turned on', '* Another program is using the device', ...
            '* Device needs rebooting', '', 'TO REBOOT:', '1. DISCONNECT THE BIRD!!!', '2. Power cycle', '3. Reconnect bird.'});
        error('plexon:init', 'Plexon: no devices available.  Is the blue box on?  Is other software accessing it?');
    otherwise
        hardware.plexon.open = true;
end

%Check to confirm there is only one Plexon attached
nstim = PS_GetNStim;
if nstim > 1
    err = PS_CloseAllStim;
    hardware.plexon.open = false;
    error('plexon:init', 'Plexon: %d devices available, but only planned for one!', nstim);
    return;
end

%Check the number of stim channels on the Plexon
[nchan, err] = PS_GetNChannels(hardware.plexon.id);
if err
    ME = MException('plexon:init', 'Plexon: invalid stimulator number "%d".', hardware.plexon.id);
    throw(ME);
else
    %disp(sprintf('Plexon show_device %d has %d channels.', hardware.plexon.id, nchan));
end
if nchan ~= 16
    ME = MException('plexon:init', 'Ben assumed that there would always be 16 channels, but there are in fact %d', nchan);
    throw(ME);
end

%Set the stimulation mode (0-software; 1-TTL pulse)
err = PS_SetTriggerMode(hardware.plexon.id, 1);
if err
    ME = MException('plexon:stimulate', 'Could not set trigger mode on stimbox %d', hardware.plexon.id);
    throw(ME);
end

%Update status
plexon_newly_initialised = true;
set(handles.text_message, 'String', 'Plexon stimulator connected and initialized')

guidata(hObject, handles);

function init_ni(hObject, handles)
global hardware;
%global ni_response_channels;
%global ni_trigger_index;
%global recording_time;

%Clear a previous open session
if ~isempty(hardware.ni.session)
    stop(hardware.ni.session);
    release(hardware.ni.session);
    hardware.ni.session = [];
end

% Open NI acquisition board session
dev='Dev3'; % location of input device
hardware.ni.session = daq.createSession('ni');

%Specify TTL Read/Write
addDigitalChannel(hardware.ni.session, dev, 'Port0/Line0:7', 'InputOnly');
addDigitalChannel(hardware.ni.session, dev, 'Port1/Line0:7', 'OutputOnly');

% ni_trigger_index = size(hardware.ni.session.Channels, 2);
% hardware.ni.session.Channels(ni_trigger_index).Name = channel_labels{ni_trigger_index};
% nscans = round(hardware.ni.session.Rate * hardware.ni.session.DurationInSeconds);
% tc = addTriggerConnection(hardware.ni.session, sprintf('%s/PFI0', dev), 'external', 'StartTrigger');
%ch = hardware.ni.session.addDigitalChannel(dev, 'Port0/Line1', 'OutputOnly');
%ch = hardware.ni.session.addCounterOutputChannel(dev, 'ctr0', 'PulseGeneration');
%disp(sprintf('Output trigger channel is ctr0 = %s', ch.Terminal));
%ch.Frequency = 0.1;
%ch.InitialDelay = 0;
%ch.DutyCycle = 0.01;
%pulseme = zeros(nscans, 1);
%pulseme(1:1000) = ones(1000, 1);
%hardware.ni.session.queueOutputData(pulseme);

prepare(hardware.ni.session);

function updatePlexonMain(handles)
%The script handles updating the Plexon main parameters -- these won't
%typically change from trial to trial
global hardware;
global plexon_newly_initialised;

%Check status
if ~hardware.plexon.open || ~plexon_newly_initialised
    set(handles.text_messages, 'String', 'Plexon doesn''t seem to be online... check that shit and try again.') 
    return
end

%Stimulator wide settings
%Monitor channel
pntr = get(handles.popup_monChan, 'Value');
errMonSel = PS_SetMonitorChannel(hardware.plexon.id, hardware.plexon.chans(pntr));

%Monitor scaling
errMonScale = PS_SetVmonScaling(hardware.plexon.id, 0.25);

%Error checking
if errMonSel || errMonScale
    set(handles.text_message, 'String', ['Error updating monitor state.'])
    return
end
   
for i = handles.chans
   %Set channel pattern type (0-rectangular; 1-arbitrary)
   errType = PS_SetPatternType (hardware.plexon.id, hardware.plexon.chans(i), 0);

   %Single biphasic pulses
   errReps = PS_SetRepetitions (hardware.plexon.id, hardware.plexon.chans(i), 1);
   
   if errType || errReps
       set(handles.text_message, 'String', ['Error updating Chan' num2str(i) ' state.'])
       return
   else
       set(handles.text_message, 'String', ['Completed initializing main Plexon parameters.'])
   end
end

function updatePlexonFromGUI(handles)
%The script handles updating the Plexon channel parameters from the GUI
%valuesAh ha!
global hardware;
global plexon_newly_initialised;

%Check status
if ~hardware.plexon.open || ~plexon_newly_initialised
    set(handles.text_messages, 'String', 'Plexon doesn''t seem to be online... check that shit and try again.') 
    return
end

for i = handles.chans
   %Scrape and package parameters 
   params = scrapeChan(handles, i);
   
   %Set channel pattern params
   errParam = PS_SetRectParam (hardware.plexon.id, hardware.plexon.chans(i), params);
   
   %Load channel settings
   errLoad = PS_LoadChannel(hardware.plexon.id, hardware.plexon.chans(i));
   
   if errParam || errLoad
       set(handles.text_message, 'String', ['Error updating Chan' num2str(i) ' state.'])
       return
   else
       set(handles.text_message, 'String', ['Completed updating Chan' num2str(i) ' state.'])
   end
end

function statusMsg = getPlexonStatus(handles)
%The script polls the current Plexon configuration and formats the data in
%a readable message

global hardware;
global plexon_newly_initialised;

%Check status
if ~hardware.plexon.open || ~plexon_newly_initialised
    set(handles.text_messages, 'String', 'Plexon doesn''t seem to be online... check that shit and try again.')
    return
end

%Grab stimulator settings
[Mode, errTrig] = PS_GetTriggerMode (hardware.plexon.id);
[MonChan, errMon] = PS_GetMonitorChannel(hardware.plexon.id);
[Scaling, errScal] = PS_GetVmonScaling (hardware.plexon.id);

if errTrig || errMon || errScal
    set(handles.text_message, 'String', ['Error retrieving stimulator settings.'])
    return
end

%Grab channel-specific settings
Rep = [];
Param = [];
for i = handles.chans
    [Rep(i), errRep] = PS_GetRepetitions (hardware.plexon.id, hardware.plexon.chans(i));
    [Param(i,:), errParam] = PS_GetRectParam (hardware.plexon.id, hardware.plexon.chans(i));
    
    if errRep || errParam
        set(handles.text_message, 'String', ['Error retrieving Chan' num2str(i) ' settings.'])
        return
    end
end

%Pre-format inputs
timestamp = ['Plexon @ ' datestr(datetime)];

if Mode == 1
    modeStr = 'Triggers from TTL';
else
    modeStr = 'Trigger is mis-configured';
end
md = ['Monitor: Chan ' num2str(handles.chans(hardware.plexon.chans == MonChan))];
vmSc = ['VMon Scaling: ' num2str(Scaling) ' V/V'];
P = cellstr(num2str(Param));
p1 = ['Chan 1: ' P{1}];
p2 = ['Chan 2: ' P{2}];
p3 = ['Chan 3: ' P{3}];
p4 = ['Chan 4: ' P{4}];
p5 = ['Chan 5: ' P{5}];
p6 = ['Chan 6: ' P{6}];

statusMsg=[timestamp ...
    sprintf('\n') modeStr ...
    sprintf('\n') md ...
    sprintf('\n') vmSc ...
    sprintf('\n') ...
    sprintf('\n') 'Pattern Settings' ...
    sprintf('\n') p1 ...
    sprintf('\n') p2 ...
    sprintf('\n') p3 ...
    sprintf('\n') p4 ...
    sprintf('\n') p5 ...
    sprintf('\n') p6];

set(handles.text_plexonStats, 'String', statusMsg, 'FontSize', 10);

function [MonChan, Scaling, Rep, Param] = getStimParams(handles)
%The script polls the current Plexon configuration and formats the data in
%a readable message

global hardware;
global plexon_newly_initialised;

%Check status
if ~hardware.plexon.open || ~plexon_newly_initialised
    set(handles.text_messages, 'String', 'Plexon doesn''t seem to be online... check that shit and try again.')
    return
end

%Grab stimulator settings
[MonChan, errMon] = PS_GetMonitorChannel(hardware.plexon.id);
[Scaling, errScal] = PS_GetVmonScaling (hardware.plexon.id);

if errMon || errScal
    set(handles.text_message, 'String', ['Error retrieving stimulator settings.'])
    return
end

%Grab channel-specific settings
Rep = [];
Param = [];
for i = handles.chans
    [Rep(i), errRep] = PS_GetRepetitions (hardware.plexon.id, hardware.plexon.chans(i));
    [Param(i,:), errParam] = PS_GetRectParam (hardware.plexon.id, hardware.plexon.chans(i));
    
    if errRep || errParam
        set(handles.text_message, 'String', ['Error retrieving Chan' num2str(i) ' settings.'])
        return
    end
end

function elements = saveStimFile(handles)
% Save the currently StimFile to disk
if(isfield(handles,'stimFilename'))
    aaSaveHashtable(handles.stimFilename, handles.stimFile);

    %Copy out the two main structures of the annotation file
    annotation = aaLoadHashtable(handles.stimFilename);
    
    if isempty(annotation)
        elements = [];
    else
        elements = annotation.elements;
    end
else
    elements = [];
    warndlg('You have not yet created a stimFile.');
    uiwait;
end

function handles = updateDestinations(handles)
%If "Guess Directory" mode is on, this function updates (and creates) new 

%Derive date folder from the current computer clock
formatOut = 'yyyy-mm-dd';
dFolder = datestr(now,formatOut);

%Assemble the monitor path
monPath = [handles.motherDir, filesep, handles.birdname, filesep, dFolder];

%If the folder doesn't yet exist, create it -- but notify the user
if ~exist(monPath, 'dir')
    mkdir(monPath);
    set(handles.text_message, 'String', ['Created a new folder: ' monPath])
end

%Update monitor folders
handles.dirname = monPath;
set(handles.text_monDir, 'String', handles.dirname)
                  
%Assemble the new stimFile path and name
dateSnip = [dFolder(3:4), dFolder(6:7), dFolder(9:10)];
stimName = [handles.birdname, '_' dateSnip '_stmFile.mat'];
stimFile = [handles.motherDir, filesep, handles.birdname, filesep, stimName];

%if one is loaded, take care of the old stimFile
if handles.stimFileLoaded
    %Save one last time
    saveStimFile(handles);
    
    %Delete the currently loaded hashtable
    handles.stimFile.delete;
end

if ~exist(stimFile, 'file')
    %Create new stimFile hashtable
    handles.stimFile = mhashtable;
    
else
    %Load the one on the disk
    handles.stimFile = aaLoadHashtable(stimFile);
    
end

%Update stimFile locations
handles.stimFilename = stimFile;
handles.stimFileLoaded = true;
set(handles.text_stimFile, 'String', handles.stimFilename);

%Save
saveStimFile(handles);

function handles = monitorLoop(handles)
%This is the main looping function for this program. Aims are to monitor
%the progress of the LabVIEW tCNS and A&N Recorder programs and log
%progress
global hardware
global loop_stopped
global in_stim_loop

%Create time object to track changes in day (for housekeeping)
loopDay = datetime('today');
updateTrigger = 60*60*2; %Approx 1 hour @ ~500us cycles
loopCounter = 0;

%Array for computing running average of loop time
mLoopTime = zeros(1,10);

filelist = [];
a = true;
while ~loop_stopped 
   %Estimate loop time
   tic;
   
   %Scan for digital inputs at the start of the cycle
   dataIn = inputSingleScan(hardware.ni.session);
   dataOut = false(1,8); %by default, set all outputs low
   
   %If folder/stimFile guessing is enabled, periodically check what day it
   %is and update the fields as required
   if get(handles.check_guessDir, 'Value')
       %Increment the loop counter
       loopCounter = loopCounter + 1;
       
       %Count until the trigger
       if loopCounter > updateTrigger
           %If the day has actually changed since the last check, update
           %the monitor directory and create new stimFile
           if days(datetime('today') - loopDay) ~= 0
               handles = updateDestinations(handles);
               loopDay = datetime('today');
           end
           
           %Reset the loop counter
           loopCounter = 0;
       end
   end
   
   %Look for a signal to indicate a change in system status
   if get(handles.check_runTTL, 'Value')
       fileReady = dataIn(2); %|| get(handles.push_simFileReady, 'Value');
       if fileReady
           %filelist = dir([handles.dirname filesep '*.stm']); %In the future, replace with TCP read
           %curFile = filelist(end).name;
           curFile = exchangeData(hardware.tcp.port, hardware.tcp.host, hardware.tcp.numRetries, hardware.tcp.sr, []);
       end
   else
       [filelist, fileReady] = scanFolder(handles.dirname, filelist);
       curFile = filelist(end).name;
   end
   
   if (fileReady || get(handles.check_forceRun, 'Value')) && ~isempty(curFile)
       %Parse the filename of the current file
       sp = regexp(curFile, '_', 'split');
       
       %Extract the useful stuff from the file
       [stimTime, targetTime, Vmon, Imon] = getSTMInfo([handles.dirname filesep curFile]);
       
       %Gather (from Plexon) specified stim parameters for this file
       [MonChan, Scaling, Rep, Param] = getStimParams(handles);
       if ~isempty(stimTime)
           for i = 1:numel(stimTime)
               ChanPattern(i,:,:) = Param;
               ChanRep(i,:) = Rep;
           end
       else
           ChanPattern = [];
           ChanRep = [];
       end
       
       %Generate an record for the experiment hash table
       curStimInfo.birdname = sp{1};
       curStimInfo.filenum = str2double(sp{2});
       curStimInfo.timestamp = datetime;
       curStimInfo.stimTimes = stimTime;
       curStimInfo.targetTimes = targetTime; %Should adding a record be contingent on having a target in the file?
       curStimInfo.monChan = MonChan;
       curStimInfo.Vmon = Vmon;
       curStimInfo.VmonScale = Scaling;
       curStimInfo.Imon = Imon;
       curStimInfo.chanPattern = ChanPattern;
       curStimInfo.chanRep = ChanRep;
       
       %Add the record to the hashtable and save to file
       handles.stimFile.put(curFile, curStimInfo);
       handles.elements = saveStimFile(handles);
       
       %Update the stimulation parameters (according to some rule)
       

       
       %Anything else?
       if ~isempty(ChanPattern)
           
       end
       
       %Notify of new record
       set(handles.text_lastFile, 'String', curFile);
       
   end
   %Pause the loop to limit cycle spinup -- this will need to be adjusted
   %as the "change signal" is improved
   pause(.5)
   
   %Send Ack to LV; reset Ack if input is low
   if dataIn(2)
       dataOut(1) = true;
   else
       dataOut(1) = false;
   end
   
   %Update DIO output
   outputSingleScan(hardware.ni.session,dataOut)
   
   %Flash the indicators
   %LabView Heartbeat
   if dataIn(1)
       set(handles.text_LVbeat, 'String', 'LV Heartbeat', 'BackgroundColor', [0,1,0])
   else
       set(handles.text_LVbeat, 'String', '', 'BackgroundColor', [.94,.94,.94])
   end
   
   %Loop Running
   mLoopTime = [round(toc*1000) mLoopTime(1:end-1)];
   if a
       set(handles.text_loopInd, 'String', ['LT: ' num2str(mean(mLoopTime)) 'ms'], 'BackgroundColor', [0,1,0])
       a = false;
   else
       set(handles.text_loopInd, 'String', ['LT: ' num2str(mean(mLoopTime)) 'ms'], 'BackgroundColor', [.94,.94,.94])
       a = true;
   end

end

%Clear global flag
in_stim_loop = false;

function [stimTimeIdx, targetTimeIdx, Vmon, Imon] = getSTMInfo(stmFile)
%This function loads the STM file from disk, detects the recorded trigger
%pulses, and returns the relative timing (i.e., wrt the start of the file) 
%of those pulses. May grab some other shit too

%Predefine
targetTimeIdx = [];
Vmon = [];
Imon = [];
snipLength = 50; %in ms

%Load file from disk and parse channels
[chans, fs] = getChannels(stmFile);

%Extract the times for the stim TTL
stimTimeIdx = getPulseTimes(chans(4,:));

%If there is a stimulation impulse present, continue to extract the other info...
if ~isempty(stimTimeIdx)
    %Targeting pulses
    targetTimeIdx = getPulseTimes(chans(5,:));
    
    for i = 1:numel(stimTimeIdx)
        snipEnd = (stimTimeIdx(i) + floor((snipLength/1000)*fs));
        snip = stimTimeIdx(i):min([snipEnd,size(chans,2)]);
        if snipEnd > snip(end) %must pad to fit
            numPad = snipEnd - snip(end);
            VmonTemp = chans(2,snip);
            ImonTemp = chans(3,snip);
            
            Vmon(i,:) = [VmonTemp, NaN(1,numPad)];
            Imon(i,:) = [ImonTemp, NaN(1,numPad)];
        else
            
            Vmon(i,:) = chans(2,snip);
            Imon(i,:) = chans(3,snip);
        end
    end
end

%Here would probably be a good place to put in any additional analysis of
%song/behavior/etc
%

%Convert sample time to filetime
stimTimeIdx = stimTimeIdx/fs;
targetTimeIdx = targetTimeIdx/fs;

function ttlTimes = getPulseTimes(ts)
%Set threshold for TTL detection
thresh = 2.5; %5VDC signal, so something lower should work fine, I think...

%Threshold crossings (only look at rising edge)
p = 2:numel(ts);
ttlTimes = find(ts(p)>thresh & ts(p-1)<thresh);
    
function [newFiles, doIt] = scanFolder(dirLoc, filelist)
%This script compares the current folder with past contents; if they are
%different, it sets a boolean flag.
doIt = false;

%Read current folder contents
newFiles = dir([dirLoc filesep '*.stm']);

%Compare length of contents (faster than comparing names... may need to change if problems arise)
if ~isempty(filelist)
    doIt = numel(filelist) ~= numel(newFiles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_loopInd_Callback(hObject, eventdata, handles)
%This is the "in-loop" indicator light -- no coding needed.

function push_startLoop_Callback(hObject, eventdata, handles)
%Button begins the monitoring/formatting loop
global hardware
global plexon_newly_initialised
global loop_stopped
global in_stim_loop

%Check the preconditions for starting loop
if ~hardware.plexon.open || ~plexon_newly_initialised
    set(handles.text_messages, 'String', 'Plexon doesn''t seem to be online... check that shit and try again.') 
    return
end

if in_stim_loop
    set(handles.text_messages, 'String', 'System is already in the monitoring loop.') 
    return
end

%Change button state
set(handles.push_startLoop, 'String', 'Looping!', 'BackgroundColor', 'g')
set(handles.push_stopLoop, 'String', 'Stop Looping!', 'BackgroundColor', 'r')

%Main Loop call
in_stim_loop = true;
loop_stopped = false;

%Update timestamp
set(handles.text_startTime, 'String', datestr(datetime));

%Start the loop
handles = monitorLoop(handles);

%Do whatever you gotta do on exit
in_stim_loop = false;

% Update handles structure
guidata(hObject, handles);

function push_stopLoop_Callback(hObject, eventdata, handles)
%Button stops the monitoring/formatting loop
global loop_stopped
global in_stim_loop

%Kill the loop
loop_stopped = true;
in_stim_loop = false;

%Change button state
set(handles.push_startLoop, 'String', 'Start Loop', 'BackgroundColor', [0.9, 0.9, 1])
set(handles.push_stopLoop, 'String', 'Stopped', 'BackgroundColor', [0.9, 0.9, 1])

% Update handles structure
guidata(hObject, handles);

function push_monPick_Callback(hObject, eventdata, handles)
%This script updates the monitoring folder


%Shift what the button does, based on the checkbox
if get(handles.check_guessDir, 'Value')
    %Get the birdname and confirm mother directory
    prompt = {'Birdname:', 'Mother Location:'};
    dlg_title = 'Gimme Somethin''';
    defAns = {'', 'D:\Dropbox'};
    answer = inputdlg(prompt,dlg_title, 1, defAns);

    %Copy out to handles structure
    handles.birdname = answer{1};
    handles.motherDir = answer{2};
    
    %Update destinations
    handles = updateDestinations(handles);
    
%     %Derive date folder from the current computer clock
%     formatOut = 'yyyy-mm-dd';
%     dFolder = datestr(now,formatOut);
%     
%     %Assemble the monitor path
%     monPath = [answer{2}, filesep, answer{1}, filesep, dFolder];
%     
%     %Confirm its a real location on the disk
%     if exist(monPath, 'dir')
%         handles.dirname = monPath;
%     else
%         %Throw an notice if it's not a real place
%         set(handles.text_message, 'String', ['Could not find derived path ' monPath '. Check and try again.']);
%         return
%     end
    
else
    %Ask user which directory to monitor
    directory_name = uigetdir('D:\Dropbox\');
    %directory_name = 'D:\Dropbox\LR73RY133\2017-12-15';
    
    %Confirm its non-empty
    if ~isstr(directory_name)
        return
    else
        handles.dirname = directory_name;
        %cd(handles.dirname)
    end
    
end

%Update GUI display
set(handles.text_monDir, 'String', handles.dirname)

% Update handles structure
guidata(hObject, handles);

function pushbutton1_Callback(hObject, eventdata, handles)

init_plexon(hObject, handles);
disp('...done');

% Update handles structure
guidata(hObject, handles);

function pushbutton2_Callback(hObject, eventdata, handles)

updatePlexonMain(handles)
updatePlexonFromGUI(handles)

% Update handles structure
guidata(hObject, handles);

function pushbutton3_Callback(hObject, eventdata, handles)

msg = getPlexonStatus(handles);

% Update handles structure
guidata(hObject, handles);

function push_simFileReady_Callback(hObject, eventdata, handles)


% Update handles structure
guidata(hObject, handles);

function push_simAck_Callback(hObject, eventdata, handles)

function push_stimFile_Callback(hObject, eventdata, handles)
%If StimFile is already open/active check whether to append or clear.
if isfield(handles, 'stimFilename')
    ansButton = questdlg('Would you like to save the current StimFile and begin a new one, or append to the existing StimFile?','StimFile Already Loaded','Save & Clear','Append','Append');
    
    switch ansButton
        case 'Save & Clear'
            %If there is already an annotation object loaded, then save current annotation.
            saveStimFile(handles);
            
            %Delete the currently loaded hashtable
            handles.stimFile.delete;
        case 'Append'
            %Do nothing and carry on.
            set(handles.text_message, 'String', 'Ok, nothing changed. Keep on trucking.')
            return
    end
end

%Get the new location for the StiFile
[file,path] = uiputfile('*.mat', 'Create or select a .mat file for the StimFile:');
if isempty(file) || isempty(path) %if it's a not real location, exit
    set(handles.text_message, 'String', 'You gotta select a real StimFile location before this is gonna work.')
    return
end

%Does the file already exist?
if ~exist([file,path], 'file')
    %Update new stimFile path
    handles.stimFilename = [path,file];
    
    %Create new annotation hashtable
    handles.stimFile = mhashtable;
    handles.stimFileLoaded = true;
    
else
    ansButton = questdlg('The selected StimFile already exists on disk. What should we do with it?','StimFile Already Exists','Load and Append','Delete and Replace','Load and Append');
    switch ansButton
        case 'Load and Append'
            %Update new stimFile path
            handles.stimFilename = [path,file];
            
            %Load StimFile into memory
            handles.stimFile = aaLoadHashtable(handles.stimFilename);
            handles.stimFileLoaded = true;
            
        case 'Delete and Replace'
            %Update new stimFile path
            handles.stimFilename = [path,file];
            
            %Create new annotation hashtable
            handles.stimFile = mhashtable;
            handles.stimFileLoaded = true;
    end
    
end

%Update stimFile path display
set(handles.text_stimFile, 'String', handles.stimFilename);

%Save StimFile
saveStimFile(handles);

guidata(hObject, handles);

function check_forceRun_Callback(hObject, eventdata, handles)

function check_runTTL_Callback(hObject, eventdata, handles)

function check_guessDir_Callback(hObject, eventdata, handles)
%This check box determines how the minotring directory is determined
if get(handles.check_guessDir, 'Value')
    %We be guessing
    set(handles.push_monPick, 'String', 'Set Bird')
else
    %Direct specify
    set(handles.push_monPick, 'String', 'Set Directory')
end

% Update handles structure
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unused Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_ch1_I1_Callback(hObject, eventdata, handles)

function edit_ch1_I1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch1_T1_Callback(hObject, eventdata, handles)

function edit_ch1_T1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch1_TI_Callback(hObject, eventdata, handles)

function edit_ch1_TI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch1_I2_Callback(hObject, eventdata, handles)

function edit_ch1_I2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch1_T2_Callback(hObject, eventdata, handles)

function edit_ch1_T2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch2_I1_Callback(hObject, eventdata, handles)

function edit_ch2_I1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch2_T1_Callback(hObject, eventdata, handles)

function edit_ch2_T1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch2_TI_Callback(hObject, eventdata, handles)

function edit_ch2_TI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch2_I2_Callback(hObject, eventdata, handles)

function edit_ch2_I2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch2_T2_Callback(hObject, eventdata, handles)

function edit_ch2_T2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch4_I1_Callback(hObject, eventdata, handles)

function edit_ch4_I1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch4_T1_Callback(hObject, eventdata, handles)

function edit_ch4_T1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch4_TI_Callback(hObject, eventdata, handles)

function edit_ch4_TI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch4_I2_Callback(hObject, eventdata, handles)

function edit_ch4_I2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch4_T2_Callback(hObject, eventdata, handles)

function edit_ch4_T2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch3_I1_Callback(hObject, eventdata, handles)

function edit_ch3_I1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch3_T1_Callback(hObject, eventdata, handles)

function edit_ch3_T1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch3_TI_Callback(hObject, eventdata, handles)

function edit_ch3_TI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch3_I2_Callback(hObject, eventdata, handles)

function edit_ch3_I2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch3_T2_Callback(hObject, eventdata, handles)

function edit_ch3_T2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch5_I1_Callback(hObject, eventdata, handles)

function edit_ch5_I1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch5_T1_Callback(hObject, eventdata, handles)

function edit_ch5_T1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch5_TI_Callback(hObject, eventdata, handles)

function edit_ch5_TI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch5_I2_Callback(hObject, eventdata, handles)

function edit_ch5_I2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch5_T2_Callback(hObject, eventdata, handles)

function edit_ch5_T2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch6_I1_Callback(hObject, eventdata, handles)

function edit_ch6_I1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch6_T1_Callback(hObject, eventdata, handles)

function edit_ch6_T1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch6_TI_Callback(hObject, eventdata, handles)

function edit_ch6_TI_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch6_I2_Callback(hObject, eventdata, handles)

function edit_ch6_I2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ch6_T2_Callback(hObject, eventdata, handles)

function edit_ch6_T2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_monChan_Callback(hObject, eventdata, handles)

function popup_monChan_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function push_LVbeat_Callback(hObject, eventdata, handles)
