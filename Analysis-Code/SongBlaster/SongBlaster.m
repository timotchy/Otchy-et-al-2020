function varargout = SongBlaster(varargin)
% SONGBLASTER M-file for SongBlaster.fig
%      SONGBLASTER, by itself, creates a new SONGBLASTER or raises the existing
%      singleton*.
%
%      H = SONGBLASTER returns the handle to a new SONGBLASTER or the handle to
%      the existing singleton*.
%
%      SONGBLASTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SONGBLASTER.M with the given input arguments.
%
%      SONGBLASTER('Property','Value',...) creates a new SONGBLASTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SongBlaster_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SongBlaster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SongBlaster

% Last Modified by GUIDE v2.5 07-Dec-2017 16:37:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SongBlaster_OpeningFcn, ...
                   'gui_OutputFcn',  @SongBlaster_OutputFcn, ...
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

% --- Executes just before SongBlaster is made visible.
function SongBlaster_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for SongBlaster
handles.output = hObject;

load('Default PreProcess.mat');
handles.PreProcess = PP;

set(handles.axes_curDisp,'XTick',[],'YTick',[]);
set(handles.axes_spec,'XTick',[],'YTick',[]);
set(handles.axes_mainDisp,'XTick',[],'YTick',[]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SongBlaster wait for user response (see UIRESUME)
% uiwait(handles.figure1);
setappdata(0  , 'hSongBlaster'    , gcf);
setappdata(gcf,   'PreProcess'    , handles.PreProcess);
setappdata(gcf, 'fhUpdateDataAxes', @updateDataAxes);

function varargout = SongBlaster_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Independent Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = loadPage(handles)
%Loads the 10 files and display on the MainDisplay

%Determine which page to load
pageNum = str2double(get(handles.edit_curPageNum,'String'));

%Create the index of which files are on this page
startIdx = (pageNum*10)+1;
endIdx = min(startIdx+9,length(handles.filelist));
handles.curPageRecIdx = startIdx:endIdx;

%Load from file each of the recordings on this page and save to pageData
%structure
fs = 44150;
HP_fNorm = 1000/(fs/2);
LP_fNorm = 6500/(fs/2);
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);

handles.pageData = [];
for i = 1:length(handles.curPageRecIdx)
    if strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'wav')
        %[handles.pageData{i}, ~] = wavread(handles.filelist{handles.curPageRecIdx(i)});
        %         [handles.pageData{i}, ~] = audioread(handles.filelist{handles.curPageRecIdx(i)});
        [t, ~] = audioread(handles.filelist{handles.curPageRecIdx(i)});
        handles.pageData{i} = filtfilt(BP_b,BP_a, t.*5);
    end
    
    if strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'dat') || strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'stm')
        [chanData, ~] = getChannels(handles.filelist{handles.curPageRecIdx(i)});   
        
        % PreProcess data, if selected
        if strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'dat') && get(handles.check_PreProcess,'Value')
            PPdata = [];
            PPdata(1,:) = chanData(1,:);
            for j = 2:size(chanData,1)
                PPdata(j,:) = PreProcess(chanData(j,:), handles.PreProcess(j-1));
            end
            handles.pageData{i} = PPdata;
        else
            handles.pageData{i} = chanData;
        end
        
    end
end

%Display the page data on the Main Axes
handles = updateMainDisp(handles);

%Update the deletion checkboxes
handles = updateChecks(handles);


function handles = loadTrashPage(handles)
%Loads the 10 files and display on the MainDisplay

%Determine which page to load
pageNum = str2double(get(handles.edit_curPageNum,'String'));

%Create the index of which files are on this page
startIdx = (pageNum*10)+1;
%endIdx = min(startIdx+9,length(handles.filelist));
endIdx = min(startIdx+9,length(handles.trashList));
handles.curPageRecIdx = startIdx:endIdx;

%Parse the data to include only those with in the delete index
%handles.trashList = handles.filelist(handles.deleteIndx);

%Load from file each of the recordings on this page and save to pageData
%structure
handles.pageData = [];
fs = 44150;
HP_fNorm = 300/(fs/2);
LP_fNorm = 6500/(fs/2);
[BP_b,BP_a] = butter(2,[HP_fNorm LP_fNorm]);
for i = 1:length(handles.curPageRecIdx)
    if strcmp(handles.trashList{handles.curPageRecIdx(i)}(end-2:end),'wav')
%         [handles.pageData{i}, ~] = wavread(handles.trashList{handles.curPageRecIdx(i)});
%             [handles.pageData{i}, ~] = audioread(handles.trashList{handles.curPageRecIdx(i)});
            [t, ~] = audioread(handles.trashList{handles.curPageRecIdx(i)});
            handles.pageData{i} = filtfilt(BP_b,BP_a, t);
    end
    
    if strcmp(handles.trashList{handles.curPageRecIdx(i)}(end-2:end),'dat') || strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'stm')
        [chanData, ~] = getChannels(handles.trashList{handles.curPageRecIdx(i)});

        % PreProcess data, if selected
        if strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'dat') && get(handles.check_PreProcess,'Value')
            PPdata = [];
            PPdata(1,:) = chanData(1,:);
            for j = 2:size(chanData,1)
                PPdata(j,:) = PreProcess(chanData(j,:), handles.PreProcess(j-1));
            end
            handles.pageData{i} = PPdata;
        else
            handles.pageData{i} = chanData;
        end

    end
end

%Display the page data on the Main Axes
handles = updateMainDisp(handles);

%Update the deletion checkboxes
handles = updateTrashChecks(handles);

function handles = updateMainDisp(handles)
hSongBlaster = getappdata(0, 'hSongBlaster');
handles.PreProcess = getappdata(hSongBlaster, 'PreProcess');

%Get scaling value from slifer bar
scaleVal = 10^get(handles.slider_scale,'Value');

%Select working axes
axes(handles.axes_mainDisp);
cla
hold on
ylim([0,11])
%Step through each record in pageData and plot it to the axes
for i = 1:length(handles.curPageRecIdx)
    
    if strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'wav')
        len = length(handles.pageData{i});
        if get(handles.popup_chan,'Value') == 1
            %If the sudio channel is selected, plot audio...
            plot((scaleVal*handles.pageData{i})+i);
        else
            %...else plot zeros
            plot(zeros(1,len)+i)
        end 
    end
    
    if strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'dat')  || strcmp(handles.filelist{handles.curPageRecIdx(i)}(end-2:end),'stm')
        %Plot the selected channel
        plot((scaleVal*handles.pageData{i}(get(handles.popup_chan,'Value'),:)+i));
    end

end

function handles = updateSpotlight(handles)
hSongBlaster = getappdata(0, 'hSongBlaster');
handles.PreProcess = getappdata(hSongBlaster, 'PreProcess');

%Get record number (relative to page start!)
recNum = str2double(get(handles.edit_curRecNum,'String'));

%Get scaling value from slifer bar
scaleVal = 10^get(handles.slider_spotScale,'Value');

%If the current recording is a Wav file
if strcmp(handles.filelist{handles.curPageRecIdx(recNum)}(end-2:end),'wav')
    %Plot Spectrogram
    axes(handles.axes_spec);
    cla
    displaySpecgramQuick(handles.pageData{recNum}, 44150,[300,8000],[],0);
    axis xy; axis tight
    set(gca,'XTick',[],'YTick',[])
    xlabel([]); ylabel([])
    
    %Plot sound pressure
    axes(handles.axes_curDisp)
    cla
    plot((1:length(handles.pageData{recNum}))/44150,scaleVal*handles.pageData{recNum});
    axis tight
    ylim([-0.5, 0.5])
% end
% 
% if strcmp(handles.filelist{handles.curPageRecIdx(recNum)}(end-2:end),'dat')
else
    %Plot Spectrogram
    axes(handles.axes_spec);
    cla
    displaySpecgramQuick(handles.pageData{recNum}(1,:), 44150,[300,8000],[],0);
    axis xy; axis tight
    set(gca,'XTick',[],'YTick',[])
    xlabel([]); ylabel([])
    
    %Plot all channels
    axes(handles.axes_curDisp)
    cla; hold on
    col = {'b','r','g','k','m','c','y'};
    %plot((1:length(handles.pageData{recNum}(1,:)))/44150,handles.pageData{recNum}(1,:)+(10),col{1}); %audio
    plot((1:length(handles.pageData{recNum}(1,:)))/44150,handles.pageData{recNum}(1,:)+(14),col{1});
    chanNum = size(handles.pageData{recNum},1);
    if get(handles.check_vref, 'Value')
        refChan = get(handles.popup_vref,'Value')+1;
        vref = handles.pageData{recNum}(refChan,:);
        for i = 2:chanNum
            tmp = handles.pageData{recNum}(i,:) - vref;
            plot((1:length(handles.pageData{recNum}(i,:)))/44150,scaleVal*tmp+(16-2*i),col{i});
        end
        
    elseif get(handles.check_CMS, 'Value')
        for i = 2:chanNum
            ref = commonMode(handles.pageData{recNum}(2:end,:),i);
            tmp = handles.pageData{recNum}(i,:) - ref;
            plot((1:length(handles.pageData{recNum}(i,:)))/44150,scaleVal*tmp+(16-2*i),col{i});
        end
    else
        for i = 2:chanNum
            plot((1:length(handles.pageData{recNum}(i,:)))/44150,scaleVal*handles.pageData{recNum}(i,:)+(16-2*i),col{i});
        end
    end
    
    axis tight
    ylim([0.5, 17.5])
end

if ~get(handles.check_showTrash,'Value' )
    set(handles.text_curFilename,'String',handles.filelist{handles.curPageRecIdx(recNum)})
else
    fIndx = handles.t2fMap(handles.curPageRecIdx,2);
    set(handles.text_curFilename,'String',handles.filelist{fIndx(recNum)})
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
CM = mean(data(commonChans,:),1);

function updateTimeDistr(handles,hourVect,minVect)

binWidth = 0.5; %width of bins in hours
ToD = hourVect+minVect/60;
hs = histc(ToD,0:binWidth:24);

axes(handles.axes_timeDist)
h = bar(.5:binWidth:23.5,hs(2:end-1));
xlim([8 22]);
set(h,'facecolor',[.5 .5 .5])
box off;
xlabel('Time of day');
ylabel('');
title(['Distribution of ' num2str(length(hourVect)) ' files'])

function handles = updateChecks(handles)
%This function updates the checkboxes based on the stored deletions
%structure
%Get the boolean pattern for this page's checkboxes
on_offs = handles.deleteIndx(handles.curPageRecIdx);

%Set each in a loop
for i = 1:length(handles.curPageRecIdx)
    set(eval(['handles.check_' num2str(i)]),'Value',on_offs(i));
end

function handles = updateTrashChecks(handles)
%This function updates the checkboxes based on the stored deletions
%structure
%Get the boolean pattern for this page's checkboxes
% on_offs = true(1,10);

%Lookup in the t2fMap which files correspond to these trash items
fIndx = handles.t2fMap(handles.curPageRecIdx,2);
on_offs = handles.deleteIndx(fIndx);

%Set each in a loop
for i = 1:length(handles.curPageRecIdx)
    set(eval(['handles.check_' num2str(i)]),'Value',on_offs(i));
end

function handles = writeChecks(handles)
%This function writes the currently selected checkboxes to file; it should
%be executed before changing the page    

%Set each in a loop
for i = 1:length(handles.curPageRecIdx)
    handles.deleteIndx(handles.curPageRecIdx(i)) = get(eval(['handles.check_' num2str(i)]),'Value');
end  

function handles = writeTrashChecks(handles)
%This function writes the currently selected checkboxes to file; it should
%be executed before changing the page    
fIndx = handles.t2fMap(handles.curPageRecIdx,2);

%Set each in a loop
for i = 1:length(handles.curPageRecIdx)
    handles.deleteIndx(fIndx(i)) = get(eval(['handles.check_' num2str(i)]),'Value');
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI Controls...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function push_delete_Callback(hObject, eventdata, handles)
%Turn delete files to the recycle bin for possible recovery (if needed)
if sum(handles.deleteIndx)==0
    set(handles.text_message,'String','There are no quarentined files to delete')
    return
else
    ansButton = questdlg('This will delete all quarentined files.  Are you sure?','Confirm Delete','Yep!','Shit! No!','Shit! No!');
end

switch ansButton
case {'Shit! No!'}
	return
case 'Yep!'
    quarFiles = handles.filelist(handles.deleteIndx);
    for i=1:length(quarFiles)
        delete(quarFiles{i});
    end
end

%Load and sort from folder
[handles.filelist, handles.listSize] = folderLoad(handles);
handles.deleteIndx = false(1,length(handles.filelist));
if isempty(handles.filelist)
    set(handles.text_message, 'String','Folder is empty');
    return
end

%Extract time values and update histogram
[~,handles.hourStamp] = getKeysSubset(handles.filelist,'h');
[~,handles.minStamp] = getKeysSubset(handles.filelist,'min');
updateTimeDistr(handles,handles.hourStamp,handles.minStamp);

% Populate display boxes
set(handles.text_curFolder,'String',handles.dirname)
set(handles.edit_curRecNum,'String','1');
set(handles.slider_curRec,'Max',10,'Min',1);
set(handles.slider_curRec,'Value',1);
set(handles.slider_curRec,'SliderStep',[1/9, 1/9]);

set(handles.edit_curPageNum,'String','0');
set(handles.slider_curPage,'Value',0);
pages = ceil(length(handles.filelist)/10)-1;
set(handles.slider_curPage,'Max',pages);
if pages ~= 0
    set(handles.slider_curPage,'SliderStep',[1/pages, 1/pages]);
    set(handles.text_page,'String',['Page of ' num2str(pages)]);
else
    set(handles.slider_curPage,'SliderStep',[1, 1]);
    set(handles.text_page,'String',['Page of 0']);
end

%Load/display first page of recordings
handles = loadPage(handles);

set(handles.text_message,'String','Checked files deleted and fileset updated.')

guidata(hObject, handles);

function [filelist, listSize] = folderLoad(handles)
%Load and sort names from an input folder
cd(handles.dirname)
list1 = [];
list2 = [];
list3 = [];
list1Names = [];
list2Names = [];
list3Names = [];
list1Size = [];
list2Size = [];
list3Size = [];
if get(handles.check_loadWAV,'Value')
    list1 = dir('*.wav');
    for i = 1:length(list1)
        list1Names{i} = list1(i).name;
        list1Size(i) = list1(i).bytes;
    end
end

if get(handles.check_loadDAT,'Value')
    list2 = dir('*.dat');
    for i = 1:length(list2)
        list2Names{i} = list2(i).name;
        
        %Switch properties based on number of channels
        if i == 1
            %check whats in the first dat file
            [data, ~] = getChannels(list2(i).name);
            numChans = size(data,1);
            
            %Build popup list options based on that info
            if numChans == 5
                set(handles.popup_chan, 'Value', 1, 'String',{'Audio';'1';'2';'3';'4'})
            elseif numChans == 7
                set(handles.popup_chan, 'Value', 1, 'String',{'Audio';'1';'2';'3';'4';'5';'6'})
            end
        end
        
        %Calculate the length of the recording
        list2Size(i) = list2(i).bytes/numChans; %Based on 1 audio channel and X neural channels
    end
end

if get(handles.check_loadSTM,'Value')
    list3 = dir('*.stm');
    for i = 1:length(list3)
        list3Names{i} = list3(i).name;
        
        %Switch properties based on number of channels
        if i == 1
            %check whats in the first dat file
            [data, ~] = getChannels(list3(i).name);
            numChans = size(data,1);
            
            %Build popup list options based on that info
            if numChans == 5
                set(handles.popup_chan, 'Value', 1, 'String',{'Audio';'1';'2';'3';'4'})
            elseif numChans == 7
                set(handles.popup_chan, 'Value', 1, 'String',{'Audio';'1';'2';'3';'4';'5';'6'})
            end
        end
        
        %Calculate the length of the recording
        list3Size(i) = list3(i).bytes/numChans; %Based on 1 audio channel and X neural channels
    end
end

%Assemble components
[filelist, indx] = sort([list1Names, list2Names, list3Names]);
listSize = [list1Size, list2Size, list3Size];
listSize = listSize(indx);

function push_load_Callback(hObject, eventdata, handles)
%Ask user to select directory to load
directory_name = uigetdir;

%If its a real directory, load files from the folder
if ~isstr(directory_name)
    return
else
    %Load and sort from folder
    handles.dirname = directory_name;
    [handles.filelist, handles.listSize] = folderLoad(handles);
 
    handles.deleteIndx = false(1,length(handles.filelist));
    if isempty(handles.filelist)
        set(handles.text_message, 'String',[directory_name ' is empty']);
        return
    end
    
    %Extract time values and update histogram
    [~,handles.hourStamp] = getKeysSubset(handles.filelist,'h');
    [~,handles.minStamp] = getKeysSubset(handles.filelist,'min');
    updateTimeDistr(handles,handles.hourStamp,handles.minStamp);
    
    % Populate display boxes
    set(handles.text_curFolder,'String',handles.dirname)
    set(handles.edit_curRecNum,'String','1');
    set(handles.slider_curRec,'Max',10,'Min',1);
    set(handles.slider_curRec,'Value',1);
    set(handles.slider_curRec,'SliderStep',[1/9, 1/9]);
    %set(handles.edit_edit_minSongs,'String','')
    
    set(handles.edit_curPageNum,'String','0');
    set(handles.slider_curPage,'Value',0);
    pages = ceil(length(handles.filelist)/10)-1;
    set(handles.slider_curPage,'Max',pages);
    if pages ~= 0
        set(handles.slider_curPage,'SliderStep',[1/pages, 1/pages]);
        set(handles.text_page,'String',['Page of ' num2str(pages)]);
    else
        set(handles.slider_curPage,'SliderStep',[1, 1]);
        set(handles.text_page,'String',['Page of 0']);
    end
    
    %Load/display first page of recordings
    handles = loadPage(handles);
    
    %Update text message display
    set(handles.text_message,'String',['Files loaded from ' directory_name])
end
guidata(hObject, handles);

function check_loadWAV_Callback(hObject, eventdata, handles)

function check_loadDAT_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Display controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% Get workspace name
hSongBlaster = getappdata(0, 'hSongBlaster');
setappdata(hSongBlaster, 'PreProcess', handles.PreProcess);

TweetVisProcOptions(handles.PreProcess,1,'hSongBlaster');

guidata(hObject, handles);

function check_PreProcess_Callback(hObject, eventdata, handles)

%Update Page
if ~get(handles.check_showTrash,'Value')
    handles = loadPage(handles);
else
    handles = loadTrashPage(handles);
end

%Update spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);

function slider_scale_Callback(hObject, eventdata, handles)
handles = updateMainDisp(handles);

%Reset Record control
set(handles.edit_curRecNum,'String','1');
set(handles.slider_curRec,'Max',10,'Min',1);
set(handles.slider_curRec,'Value',1);
set(handles.slider_curRec,'SliderStep',[1/9, 1/9]);

guidata(hObject, handles);

function slider_scale_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function popup_chan_Callback(hObject, eventdata, handles)
%Update displaye according to the selected channel
handles = updateMainDisp(handles);

guidata(hObject, handles);  

function popup_chan_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_CMS_Callback(hObject, eventdata, handles)
%At most, one box can be checked
if get(handles.check_vref, 'Value')
    set(handles.check_vref, 'Value', 0)
end

handles = updateSpotlight(handles);

guidata(hObject, handles);

function check_vref_Callback(hObject, eventdata, handles)
%At most, one box can be checked
if get(handles.check_CMS, 'Value')
    set(handles.check_CMS, 'Value', 0)
end

handles = updateSpotlight(handles);

guidata(hObject, handles);

function popup_vref_Callback(hObject, eventdata, handles)

handles = updateSpotlight(handles);

guidata(hObject, handles);

function popup_vref_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Navigation controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_curPageNum_Callback(hObject, eventdata, handles)
%Update the deletionIndx
if ~get(handles.check_showTrash,'Value')
    handles = writeChecks(handles);
else
    handles = writeTrashChecks(handles);
end

%Advance the display to the new page selection
newPage = str2double(get(handles.edit_curPageNum,'String'));
set(handles.slider_curPage,'Value',newPage)

%Update page display
if ~get(handles.check_showTrash,'Value')
    handles = loadPage(handles);
else
    handles = loadTrashPage(handles);
end

%Reset Record control
set(handles.edit_curRecNum,'String','1');
set(handles.slider_curRec,'Max',10,'Min',1);
set(handles.slider_curRec,'Value',1);
set(handles.slider_curRec,'SliderStep',[1/9, 1/9]);

%Update spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);    

function edit_curPageNum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_curRecNum_Callback(hObject, eventdata, handles)
%Advance the display to the new page selection
newRec = str2double(get(handles.edit_curRecNum,'String'));
set(handles.slider_curRec,'Value',newRec)

%Update spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);  

function edit_curRecNum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function slider_curRec_Callback(hObject, eventdata, handles)
%Advance the display to the new record selection
newRec = get(handles.slider_curRec,'Value');
set(handles.edit_curRecNum,'String',num2str(newRec))

%Update spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);

function slider_curRec_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider_curPage_Callback(hObject, eventdata, handles)
%Update the deletionIndx
if ~get(handles.check_showTrash,'Value')
    handles = writeChecks(handles);
else
    handles = writeTrashChecks(handles);
end
    
%Advance the display to the new page selection
newPage = round(get(handles.slider_curPage,'Value'));
set(handles.edit_curPageNum,'String',num2str(newPage))

%Update Page
if ~get(handles.check_showTrash,'Value')
    handles = loadPage(handles);
else
    handles = loadTrashPage(handles);
end

%Reset Record control
set(handles.edit_curRecNum,'String','1');
set(handles.slider_curRec,'Max',10,'Min',1);
set(handles.slider_curRec,'Value',1);
set(handles.slider_curRec,'SliderStep',[1/9, 1/9]);

%Update spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);

function slider_curPage_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider_spotScale_Callback(hObject, eventdata, handles)
%Update spotlight
handles = updateSpotlight(handles);

guidata(hObject, handles);

function slider_spotScale_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stay/Go Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_selectAll_Callback(hObject, eventdata, handles)
%This function trues the checkboxes

%Set each in a loop
for i = 1:10
    set(eval(['handles.check_' num2str(i)]),'Value',true);
end

function push_unselectAll_Callback(hObject, eventdata, handles)
%This function clears the checkboxes

%Set each in a loop
for i = 1:10
    set(eval(['handles.check_' num2str(i)]),'Value',false);
end

function check_10_Callback(hObject, eventdata, handles)

function check_9_Callback(hObject, eventdata, handles)

function check_8_Callback(hObject, eventdata, handles)

function check_7_Callback(hObject, eventdata, handles)

function check_6_Callback(hObject, eventdata, handles)

function check_5_Callback(hObject, eventdata, handles)

function check_4_Callback(hObject, eventdata, handles)

function check_3_Callback(hObject, eventdata, handles)

function check_2_Callback(hObject, eventdata, handles)

function check_1_Callback(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AutoTrash Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function push_autoTrash_Callback(hObject, eventdata, handles)
%Simple script to automatically scan folder for trash

%Calculate the sizes of the minSong files (in KB for now)
mintext = get(handles.edit_minSongs,'String');
eval(['minIndx = [' mintext '];'])
minLength = min(handles.listSize(minIndx'));

%Calculate a threshold with a little safety factor
safetyBuffer = 10000; %10kB
trashThreshold = minLength-safetyBuffer;

%Generate trash index using simple threshold
handles.deleteIndx = (handles.listSize < trashThreshold);

%Reload page to update view
handles = loadPage(handles);

guidata(hObject, handles);

function edit_minSongs_Callback(hObject, eventdata, handles)

function edit_minSongs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function check_showTrash_Callback(hObject, eventdata, handles)
%Update display mode basedon user selection
if ~get(handles.check_showTrash,'Value')
    %All files view
    pages = ceil(length(handles.filelist)/10)-1;
    set(handles.edit_curPageNum,'String','0')
    handles = loadPage(handles);
else
    %Trash Bin view
    handles.trashList = handles.filelist(handles.deleteIndx);
    
    %Generate a mapping from files in the trashlist to entries in the fillelist
    tIndx = 1:length(handles.trashList);
    fIndx = find(handles.deleteIndx);
    handles.t2fMap = [tIndx', fIndx'];
    
    pages = ceil(length(handles.trashList)/10)-1;
    set(handles.edit_curPageNum,'String','0')
    handles = loadTrashPage(handles);
end

set(handles.slider_curPage,'Max',pages);
set(handles.slider_curPage,'Value',1);
if pages ~= 0
    set(handles.slider_curPage,'SliderStep',[1/pages, 1/pages]);
    set(handles.text_page,'String',['Page of ' num2str(pages)]);
else
    set(handles.slider_curPage,'SliderStep',[1, 1]);
    set(handles.text_page,'String',['Page of 0']);
end

guidata(hObject, handles);

function push_resetDelete_Callback(hObject, eventdata, handles)
%Resets the entire deletion index to zeros (whether in auto or manual trash modes)

%Generate all false trash index
handles.deleteIndx = false(1,length(handles.filelist));

%Display all files
handles = loadPage(handles);
set(handles.check_showTrash,'Value',0)

guidata(hObject, handles);

function check_loadSTM_Callback(hObject, eventdata, handles)
