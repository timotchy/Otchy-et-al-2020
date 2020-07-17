function varargout = TweetVisProcOptions(varargin)
% TWEETVISPROCOPTIONS M-file for TweetVisProcOptions.fig
%      TWEETVISPROCOPTIONS, by itself, creates a new TWEETVISPROCOPTIONS or raises the existing
%      singleton*.
%
%      H = TWEETVISPROCOPTIONS returns the handle to a new TWEETVISPROCOPTIONS or the handle to
%      the existing singleton*.
%
%      TWEETVISPROCOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWEETVISPROCOPTIONS.M with the given input arguments.
%
%      TWEETVISPROCOPTIONS('Property','Value',...) creates a new TWEETVISPROCOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TweetVisProcOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TweetVisProcOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TweetVisProcOptions

% Last Modified by GUIDE v2.5 20-Mar-2012 16:13:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TweetVisProcOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @TweetVisProcOptions_OutputFcn, ...
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


% --- Executes just before TweetVisProcOptions is made visible.
function TweetVisProcOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TweetVisProcOptions (see VARARGIN)

% Choose default command line output for TweetVisProcOptions
handles.output = hObject;

% Grab the config data that is passed in and store in the handles structure
if length(varargin) == 3
    handles.PreProcess = varargin{1};
    handles.curChan = varargin{2};
    handles.parentName = varargin{3};
elseif length(varargin) == 2
    handles.PreProcess = varargin{1};
    handles.curChan = varargin{2};
    handles.parentName = gcf;
elseif length(varargin) == 1
    handles.PreProcess = varargin{1};
    handles.curChan = 1;
    handles.parentName = gcf;
else
    handles.PreProcess = createPreProcess; %use to generate default
    handles.curChan = 1;
    handles.parentName = gcf;
end

% Call the routine to populate all of the GUI boxes before window opens
handles = PopulateGUI(handles);

%hTweetVision = getappdata(0, 'hTweetVision');
%handles.PreProcess = getappdata(hTweetVision, 'PreProcess');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TweetVisProcOptions wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TweetVisProcOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function handles = PopulateGUI(handles)
% Takes the data the program is passed and uses it to populate the default
% values of the GUI.  Should run only on launch, when the channel is switched, 
% or when values are loaded from file.

% Set index
i = handles.curChan;

% Determine the number of channels in the config file and set the menu
chanNum = length(handles.PreProcess);
set(handles.popup_channel,'String',int2str((1:chanNum)'));
set(handles.popup_channel,'Value',i);

% Pop Processing Steps popups
popStrings = {'None'; 'Local Detrend'; 'Denoise'; 'High Pass Filter'; 'Low Pass Filter'; 'Amplify'; 'Smooth'; 'Down Sample'};
for j=1:7
    set(eval(['handles.popup_Step' num2str(j)]),'String',popStrings');
    set(eval(['handles.popup_Step' num2str(j)]),'Value',handles.PreProcess(i).Steps(j));
end

% Pop Gen Params boxes
set(handles.edit_genFs,'String',num2str(handles.PreProcess(i).params.Fs));
set(handles.edit_genTW,'String',num2str(handles.PreProcess(i).params.tapers(1)));
set(handles.edit_genK,'String',num2str(handles.PreProcess(i).params.tapers(2)));
set(handles.edit_genZeroPad,'String',num2str(handles.PreProcess(i).params.pad));
set(handles.edit_genLoLim,'String',num2str(handles.PreProcess(i).params.fpass(1)));
set(handles.edit_genHiLim,'String',num2str(handles.PreProcess(i).params.fpass(2)));
set(handles.edit_genErrorType,'String',num2str(handles.PreProcess(i).params.err(1)));
set(handles.edit_genErrorP,'String',num2str(handles.PreProcess(i).params.err(2)));
set(handles.edit_genTrialAve,'String',num2str(handles.PreProcess(i).params.trialave));

% Pop Detrending boxes
set(handles.edit_detrendWindow,'String',num2str(handles.PreProcess(i).Detrend.window));
set(handles.edit_detrendWinstep,'String',num2str(handles.PreProcess(i).Detrend.winstep));

% Pop Denoise boxes
set(handles.edit_denoiseLP,'String',num2str(handles.PreProcess(i).Denoise.LP));
set(handles.edit_denoiseHP,'String',num2str(handles.PreProcess(i).Denoise.HP));
set(handles.edit_denoiseP,'String',num2str(handles.PreProcess(i).Denoise.p));
set(handles.edit_denoiseF0,'String',num2str(handles.PreProcess(i).Denoise.f0));

% Pop HP filter boxes
popStrings = {'butter'; 'cheby1'; 'cheby2'; 'ellip'};
set(handles.popup_hpFilter,'String',popStrings');
for j=1:length(popStrings)
    if strcmp(popStrings{j},handles.PreProcess(i).HP.type)
        ind = j;
    end
end
set(handles.popup_hpFilter,'Value',ind);
set(handles.edit_hpCutoff,'String',num2str(handles.PreProcess(i).HP.cutoff));
set(handles.edit_hpOrder,'String',num2str(handles.PreProcess(i).HP.order));
set(handles.edit_hpPB,'String',num2str(handles.PreProcess(i).HP.rpass_ripple));
set(handles.edit_hpSB,'String',num2str(handles.PreProcess(i).HP.rstop_ripple));

% Pop LP filter boxes
popStrings = {'butter'; 'cheby1'; 'cheby2'; 'ellip'};
set(handles.popup_lpFilter,'String',popStrings');
for j=1:length(popStrings)
    if strcmp(popStrings{j},handles.PreProcess(i).LP.type)
        ind = j;
    end
end
set(handles.popup_lpFilter,'Value',ind);
set(handles.edit_lpCutoff,'String',num2str(handles.PreProcess(i).LP.cutoff));
set(handles.edit_lpOrder,'String',num2str(handles.PreProcess(i).LP.order));
set(handles.edit_lpPB,'String',num2str(handles.PreProcess(i).LP.rpass_ripple));
set(handles.edit_lpSB,'String',num2str(handles.PreProcess(i).LP.rstop_ripple));

% Pop Amplifier boxes
set(handles.edit_ampGain,'String',num2str(handles.PreProcess(i).Amp.gain));

% Pop Down Sample boxes
factor = floor(handles.PreProcess(i).params.Fs/handles.PreProcess(i).DSample.target);
trueRate = floor(handles.PreProcess(i).params.Fs*100/factor)/100;
set(handles.edit_dsampleTarget,'String',num2str(handles.PreProcess(i).DSample.target));
set(handles.text_dsampleRate,'String',num2str(trueRate));

% Pop Smoothing boxes
popStrings = {'moving'; 'lowess'; 'loess'; 'sgolay'; 'rlowess'; 'rloess'};
set(handles.popup_smoothType,'String',popStrings');
for j=1:length(popStrings)
    if strcmp(popStrings{j},handles.PreProcess(i).Smooth.type)
        ind = j;
    end
end
set(handles.popup_smoothType,'Value',ind);
set(handles.edit_smoothWindow,'String',num2str(handles.PreProcess(i).Smooth.window));
set(handles.edit_smoothDegree,'String',num2str(handles.PreProcess(i).Smooth.degree));

function edit_detrendWindow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_detrendWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_detrendWindow as text
%        str2double(get(hObject,'String')) returns contents of edit_detrendWindow as a double

handles.PreProcess(handles.curChan).Detrend.window = str2double(get(handles.edit_detrendWindow,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_detrendWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_detrendWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_detrendWinstep_Callback(hObject, eventdata, handles)
% hObject    handle to edit_detrendWinstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_detrendWinstep as text
%        str2double(get(hObject,'String')) returns contents of edit_detrendWinstep as a double
handles.PreProcess(handles.curChan).Detrend.winstep = str2double(get(handles.edit_detrendWinstep,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_detrendWinstep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_detrendWinstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_pushout.
function push_pushout_Callback(hObject, eventdata, handles)
% hObject    handle to push_pushout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the current position of the GUI from the handles structure
% to pass to the modal dialog.

parentName = getappdata(0, handles.parentName);
fhUpdateAxes = getappdata(parentName, 'fhUpdateDataAxes');

setappdata(parentName, 'PreProcess', handles.PreProcess);
%feval(fhUpdateAxes);
guidata(hObject, handles);


% --- Executes on selection change in popup_Step1.
function popup_Step1_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Step1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Step1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Step1
handles.PreProcess(handles.curChan).Steps(1) = get(handles.popup_Step1,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_Step1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Step1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_Step2.
function popup_Step2_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Step2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Step2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Step2
handles.PreProcess(handles.curChan).Steps(2) = get(handles.popup_Step2,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_Step2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Step2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_Step3.
function popup_Step3_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Step3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Step3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Step3
handles.PreProcess(handles.curChan).Steps(3) = get(handles.popup_Step3,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_Step3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Step3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_Step4.
function popup_Step4_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Step4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Step4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Step4
handles.PreProcess(handles.curChan).Steps(4) = get(handles.popup_Step4,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_Step4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Step4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_channel.
function popup_channel_Callback(hObject, eventdata, handles)
% hObject    handle to popup_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_channel
handles.curChan = get(handles.popup_channel,'Value');
handles = PopulateGUI(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popup_Step5.
function popup_Step5_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Step5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Step5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Step5
handles.PreProcess(handles.curChan).Steps(5) = get(handles.popup_Step5,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_Step5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Step5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_Step6.
function popup_Step6_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Step6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Step6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Step6
handles.PreProcess(handles.curChan).Steps(6) = get(handles.popup_Step6,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_Step6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Step6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_Step7.
function popup_Step7_Callback(hObject, eventdata, handles)
% hObject    handle to popup_Step7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_Step7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_Step7

handles.PreProcess(handles.curChan).Steps(7) = get(handles.popup_Step7,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_Step7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_Step7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_denoiseLP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_denoiseLP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_denoiseLP as text
%        str2double(get(hObject,'String')) returns contents of edit_denoiseLP as a double
handles.PreProcess(handles.curChan).Denoise.LP = str2double(get(handles.edit_denoiseLP,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_denoiseLP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_denoiseLP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_denoiseHP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_denoiseHP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_denoiseHP as text
%        str2double(get(hObject,'String')) returns contents of edit_denoiseHP as a double
handles.PreProcess(handles.curChan).Denoise.HP = str2double(get(handles.edit_denoiseHP,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_denoiseHP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_denoiseHP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_denoiseP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_denoisep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_denoisep as text
%        str2double(get(hObject,'String')) returns contents of edit_denoisep as a double
handles.PreProcess(handles.curChan).Denoise.p = str2double(get(handles.edit_denoiseP,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_denoiseP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_denoisep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_denoiseF0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_denoisef0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_denoisef0 as text
%        str2double(get(hObject,'String')) returns contents of edit_denoisef0 as a double
handles.PreProcess(handles.curChan).Denoise.f0 = str2double(get(handles.edit_denoiseF0,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_denoiseF0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_denoisef0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_hpCutoff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hpCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hpCutoff as text
%        str2double(get(hObject,'String')) returns contents of edit_hpCutoff as a double
handles.PreProcess(handles.curChan).HP.cutoff = str2double(get(handles.edit_hpCutoff,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_hpCutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hpCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_hpOrder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hpOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hpOrder as text
%        str2double(get(hObject,'String')) returns contents of edit_hpOrder as a double
handles.PreProcess(handles.curChan).HP.order = str2double(get(handles.edit_hpOrder,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_hpOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hpOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_hpPB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hpPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hpPB as text
%        str2double(get(hObject,'String')) returns contents of edit_hpPB as a double
handles.PreProcess(handles.curChan).HP.rpass_ripple = str2double(get(handles.edit_hpPB,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_hpPB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hpPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_hpSB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hpSB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hpSB as text
%        str2double(get(hObject,'String')) returns contents of edit_hpSB as a double
handles.PreProcess(handles.curChan).HP.rstop_ripple = str2double(get(handles.edit_hpSB,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_hpSB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hpSB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_hpFilter.
function popup_hpFilter_Callback(hObject, eventdata, handles)
% hObject    handle to popup_hpFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_hpFilter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_hpFilter
popStrings = {'butter'; 'cheby1'; 'cheby2'; 'ellip'};
handles.PreProcess(handles.curChan).HP.type = popStrings{get(handles.popup_hpFilter,'Value')};
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_hpFilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_hpFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_lpFilter.
function popup_lpFilter_Callback(hObject, eventdata, handles)
% hObject    handle to popup_lpFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_lpFilter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_lpFilter
popStrings = {'butter'; 'cheby1'; 'cheby2'; 'ellip'};
handles.PreProcess(handles.curChan).LP.type = popStrings{get(handles.popup_lpFilter,'Value')};
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_lpFilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_lpFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lpSB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lpSB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lpSB as text
%        str2double(get(hObject,'String')) returns contents of edit_lpSB as a double
handles.PreProcess(handles.curChan).LP.rstop_ripple = str2double(get(handles.edit_lpSB,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_lpSB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lpSB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_lpPB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lpPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lpPB as text
%        str2double(get(hObject,'String')) returns contents of edit_lpPB as a double
handles.PreProcess(handles.curChan).LP.rpass_ripple = str2double(get(handles.edit_lpPB,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_lpPB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lpPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_lpOrder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lpOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lpOrder as text
%        str2double(get(hObject,'String')) returns contents of edit_lpOrder as a double
handles.PreProcess(handles.curChan).LP.order = str2double(get(handles.edit_lpOrder,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_lpOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lpOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_lpCutoff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lpCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lpCutoff as text
%        str2double(get(hObject,'String')) returns contents of edit_lpCutoff as a double
handles.PreProcess(handles.curChan).LP.cutoff = str2double(get(handles.edit_lpCutoff,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_lpCutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lpCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dsampleTarget_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dsampleTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dsampleTarget as text
%        str2double(get(hObject,'String')) returns contents of edit_dsampleTarget as a double
handles.PreProcess(handles.curChan).DSample.target = str2double(get(handles.edit_dsampleTarget,'String'));

factor = floor(handles.PreProcess(handles.curChan).params.Fs/handles.PreProcess(handles.curChan).DSample.target);
trueRate = floor(handles.PreProcess(handles.curChan).params.Fs*100/factor)/100;
set(handles.text_dsampleRate,'String',num2str(trueRate));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_dsampleTarget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dsampleTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text_dsampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_dsampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function edit_ampGain_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ampGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ampGain as text
%        str2double(get(hObject,'String')) returns contents of edit_ampGain as a double
handles.PreProcess(handles.curChan).Amp.gain = str2double(get(handles.edit_ampGain,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_ampGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ampGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_smoothWindow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_smoothWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_smoothWindow as text
%        str2double(get(hObject,'String')) returns contents of edit_smoothWindow as a double
handles.PreProcess(handles.curChan).Smooth.window = str2double(get(handles.edit_smoothWindow,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_smoothWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_smoothWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_smoothDegree_Callback(hObject, eventdata, handles)
% hObject    handle to edit_smoothDegree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_smoothDegree as text
%        str2double(get(hObject,'String')) returns contents of edit_smoothDegree as a double
handles.PreProcess(handles.curChan).Smooth.degree = str2double(get(handles.edit_smoothDegree,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_smoothDegree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_smoothDegree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_smoothType.
function popup_smoothType_Callback(hObject, eventdata, handles)
% hObject    handle to popup_smoothType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_smoothType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_smoothType
popStrings = {'moving'; 'lowess'; 'loess'; 'sgolay'; 'rlowess'; 'rloess'};
handles.PreProcess(handles.curChan).Smooth.type = popStrings{get(handles.popup_smoothType,'Value')};
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_smoothType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_smoothType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text69_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in push_save.
function push_save_Callback(hObject, eventdata, handles)
% hObject    handle to push_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PP = handles.PreProcess;
uisave('PP');
guidata(hObject, handles);


% --- Executes on button press in push_import.
function push_import_Callback(hObject, eventdata, handles)
% hObject    handle to push_import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,FilterIndex] = uigetfile('*.mat');
load([PathName FileName]);
handles.PreProcess = PP;
handles.curChan = 1;
handles = PopulateGUI(handles);
guidata(hObject, handles);


% --- Executes on button press in push_copyall.
function push_copyall_Callback(hObject, eventdata, handles)
% hObject    handle to push_copyall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for i = 1:length(handles.PreProcess)
   if i~=handles.curChan
      handles.PreProcess(i) = handles.PreProcess(handles.curChan);
   end
end

guidata(hObject, handles);


% --- Executes on button press in push_addChan.
function push_addChan_Callback(hObject, eventdata, handles)
% hObject    handle to push_addChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PreProcess(length(handles.PreProcess)+1) = createPreProcess;
handles.curChan = length(handles.PreProcess)+1;
handles = PopulateGUI(handles);
guidata(hObject, handles);


function edit_genFs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genFs as text
%        str2double(get(hObject,'String')) returns contents of edit_genFs as a double
handles.PreProcess(handles.curChan).params.Fs = str2double(get(handles.edit_genFs,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genFs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_genZeroPad_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genZeroPad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genZeroPad as text
%        str2double(get(hObject,'String')) returns contents of edit_genZeroPad as a double
handles.PreProcess(handles.curChan).params.pad = str2double(get(handles.edit_genZeroPad,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genZeroPad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genZeroPad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_genTrialAve_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genTrialAve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genTrialAve as text
%        str2double(get(hObject,'String')) returns contents of edit_genTrialAve as a double
handles.PreProcess(handles.curChan).params.trialave = str2double(get(handles.edit_genTrialAve,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genTrialAve_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genTrialAve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_genTW_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genTW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genTW as text
%        str2double(get(hObject,'String')) returns contents of edit_genTW as a double
handles.PreProcess(handles.curChan).params.tapers(1) = str2double(get(handles.edit_genTW,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genTW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genTW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_genK_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genK as text
%        str2double(get(hObject,'String')) returns contents of edit_genK as a double
handles.PreProcess(handles.curChan).params.tapers(2) = str2double(get(handles.edit_genK,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genK_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_genLoLim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genLoLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genLoLim as text
%        str2double(get(hObject,'String')) returns contents of edit_genLoLim as a double
handles.PreProcess(handles.curChan).params.fpass(1) = str2double(get(handles.edit_genLoLim,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genLoLim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genLoLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_genHiLim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genHiLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genHiLim as text
%        str2double(get(hObject,'String')) returns contents of edit_genHiLim as a double
handles.PreProcess(handles.curChan).params.fpass(2) = str2double(get(handles.edit_genHiLim,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genHiLim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genHiLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_genErrorType_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genErrorType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genErrorType as text
%        str2double(get(hObject,'String')) returns contents of edit_genErrorType as a double
handles.PreProcess(handles.curChan).params.error(1) = str2double(get(handles.edit_genErrorType,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genErrorType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genErrorType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_genErrorP_Callback(hObject, eventdata, handles)
% hObject    handle to edit_genErrorP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_genErrorP as text
%        str2double(get(hObject,'String')) returns contents of edit_genErrorP as a double
handles.PreProcess(handles.curChan).params.error(2) = str2double(get(handles.edit_genErrorP,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_genErrorP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_genErrorP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
