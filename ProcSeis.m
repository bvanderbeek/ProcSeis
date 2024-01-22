function varargout = ProcSeis(varargin)
% PROCSEIS MATLAB code for ProcSeis.fig
%      PROCSEIS, by itself, creates a new PROCSEIS or raises the existing
%      singleton*.
%
%      H = PROCSEIS returns the handle to a new PROCSEIS or the handle to
%      the existing singleton*.
%
%      PROCSEIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROCSEIS.M with the given input arguments.
%
%      PROCSEIS('Property','Value',...) creates a new PROCSEIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProcSeis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProcSeis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProcSeis

% Last Modified by GUIDE v2.5 02-May-2022 08:12:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProcSeis_OpeningFcn, ...
                   'gui_OutputFcn',  @ProcSeis_OutputFcn, ...
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- I N I T I A L I Z A T I O N ----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before ProcSeis is made visible.
function ProcSeis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProcSeis (see VARARGIN)

% Choose default command line output for ProcSeis
handles.output = hObject;

% Update fields in handles
handles = ps_reset_handles(handles);

% Make plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ProcSeis wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = ProcSeis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------- D A T A   P A N E L --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. DYNAMIC TEXT

function DTxtDataDir_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtDataDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtDataDir as text
%        str2double(get(hObject,'String')) returns contents of DTxtDataDir as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtDataDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtDataDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtPickDir_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtPickDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtPickDir as text
%        str2double(get(hObject,'String')) returns contents of DTxtPickDir as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtPickDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtPickDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtChan_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtChan as text
%        str2double(get(hObject,'String')) returns contents of DTxtChan as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtChan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtExt_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtExt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtExt as text
%        str2double(get(hObject,'String')) returns contents of DTxtExt as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtExt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtExt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtRefMod_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtRefMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtRefMod as text
%        str2double(get(hObject,'String')) returns contents of DTxtRefMod as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtRefMod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtRefMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtPhase_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtPhase as text
%        str2double(get(hObject,'String')) returns contents of DTxtPhase as a double
% 
% % Update handles structure
% guidata(hObject, handles);

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtPhase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtMinTime_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtMinTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtMinTime as text
%        str2double(get(hObject,'String')) returns contents of DTxtMinTime as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtMinTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtMinTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtMaxTime_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtMaxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtMaxTime as text
%        str2double(get(hObject,'String')) returns contents of DTxtMaxTime as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtMaxTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtMaxTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% 2. PUSH BUTTONS

% --- Executes on button press in ButtLoadData.
function ButtLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ProcSeis_OpeningFcn(hObject,eventdata,handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------- E V E N T   P A N E L -------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. DYNAMIC TEXT

function DTxtEind_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtEind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtEind as text
%        str2double(get(hObject,'String')) returns contents of DTxtEind as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Idiot-proof user-input event index (defaults to be integer within limits
% of number of events)
ievt = round(str2double(handles.DTxtEind.String));
ievt = max(ievt,1);
ievt = min(ievt,length(handles.event.id));
handles.DTxtEind.String = num2str(ievt);

% Re-initialize
ProcSeis_OpeningFcn(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function DTxtEind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtEind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% 2. PUSH BUTTONS

% --- Executes on button press in ButtPreEvt.
function ButtPreEvt_Callback(hObject, eventdata, handles)
% hObject    handle to ButtPreEvt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update event index
ievt = str2double(handles.DTxtEind.String);
ievt = max(ievt - 1,1);
handles.DTxtEind.String = num2str(ievt);

% Re-initialize
ProcSeis_OpeningFcn(hObject,eventdata,handles);



% --- Executes on button press in ButtNxtEvt.
function ButtNxtEvt_Callback(hObject, eventdata, handles)
% hObject    handle to ButtNxtEvt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update event index
ievt = str2double(handles.DTxtEind.String);
ievt = min(ievt + 1,length(handles.event.id));
handles.DTxtEind.String = num2str(ievt);

% Re-initialize
ProcSeis_OpeningFcn(hObject,eventdata,handles);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ F I L T E R   P A N E L ------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. DYNAMIC TEXT

function DTxtMinFreq_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtMinFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtMinFreq as text
%        str2double(get(hObject,'String')) returns contents of DTxtMinFreq as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtMinFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtMinFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtMaxFreq_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtMaxFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtMaxFreq as text
%        str2double(get(hObject,'String')) returns contents of DTxtMaxFreq as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtMaxFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtMaxFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtOrder_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtOrder as text
%        str2double(get(hObject,'String')) returns contents of DTxtOrder as a double

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% 2. PUSH BUTTONS

% --- Executes on button press in ButtFilter.
function ButtFilter_Callback(hObject, eventdata, handles)
% hObject    handle to ButtFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If a new filter is desired, we must re-initialize everything from the
% beginning.
ProcSeis_OpeningFcn(hObject,eventdata,handles);



% 3. CHECK BOXES

% --- Executes on button press in ChkZeroPhase.
function ChkZeroPhase_Callback(hObject, eventdata, handles)
% hObject    handle to ChkZeroPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ChkZeroPhase

% Toggle plot flags to out-of-date
ps_switch_axis_stat(handles.AxStat1,1);
ps_switch_axis_stat(handles.AxStat2,1);

% Update handles structure
guidata(hObject, handles);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------- D I S P L A Y  P A N E L ------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. DYNAMIC TEXT

function DTxtXmin_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtXmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtXmin as text
%        str2double(get(hObject,'String')) returns contents of DTxtXmin as a double

% Update time axis limits
handles.AxStack.XLim(1)  = str2double(handles.DTxtXmin.String);
handles.AxTraces.XLim(1) = str2double(handles.DTxtXmin.String);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtXmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtXmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtXmax_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtXmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtXmax as text
%        str2double(get(hObject,'String')) returns contents of DTxtXmax as a double

% Update time axis limits
handles.AxStack.XLim(2)  = str2double(handles.DTxtXmax.String);
handles.AxTraces.XLim(2) = str2double(handles.DTxtXmax.String);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtXmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtXmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function DTxtGain_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtGain as text
%        str2double(get(hObject,'String')) returns contents of DTxtGain as a double

% Make plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtGain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtGain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtNTrace_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtNTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtNTrace as text
%        str2double(get(hObject,'String')) returns contents of DTxtNTrace as a double

% Get trace indices currently plotted
itrace = strsplit(handles.DTxtNTrace.String,':');
jtrace = str2double(itrace{end});
itrace = str2double(itrace{1});

% Make sure traces are in bounds
itrace = min(max(itrace,1),handles.SeisDat.ntrace);
jtrace = max(min(jtrace,handles.SeisDat.ntrace),1);

% Update trace index text field
handles.DTxtNTrace.String = [num2str(itrace),':',num2str(jtrace)];

% Update plot
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtNTrace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtNTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% 2. PUSH BUTTONS

% --- Executes on button press in ButtPreStat.
function ButtPreStat_Callback(hObject, eventdata, handles)
% hObject    handle to ButtPreStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get trace indices currently plotted
itrace = strsplit(handles.DTxtNTrace.String,':');
jtrace = str2double(itrace{end});
itrace = str2double(itrace{1});
ntrace = 1 + jtrace - itrace;

% Shift trace indices
itrace = max(itrace - ntrace,1);
jtrace = min(itrace + ntrace - 1,handles.SeisDat.ntrace);

% Update trace index text field
handles.DTxtNTrace.String = [num2str(itrace),':',num2str(jtrace)];

% Update plot
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ButtNxtStat.
function ButtNxtStat_Callback(hObject, eventdata, handles)
% hObject    handle to ButtNxtStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get trace indices currently plotted
itrace = strsplit(handles.DTxtNTrace.String,':');
jtrace = str2double(itrace{end});
itrace = str2double(itrace{1});
ntrace = 1 + jtrace - itrace;

% Shift trace indices
jtrace = min(jtrace + ntrace,handles.SeisDat.ntrace);
itrace = max(jtrace - ntrace + 1,1);

% Update trace index text field
handles.DTxtNTrace.String = [num2str(itrace),':',num2str(jtrace)];

% Update plot
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ButtPlotAllTraces.
function ButtPlotAllTraces_Callback(hObject, eventdata, handles)
% hObject    handle to ButtPlotAllTraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get parameters
ichn    = handles.MenuChan.Value;
wlength = str2double(handles.DTxtWLen.String);
tbuff   = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
xmin    = str2double(handles.DTxtXmin.String);
xmax    = str2double(handles.DTxtXmax.String);

% Stacked phase arrival time vector
t = handles.SeisDat.t - handles.picks.t_stack;

% Normalize trace amplitudes
if (wlength > 0)
    A = get_rms_amp(t,squeeze(handles.SeisDat.seis(:,ichn,:)),handles.picks.DT-tbuff(1),wlength+sum(tbuff));
else
    A = ones(handles.SeisDat.ntrace,1);
end
A = A(~handles.picks.tf_bad);

% Normalized seismic traces to plot
seis = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,ichn,:));
seis = seis./repmat(A,1,handles.SeisDat.nsamp);
seis = seis'; % Plot function works on columns

% Phase arrival time array
t = repmat(t(:),1,sum(~handles.picks.tf_bad));
t = t - repmat(handles.picks.DT(~handles.picks.tf_bad)',handles.SeisDat.nsamp,1);

% Plot Traces
H = figure; hold on;
plot(t,seis,'linewidth',1);
box on; grid on; axis tight;
xlim([xmin,xmax]);
title('All Traces');
pbaspect(handles.AxStack.PlotBoxAspectRatio);

% Plot Correlation Window
if wlength > 0
    xbox = [-tbuff(1),wlength+tbuff(2),wlength+tbuff(2),-tbuff(1),-tbuff(1)];
    ybox = [H.Children.YLim(1),H.Children.YLim(1),H.Children.YLim(2),...
        H.Children.YLim(2),H.Children.YLim(1)];
    fill(xbox,ybox,[0.5,0.5,0.5],'EdgeColor','None','FaceAlpha',0.25);
end



% --- Executes on button press in ButtLoadPicks.
function ButtLoadPicks_Callback(hObject, eventdata, handles)
% hObject    handle to ButtLoadPicks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get parameters
EID      = num2str(handles.SeisDat.event);
aPhase   = handles.DTxtPhase.String;
aChannel = handles.MenuChan.String{handles.MenuChan.Value};
T        = num2str(round(1000*str2double(handles.DTxtWLen.String)));
pickDir  = handles.DTxtPickDir.String;

% Construct file name
theFile = ['ProcSeisPicks_evtid_',EID,'_phase_',aPhase,...
    '_channel_',aChannel,'_period_',T,'ms.mat'];

% % Ad-hoc file name for kludges
% theFile = ['ProcSeisPicks_evtid_',EID,'_phase_',aPhase,...
%     '_channel_MDC','_period_',T,'ms.mat'];

% Load pick file and update handles
if isfile([pickDir,'/',theFile])
    fprintf(['\n Loading pick file, ''',theFile,'''.\n']);
    load([pickDir,'/',theFile],'picks');
    
    % Pick indexing
    [~,itrace] = ismember(picks.station,handles.SeisDat.station);
    
    % Reset fields
    handles.picks.t_stack = 0;
    handles.picks.DT      = nan(handles.SeisDat.ntrace,1);
    handles.picks.ddt     = nan(handles.SeisDat.ntrace,1);
    handles.picks.CCF     = zeros(handles.SeisDat.ntrace,1);
    handles.picks.SI      = [];
    % Update pick structure values
    handles.picks.DT(itrace)  = picks.traveltime - handles.SeisDat.tt1D(itrace);
    handles.picks.ddt(itrace) = picks.ddt;
    handles.picks.t_stack     = mean(handles.picks.DT(itrace));
    handles.picks.DT(itrace)  = handles.picks.DT(itrace) - handles.picks.t_stack;
    handles.picks.tf_bad      = isnan(handles.picks.DT);
    handles.picks.tf_mcc      = picks.tf_mcc;
    % Reset bad pick delays
    handles.picks.DT(handles.picks.tf_bad) = 0;
    
    % Make plots
    handles = ps_plot_stack(handles);
    handles = ps_plot_traces(handles);
    
    % Update handles structure
    guidata(hObject, handles);
else
    fprintf(2,['\n Could not find pick file, ''',theFile,'''.\n']);
end



% --- Executes on button press in WienerFilter.
function WienerFilter_Callback(hObject, eventdata, handles)
% hObject    handle to WienerFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Apply Wiener Filter to all traces
handles = ps_WienerFilter(handles);

% Re-make plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);


% 3. CHECK BOXES

% --- Executes on button press in ChkNormAmp.
function ChkNormAmp_Callback(hObject, eventdata, handles)
% hObject    handle to ChkNormAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ChkNormAmp

% Make plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ChkFillSeis.
function ChkFillSeis_Callback(hObject, eventdata, handles)
% hObject    handle to ChkFillSeis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ChkFillSeis

% Make plots
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);



% 4. POP-UP MENUES

% --- Executes on selection change in MenuChan.
function MenuChan_Callback(hObject, eventdata, handles)
% hObject    handle to MenuChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MenuChan contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MenuChan

% Update plot
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MenuChan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MenuChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in MenuCoord.
function MenuCoord_Callback(hObject, eventdata, handles)
% hObject    handle to MenuCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MenuCoord contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MenuCoord

% Can only do rotations if 2 or 3 channels have been loaded
if (handles.SeisDat.nchan > 1) && (handles.SeisDat.nchan < 4)
    % Toggle axes flag to working
    ps_switch_axis_stat(handles.AxStat1,2);
    ps_switch_axis_stat(handles.AxStat2,2);
    
    % Requested coordinate system
    xyz = handles.MenuCoord.String{handles.MenuCoord.Value};
    
    % Undo prior rotations
    for ii = 1:handles.SeisDat.ntrace
%         % Reverse any sign changes introduced after splitting analysis
%         if strcmp(xyz,'FSZ') && isfield(handles.picks,'dt_split')
%             if handles.picks.dt_split(ii) > 0
%                 handles.SeisDat.seis(ii,1,:) = -handles.SeisDat.seis(ii,1,:);
%             end
%         end
        
        % Prior rotation matrix
        R = reshape(handles.SeisDat.Rmat(ii,:),handles.SeisDat.nchan,handles.SeisDat.nchan);
        % Reverse prior rotation
        handles.SeisDat.seis(ii,:,:) = (R')*squeeze(handles.SeisDat.seis(ii,:,:));
    end
    
    % Rotate from ENZ into desired coordinate system
    switch xyz
        case 'ENZ'
            % This is always the starting coordinate system so no further
            % rotation is necessary.
            
            % Reset rotation matrices
            R = eye(handles.SeisDat.nchan);
            handles.SeisDat.Rmat = repmat(R(:)',handles.SeisDat.ntrace,1);
            
            % Update channel menu
            handles.MenuChan.String = strsplit(handles.DTxtChan.String,',')';
            
        case {'TRZ','PAZ','FSZ'}
            
            % Rotate
            [handles.SeisDat.seis,handles.SeisDat.Rmat] = ...
                rotTRZ(handles.SeisDat.seis,handles.SeisDat.baz);
            
            % Update channel menu
            if handles.SeisDat.nchan == 3
                handles.MenuChan.String = {'T';'R';'Z'};
            else
                handles.MenuChan.String = {'T';'R'};
            end
            
            % Define default polarization azimuth in TRZ coordinates which
            % is assumed to parallel the R-channel.
            handles.SeisDat.PAZ = (pi/2);
            
            % Option to further rotate into polarization coordinates
            if strcmp(xyz,'PAZ') || strcmp(xyz,'FSZ')
                % Stack components
                % Get parameters relevant for polarization analysis
                ichn          = handles.MenuChan.Value;
                tf_norm       = handles.ChkNormAmp.Value;
                DT            = handles.picks.DT;
                wlength       = str2double(handles.DTxtWLen.String);
                
                % Define default window for polarization analysis
                if wlength <= 0
                    % Default window is first 100 samples
                    wlength = 100/handles.SeisDat.fs;
                end
                
                % Phase arrival time vector
                t = handles.SeisDat.t - handles.picks.t_stack;
                
                % Define window for polarization analysis
                winpaz = (t >= 0) & (t <= wlength);
                
                % Define amplitude scaling
                if tf_norm
                    % Get RMS trace amplitudes in phase window
                    A = get_rms_amp(t,squeeze(handles.SeisDat.seis(:,ichn,:)),DT,wlength);
                else
                    A = ones(handles.SeisDat.ntrace,1);
                end
                A = A(~handles.picks.tf_bad);
                A = repmat(A,1,handles.SeisDat.nsamp);
                
                % Select good traces, define normalized stack, subset trace
                % to polarization window
                DT = DT(~handles.picks.tf_bad);
                ut = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,1,:));
                ut = stack_seis(t,ut./A,DT);
                ut = ut(winpaz);
                ur = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,2,:));
                ur = stack_seis(t,ur./A,DT);
                ur = ur(winpaz);
                if handles.SeisDat.nchan == 3
                    uz = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,3,:));
                    uz = stack_seis(t,uz./A,DT);
                    uz = uz(winpaz);
                else
                    uz = zeros(1,sum(winpaz));
                end
                
                % Compute polarization
                [azm,elv,D] = compute_polarization(ut,ur,uz);
                
                % Display result
                fprintf(['\n Polarization Azimuth: ',num2str(round(rad2deg(azm))),' deg. \n']);
                fprintf([' Polarization Elevation: ',num2str(round(rad2deg(elv))),' deg. \n']);
                fprintf([' Minimum Flattening Factor: ',num2str(round(100*(D(3)-D(2))/D(3))),'%% \n']);
                
                % The eigen-vector matrix will rotate the seismograms into
                % the principle coordinate system. However, we just want to
                % rotate channel 1 into the polarization direction such
                % channel 2 remains in the plane of polarization.
                V = [cos(elv), 0,  sin(elv);...
                     0,        1,  0       ;...
                    -sin(elv), 0,  cos(elv)];
                V = V*[cos(-azm), -sin(-azm), 0;...
                       sin(-azm),  cos(-azm), 0;...
                       0,          0,         1];
                V = V';
                
%                 % !!!!!!TRACE PAZ!!!!!!!!
%                 [azi,eli,~] = TracePAZ(handles);
%                 % !!!!!!TRACE PAZ!!!!!!!!
                
                % Rotate
                for ii = 1:handles.SeisDat.ntrace
                    
%                     % !!!!!!TRACE PAZ!!!!!!!!
%                     V = [cos(eli(ii)), 0,  sin(eli(ii));...
%                         0,        1,  0       ;...
%                         -sin(eli(ii)), 0,  cos(eli(ii))];
%                     V = V*[cos(-azi(ii)), -sin(-azi(ii)), 0;...
%                         sin(-azi(ii)),  cos(-azi(ii)), 0;...
%                         0,          0,         1];
%                     V = V';
%                     % !!!!!!TRACE PAZ!!!!!!!!
                    
                    % Apply rotation
                    handles.SeisDat.seis(ii,:,:) = (V')*squeeze(handles.SeisDat.seis(ii,:,:));
                    % Store rotation matrix
                    R = reshape(handles.SeisDat.Rmat(ii,:),handles.SeisDat.nchan,handles.SeisDat.nchan);
                    R = (V')*R;
                    handles.SeisDat.Rmat(ii,:) = R(:)';
                end
                
                % Update channel menu
                if handles.SeisDat.nchan == 3
                    handles.MenuChan.String = {'E1';'E2';'E3'};
                else
                    handles.MenuChan.String = {'E1';'E2'};
                end
                
                % Reset channel to PAZ
                handles.MenuChan.Value = 1;
                
                % Store polarization azimuth in PAZ coordinate system where
                % channel 1 is parallel to polarization hence the 90-deg.
                % phase shift with respect to TRZ. So, the PAZ is the
                % CCW-positive angle measured from the radial direction.
                handles.SeisDat.PAZ = azm - (pi/2);
                
%                 % !!!!!!TRACE PAZ!!!!!!!!
%                 handles.SeisDat.PAZ = azi - (pi/2);
%                 handles.SeisDat.PEL = eli;
%                 % !!!!!!TRACE PAZ!!!!!!!!
                
                if strcmp(xyz,'FSZ')
                    handles = RotCorrFSZ(handles);
                    handles.MenuChan.String = {'F';'S';'Z'};
                end
            end
            
        case 'TQL'
            if handles.SeisDat.nchan == 3
                % Rotate
                [handles.SeisDat.seis,handles.SeisDat.Rmat] = ...
                    rotTQL(handles.SeisDat.seis,handles.SeisDat.baz,handles.SeisDat.incidence);
                
                % Update channel menu
                handles.MenuChan.String = {'T';'Q';'L'};
                
                % Define default polarization azimuth in TQL coordinates which
                % is assumed to parallel the Q-channel.
                handles.SeisDat.PAZ = (pi/2);
                
                
                % FORCE PAZ ESTIMATE AND ROTATION
                warning('ROTATING TQL TO FIXED PNL ESTIMATED ON TQ STACK!');
                % Stack components
                % Get parameters relevant for polarization analysis
                ichn          = handles.MenuChan.Value;
                tf_norm       = handles.ChkNormAmp.Value;
                DT            = handles.picks.DT;
                wlength       = str2double(handles.DTxtWLen.String);
                
                % Define default window for polarization analysis
                if wlength <= 0
                    % Default window is first 100 samples
                    wlength = 100/handles.SeisDat.fs;
                end
                
                % Phase arrival time vector
                t = handles.SeisDat.t - handles.picks.t_stack;
                
                % Define window for polarization analysis
                winpaz = (t >= 0) & (t <= wlength);
                
                % Define amplitude scaling
                if tf_norm
                    % Get RMS trace amplitudes in phase window
                    A = get_rms_amp(t,squeeze(handles.SeisDat.seis(:,ichn,:)),DT,wlength);
                else
                    A = ones(handles.SeisDat.ntrace,1);
                end
                A = A(~handles.picks.tf_bad);
                A = repmat(A,1,handles.SeisDat.nsamp);
                
                % Select good traces, define normalized stack, subset trace
                % to polarization window
                % Delays
                DT = DT(~handles.picks.tf_bad);
                % Transverse stack
                ut = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,1,:));
                ut = stack_seis(t,ut./A,DT);
                ut = ut(winpaz);
                % Radial stack
                ur = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,2,:));
                ur = stack_seis(t,ur./A,DT);
                ur = ur(winpaz);
                % Longitudinal stack
                uz = squeeze(handles.SeisDat.seis(~handles.picks.tf_bad,3,:));
                uz = stack_seis(t,uz./A,DT);
                uz = uz(winpaz);
                % Assume zeros
                % uz = zeros(1,sum(winpaz));
                
                % Compute polarization
                [azm,elv] = compute_polarization(ut,ur,uz);
                azm = atan(tan(azm));
                fprintf(['\n Polarization Azimuth  : ',num2str(round(rad2deg(azm))),' deg. \n']);
                fprintf(['\n Polarization Elevation: ',num2str(round(rad2deg(elv))),' deg. (not included in rotation) \n']);
                
                
                % ROTATION INTO CONSTANT POLARISATION DIRECTION
                W  = rad2deg(azm); % CCW-positive from T
                Rp = [cosd(-W), -sind(-W), 0;...
                      sind(-W),  cosd(-W), 0;...
                      0,         0,        1];
                for ii = 1:handles.SeisDat.ntrace
                    % Apply rotation
                    handles.SeisDat.seis(ii,:,:) = Rp*squeeze(handles.SeisDat.seis(ii,:,:));
                    % Store rotation matrix
                    R = reshape(handles.SeisDat.Rmat(ii,:),handles.SeisDat.nchan,handles.SeisDat.nchan);
                    R = Rp*R;
                    handles.SeisDat.Rmat(ii,:) = R(:)';
                end
                handles.SeisDat.PAZ = deg2rad(W - 90);
                % AD-HOC ROTATION INTO CONSTANT POLARISATION DIRECTION
                
                
                
%                 % AD-HOC S-WAVE POLARISATION FILTER
%                 % NOT REVERSIBLE WITHOUT RELOADING DATA
%                 % Note that it doesn't mater if you compute and apply this
%                 % filter in TQL or polarisation-aligned coordinates. The
%                 % effect is the same.
%                 warning('APPLYING S-WAVE POLARISATION FILTER!');
%                 twin = 15;
%                 fs   = handles.SeisDat.fs;
%                 for ii = 1:handles.SeisDat.ntrace
%                     % The i'th seismogram
%                     S = squeeze(handles.SeisDat.seis(ii,:,:))';
%                     
%                     % Compute filter
%                     [R,AZM,ELV,Pfilt,Sfilt] = PolarisationFilter(fs,twin,S);
%                     
%                     % Apply filter
%                     % S-wave polarisation filter
%                     S = S.*repmat(Sfilt,1,3);
%                     % Filter out energy in out-of-polarisation direction?
%                     S = S.*repmat(abs(cos(AZM)),1,3);
%                     
%                     % Store result
%                     handles.SeisDat.seis(ii,1,:) = S(:,1)';
%                     handles.SeisDat.seis(ii,2,:) = S(:,2)';
%                     handles.SeisDat.seis(ii,3,:) = S(:,3)';
%                 end
%                 % AD-HOC POLARISATION FILTER
                
            else
                fprintf(2,'\n No rotation applied. Rotation to TQL requires 3 channels.\n');
            end
            
%         case 'FSZ'
%             handles = RotCorrFSZ(handles);
%             
% %             % Window parameters
% %             tstack  = handles.picks.t_stack;
% %             DT      = handles.picks.DT;
% %             wlength = str2double(handles.DTxtWLen.String);
% %             tbuff   = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
% %             
% %             % Rotate into Fast-Slow-Vertical coordinates as suggested by
% %             % Sieminski et al. (2007).
% %             for ii = 1:handles.SeisDat.ntrace
% %                 % Index station
% %                 [~,ista] = ismember(handles.SeisDat.station(ii),handles.station.id);
% %                 % Get stored fast axis orientation
% %                 alpha = handles.station.rotation.FSZ(ista);
% %                 % Rotate fast direction to channel 1 from ENZ
% %                 R = [cos(alpha),-sin(alpha), 0;...
% %                      sin(alpha), cos(alpha), 0;...
% %                      0,          0,          1];
% %                 R = R(1:handles.SeisDat.nchan,1:handles.SeisDat.nchan);
% %                 % We are undoing the rotation hence the transpose
% %                 R = R';
% %                 % Apply rotation
% %                 handles.SeisDat.seis(ii,:,:) = R*squeeze(handles.SeisDat.seis(ii,:,:));
% %                 % Store rotation matrix
% %                 handles.SeisDat.Rmat(ii,:) = R(:)';
% %                 
% %                 % Window on phase arrival
% %                 t    = handles.SeisDat.t - (tstack + DT(ii));
% %                 lwin = (t >= -tbuff(1)) & (t <= wlength + tbuff(2));
% %                 
% %                 % Check correlation between channels
% %                 maxlag = round(handles.SeisDat.fs*wlength/2);
% %                 r12 = xcorr(squeeze(handles.SeisDat.seis(ii,1,lwin)),...
% %                     squeeze(handles.SeisDat.seis(ii,2,lwin)),maxlag,'coeff');
% %                 ir  = abs(r12) == max(abs(r12));
% %                 pm  = unique(sign(r12(ir)));
% %                 pm(pm == 0) = 1;
% %                 % Force positive correlation
% %                 if length(pm) == 1
% %                     handles.SeisDat.seis(ii,2,:) = pm*handles.SeisDat.seis(ii,2,:);
% %                 else
% %                     error('Positive and negative correlations.');
% %                 end
% %             end
%             
%             % Update channel menu
%             handles.MenuChan.String = {'F';'S';'Z'};
%             
        otherwise
            fprintf(2,['\n No rotation applied. Unrecognized coordinate system ',xyz,'.\n']);
    end
else
    fprintf(2,'No rotation applied. Seismogram rotation requires 2 or 3 channels to be loaded.')
end

% Update plot
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MenuCoord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MenuCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in MenuSort.
function MenuSort_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MenuSort contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MenuSort

% Update plot
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function MenuSort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MenuSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ P H A S E   W I N D O W ------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DTxtWLen_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtWLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtWLen as text
%        str2double(get(hObject,'String')) returns contents of DTxtWLen as a double

% Update plot
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtWLen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtWLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtWTap_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtWTap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtWTap as text
%        str2double(get(hObject,'String')) returns contents of DTxtWTap as a double

% Update plot
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtWTap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtWTap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtPreBuff_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtPreBuff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtPreBuff as text
%        str2double(get(hObject,'String')) returns contents of DTxtPreBuff as a double

% Update plot
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtPreBuff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtPreBuff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTxtPostBuff_Callback(hObject, eventdata, handles)
% hObject    handle to DTxtPostBuff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTxtPostBuff as text
%        str2double(get(hObject,'String')) returns contents of DTxtPostBuff as a double

% Update plot
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function DTxtPostBuff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTxtPostBuff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in ButtPickStack.
function ButtPickStack_Callback(hObject, eventdata, handles)
% hObject    handle to ButtPickStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Activate stack axes
axes(handles.AxStack);
hold(handles.AxStack,'on');

% Picking Loop
fprintf('\n Pick stacked waveform. Press ''q'' to exit picking.\n');
input = 0;
qid   = double('q');
h     = plot(handles.AxStack,[0,0],handles.AxStack.YLim,'-r','linewidth',2);
t0    = 0;
t1    = 0;
while ~(input == qid)
    % Update time
    t1 = t0;
    
    % Update Plotted Pick
    h.XData = [t1,t1];
    drawnow;
    
    % Make New Pick
    [t0,~,input] = ginput(1);
end

% Update stacked arrival time
handles.picks.t_stack = handles.picks.t_stack + t1;

% Update plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ButtPickTrace.
function ButtPickTrace_Callback(hObject, eventdata, handles)
% hObject    handle to ButtPickTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Picking codes
qid = double('q'); % Quit picking
uid = double('u'); % Undo last pick

% Trace indexing
sort_var = handles.MenuSort.String{handles.MenuSort.Value};
if strcmp(sort_var,'traceid')
    sort_var = (1:handles.SeisDat.ntrace)';
else
    sort_var = handles.SeisDat.(sort_var);
end
% Fixed index in SeisDat structure
[~,ind] = sort(sort_var,'ascend');
% Sorted trace index range
jmin = strsplit(handles.DTxtNTrace.String,':');
jmax = str2double(jmin{end});
jmin = str2double(jmin{1});

% Picking fields
h  = cell(jmax-jmin+1,1);

% Activate traces axes
axes(handles.AxTraces);
hold(handles.AxTraces,'on');

% Picking Loop
fprintf('\n Pick individual waveforms. Press ''u'' to undo last pick and ''q'' to exit picking.\n');
input = 0;
while ~(input == qid)
    % Make Pick
    [t,jtrace,input] = ginput(1);
    
    % Consider input
    switch input
        case 1 % Pick Made
            % Get nearest trace index. Must shift trace index to account
            % for the fact that the first trace is plotted at the top of
            % the figure window.
            jtrace = jmax - round(jtrace) + jmin;
            jtrace = max(jtrace,jmin);
            jtrace = min(jtrace,jmax);
            
            % Store prior pick for undo option
            t0      = handles.picks.DT(ind(jtrace));
            jtrace0 = jtrace;
            
            % Store pick
            handles.picks.DT(ind(jtrace)) = handles.picks.DT(ind(jtrace)) + t;
            
            % Plot pick
            delete(h{jmax-jtrace+1});
            h{jmax-jtrace+1} = plot(handles.AxTraces,[t,t],...
                jmax - jtrace + jmin + [-0.4,0.4],'-r','linewidth',2);
            drawnow;
            
        case uid % Undo Last pick
            % Reset pick to prior value
            handles.picks.DT(ind(jtrace0)) = t0;
            
            % Plot pick
            delete(h{jmax-jtrace0+1});
            h{jmax-jtrace0+1} = plot(handles.AxTraces,[t0,t0],...
                jmax - jtrace0 + jmin + [-0.4,0.4],'-r','linewidth',2);
            drawnow;
        otherwise
            % Nothing
    end
    
end

% Update plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ButtExTrace.
function ButtExTrace_Callback(hObject, eventdata, handles)
% hObject    handle to ButtExTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Picking codes
qid = double('q'); % Quit picking

% Trace indexing
sort_var = handles.MenuSort.String{handles.MenuSort.Value};
if strcmp(sort_var,'traceid')
    sort_var = (1:handles.SeisDat.ntrace)';
else
    sort_var = handles.SeisDat.(sort_var);
end
% Fixed index in SeisDat structure
[~,ind] = sort(sort_var);
% Sorted trace index range
jmin = strsplit(handles.DTxtNTrace.String,':');
jmax = str2double(jmin{end});
jmin = str2double(jmin{1});

% Activate traces axes
axes(handles.AxTraces);
hold(handles.AxTraces,'on');

% Picking Loop
fprintf('\n Select traces to exclude/include. Press ''q'' to exit.\n');
input = 0;
while ~(input == qid)
    % Make Pick
    [~,jtrace,input] = ginput(1);
    
    % Consider input
    if input == 1 % Pick Made
        % Get nearest trace index. Must shift trace index to account
        % for the fact that the first trace is plotted at the top of
        % the figure window.
        jtrace = jmax - round(jtrace) + jmin;
        jtrace = max(jtrace,jmin);
        jtrace = min(jtrace,jmax);
        
        % Reverse trace status
        handles.picks.tf_bad(ind(jtrace)) = ~handles.picks.tf_bad(ind(jtrace));
        if handles.picks.tf_bad(ind(jtrace))
            lspec = '-r';
        else
            lspec = '-g';
        end
        
        % Indicate on plot
        plot(handles.AxTraces,handles.AxTraces.XLim,...
            (jmax - jtrace + jmin)*[1,1],lspec,'linewidth',2);
    end
    
end

% Update plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ButtRunSCC.
function ButtRunSCC_Callback(hObject, eventdata, handles)
% hObject    handle to ButtRunSCC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Toggle axes flag to working
ps_switch_axis_stat(handles.AxStat1,2);
ps_switch_axis_stat(handles.AxStat2,2);

% Hard-coded quality controls
ccf_min = 0.6; % Minimum acceptable normalized cross-correlation coefficient
dt_max  = 6; % Maximum allowed delay time (s)

% Define SCC parameters
ichan   = handles.MenuChan.Value;
wlength = str2double(handles.DTxtWLen.String);
tbuff   = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
ttaper  = str2double(handles.DTxtWTap.String);
t       = handles.SeisDat.t - handles.picks.t_stack;

fprintf('\n Running SCC...\n');
[DDT,eddt,CCF] = SCC(t,handles.SeisDat.fs,squeeze(handles.SeisDat.seis(:,ichan,:)),...
    handles.picks.DT,handles.picks.tf_bad,wlength,tbuff,ttaper,handles.ChkNormAmp.Value);
fprintf(' ...Finished SCC.\n');

% Update Delays
handles.picks.DT     = handles.picks.DT + DDT;
handles.picks.ddt    = eddt;
handles.picks.tf_bad = handles.picks.tf_bad | (CCF < ccf_min) | (abs(handles.picks.DT) > dt_max);

% Print some info
fprintf([' Excluded ',num2str(sum(handles.picks.tf_bad)),...
    ' traces with CCF < ',num2str(ccf_min),...
    ' or |DT| > ',num2str(round(1000*dt_max)),' ms.\n']);
fprintf([' Maximum time adjustment is ',num2str(round(1000*max(abs(DDT)))),' ms. \n']);
fprintf([' The estimated RMS error is ',num2str(round(1000*rms(eddt))),' ms. \n']);

% Update plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Reset pick method flag
handles.picks.tf_mcc = false;

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ButtRunMCC.
function ButtRunMCC_Callback(hObject, eventdata, handles)
% hObject    handle to ButtRunMCC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Toggle axes flag to working
ps_switch_axis_stat(handles.AxStat1,2);
ps_switch_axis_stat(handles.AxStat2,2);

% Define MCC parameters
ichan   = handles.MenuChan.Value;
wlength = str2double(handles.DTxtWLen.String);
tbuff   = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
ttaper  = str2double(handles.DTxtWTap.String);
t       = handles.SeisDat.t - handles.picks.t_stack;

fprintf('\n Running MCC...\n');
% Only pass good traces
seis     = squeeze(handles.SeisDat.seis(:,ichan,:));
seis     = seis(~handles.picks.tf_bad,:);
dt_trace = handles.picks.DT(~handles.picks.tf_bad);
% warning('USING COVARIANCE MATRIX DETERMINANT MINIMISATION!');
% tic
% [DT,ddt] = MCC_MDC(t,handles.SeisDat.fs,seis,dt_trace,wlength,tbuff,ttaper,handles.ChkNormAmp.Value);
% toc
tic
[DT,ddt] = MCC(t,handles.SeisDat.fs,seis,dt_trace,wlength,tbuff,ttaper,handles.ChkNormAmp.Value);
toc
fprintf(' ...Finished MCC.\n');
fprintf([' Maximum delay measured is ',num2str(round(1000*max(abs(DT)))),' ms. \n']);
fprintf([' The estimated RMS error is ',num2str(round(1000*rms(ddt(:,1)))),...
    '/',num2str(round(1000*rms(ddt(:,2)))),' ms. \n']);

% Update delays
handles.picks.DT(~handles.picks.tf_bad)  = DT;
handles.picks.ddt(~handles.picks.tf_bad) = ddt(:,2);

% Update plots
handles = ps_plot_stack(handles);
handles = ps_plot_traces(handles);

% Reset pick method flag
handles.picks.tf_mcc = true;

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in ButtPlotDelays.
function ButtPlotDelays_Callback(hObject, eventdata, handles)
% hObject    handle to ButtPlotDelays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hard-coded color axis
clims = [-2,2];

% Index traces to stations
[~,ista] = ismember(handles.SeisDat.station,handles.station.id);
ibad     = handles.picks.tf_bad;

% Make plot
figure; hold on;
scatter(handles.station.longitude(ista(~ibad)),handles.station.latitude(ista(~ibad)),100,...
    handles.picks.DT(~ibad) - mean(handles.picks.DT(~ibad)),'filled');
scatter(handles.station.longitude(ista(ibad)),handles.station.latitude(ista(ibad)),...
    100,[0.5,0.5,0.5],'filled');
axis image;
box on; grid on;
xlabel('longitude (deg.)');
ylabel('latitude (deg.)');
title('Demeaned Station Delays (s)');
colormap(jet); caxis(clims);
colorbar;



% --- Executes on button press in ButtSavePicks.
function ButtSavePicks_Callback(hObject, eventdata, handles)
% hObject    handle to ButtSavePicks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Extract parameters
fmin    = str2double(handles.DTxtMinFreq.String);
fmax    = str2double(handles.DTxtMaxFreq.String);
order   = str2double(handles.DTxtOrder.String);
tf_zero = handles.ChkZeroPhase.Value;
wlength = str2double(handles.DTxtWLen.String);
tbuff   = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];
ttaper  = str2double(handles.DTxtWTap.String);

% Define picks structure with relevant parameters
picks.eventid    = handles.SeisDat.event; % Event ID picked
picks.phase      = handles.DTxtPhase.String; % Phase Picked
picks.channel    = handles.MenuChan.String{handles.MenuChan.Value}; % Channel Picked
picks.station    = handles.SeisDat.station(~handles.picks.tf_bad); % Station ID's picked
picks.traveltime = handles.picks.t_stack + handles.picks.DT(~handles.picks.tf_bad) + handles.SeisDat.tt1D(~handles.picks.tf_bad); % Total travel-times
picks.ddt        = handles.picks.ddt(~handles.picks.tf_bad); % Estimated Errors
if isfield(handles.SeisDat,'Rmat')
    picks.Rmat       = handles.SeisDat.Rmat(~handles.picks.tf_bad,:); % Rotation matrices
end
picks.filtopts   = [fmin,fmax,order,tf_zero]; % Filter parameters
picks.tf_mcc     = handles.picks.tf_mcc; % If true, delays measured via MCC. If false, delays measured via SCC
picks.tf_norm    = handles.ChkNormAmp.Value; % If true, trace normalization was used
picks.wlength    = wlength; % Length of phase window
picks.tbuff      = tbuff; % Buffer around phase window
picks.ttaper     = ttaper; % Taper at ends of phase window
picks.DataDir    = handles.DTxtDataDir.String; % Location of seismic data picked
% Supplementary fields
if isfield(handles.SeisDat,'PAZ')
    picks.PAZ = handles.SeisDat.PAZ;
else
    picks.PAZ = 0;
end

% !!!!!!TRACE PAZ!!!!!!!!
if isfield(handles.SeisDat,'PEL')
    picks.PEL = handles.SeisDat.PEL;
else
    picks.PEL = 0;
end
% !!!!!!TRACE PAZ!!!!!!!!

if isfield(handles.picks,'SI')
    picks.SI = handles.picks.SI;
end

% Construct file name
theFile = ['ProcSeisPicks_evtid_',num2str(picks.eventid),'_phase_',picks.phase,...
    '_channel_',picks.channel,'_period_',num2str(round(1000*picks.wlength)),'ms.mat'];

% !!!!!!TRACE PAZ!!!!!!!!
if isfield(handles.SeisDat,'PEL')
    theFile = ['ProcSeisPicks_evtid_',num2str(picks.eventid),'_phase_',picks.phase,...
        '_channel_',picks.channel,'i_period_',num2str(round(1000*picks.wlength)),'ms.mat'];
end
% !!!!!!TRACE PAZ!!!!!!!!

% theFile = ['ProcSeisPicks_evtid_',num2str(picks.eventid),'_phase_',picks.phase,...
%     '_channel_MDC_','_period_',num2str(round(1000*picks.wlength)),'ms.mat'];

% Save
save([handles.DTxtPickDir.String,'/',theFile],'picks');

fprintf(['\n Saved pick file ''',theFile,'''.\n']);


% --- Executes on button press in ButtRunSplits.
function ButtRunSplits_Callback(hObject, eventdata, handles)
% hObject    handle to ButtRunSplits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf('\n Running splitting analysis... \n');

% Get Parameters
xyz           = handles.MenuCoord.String{handles.MenuCoord.Value};
dt            = handles.picks.DT;
wlength       = str2double(handles.DTxtWLen.String);
tbuff         = [str2double(handles.DTxtPreBuff.String),str2double(handles.DTxtPostBuff.String)];

% Check for valid window
if (wlength + sum(tbuff)) <= (2/handles.SeisDat.fs)
    error('Window too short for splitting analysis.');
end

% Phase arrival time vector
t = handles.SeisDat.t - handles.picks.t_stack;

% Define shear wave components
if strcmp(xyz,'TRZ') || strcmp(xyz,'TQL')
    Rij = squeeze(handles.SeisDat.seis(:,2,:))';
    Tij = squeeze(handles.SeisDat.seis(:,1,:))';
elseif strcmp(xyz,'PAZ')
    % In polarization coordinate system, channel 1 should be taken as the
    % radial while channel 2 is the transverse equivalent
    Rij = squeeze(handles.SeisDat.seis(:,1,:))';
    Tij = squeeze(handles.SeisDat.seis(:,2,:))';
else
    error(['Cannot do splitting analysis in ',xyz,'-coordinate system.']);
end

% % Rotation Correlation Splitting
% keyboard
% maxlag = wlength/2;
% dphi   = 1;
% phirc  = zeros(handles.SeisDat.ntrace,1);
% dtrc   = zeros(handles.SeisDat.ntrace,1);
% for jj = 1:handles.SeisDat.ntrace
%     % Define phase window indices
%     iwin = (t >= (dt(jj) - tbuff(1))) & (t <= (dt(jj) + wlength + tbuff(2)));
%     Rij(~iwin,jj) = 0;
%     Tij(~iwin,jj) = 0;
%     iwin = find(iwin);
%     % Use SplitLab function to compute fast directions and delays via
%     % the rotation correlation method
%     [phirc(jj),dtrc(jj)] = splitRotCorr(Rij(:,jj),Tij(:,jj),iwin(:),...
%         maxlag,1/handles.SeisDat.fs,false,dphi,false);
% end
% % Fast axes are measured with respect to polarization direction (e.g. the 
% % radial axis for SKS). Convert to geographic orientation.
% phirc = deg2rad(phirc) + handles.SeisDat.PAZ + deg2rad(90 - handles.SeisDat.baz);

% Interpolate traces to align phases
Rij = interp2(1:handles.SeisDat.ntrace,t(:),Rij,repmat(1:handles.SeisDat.ntrace,handles.SeisDat.nsamp,1),...
    repmat(t,1,handles.SeisDat.ntrace) + repmat(dt(:)',handles.SeisDat.nsamp,1),'linear',0);
Tij = interp2(1:handles.SeisDat.ntrace,t(:),Tij,repmat(1:handles.SeisDat.ntrace,handles.SeisDat.nsamp,1),...
    repmat(t,1,handles.SeisDat.ntrace) + repmat(dt(:)',handles.SeisDat.nsamp,1),'linear',0);

% Select phase window
lwin = (t >= -tbuff(1)) & (t <= (wlength + tbuff(2)));
Rij  = Rij(lwin(:),:);
Tij  = Tij(lwin(:),:);

% Radial channel derivative
dRijdt = handles.SeisDat.fs*[2*(Rij(2,:) - Rij(1,:));...
    (Rij(3:end,:) - Rij(1:end-2,:)); 2*(Rij(end,:) - Rij(end-1,:))]./2;

% Demean--uneccesary
Tij    = Tij - repmat(mean(Tij,1),sum(lwin),1);
dRijdt = dRijdt - repmat(mean(dRijdt,1),sum(lwin),1);

% Splitting intensity
handles.picks.SI = -2*dot(Tij,dRijdt)./dot(dRijdt,dRijdt);
handles.picks.SI = handles.picks.SI(:);

% Update handles structure
guidata(hObject, handles);

fprintf(' ...Finished splitting analysis... \n');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------- Q U I T   B U T T O N -------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in ButtQuit.
function ButtQuit_Callback(hObject, eventdata, handles)
% hObject    handle to ButtQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure1);
