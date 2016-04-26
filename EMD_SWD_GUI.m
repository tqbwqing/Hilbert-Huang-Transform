function varargout = EMD_SWD_GUI(varargin)
% EMD_SWD_GUI MATLAB code for EMD_SWD_GUI.fig
%      EMD_SWD_GUI, by itself, creates a new EMD_SWD_GUI or raises the existing
%      singleton*.
%
%      H = EMD_SWD_GUI returns the handle to a new EMD_SWD_GUI or the handle to
%      the existing singleton*.
%
%      EMD_SWD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMD_SWD_GUI.M with the given input arguments.
%
%      EMD_SWD_GUI('Property','Value',...) creates a new EMD_SWD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EMD_SWD_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EMD_SWD_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EMD_SWD_GUI

% Last Modified by GUIDE v2.5 23-Apr-2016 21:08:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EMD_SWD_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @EMD_SWD_GUI_OutputFcn, ...
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


% --- Executes just before EMD_SWD_GUI is made visible.
function EMD_SWD_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EMD_SWD_GUI (see VARARGIN)

% Choose default command line output for EMD_SWD_GUI
handles.output = hObject;
handles.data = [];
handles.imfs = [];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EMD_SWD_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EMD_SWD_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_pushbutton.
function load_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName, PathName] = uigetfile('*.mat','Select a .mat file with your signal');
handles.data = importdata([PathName, FileName]);

axes(handles.time_domain_axis);
plot(handles.data);

guidata(hObject, handles);



function Nstd_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Nstd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nstd_edit as text
%        str2double(get(hObject,'String')) returns contents of Nstd_edit as a double


% --- Executes during object creation, after setting all properties.
function Nstd_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nstd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NE_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NE_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NE_edit as text
%        str2double(get(hObject,'String')) returns contents of NE_edit as a double


% --- Executes during object creation, after setting all properties.
function NE_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NE_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calc_EEMD_pushbutton.
function calc_EEMD_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to calc_EEMD_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = handles.data;
Nstd = str2double(get(handles.Nstd_edit, 'String'));
NE = str2double(get(handles.NE_edit, 'String'));

allmode = eemd(data, Nstd, NE);
imfs = allmode(:, 2:end);

thresh = sum(data.^2);
A = sum(imfs.^2, 1) > 0.1 * thresh ;
idx = find(A > 0);

figure;
for i = 1:1:length(idx)
    subplot(length(idx)+1, 1, i); plot(imfs(:, i));
end
subplot(length(idx)+1, 1, length(idx)+1); plot(sum(imfs(:, idx+1:end), 2));

handles.imfs = imfs;

guidata(hObject, handles);



function freq_bins_edit_Callback(hObject, eventdata, handles)
% hObject    handle to freq_bins_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_bins_edit as text
%        str2double(get(hObject,'String')) returns contents of freq_bins_edit as a double


% --- Executes during object creation, after setting all properties.
function freq_bins_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_bins_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in HS_EEMD_checkbox.
function HS_EEMD_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to HS_EEMD_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HS_EEMD_checkbox
set(handles.HS_pushbutton, 'Visible', 'on');



% --- Executes on button press in HS_pushbutton.
function HS_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to HS_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

isEEMD = get(handles.HS_EEMD_checkbox, 'Value');
isSWD = get(handles.HS_SWD_checkbox, 'Value');

if isEEMD
    
    input = handles.imfs;
    frequency_bins = str2double(get(handles.freq_bins_edit, 'String'));
    freq_range(1) = str2double(get(handles.freq_range_1_edit, 'String'));
    freq_range(2) = str2double(get(handles.freq_range_2_edit, 'String'));
    [eemd_HS, eemd_frequency_axis, eemd_time_axis] = HilbertSpectrum(input, frequency_bins, freq_range);
    
    figure;
    imagesc('XData', eemd_time_axis, 'YData', eemd_frequency_axis, 'CData', eemd_HS);
    axis([eemd_time_axis(1) eemd_time_axis(end) eemd_frequency_axis(1) eemd_frequency_axis(end)]);
    
end

if isSWD
    
    input = handles.swd_components;
    frequency_bins = str2double(get(handles.freq_bins_edit, 'String'));
    freq_range(1) = str2double(get(handles.freq_range_1_edit, 'String'));
    freq_range(2) = str2double(get(handles.freq_range_2_edit, 'String'));
    [eemd_HS, eemd_frequency_axis, eemd_time_axis] = HilbertSpectrum(input, frequency_bins, freq_range);
    
    figure;
    imagesc('XData', eemd_time_axis, 'YData', eemd_frequency_axis, 'CData', eemd_HS);
    axis([eemd_time_axis(1) eemd_time_axis(end) eemd_frequency_axis(1) eemd_frequency_axis(end)]);
    
end



function freq_range_1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to freq_range_1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_range_1_edit as text
%        str2double(get(hObject,'String')) returns contents of freq_range_1_edit as a double


% --- Executes during object creation, after setting all properties.
function freq_range_1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_range_1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freq_range_2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to freq_range_2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_range_2_edit as text
%        str2double(get(hObject,'String')) returns contents of freq_range_2_edit as a double


% --- Executes during object creation, after setting all properties.
function freq_range_2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_range_2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min_peak_edit_Callback(hObject, eventdata, handles)
% hObject    handle to min_peak_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_peak_edit as text
%        str2double(get(hObject,'String')) returns contents of min_peak_edit as a double


% --- Executes during object creation, after setting all properties.
function min_peak_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_peak_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function component_std_edit_Callback(hObject, eventdata, handles)
% hObject    handle to component_std_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of component_std_edit as text
%        str2double(get(hObject,'String')) returns contents of component_std_edit as a double


% --- Executes during object creation, after setting all properties.
function component_std_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to component_std_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calc_SWD_pushbutton.
function calc_SWD_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to calc_SWD_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = handles.data;
min_peak = str2double(get(handles.min_peak_edit, 'String'));
component_std = str2double(get(handles.component_std_edit, 'String'));
param_struct = struct('min_peak_thresh', min_peak, 'component_std', component_std, 'welch_window', 256, 'welch_noverlap', 128);

[swd_components, ~, ~] = SwarmDecomposition(data, param_struct);
number_of_components = size(swd_components, 2);
figure;
for i = 1:1:length(number_of_components)  
    subplot(length(number_of_components), 1, i); plot(swd_components(:, i));
end

handles.swd_components = swd_components;

guidata(hObject, handles);


% --- Executes on button press in HS_SWD_checkbox.
function HS_SWD_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to HS_SWD_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HS_SWD_checkbox
set(handles.HS_pushbutton, 'Visible', 'on');
