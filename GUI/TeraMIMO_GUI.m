function varargout = TeraMIMO_GUI(varargin)
% TERAMIMO_GUI MATLAB code for TeraMIMO_GUI.fig
%      TERAMIMO_GUI, by itself, creates a new TERAMIMO_GUI or raises the existing
%      singleton*.
%
%      H = TERAMIMO_GUI returns the handle to a new TERAMIMO_GUI or the handle to
%      the existing singleton*.
%
%      TERAMIMO_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TERAMIMO_GUI.M with the given input arguments.
%
%      TERAMIMO_GUI('Property','Value',...) creates a new TERAMIMO_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TeraMIMO_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TeraMIMO_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TeraMIMO_GUI

% Last Modified by GUIDE v2.5 17-May-2021 12:32:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TeraMIMO_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @TeraMIMO_GUI_OutputFcn, ...
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


% --- Executes just before TeraMIMO_GUI is made visible.
function TeraMIMO_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TeraMIMO_GUI (see VARARGIN)
path(pathdef); addpath(pwd); cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd Molecular_Absorption;addpath(genpath(pwd)); cd ..; 
cd Visualization;addpath(genpath(pwd)); cd ..;
cd GUI;

axes(handles.logo_kaust);
i = imread('kaust_logo_gui.png');
imshow(i);
axis off
axis image
% Choose default command line output for TeraMIMO_GUI
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TeraMIMO_GUI wait for user response (see UIRESUME)
% uiwait(handles.TeraMIMO_gui);


% --- Outputs from this function are returned to the command line.
function varargout = TeraMIMO_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in isl_btn.
function isl_btn_Callback(hObject, eventdata, handles)
% hObject    handle to isl_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ispc
    dos('explorer https://cemse.kaust.edu.sa/isl'); 
elseif ismac
    web https://cemse.kaust.edu.sa/isl;
elseif isunix
    web https://cemse.kaust.edu.sa/isl;
end

% --- Executes on button press in ctl_btn.
function ctl_btn_Callback(hObject, eventdata, handles)
% hObject    handle to ctl_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% dos('explorer https://cemse.kaust.edu.sa/ctl'); 
if ispc
    dos('explorer https://cemse.kaust.edu.sa/ctl'); 
elseif ismac
    web https://cemse.kaust.edu.sa/ctl;
elseif isunix
    web https://cemse.kaust.edu.sa/ctl;
end
function tx_pos_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tx_pos_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tx_pos_edit as text
%        str2double(get(hObject,'String')) returns contents of tx_pos_edit as a double


% --- Executes during object creation, after setting all properties.
function tx_pos_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_pos_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function subb_edit_Callback(hObject, eventdata, handles)
% hObject    handle to subb_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subb_edit as text
%        str2double(get(hObject,'String')) returns contents of subb_edit as a double


% --- Executes during object creation, after setting all properties.
function subb_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subb_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function subc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to subc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subc_edit as text
%        str2double(get(hObject,'String')) returns contents of subc_edit as a double


% --- Executes during object creation, after setting all properties.
function subc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bw_edit_Callback(hObject, eventdata, handles)
% hObject    handle to bw_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bw_edit as text
%        str2double(get(hObject,'String')) returns contents of bw_edit as a double


% --- Executes during object creation, after setting all properties.
function bw_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bw_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fc_edit as text
%        str2double(get(hObject,'String')) returns contents of fc_edit as a double


% --- Executes during object creation, after setting all properties.
function fc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in scenar_pop.
function scenar_pop_Callback(hObject, eventdata, handles)
% hObject    handle to scenar_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns scenar_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scenar_pop


% --- Executes during object creation, after setting all properties.
function scenar_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scenar_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tx_sar_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tx_sar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tx_sar_edit as text
%        str2double(get(hObject,'String')) returns contents of tx_sar_edit as a double


% --- Executes during object creation, after setting all properties.
function tx_sar_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_sar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tx_sac_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tx_sac_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tx_sac_edit as text
%        str2double(get(hObject,'String')) returns contents of tx_sac_edit as a double


% --- Executes during object creation, after setting all properties.
function tx_sac_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_sac_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rx_sar_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rx_sar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rx_sar_edit as text
%        str2double(get(hObject,'String')) returns contents of rx_sar_edit as a double


% --- Executes during object creation, after setting all properties.
function rx_sar_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rx_sar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rx_sac_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rx_sac_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rx_sac_edit as text
%        str2double(get(hObject,'String')) returns contents of rx_sac_edit as a double


% --- Executes during object creation, after setting all properties.
function rx_sac_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rx_sac_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Deltasa_tx_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Deltasa_tx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Deltasa_tx_edit as text
%        str2double(get(hObject,'String')) returns contents of Deltasa_tx_edit as a double


% --- Executes during object creation, after setting all properties.
function Deltasa_tx_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Deltasa_tx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mol_abs_pop.
function mol_abs_pop_Callback(hObject, eventdata, handles)
% hObject    handle to mol_abs_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mol_abs_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mol_abs_pop


% --- Executes during object creation, after setting all properties.
function mol_abs_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mol_abs_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prof1_check.
function prof1_check_Callback(hObject, eventdata, handles)
% hObject    handle to prof1_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of prof1_check


% --- Executes on button press in prof2_check.
function prof2_check_Callback(hObject, eventdata, handles)
% hObject    handle to prof2_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of prof2_check

% --- Executes on button press in mis_check.
function mis_check_Callback(hObject, eventdata, handles)
% hObject    handle to mis_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mis_check



function mis_var_edit_Callback(hObject, eventdata, handles)
% hObject    handle to mis_var_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mis_var_edit as text
%        str2double(get(hObject,'String')) returns contents of mis_var_edit as a double


% --- Executes during object creation, after setting all properties.
function mis_var_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mis_var_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in beam_sp_check.
function beam_sp_check_Callback(hObject, eventdata, handles)
% hObject    handle to beam_sp_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of beam_sp_check


% --- Executes on button press in run_btn.
function run_btn_Callback(hObject, eventdata, handles)
% hObject    handle to run_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ChType_choice = get(handles.chtype_pop,'Value');
  switch ChType_choice
      case 1 %means 'time-invariant'
        p.tv_tiv_channel = 1;
      case 2 % means 'time-variant'
        p.tv_tiv_channel = 0;
  end
if p.tv_tiv_channel

    ChScenario_choice = get(handles.scenar_pop,'Value');
    switch ChScenario_choice
      case 1
          p.channelType = 'LoS';
      case 2
          p.channelType = 'Multipath';
      case 3
          p.channelType = 'Multipath+LoS';
    end
    
    p.Fc = str2double(get(handles.fc_edit,'String'));
    
    if p.Fc < 100e9 || p.Fc > 10e12
        error('Center frequency should be in the frequency range [0.1-10] THz');
    end
    p.BW = str2double(get(handles.bw_edit,'String'));
    if p.BW > 9.9e12
        error('Bandwidth should not exceed 10 THz');
    end
    
    p.Nsub_c = str2double(get(handles.subc_edit,'String'));
    if ~CheckPowerofTwo(p.Nsub_c)
      error('Number of subcarriers should be power of 2');  
    end
    
    p.Nsub_b = str2double(get(handles.subb_edit,'String'));
    if ~CheckPowerofTwo(p.Nsub_b)
      error('Number of sub-bands should be power of 2');  
    end
    
    MolAbs_choice = get(handles.mol_abs_pop,'Value');
    switch MolAbs_choice
      case 1
        p.absorptionType = 'Hitran';
      case 2
          p.absorptionType = 'Approx1';
      case 3
          p.absorptionType = 'Approx2';
    end
        
    if get(handles.prof1_check,'Value') && get(handles.prof2_check,'Value')
        error('You should select only one profile!');
    elseif ~get(handles.prof1_check,'Value') && ~get(handles.prof2_check,'Value')
        error('You should select at least one profile!');    
    else
        if get(handles.prof1_check,'Value')
            % Profile 1
            p.molecules = {'N2','O2','H2O','CO2','CH4'};
            p.moleculesRatio = [76.5450 20.946 1.57 0.033 0.906]/100;
        end
        
        if get(handles.prof2_check,'Value')
            % Profile 2
            p.molecules = {'N2','O2','H2O','CO2','CH4'};
            p.moleculesRatio = [75.295 20.946 2.82 0.033 0.906]/100;
        end
    end
    
    p.positionTx = str2num(get(handles.tx_pos_edit,'String'));
    p.eulerTx = str2num(get(handles.tx_rot_edit,'String'));
    
    if p.eulerTx(1) <= -pi || p.eulerTx(1) > pi
        error('First Tx rotation angle should be in the range ]-pi,pi] !');
    end
    if p.eulerTx(2) < -pi/2 || p.eulerTx(2) > pi/2
        error('Second Tx rotation angle should be in the range [-pi/2,pi/2] !');
    end
    if p.eulerTx(3) <= -pi || p.eulerTx(3) > pi
        error('Third Tx rotation angle should be in the range ]-pi,pi] !');
    end
    
    p.positionRx = str2num(get(handles.rx_pos_edit,'String'));
    p.eulerRx = str2num(get(handles.rx_rot_edit,'String'));
    
    if p.eulerRx(1) <= -pi || p.eulerRx(1) > pi
        error('First Rx rotation angle should be in the range ]-pi,pi] !');
    end
    if p.eulerRx(2) < -pi/2 || p.eulerRx(2) > pi/2
        error('Second Rx rotation angle should be in the range [-pi/2,pi/2] !');
    end
    if p.eulerRx(3) <= -pi || p.eulerRx(3) > pi
        error('Third Rx rotation angle should be in the range ]-pi,pi] !');
    end
    
    p.Mt = str2double(get(handles.tx_sar_edit,'String'));
    p.Nt = str2double(get(handles.tx_sac_edit,'String'));    
    p.Mr = str2double(get(handles.rx_sar_edit,'String'));  
    p.Nr = str2double(get(handles.rx_sac_edit,'String'));  
    
    p.DeltaMt = str2double(get(handles.Deltasa_tx_edit,'String'));
    p.DeltaNt = str2double(get(handles.Deltasa_tx_edit,'String'));
    p.DeltaMr = str2double(get(handles.Deltasa_rx_edit,'String'));
    p.DeltaNr = str2double(get(handles.Deltasa_rx_edit,'String'));
    
    p.Mat = str2double(get(handles.tx_aer_edit,'String'));
    p.Nat = str2double(get(handles.tx_aec_edit,'String'));
    p.Mar = str2double(get(handles.rx_aer_edit,'String'));
    p.Nar = str2double(get(handles.rx_aec_edit,'String'));
     
    p.deltaMt = str2double(get(handles.deltaae_tx_edit,'String'));
    p.deltaNt = str2double(get(handles.deltaae_tx_edit,'String'));
    p.deltaMr = str2double(get(handles.deltaae_rx_edit,'String'));
    p.deltaNr = str2double(get(handles.deltaae_rx_edit,'String'));
     
    p.clusterArrivalRate = str2double(get(handles.clus_ar_edit,'String'))/1e-9;
    p.rayArrivalRate = str2double(get(handles.ray_ar_edit,'String'))/1e-9;
    
    p.clusterDecayFactor = str2double(get(handles.clus_df_edit,'String'))*1e-9;
    p.rayDecayFactor = str2double(get(handles.ray_df_edit,'String'))*1e-9;
    
    WM_choice = get(handles.wpm_pop,'Value');
    switch WM_choice
      case 1
          p.WaveModelSA = 'Plane';
      case 2
          p.WaveModelSA = 'Sphere';
    end
    
    if get(handles.beam_sp_check,'Value')
        p.BeamSplitEffect = 'On';
    end
    
    if get(handles.mis_check,'Value') 
        p.Misalignment = 'On';
    end
    
    p.sigma_s = str2double(get(handles.mis_var_edit,'String'));
    
    p.GaindBi = str2double(get(handles.ant_g_edit,'String'));
     
    % Update Parameters
    p_ch = update_channel_param_TIV(p); 
    % Calculation of Absorption Coefficient
    K_abs = compute_Abs_Coef(p_ch);
    % Call Channel
    [CH_Response, CH_Info] = channel_TIV(p_ch, K_abs);
else
    ChScenario_choice = get(handles.scenar_pop,'Value');
    switch ChScenario_choice
      case 1
          p.channelType = 'LoS';
      case 2
          p.channelType = 'Multipath';
      case 3
          p.channelType = 'Multipath+LoS';
    end
    
    p.Fc = str2double(get(handles.fc_edit,'String'));
    
    if p.Fc < 100e9 || p.Fc > 10e12
        error('Center frequency should be in the frequency range [0.1-10] THz');
    end
    p.BW = str2double(get(handles.bw_edit,'String'));
    if p.BW > 9.9e12
        error('Bandwidth should not exceed 10 THz');
    end
    
    p.Nsub_c = str2double(get(handles.subc_edit,'String'));
    if ~CheckPowerofTwo(p.Nsub_c)
      error('Number of subcarriers should be power of 2');  
    end
    
    p.Nsub_b = str2double(get(handles.subb_edit,'String'));
    if ~CheckPowerofTwo(p.Nsub_b)
      error('Number of sub-bands should be power of 2');  
    end
    
    MolAbs_choice = get(handles.mol_abs_pop,'Value');
    switch MolAbs_choice
      case 1
        p.absorptionType = 'Hitran';
      case 2
          p.absorptionType = 'Approx1';
      case 3
          p.absorptionType = 'Approx2';
    end
        
    if get(handles.prof1_check,'Value') && get(handles.prof2_check,'Value')
        error('You should select only one profile!');
    elseif ~get(handles.prof1_check,'Value') && ~get(handles.prof2_check,'Value')
        error('You should select at least one profile!');    
    else
        if get(handles.prof1_check,'Value')
            % Profile 1
            p.molecules = {'N2','O2','H2O','CO2','CH4'};
            p.moleculesRatio = [76.5450 20.946 1.57 0.033 0.906]/100;
        end
        
        if get(handles.prof2_check,'Value')
            % Profile 2
            p.molecules = {'N2','O2','H2O','CO2','CH4'};
            p.moleculesRatio = [75.295 20.946 2.82 0.033 0.906]/100;
        end
    end
    
    p.positionTx = str2num(get(handles.tx_pos_edit,'String'));
    p.eulerTx = str2num(get(handles.tx_rot_edit,'String'));
    
    if p.eulerTx(1) <= -pi || p.eulerTx(1) > pi
        error('First Tx rotation angle should be in the range ]-pi,pi] !');
    end
    if p.eulerTx(2) < -pi/2 || p.eulerTx(2) > pi/2
        error('Second Tx rotation angle should be in the range [-pi/2,pi/2] !');
    end
    if p.eulerTx(3) <= -pi || p.eulerTx(3) > pi
        error('Third Tx rotation angle should be in the range ]-pi,pi] !');
    end
    
    p.positionRx = str2num(get(handles.rx_pos_edit,'String'));
    p.eulerRx = str2num(get(handles.rx_rot_edit,'String'));
    
    if p.eulerRx(1) <= -pi || p.eulerRx(1) > pi
        error('First Rx rotation angle should be in the range ]-pi,pi] !');
    end
    if p.eulerRx(2) < -pi/2 || p.eulerRx(2) > pi/2
        error('Second Rx rotation angle should be in the range [-pi/2,pi/2] !');
    end
    if p.eulerRx(3) <= -pi || p.eulerRx(3) > pi
        error('Third Rx rotation angle should be in the range ]-pi,pi] !');
    end
    
    p.Mt = str2double(get(handles.tx_sar_edit,'String'));
    p.Nt = str2double(get(handles.tx_sac_edit,'String'));    
    p.Mr = str2double(get(handles.rx_sar_edit,'String'));  
    p.Nr = str2double(get(handles.rx_sac_edit,'String'));  
    
    p.DeltaMt = str2double(get(handles.Deltasa_tx_edit,'String'));
    p.DeltaNt = str2double(get(handles.Deltasa_tx_edit,'String'));
    p.DeltaMr = str2double(get(handles.Deltasa_rx_edit,'String'));
    p.DeltaNr = str2double(get(handles.Deltasa_rx_edit,'String'));
    
    p.Mat = str2double(get(handles.tx_aer_edit,'String'));
    p.Nat = str2double(get(handles.tx_aec_edit,'String'));
    p.Mar = str2double(get(handles.rx_aer_edit,'String'));
    p.Nar = str2double(get(handles.rx_aec_edit,'String'));
     
    p.deltaMt = str2double(get(handles.deltaae_tx_edit,'String'));
    p.deltaNt = str2double(get(handles.deltaae_tx_edit,'String'));
    p.deltaMr = str2double(get(handles.deltaae_rx_edit,'String'));
    p.deltaNr = str2double(get(handles.deltaae_rx_edit,'String'));
     
    p.clusterArrivalRate = str2double(get(handles.clus_ar_edit,'String'))/1e-9;
    p.rayArrivalRate = str2double(get(handles.ray_ar_edit,'String'))/1e-9;
    
    p.clusterDecayFactor = str2double(get(handles.clus_df_edit,'String'))*1e-9;
    p.rayDecayFactor = str2double(get(handles.ray_df_edit,'String'))*1e-9;
    
    WM_choice = get(handles.wpm_pop,'Value');
    switch WM_choice
      case 1
          p.WaveModelSA = 'Plane';
      case 2
          p.WaveModelSA = 'Sphere';
    end
    
    if get(handles.beam_sp_check,'Value')
        p.BeamSplitEffect = 'On';
    end
    
    if get(handles.mis_check,'Value') 
        p.Misalignment = 'On';
    end
    
    p.sigma_s = str2double(get(handles.mis_var_edit,'String'));
    
    p.GaindBi = str2double(get(handles.ant_g_edit,'String'));
    
    DoppShape_choice = get(handles.dop_spec_pop,'Value');
    switch DoppShape_choice
      case 1
          p.DopplerSpecShape = 'Jakes';
      case 2
          p.DopplerSpecShape = 'Flat';
    end
    
    p.nSamplesperFrame = str2double(get(handles.frm_dur_edit,'String'));
    
    % Update Parameters
    p_ch = update_channel_param_TV(p); 
    % Calculation of Absorption Coefficient
    K_abs = compute_Abs_Coef(p_ch);
    % Call Channel
    [CH_Response, CH_Info] = channel_TV(p_ch, K_abs);
end
handles.p_ch = p_ch;
handles.K_abs = K_abs;
handles.CH_Response = CH_Response;
handles.CH_Info = CH_Info;
assignin('base', 'p_ch', p_ch);
assignin('base', 'K_abs', K_abs);
assignin('base', 'CH_Response', CH_Response);
assignin('base', 'CH_Info', CH_Info);
guidata(hObject,handles);
msgbox('Simulation Completed !!');
% --- Executes on selection change in chtype_pop.
function chtype_pop_Callback(hObject, eventdata, handles)
% hObject    handle to chtype_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chtype_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chtype_pop


% --- Executes during object creation, after setting all properties.
function chtype_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chtype_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rx_pos_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rx_pos_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rx_pos_edit as text
%        str2double(get(hObject,'String')) returns contents of rx_pos_edit as a double


% --- Executes during object creation, after setting all properties.
function rx_pos_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rx_pos_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tx_rot_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tx_rot_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tx_rot_edit as text
%        str2double(get(hObject,'String')) returns contents of tx_rot_edit as a double


% --- Executes during object creation, after setting all properties.
function tx_rot_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_rot_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rx_rot_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rx_rot_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rx_rot_edit as text
%        str2double(get(hObject,'String')) returns contents of rx_rot_edit as a double


% --- Executes during object creation, after setting all properties.
function rx_rot_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rx_rot_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function deltaae_tx_edit_Callback(hObject, eventdata, handles)
% hObject    handle to deltaae_tx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deltaae_tx_edit as text
%        str2double(get(hObject,'String')) returns contents of deltaae_tx_edit as a double


% --- Executes during object creation, after setting all properties.
function deltaae_tx_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltaae_tx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rx_aec_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rx_aec_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rx_aec_edit as text
%        str2double(get(hObject,'String')) returns contents of rx_aec_edit as a double


% --- Executes during object creation, after setting all properties.
function rx_aec_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rx_aec_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rx_aer_edit_Callback(hObject, eventdata, handles)
% hObject    handle to rx_aer_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rx_aer_edit as text
%        str2double(get(hObject,'String')) returns contents of rx_aer_edit as a double


% --- Executes during object creation, after setting all properties.
function rx_aer_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rx_aer_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tx_aec_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tx_aec_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tx_aec_edit as text
%        str2double(get(hObject,'String')) returns contents of tx_aec_edit as a double


% --- Executes during object creation, after setting all properties.
function tx_aec_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_aec_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tx_aer_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tx_aer_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tx_aer_edit as text
%        str2double(get(hObject,'String')) returns contents of tx_aer_edit as a double


% --- Executes during object creation, after setting all properties.
function tx_aer_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_aer_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in wpm_pop.
function wpm_pop_Callback(hObject, eventdata, handles)
% hObject    handle to wpm_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns wpm_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wpm_pop


% --- Executes during object creation, after setting all properties.
function wpm_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wpm_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function clus_ar_edit_Callback(hObject, eventdata, handles)
% hObject    handle to clus_ar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clus_ar_edit as text
%        str2double(get(hObject,'String')) returns contents of clus_ar_edit as a double


% --- Executes during object creation, after setting all properties.
function clus_ar_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clus_ar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ray_ar_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ray_ar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ray_ar_edit as text
%        str2double(get(hObject,'String')) returns contents of ray_ar_edit as a double


% --- Executes during object creation, after setting all properties.
function ray_ar_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ray_ar_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function clus_df_edit_Callback(hObject, eventdata, handles)
% hObject    handle to clus_df_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clus_df_edit as text
%        str2double(get(hObject,'String')) returns contents of clus_df_edit as a double


% --- Executes during object creation, after setting all properties.
function clus_df_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clus_df_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ray_df_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ray_df_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ray_df_edit as text
%        str2double(get(hObject,'String')) returns contents of ray_df_edit as a double


% --- Executes during object creation, after setting all properties.
function ray_df_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ray_df_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Deltasa_rx_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Deltasa_rx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Deltasa_rx_edit as text
%        str2double(get(hObject,'String')) returns contents of Deltasa_rx_edit as a double


% --- Executes during object creation, after setting all properties.
function Deltasa_rx_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Deltasa_rx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function deltaae_rx_edit_Callback(hObject, eventdata, handles)
% hObject    handle to deltaae_rx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deltaae_rx_edit as text
%        str2double(get(hObject,'String')) returns contents of deltaae_rx_edit as a double


% --- Executes during object creation, after setting all properties.
function deltaae_rx_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltaae_rx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frm_dur_edit_Callback(hObject, eventdata, handles)
% hObject    handle to frm_dur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frm_dur_edit as text
%        str2double(get(hObject,'String')) returns contents of frm_dur_edit as a double


% --- Executes during object creation, after setting all properties.
function frm_dur_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frm_dur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in dop_spec_pop.
function dop_spec_pop_Callback(hObject, eventdata, handles)
% hObject    handle to dop_spec_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dop_spec_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dop_spec_pop


% --- Executes during object creation, after setting all properties.
function dop_spec_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dop_spec_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ant_g_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ant_g_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ant_g_edit as text
%        str2double(get(hObject,'String')) returns contents of ant_g_edit as a double


% --- Executes during object creation, after setting all properties.
function ant_g_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ant_g_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
