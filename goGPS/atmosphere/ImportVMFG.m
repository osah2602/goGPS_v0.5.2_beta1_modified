function varargout = ImportVMFG(varargin)
% IMPORTVMFG MATLAB code for ImportVMFG.fig
%      IMPORTVMFG, by itself, creates a new IMPORTVMFG or raises the existing
%      singleton*.
%
%      H = IMPORTVMFG returns the handle to a new IMPORTVMFG or the handle to
%      the existing singleton*.
%
%      IMPORTVMFG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPORTVMFG.M with the given input arguments.
%
%      IMPORTVMFG('Property','Value',...) creates a new IMPORTVMFG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImportVMFG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImportVMFG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImportVMFG

% Last Modified by GUIDE v2.5 24-Sep-2019 14:48:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImportVMFG_OpeningFcn, ...
                   'gui_OutputFcn',  @ImportVMFG_OutputFcn, ...
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


% --- Executes just before ImportVMFG is made visible.
function ImportVMFG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImportVMFG (see VARARGIN)

% Choose default command line output for ImportVMFG
handles.output = hObject;

%******MOVE FIGURE TO THE CENTER ON SCREEN
movegui(gcf,'center')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImportVMFG wait for user response (see UIRESUME)
% uiwait(handles.fig_importVMFG);


% --- Outputs from this function are returned to the command line.
function varargout = ImportVMFG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb_cancel.
function pb_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(gcf)


% --- Executes on button press in pb_importVMFG_files.
function pb_importVMFG_files_Callback(hObject, eventdata, handles)
% hObject    handle to pb_importVMFG_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
[filename, pathname, filterindex] = uigetfile( ...
{  '*.mat','MAT-files (*.mat)'; ...
   '*.slx;*.mdl','Models (*.slx, *.mdl)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file', ...
   'MultiSelect', 'on')
%Close Figure
% delete(gcf)

%CALLING THE "ImportVMFgrid1" function
% [H0files,H6files,H12files,H18files,orographyFILE] = ImportVMFgrid1() ;

% --- Executes on button press in pb_importnothing.
function pb_importnothing_Callback(hObject, eventdata, handles)
% hObject    handle to pb_importnothing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Close Figure
delete(gcf)

% --- Executes on button press in pb_importVMFG_folders.
function pb_importVMFG_folders_Callback(hObject, eventdata, handles)
% hObject    handle to pb_importVMFG_folders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+
%*************GET CHECKBOX STATE
cb=get(handles.cb_importOROGRAPHY,'value');

%SAVE import OROGRAPHY checkbox state
setappdata(0,'orography',cb)

%Close Figure
delete(gcf)

%CALLING THE "ImportVMFgrid" function
 ImportVMFgrid();




% TRYreadVMF1grid1([2018 01 01],UserXYZ,SatXYZ,H0files,H6files,H12files)
