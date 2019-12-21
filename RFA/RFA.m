function varargout = RFA(varargin)
% RFA MATLAB code for RFA.fig
%      RFA, by itself, creates a new RFA or raises the existing
%      singleton*.
%
%      H = RFA returns the handle to a new RFA or the handle to
%      the existing singleton*.
%
%      RFA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RFA.M with the given input arguments.
%
%      RFA('Property','Value',...) creates a new RFA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RFA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RFA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RFA

% Last Modified by GUIDE v2.5 11-Dec-2019 15:00:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RFA_OpeningFcn, ...
                   'gui_OutputFcn',  @RFA_OutputFcn, ...
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


% --- Executes just before RFA is made visible.
function RFA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RFA (see VARARGIN)

% Choose default command line output for RFA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RFA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RFA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% function gambar=filterBinary(gambar)
% gambar=bwareaopen(gambar,10);

function grayImage = grayscale(rgbImage)
try
  [rows, columns, numberOfColorChannels ] = size(rgbImage);
  if numberOfColorChannels  == 3
      
      redChannel = rgbImage(:, :, 1);
      greenChannel = rgbImage(:, :, 2);
      blueChannel = rgbImage(:, :, 3);
      
      grayImage = .299*double(redChannel) + ...
                  .587*double(greenChannel) + ...
                  .114*double(blueChannel);
      
      grayImage = uint8(grayImage);
  else
      
      grayImage = rgbImage;  
  end
catch ME
  errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
    ME.stack(1).name, ME.stack(1).line, ME.message);
  fprintf(1, '%s\n', errorMessage);
  uiwait(warndlg(errorMessage));
end

function gambar=doHPF(gambar)
kernelFilter=[0 -1/4  0;-1/4 2 -1/4;0 -1/4 0];
gambar=imfilter(gambar,kernelFilter,'conv');

function gambar=doLPF(gambar)
kernelFilter=[ 1/9 1/9 1/9;1/9 1/9 1/9;1/9 1/9 1/9; ];
gambar=imfilter(gambar,kernelFilter,'conv');

function preProcessing(hObject, eventdata, handles,gambar)
% hObject    handle to selectFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% gambar     dalam bentuk matrix image

% Noise reduction
%# Create the gaussian filter with hsize = [10 10] and sigma = 2
% https://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm

global centerxtot;
global centerytot;
gambar=grayscale(gambar);
figure(2);
% axes(handles.real);
imshow(gambar);
AInv = imcomplement(gambar);
BInv = imreducehaze(AInv, 'Method','approx','ContrastEnhancement','boost');
enhanc = imcomplement(BInv);
figure(3);
imshow(enhanc)
G = fspecial('gaussian',[5 5],5);
gambar=imfilter(enhanc,G,'same');
figure(4);
imshow(gambar)

% do High pass filter
gambar=doHPF(gambar);
C=double(gambar);
figure(5);
imshow(gambar);
for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Sobel mask for x-direction:
        Gx=((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        %Sobel mask for y-direction:
        Gy=((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
      
        %The gradient of the image
        %B(i,j)=abs(Gx)+abs(Gy);
        gambar(i,j)=sqrt(Gx.^2+Gy.^2);
      
    end
end
figure(6);
imshow(gambar);
% filtering implement low pass filter 
gambar=doLPF(gambar);
figure(7);
imshow(gambar);
%Template Matching
template=imread('Point.png');
% figure(2);
% imshow(template);
template=grayscale(template);
AInvtemplate = imcomplement(template);
BInvtemplate = imreducehaze(AInvtemplate, 'Method','approx','ContrastEnhancement','boost');
enhanctemplate = imcomplement(BInvtemplate);
G = fspecial('gaussian',[5 5],5);
template=imfilter(enhanctemplate,G,'same');

% do High pass filter
template=doHPF(template);
C=double(template);

for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Sobel mask for x-direction:
        Gx=((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        %Sobel mask for y-direction:
        Gy=((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
      
        %The gradient of the image
        %B(i,j)=abs(Gx)+abs(Gy);
        template(i,j)=sqrt(Gx.^2+Gy.^2);
      
    end
end

% filtering implement low pass filter 
template=doLPF(template);

R = normxcorr2(template,gambar);
% figure(1);
% surf(R),shading flat

[ypeak, xpeak] = find(R>0.67);

y=ypeak;
x=xpeak;
    
for n = 1:length(y)-1
    if(abs(ypeak(n)-ypeak(n+1))<7)
        y(n)=0;
        x(n)=0;
    end
end

x=nonzeros(x);
y=nonzeros(y);

yoffset=y-size(template,1);
xoffset=x-size(template,2);



axes(handles.Proc);
imshow(gambar);
imtool(gambar);
centerxtot=[];
centerytot=[];
for m=1:2
%     disp(length(y));
    a=imrect(gca, [xoffset(m)+1, yoffset(m)+1, size(template,2), size(template,1)]);

    x=xoffset(m)+1;
    y=yoffset(m)+1;
    width=(size(template,2));
    height=(size(template,1));

%     disp(x);
%     disp(y);
%     disp(width);
%     disp(height);
%     disp('====');

    centerx=(round(plus(x,(width/2))));
    centerxtot = [centerxtot, centerx];
    centery=(round(plus(y,(height/2))));
    centerytot = [centerytot, centery];
%     disp(centerx);
%     disp(centery);
%     disp('====');
%     disp(centerxtot);
%     disp(centerytot);

end


% --- Executes on button press in BrwImg.
function BrwImg_Callback(hObject, eventdata, handles)
% hObject    handle to BrwImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fullName;

[FileName,Pathname]=uigetfile('*.jpg',sprintf('Pilih sampel untuk scan'));
if FileName==0 
    return
end

fullName=fullfile(Pathname,FileName);

disp(fullName);

imdat=imread(fullName);
% imdat=imresize(imdat,[1080 1920]);
% figure(1);
axes(handles.real);
imshow(imdat);

preProcessing(hObject, eventdata, handles,imdat);

global indikator;



% --- Executes on button press in ClcDat.
function ClcDat_Callback(hObject, eventdata, handles)
% hObject    handle to ClcDat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.real,'reset');
cla(handles.Proc,'reset');
set(handles.grad,'String','0');
set(handles.RFA,'String','0');
set(handles.text8,'String','None');

clear all;
clc;


% --- Executes on button press in Calc.
function Calc_Callback(hObject, eventdata, handles)
% hObject    handle to Calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global centerxtot;
global centerytot;
% disp('==calculateRFA==')

% disp(centerxtot);
% disp(centerytot);

linemidcalc=line(centerxtot,centerytot,'LineWidth',2);
meangradient_calca=mean(gradient([centerxtot],[centerytot]));

if meangradient_calca < 0
   meangradient_calca = meangradient_calca * -1;
else
    meangradient_calca = meangradient_calca;
end

% disp(meangradient_calca);

calca_angle=atan(meangradient_calca)*57.2957795131;
% if calca_angle < 0
%    calca_angle = calca_angle * -1;
% else
%     calca_angle = calca_angle;
% end
% disp(calca_angle);

set(handles.grad,'String',sprintf(num2str(meangradient_calca)));
set(handles.RFA,'String',sprintf(num2str(calca_angle)));

if calca_angle < 4
    prediksi = "Normal";
else
    prediksi = "Pronasi";
end

set(handles.text8,'String',sprintf(num2str(prediksi)));
