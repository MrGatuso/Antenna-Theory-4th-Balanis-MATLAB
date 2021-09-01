function SMART
% ********************************************************************
%    Beamforming
%*********************************************************************
%
%    THIS IS A MATLAB BASED PROGRAM THAT COMPUTES THE WEIGHTS AND
%    BEAMFORMED PATTERN OF:
%
%     I.   LINEAR ANTENNA
%     II.  RECTANGULAR PLANAR ANTENNA
%
%     THE LINEAR ARRAY HAS M ELEMENTS PLACED EQUIDISTANTLY ALONG
%     THE X-AXIS. 
%
%     THE RECTANGULAR PLANAR ARRAY HAS M x N ELEMENTS PLACED EQUIDISTANTLY
%     ALONG THE X AND Y AXES.
%     
%     OPTION I.  LINEAR ARRAY
%
%       ** INPUT PARAMETERS:
%
%          1. NUMBER OF ELEMENTS
%          2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%          3. AMPLITUDE AND DIRECTION (THETA)OF SIGNAL-OF-INTEREST(SOI)
%          4. AMPLITUDE AND DIRECTION (THETA)OF SIGNALS-NOT-OF-INTEREST
%             (SNOIs)
%          5. NOISE MEAN AND VARIANCE
%          6. MU OF LMS ALGORITHM
%
%       ** PROGRAM OUTPUT:
%
%          1. ARRAY WEIGHTS (AMPLITUDE AND PHASE)
%          2. BEAMFORMED PATTERN 
%
%  OPTION II.  RECTANGULAR PLANAR ARRAY
%
%       ** INPUT PARAMETERS:
%
%          1. NUMBER OF ARRAY ELEMENTS IN X-DIRECTION
%          2. NUMBER OF ARRAY ELEMENTS IN Y-DIRECTION
%          3. SPACING BETWEEN THE ELEMENTS IN X-DIRECTION (IN WAVELENGTHS)
%          4. SPACING BETWEEN THE ELEMENTS IN Y-DIRECTION (IN WAVELENGTHS)
%          5. AMPLITUDE AND DIRECTION (THETA AND PHI) OF SIGNAL-OF-INTEREST
%             (SOI)
%          6. AMPLITUDE AND DIRECTION (THETA AND PHI) OF SIGNALS-NOT-OF-
%             -INTEREST (SNOIs) 
%          7. NOISE MEAN AND VARIANCE
%          8. MU OF LMS ALGORITHM
%
%       ** PROGRAM OUTPUT:
%
%          1. ARRAY WEIGHTS (AMPLITUDE AND PHASE)
%          2. BEAMFORMED PATTERN 
%
%     NOTES: ALL THE INPUT PARAMETERS ARE IN TERMS OF THE WAVELENGTH.
%************************************************************************
%     Written by: Salvatore Bellofiore, 2002
%     Revised by: Zhiyong Huang, Arizona State University,  10/2/2004
%************************************************************************

% Clean-up variables, figures
clear all;
close all;
clc;
fprintf('\n *** NOTICE: PLEASE CONTACT DR C.A BALANIS IF YOU HAVE ANY QUESTION OR COMMENT. THANKS.\n\n');

%---Choice of array, linear or planar---------------------------------------------
fprintf('Smart antenna structure option: \n\tOption (1): Linear\n\tOption (2): Planar \n');
ERR = 1;
while(ERR ~= 0)
   array = str2num(input('\n Array structure = ','s'));
   if(array == 1)
       fprintf('\n Linear Array is chosen.\n');
       ERR = 0;
       beamforming_linear;
   elseif(array == 2)
       fprintf('\n Planar Array is chosen.\n');
       ERR = 0;
       beamforming_planar;
   else
      fprintf('\n Antenna structure number should be either 1 or 2\n');
   end
end
end
function beamforming_linear
%***********************************************************************
%	Beamforming_linear.m
%***********************************************************************
%	It is a MATLAB function that simulates beamforming for linear arrays.
%***********************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

% Set output data format
format short e

% Generate output filename based on current time
cur_time = datevec(datestr(now));

ERR = 1;
save_in_file = input('\n Do you want to save reusults in file? ([n]/y):','s');
while(ERR ~= 0)
   if (isempty(save_in_file) | save_in_file == 'n')
       fprintf('------------Results will not be saved in file.-----------\n\n');
       ERR=0;
   elseif(save_in_file =='y')
       file_string = input('Input the desired output filename: ','s');
       ERR=0;
       fprintf('\n');
   else
       save_in_file=input('Do you want to save the results in file or not? ([n]/y):','s');
   end;
end



% Start timer
tic;

% Start recording
if (save_in_file=='y')
    out_file = sprintf('%s.txt',file_string);
    diary(out_file);
end
    

% User input
[N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = linear_data_entry;

% Parameters initialization
FIG           = 'figure(1)';             % Figure to record
SKIP_STEP     = 40;                      % Plot every SKIP_STEP iterations
w             = zeros(N,1);              % iteration initialization

% Generate signals
[dd, X, fm] = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern);

for i = 1 : length(dd)
    w0 = w;
    [w, err(i)] = LMS(w,Mu,X(:,i),dd(i));
    mse(i) = sum(abs(err(i))^2);
    w_err(i) = norm(w0 - w);
    if i>1
        array_factor = linear_AF(N,d,w,sig(1:size(sig,1),2),E_pattern);
        if (abs(w_err(i) - w_err(i-1)) < eps) | (array_factor <= AF_thresh)
            linear_plot_pattern(sig,w,N,d,E_pattern,AF_thresh,'half',4,'-');
            break;
        end;
    end;
    if rem(i,SKIP_STEP) == 0             % Plot every SKIP_STEP iterations
        linear_plot_pattern(sig,w,N,d,E_pattern,AF_thresh,'half',4,'-');
        if i == SKIP_STEP
            FIG_HANDLE = eval(FIG);
        if (isunix)
           pause;
        end;
        end;
    end;
end;
    
% Final weights and betas
W    = abs(w);
beta = angle(w);
iterationnumber=i;

%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n');
disp('*************************************************************************');
fprintf('\n');
disp('=== Signal Infomation ===');
string = sprintf('       Amplitude   Theta(degrees)'); disp(string);

[rowofsig,colofsig]=size(sig);
fprintf('SOI     %2.4f        %3d\n',sig(1,1),sig(1,2));
for i=2:rowofsig
    fprintf('SNOI%d   %2.4f        %3d\n',i-1,sig(i,1),sig(i,2));
end

if(isempty(noise))
    fprintf('\n');
    disp('NO noise.');
else
    fprintf('\n');
    disp('=== Noise Information ===');
    fprintf('Mean:  %f,  Variance:  %f\n', noise(1), noise(2));
end

fprintf('\n');
disp('=== Iterations Number of Beamforming ===');              % Displays number of iterations
fprintf('        %d\n',iterationnumber);

fprintf('\n');
disp('=== Weights (Amplitude) of Each Element ===');                 % Displays computed weights
disp('   Exact Value        Normalized');
for i=1:length(W)
    fprintf('%d   %f           %f\n',i,W(i),W(i)/W(1));
end

fprintf('\n');
disp('=== Beta (Phase in degrees) of Each Element ===');             % Displays computed beta [Ensemble]
disp('      Exact Value          Normalized');
nbeta=unwrap(beta)*180/pi;
for i=1:length(beta)
    fprintf('%d   %12f         %12f\n',i,nbeta(i), mod((nbeta(i)-nbeta(1)),360));
end

fprintf('\nThe end.\n');
disp('*************************************************************************');
fprintf('\n');

% Stop recording

if (save_in_file=='y')
    diary off;
    fprintf('Data saved in file %s.txt\n\n', file_string);
end


% Plot results
figure; plot(W/W(1),'*-'); grid;
xlabel('Antenna Element (n)');
ylabel('Excitation of Antenna Element');
title('Normalized Magnitude Distribution');

figure; plot(unwrap(beta)*180/pi,'*-'); grid;
xlabel('Antenna Element (n)');
ylabel('Phase of Antenna Element (degrees)');
title('Phase Distribution');

figure; plot(10*log10(mse)); grid;
xlabel('Iteration Number');
ylabel('MSE in dB');
title('Learning Curve');
    
figure; plot(10*log10(w_err)); grid;
xlabel('Iteration Number');
ylabel('Weight Error in dB');
title('Weight Estimation Error');

% Stop timer
disp('Elapsed Time (min):'); disp(toc/60);


end
function beamforming_planar
%************************************************************************
%	Beamforming_planar
%************************************************************************
%   It is the Matlab function that simulates beamforming of planar array.
%*************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

% Set output format data
format short e

% Generate output filename based on current time
cur_time = datevec(datestr(now));

ERR = 1;

save_in_file = input('\nDo you want to save reusults in file? ([n]/y):','s');
while(ERR ~= 0)
   if (isempty(save_in_file) | save_in_file == 'n')
       fprintf('------------Results will not be saved in file.-----------\n\n');
       ERR=0;
   elseif(save_in_file =='y')
       file_string = input('Input the desired output filename: ','s');
       ERR=0;
       fprintf('\n');
   else
       save_in_file=input('Do you want to save the results in file or not? ([n]/y):','s');
   end;
end

% Start timer
tic;

% Start recording
if (save_in_file=='y')
    out_file = sprintf('%s.txt',file_string);
    diary(out_file);
end

% User input
[N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = planar_data_entry;

% Parameters initialization
FIG           = 'figure(1)';             % Figure to record
w             = zeros(N(1)*N(2),1);

% Generate signals
[dd, X, fm] = planar_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern);

for i = 1 : length(dd),
    w0 = w;
    [w, err(i)] = LMS(w,Mu,X(:,i),dd(i));
    mse(i)   = sum(abs(err(i))^2);
    w_err(i) = norm(w0 - w);
    if i>1                           
        if 10*log10(mse(i)) < -1000    
            break;
        end;
    end;
end;

% Final weights and betas
W    = abs(w);
beta = angle(w);
iterationnumber=i;    

%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n\n');
disp('*************************************************************************');
fprintf('\n');
disp('=== Signal Infomation ===');
disp('       Amplitude   Theta(degrees)    Phi(degrees)'); 
[rowofsig,colofsig]=size(sig);
fprintf('SOI     %2.4f        %3d               %3d\n',sig(1,1),sig(1,2),sig(1,3));
for i=2:rowofsig
    fprintf('SNOI%d   %2.4f        %3d               %3d\n',i-1,sig(i,1),sig(i,2),sig(i,3));
end

if(isempty(noise))
    fprintf('\n');
    disp('NO noise in the system.');
else
    fprintf('\n');
    disp('=== Noise Information ===');
    fprintf('Mean:  %f,  Variance:  %f\n', noise(1), noise(2));
end

fprintf('\n');
disp('=== Iterations Number of the Beamforming ===');              % Displays number of iterations
fprintf('        %d\n',iterationnumber);

fprintf('\n');
disp('===The Exact Weight (Amplitude) of Each Element ===');                 % Displays computed weights in x*y matrix
disp('Row number represents x and colunm number represents y.');
disp('For example, the value in second row and third colunm is the weight of element in position (x,y)=(2,3).') 
fprintf('x\\y    ');
for i=1:N(2)
    fprintf('%d          ',i);
end
fprintf('\n');
for i=1:N(1)
    fprintf('%d   ', i);
    for j=1:N(2)
        fprintf('%f   ',W((i-1)*N(2)+j));
    end
    fprintf('\n');
end

fprintf('\n');
disp('===The Normalized Weight (Amplitude) of Each Element ===');                 % Displays computed weights in x*y matrix
fprintf('They are normalized with respect to amplitude of the first element, (x,y)=(1,1).\n');
disp('Row number represents x and colunm number represents y.');
disp('For example, the value in second row and third colunm is the weight of element in position (x,y)=(2,3).') 
fprintf('x\\y    ');
for i=1:N(2)
    fprintf('%d          ',i);
end
fprintf('\n');
for i=1:N(1)
    fprintf('%d   ', i);
    for j=1:N(2)
        fprintf('%f   ',W((i-1)*N(2)+j)/W(1));
    end
    fprintf('\n');
end

nbeta=unwrap(beta)*180/pi;
fprintf('\n');
disp('===The Exact Beta (Phase in Degrees) of Each Element ===');                 % Displays computed beta [Ensemble]
disp('Row number represents x and colunm number represents y.');
disp('For example, the value in second row and third colunm is the weight of element in position (x,y)=(2,3).') 
fprintf('x\\y      ');
for i=1:N(2)
    fprintf('%d             ',i);
end
fprintf('\n');
for i=1:N(1)
    fprintf('%d   ', i);
    for j=1:N(2)
        fprintf('%12f  ',nbeta((i-1)*N(2)+j));
    end
    fprintf('\n');
end

fprintf('\n');
disp('===The Normalized Beta (Phase in Degrees) of Each Element ===');                 % Displays computed beta [Ensemble]
fprintf('The phase of the first element, (x,y)=(1,1), is set to 0 degree, and other phases are normalized with repect to it.\n');
disp('Row number represents x and colunm number represents y.');
disp('For example, the value in second row and third colunm is the weight of element in position (x,y)=(2,3).') 
fprintf('x\\y    ');
for i=1:N(2)
    fprintf('%d           ',i);
end
fprintf('\n');
for i=1:N(1)
    fprintf('%d   ', i);
    for j=1:N(2)
        fprintf('%10f  ',mod((nbeta((i-1)*N(2)+j)-nbeta(1)),360));
    end
    fprintf('\n');
end

fprintf('\nThe end of the data display. Please have a look at the figures\n');
disp('*************************************************************************');
fprintf('\n');

% Stop recording
if (save_in_file=='y')
    diary off;
end

% Plot results
polar_plot(N,d,w,AF_thresh,E_pattern);
title(['Planar Array Beamforming Pattern (X*Y=', num2str(N(1)),'*',num2str(N(2)),', dx=', num2str(d(1)),'\lambda, dy=', num2str(d(2)), '\lambda, SOI \theta=', num2str(sig(1,2)),'\circ,\phi=',num2str(sig(1,3)), '\circ)']);

AF=planar_plot_pattern(N,d,w,E_pattern,AF_thresh,'planar');
title(['Planar Array Beamforming Pattern (X*Y=', num2str(N(1)),'*',num2str(N(2)),', dx=', num2str(d(1)),'\lambda, dy=', num2str(d(2)), '\lambda, SOI \theta=', num2str(sig(1,2)),'\circ,\phi=',num2str(sig(1,3)), '\circ)']);

figure; plot(W/W(1),'*-'); grid;
xlabel('Antenna Element (n)');
ylabel('Excitation of Antenna Element');
title('Normalized Magnitude Distribution');

figure; plot(unwrap(beta')*180/pi,'*-'); grid;
xlabel('Antenna Element (n)');
ylabel('Phase of Antenna Element (degrees)');
title('Phase Distribution');

figure; plot(10*log10(mse)); grid;
xlabel('Iteration Number');
ylabel('MSE in dB');
title('Learning Curve');

figure; plot(10*log10(w_err)); grid;
xlabel('Iteration Number');
ylabel('Weight Error in dB');
title('Weight Estimation Error');

% Stop timer
disp('Elapsed Time (min):'); disp(toc/60);

% zoom in to see the beamforming pattern in azimuthal plane or elevation
% plane
zoom_in = input('\nDo you want to view the beamforming pattern in azimuthal plane or elevation plane? ([n]/y): ','s');
ERR=1;
while(ERR ~= 0)
   if (isempty(zoom_in) | zoom_in == 'n')
       fprintf('OK. That is all for this simulation.\n\n');
       ERR=0;
   elseif(zoom_in =='y')
       ERR1 = 1;
       while(ERR1 ~= 0)
           fprintf('Plane option: \n\tOption (1): Azimuth\n\tOption (2): Elevation \n');
           plane1 = str2num(input('\n Choice of plane = ','s'));
           if(plane1 == 1)
               fprintf('\n Azimuthal plane is chosen.\n');
               zoom_ele=input('Please enter the value of theta in degrees (between 0 and 180):  ');
               zoom_ele0=zoom_ele;
               if zoom_ele>90    
                   zoom_ele=180-zoom_ele; %simetry of the theta=90.
               end
               figure; 
               plot(linspace(0,360,361),AF(round(zoom_ele*4+1),:));
               grid;
               title(['Beamforming pattern vs \phi at \theta=', num2str(zoom_ele0),'\circ']);
               xlabel('Azimuth (degrees)');
               ylabel('Beamforming pattern (dB)');
           elseif(plane1 == 2)
               fprintf('\n Elevation plane is chosen.\n');
               zoom_azm=input('Please enter the value of phi in degrees (between 0 and 360):  ');
               for i=1:360
                   AF_azm_flip(i)=AF(361-i,round(zoom_azm+1));
               end
               AF_azm=[AF(:,round(zoom_azm+1))',AF_azm_flip];
               figure; 
               plot(linspace(0,180,721),AF_azm);
               grid;
               title(['Beamforming pattern vs \theta at \phi=', num2str(zoom_azm),'\circ']);
               xlabel('Elevation (degrees)');
               ylabel('Beamforming pattern (dB)');
           else
               fprintf('\n Plane number should be either 1 or 2\n');    
           end
               
           zoom_again=input('\n Do you want to try again? ([n]/y): ','s');
           ERR2=1;
           while (ERR2~=0)
               if (isempty(zoom_again) | zoom_again == 'n')
                   fprintf('OK. That is all for this simulation.\n\n')
                   ERR = 0;
                   ERR1 = 0;
                   ERR2 = 0;
               elseif(zoom_again =='y')
                   ERR = 1;
                   ERR1 = 1;
                   ERR2 = 0;
               else
                   zoom_again=input('\n Do you want to try again or not? ([n]/y): ','s');
               end
           end 
       end
   else
      zoom_in = input('\n Do you want to view the beamforming pattern in a specific plane or not? ([n]/y):','s');
   end;
end

% Stop recording
if (save_in_file=='y')
    fprintf('Data saved in file %s.txt\n\n', file_string);
end
end



function Z = linear_AF(N,d,w,angles_0,E_pattern)
%************************************************************************
%	Z = Linear_AF(N,d,w,angles_0,E_pattern)
%************************************************************************
%	It computes the array factor and returns the level of the nulls (in dB)
%
%	Input Parameters Description
%	----------------------------
%	- N          number of elements
%	- d          inter-element spacing (in wavelength)
%                default value is 0.5 (Nyquist rate)
%	- w          weight coefficients (row or column vector)
%   - angle_0    direction of the nulls
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%
%	Output Parameters Description
%	-----------------------------
%   - Z          returns the level of the nulls (in dB)
%*************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%--------------------------------------------------------------------------
%---
k0 = 2*pi;

m     = linspace(0,N-1,N);
theta = linspace(-pi/2,pi/2,181);
    
AF = E_pattern .* sum(diag(w)*exp(i*(k0*d*m'*sin(theta))));
AF = 20*log10(abs(AF)./max(abs(AF)));
    
Z = AF(91+angles_0(2:length(angles_0)));
end



function [N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = linear_data_entry
%************************************************************************
%	[N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = linear_data_entry
%*************************************************************************
%	DATA_ENTRY is a MATLAB function that collects the user's input data in
%	question form.
%
%	Output Parameters Description
%	-----------------------------
%	- N          number of elements
%	- d          inter-element spacing (in wavelength)
%                default value is 0.5 (Nyquist rate)
%   - sig        signal parameters (amplitudes and directions)
%                for SOI and SNOI
%   - noise      noise parameter (mean and variance)
%   - type       sinusoidal or BPSK signal
%   - nn         number of samples
%   - NN         number of samples per cycle or symbol
%   - AF_thresh  nulls' depth
%   - Mu         convergence factor for LMS algorithm
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%--------------------------------------------------------------------------
%---
%%%%%%%%%%%%%%%%%%% Default Values %%%%%%%%%%%%%%%%%%
def_N              = 8;
def_d              = 0.5;
def_SOI            = 1;
def_q              = 1;
def_SNOI           = 1;
def_noise_mean     = 0;
def_noise_var      = 0.1;
def_type_n         = 2;
def_AF_thresh      = -40;
def_Mu             = 0.001;
def_nn             = 500;
def_E_pattern_file = 'linear_isotropic.e';
%---------------------------------------------------%

%%%%%%%%%%%%% Strings initialization %%%%%%%%%%%%%%%%
N_string           = sprintf('Enter number of elements in linear smart antenna [%d]: ',def_N);
d_string           = sprintf('Enter the spacing d (in lambda) between adjacent elements [%2.1f]: ',def_d);
SOI_string         = sprintf('Enter the Pilot signal (SOI) amplitude [%d]: ',def_SOI);
SOId_string        = sprintf('Enter the Pilot signal (SOI) direction (degrees between 0 and 90): ');
q_string           = sprintf('Enter number of interfering signals (SNOI) [%d]: ',def_q);
noise_mean_string  = sprintf('Enter the mean of the noise [%d]: ',def_noise_mean);
noise_var_string   = sprintf('Enter the variance of noise [%2.1f]: ',def_noise_var);
type               = sprintf('Type of signal:\n\t[1] sinusoid\n\t[2] BPSK\nEnter number [%d]: ',def_type_n);
nn_string          = sprintf('Enter the number of data samples [%d]: ',def_nn); 
AF_thresh_string   = sprintf('Enter AF threshold (dB) [%d]: ',def_AF_thresh);
Mu_string          = sprintf('Enter a value for Mu of LMS algorithm (0 < Mu < 1) [%5.4f]: ',def_Mu);
E_pattern_string   = sprintf('Enter element pattern filename (*.e) [%s]: ',def_E_pattern_file);
%---------------------------------------------------%
    
%%%%%%%%%%%%%%%%% Error Messages %%%%%%%%%%%%%%%%%%%%
err_1             = sprintf('\nSignal type not supported...');
%---------------------------------------------------%
    
%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%
% ------------ Number of elements? ------------- %
N = input(N_string);
if isempty(N)
   N = def_N;
end;
% ---------------------------------------------- %

% ------------ Inter-element spacing? ---------- %
d = input(d_string);
if isempty(d)
   d = def_d;
end;
% ---------------------------------------------- %
    
% -------------- SOI amplitude? ---------------- %
SOI = input(SOI_string);
if isempty(SOI)
   SOI = def_SOI;
end;
% ---------------------------------------------- %
    
% -------------- SOI direction? ---------------- %
SOId = [];
while isempty(SOId)
    SOId = input(SOId_string);
end;
% ---------------------------------------------- %
    
% ------ SNOIs amplitudes and directions? ------ %
q = input(q_string);
if isempty(q)
   q = def_q;
end;
    
sig = [SOI SOId];
    
for k = 1 : q,
    SNOI_k_string      = sprintf('Enter the amplitude of No. %d Interference signal (SNOI_%d) [%d]: ',k,k,def_SNOI);
    SNOId_k_string     = sprintf('Enter the direction of No. %d Interference signal (SNOI_%d) (degrees between 0 and 90): ',k,k);
    SNOI_k = input(SNOI_k_string);
    if isempty(SNOI_k)
       SNOI_k = def_SNOI;
    end;
    SNOId_k = [];
    while isempty(SNOId_k)
        SNOId_k = input(SNOId_k_string);
    end;
    sig = [sig; SNOI_k SNOId_k];
end;
% ---------------------------------------------- %
    
% ----------------- Noise data? ---------------- %
noise_string = input('Insert noise? ([y]/n):','s');
if (isempty(noise_string) | noise_string == 'y')
   noise_mean = input(noise_mean_string);
   if isempty(noise_mean)
      noise_mean = def_noise_mean;
   end;
   noise_var = input(noise_var_string);
   if isempty(noise_var)
      noise_var = def_noise_var;
   end;
   noise = [noise_mean noise_var];
else
   fprintf('------------No noise is inserted.-----------\n');
   noise = [];
end;
% ---------------------------------------------- %
    
% ----------------- Signal type? --------------- %
type_n = [];
if isempty(type_n)
   type_n = def_type_n;
end;
switch type_n
case 1
    type = 'sinusoid';
    def_NN = 100;
case 2
    type = 'bpsk';
    def_NN = 1;
otherwise
    error(err_1);
end;
% ---------------------------------------------- %
    
% ------------ Number of data samples? --------- %
nn = input(nn_string);
if isempty(nn)
    nn = def_nn;
end;
% ---------------------------------------------- %

% ------- Number of samples per symbol? -------- %
NN = [];
if isempty(NN)
   NN = def_NN;
end;
% ---------------------------------------------- %

% ---------------- Mu for LMS? ----------------- %
Mu = input(Mu_string);
if isempty(Mu)
   Mu = def_Mu;
end;
% ---------------------------------------------- %

% ---------------- Nulls depth? ---------------- %
if q==0
   AF_thresh = def_AF_thresh;
else
   AF_thresh = input(AF_thresh_string);
   if isempty(AF_thresh)
      AF_thresh = def_AF_thresh;
   end;
end;
% ---------------------------------------------- %

% -------------- Element pattern? -------------- %
E_pattern_file = [];
if isempty(E_pattern_file)
   E_pattern_file = def_E_pattern_file;
end;
E_pattern = load(E_pattern_file)';
% ---------------------------------------------- %

warning off;
end


function linear_plot_pattern(sig,w,N,d,E_pattern,AF_limit,TYPE,rticks,line_style)
%************************************************************************
%	Linear_plot_pattern(w,N,d,E_pattern,AF_limit,TYPE,rticks,line_style)
%************************************************************************
%	PLOT_PATTERN is a MATLAB function that plots 2-D patterns in
%	polar and rectangular coordinates where:
%
%	Input Parameters Description
%	----------------------------
%	- w          weight coefficients (row or column vector)
%	- N          number of elements
%	- d          inter-element spacing (in wavelength)
%                default value is 0.5 (Nyquist rate)
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%	- AF_limit   lower radial tick limit for polar plot (in dB)
%                default value is -40
%	- TYPE       options are:
%			     'half' plots from -90 to 90 degrees
%			     'full' plots from -180 to 180 degrees
%                default value is 'half'   
%   - rticks     is the # of radial ticks (or circles) desired.
%	             default value is 4
%   - line_style is solid (e.g., '-') or dashed (e.g., '--')
%                default value is '-'                 
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%--------------------------------------------------------------------------
%---
k0 = 2*pi;

%%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
%---- default values ----%
def_d = .5;
def_E_pattern = 1;
def_AF_limit = -40;
def_TYPE = 'half';
def_rticks = 4;
def_line_style = '-';
H_FIG = 1;
%------------------------%

switch nargin
case 2
    d = def_d; E_pattern = def_E_pattern; AF_limit = def_AF_limit;
    TYPE = def_TYPE; rticks = def_rticks; line_style = def_line_style;
case 3
    E_pattern = def_E_pattern; AF_limit = def_AF_limit;
    TYPE = def_TYPE; rticks = def_rticks; line_style = def_line_style;
case 4
    AF_limit = def_AF_limit;
    TYPE = def_TYPE; rticks = def_rticks; line_style = def_line_style;
case 5
    TYPE = def_TYPE; rticks = def_rticks; line_style = def_line_style;
case 6
    rticks = def_rticks; line_style = def_line_style;
case 7
    line_style = def_line_style;
end

m = linspace(0,N-1,N);
switch TYPE
case 'half'
    theta_1 = -pi/2; theta_2 = pi/2;
    theta = linspace(-pi/2,pi/2,181);
case 'full'
    theta_1 = -pi; theta_2 = pi;
    theta = linspace(-pi,pi,361);
otherwise
    disp('ERROR => Only two options are allowed for TYPE: half or full');
    return
end 
%-------------------------------------------------------%
    
%%%%%%%%%%%%%% Total Pattern [Element * AF] %%%%%%%%%%%%%
AF = E_pattern .* sum(diag(w)*exp(i*(k0*d*m'*sin(theta))));
AF = 20*log10(abs(AF)./max(abs(AF)));
%-------------------------------------------------------%
    
%%%%%%%%%%%%%%%%%%%%%%%% Polar Plot %%%%%%%%%%%%%%%%%%%%%
figure(H_FIG);
switch TYPE
case 'half'
    linear_semipolar_dB(theta*180/pi,AF,AF_limit,0,rticks,line_style);
case 'full'
    linear_polar_dB(theta*180/pi,AF,AF_limit,0,rticks,line_style);
end

title(['Linear Array Beamforming Pattern ( N=', num2str(N),', d=', num2str(d), '\lambda, SOI \theta=', num2str(sig(1,2)), '\circ)']);

%-------------------------------------------------------%
    
%%%%%%%%%%%%%%%%% Rectangularr Plot %%%%%%%%%%%%%%%%%%%%%
H_FIG = H_FIG + 1; figure(H_FIG);
plot(theta*180/pi,AF,line_style);
axis([theta_1*180/pi theta_2*180/pi AF_limit 0]); grid on;
xlabel('\theta (degrees)','fontsize',14);
ylabel('Magnitude (dB)','fontsize',14);
a=sprintf('\theta (degrees)','fontsize',14);
fs  = get(gca, 'FontSize'); set(gca, 'FontSize', 14);
title(['Linear Array Beamforming Pattern ( N=', num2str(N),', d=', num2str(d), '\lambda, SOI \theta=', num2str(sig(1,2)), '\circ)']);

%-------------------------------------------------------%
end


function hpol = linear_polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%-----------------------------------------------------------------------------
%       linear_polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%
%	POLAR_DB is a MATLAB function that plots 2-D patterns in
%	polar coordinates where:
%		   0      <= THETA (in degrees) <= 360
%		-infinity <  RHO   (in dB)      <  +infinity
%
%	Input Parameters Description
%	----------------------------
%	- theta (in degrees) must be a row vector from 0 to 360 degrees
%	- rho (in dB) must be a row vector
%	- rmin (in dB) sets the minimum limit of the plot (e.g., -60 dB)
%	- rmax (in dB) sets the maximum limit of the plot (e.g.,   0 dB)
%	- rticks is the # of radial ticks (or circles) desired. (e.g., 4)
%	- linestyle is solid (e.g., '-') or dashed (e.g., '--')
%
%	Output Parameters Description
%	-----------------------------
%       - hpol handle of figure
%
%	Credits:
%		S. Bellofiore
%		S. Georgakopoulos
%		A. C. Polycarpou
%		C. Wangsvick
%		C. Bishop
%       Zhiyong Huang
%
%	Tabulate your data accordingly, and call polar_dB to provide the
%	2-D polar plot
%
%	Note:  This function is different from the polar.m (provided by
%	       MATLAB) because RHO is given in dB, and it can be negative
%--------------------------------------------------------------------------
%---
% Convert degrees into radians
theta = theta * pi/180;

% Font size, font style and line width parameters
font_size  = 16;
font_name  = 'Times';
line_width = 1.5;

if nargin < 5
	error('Requires 5 or 6 input arguments.')
elseif nargin == 5 
	if isstr(rho)
		line_style = rho;
		rho = theta;
		[mr,nr] = size(rho);
		if mr == 1
			theta = 1:nr;
		else
			th = (1:mr)';
			theta = th(:,ones(1,nr));
		end
	else
		line_style = 'auto';
	end
elseif nargin == 1
	line_style = 'auto';
	rho = theta;
	[mr,nr] = size(rho);
	if mr == 1
		theta = 1:nr;
	else
		th = (1:mr)';
		theta = th(:,ones(1,nr));
	end
end
if isstr(theta) | isstr(rho)
	error('Input arguments must be numeric.');
end
if any(size(theta) ~= size(rho))
	error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
set(cax, 'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
	'DefaultTextFontName',   font_name, ...
	'DefaultTextFontSize',   font_size, ...
	'DefaultTextFontWeight', get(cax, 'FontWeight') )

% only do grids if hold is off
if ~hold_state

% make a radial grid
	hold on;
  % v returns the axis limits
  % changed the following line to let the y limits become negative
	hhh=plot([0 max(theta(:))],[min(rho(:)) max(rho(:))]);
	v = [get(cax,'xlim') get(cax,'ylim')];
	ticks = length(get(cax,'ytick'));
	delete(hhh);

% check radial limits (rticks)

 	if rticks > 5   % see if we can reduce the number
 		if rem(rticks,2) == 0
 			rticks = rticks/2;
 		elseif rem(rticks,3) == 0
 			rticks = rticks/3;
 		end
 	end

% define a circle
	th = 0:pi/50:2*pi;
	xunit = cos(th);
	yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = [1:(length(th)-1)/4:length(th)];
    xunits(inds(2:2:4)) = zeros(2,1);
    yunits(inds(1:2:5)) = zeros(3,1);

	rinc = (rmax-rmin)/rticks;

% label r
  % change the following line so that the unit circle is not multiplied
  % by a negative number.  Ditto for the text locations.
	for i=(rmin+rinc):rinc:rmax
                is = i - rmin;
		plot(xunit*is,yunit*is,'-','color',tc,'linewidth',0.5);
		text(0,is+rinc/20,['  ' num2str(i)],'verticalalignment','bottom' );
	end
% plot spokes
	th = (1:6)*2*pi/12;
	cst = cos(th); snt = sin(th);
	cs = [-cst; cst];
	sn = [-snt; snt];
	plot((rmax-rmin)*cs,(rmax-rmin)*sn,'-','color',tc,'linewidth',0.5);

% plot the ticks
	george=(rmax-rmin)/30; % Length of the ticks
        th2 = (0:36)*2*pi/72;
        cst2 = cos(th2); snt2 = sin(th2);
	cs2 = [(rmax-rmin-george)*cst2; (rmax-rmin)*cst2];
	sn2 = [(rmax-rmin-george)*snt2; (rmax-rmin)*snt2];
	plot(cs2,sn2,'-','color',tc,'linewidth',0.15); % 0.5
        plot(-cs2,-sn2,'-','color',tc,'linewidth',0.15); % 0.5


% annotate spokes in degrees
  % Changed the next line to make the spokes long enough
	rt = 1.1*(rmax-rmin);
	for i = 1:max(size(th))
		text(rt*cst(i),rt*snt(i),int2str(abs(i*30-90)),'horizontalalignment','center' );
		if i == max(size(th))
			loc = int2str(90);
		elseif i*30+90<=180
			loc = int2str(i*30+90);
                else
                        loc = int2str(180-(i*30+90-180));  
		end
		text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center' );
	end
% set viewto 2-D
	view(0,90);

% set axis limits
  % Changed the next line to scale things properly
	axis((rmax-rmin)*[-1 1 -1.1 1.1]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
	'DefaultTextFontName',   font_name, ...
	'DefaultTextFontSize',   fSize, ...
	'DefaultTextFontWeight', fWeight );

% transform data to Cartesian coordinates.
  % changed the next line so negative rho are not plotted on the other side
  
  for i = 1:length(rho)
    if (rho(i) > rmin)
      if theta(i)*180/pi >=0 & theta(i)*180/pi <=90
          xx(i) = (rho(i)-rmin)*cos(pi/2-theta(i));
          yy(i) = (rho(i)-rmin)*sin(pi/2-theta(i));
      elseif theta(i)*180/pi >=90
          xx(i) = (rho(i)-rmin)*cos(-theta(i)+pi/2);
          yy(i) = (rho(i)-rmin)*sin(-theta(i)+pi/2);
      elseif theta(i)*180/pi < 0 
          xx(i) = (rho(i)-rmin)*cos(abs(theta(i))+pi/2);
          yy(i) = (rho(i)-rmin)*sin(abs(theta(i))+pi/2);
      end
    else
      xx(i) = 0;
      yy(i) = 0;
    end
  end

% plot data on top of grid
if strcmp(line_style,'auto')
	q = plot(xx,yy);
else
	q = plot(xx,yy,line_style);
end
if nargout > 0
	hpol = q;
end
if ~hold_state
	axis('equal');axis('off');
end

% reset hold state
if ~hold_state, set(cax,'NextPlot',next); end
end


function hpol = linear_semipolar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%**********************************************************************
%	semipolar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%**********************************************************************
%	SEMIPOLAR_DB is a MATLAB function that plots 2-D patterns in
%	polar coordinates where:
%		-90        <= THETA (in degrees) <= 90
%		-infinity  <  RHO   (in dB)      <  +infinity
%
%	Input Parameters Description
%	----------------------------
%	- theta (in degrees) must be a row vector from -90 to 90 degrees
%	- rho (in dB) must be a row vector
%	- rmin (in dB) sets the minimum limit of the plot (e.g., -60 dB)
%	- rmax (in dB) sets the maximum limit of the plot (e.g.,   0 dB)
%	- rticks is the # of radial ticks (or circles) desired. (e.g., 4)
%	- linestyle is solid (e.g., '-') or dashed (e.g., '--')
%
%	Output Parameters Description
%	-----------------------------
%	- hpol handle of figure
%************************************************************************
%	Credits:
%		S. Bellofiore
%		S. Georgakopoulos
%		A. C. Polycarpou
%		C. Wangsvick
%		C. Bishop
%
%	Tabulate your data accordingly, and call semipolar_dB to provide the
%	2-D polar plot
%--------------------------------------------------------------------------
%---
% Convert degrees into radians
theta = theta * pi/180;

% Font size, font style and line width parameters
font_size  = 16;
font_name  = 'Times';
line_width = 1.5;

% Parameters initialization
count = 0;

for i=1:length(theta),
     temp(i)=theta(length(theta)+1-i);
end;
theta=temp;

if nargin < 5
	error('Requires 5 or 6 input arguments.')
elseif nargin == 5 
	if isstr(rho)
		line_style = rho;
		rho = theta;
		[mr,nr] = size(rho);
		if mr == 1
			theta = 1:nr;
		else
			th = (1:mr)';
			theta = th(:,ones(1,nr));
		end
	else
		line_style = 'auto';
	end
elseif nargin == 1
	line_style = 'auto';
	rho = theta;
	[mr,nr] = size(rho);
	if mr == 1
		theta = 1:nr;
	else
		th = (1:mr)';
		theta = th(:,ones(1,nr));
	end
end
if isstr(theta) | isstr(rho)
	error('Input arguments must be numeric.');
end
if any(size(theta) ~= size(rho))
	error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
set(cax, 'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
	'DefaultTextFontName',   font_name, ...
	'DefaultTextFontSize',   font_size, ...
	'DefaultTextFontWeight', get(cax, 'FontWeight') )

% only do grids if hold is off
if ~hold_state

% make a radial grid
	hold on;
  % v returns the axis limits
  % changed the following line to let the y limits become negative
	hhh=plot([0 max(theta(:))],[min(rho(:)) max(rho(:))]);
	v = [get(cax,'xlim') get(cax,'ylim')];
	ticks = length(get(cax,'ytick'));
	delete(hhh);

% check radial limits (rticks)

 	if rticks > 10   % see if we can reduce the number
 		if rem(rticks,2) == 0
 			rticks = rticks/2;
 		elseif rem(rticks,3) == 0
 			rticks = rticks/3;
 		end
 	end

% define a circle
	th = 0:pi/50:pi;
	xunit = cos(th);
	yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = [1:floor((length(th)-1)/4):length(th)];
    xunits(inds(2:2:4)) = zeros(2,1);
    yunits(inds(1:2:5)) = zeros(3,1);

    rinc = (rmax-rmin)/rticks;

% label r
  % change the following line so that the unit circle is not multiplied
  % by a negative number.  Ditto for the text locations.
	for i=(rmin+rinc):rinc:rmax
                is = i - rmin;
		plot(xunit*is,yunit*is,'-','color',tc,'linewidth',0.15); % 0.5
		text(is+rinc/20,-0.085*(rmax-rmin),[num2str(i)],'horizontalalignment','center' );
		text(-(is+rinc/20),-0.085*(rmax-rmin),[num2str(i)],'horizontalalignment','center' );
	end

% plot spokes
	th = (0:6)*2*pi/12;
	cst = cos(th); snt = sin(th);
	cs = [-cst*0; cst];
	sn = [-snt*0; snt];
	plot((rmax-rmin)*cs,(rmax-rmin)*sn,'-','color',tc,'linewidth',0.15); % 0.5

% plot the ticks
	tmp_sb = (rmax-rmin)/30; % Length of the ticks
        th2 = (0:36)*2*pi/72;
        cst2 = cos(th2); snt2 = sin(th2);
	cs2 = [(rmax-rmin-tmp_sb)*cst2; (rmax-rmin)*cst2];
	sn2 = [(rmax-rmin-tmp_sb)*snt2; (rmax-rmin)*snt2];
	plot(cs2,sn2,'-','color',tc,'linewidth',0.15); % 0.5

% annotate spokes in degrees
  % Changed the next line to make the spokes long enough
	rt = 1.1*(rmax-rmin);
	for i = 1:max(size(th))
		text(rt*cst(i),rt*snt(i),int2str(abs(i*30-120)),'horizontalalignment','center' );
		if i == max(size(th))
			loc = int2str(0);
		else
			loc = int2str(180+i*30);
		end
	end
% set viewto 2-D
	view(0,90);

end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle, ...
	'DefaultTextFontName',   font_name, ...
	'DefaultTextFontSize',   fSize, ...
	'DefaultTextFontWeight', fWeight );

% transform data to Cartesian coordinates.
  % changed the next line so negative rho are not plotted on the other side

  for i = 1:length(rho)
    if (rho(i) > rmin)
      count = count + 1;
      xx(i) = (rho(i)-rmin)*cos(theta(i)+pi/2);
      yy(i) = (rho(i)-rmin)*sin(theta(i)+pi/2);
    else
      xx(i) = 0;
      yy(i) = 0;
    end
  end

% Text 'dB' at the center of the plot
% ===================================
htext = text(0,-0.085*(rmax-rmin),'dB','horizontalalignment','center');
set(htext,'FontSize',font_size);

% plot data on top of grid
if strcmp(line_style,'auto')
	q = plot(xx,yy);
else
	q = plot(xx,yy,line_style,'linewidth',line_width);
end

if nargout > 0
	hpol = q;
end
if ~hold_state
	axis('equal');axis('off');
end

% reset hold state
if ~hold_state, set(cax,'NextPlot',next); end
end


function [S_SOI, s, m] = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern)
%**********************************************************************
%	[S_SOI, s, m] = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern)
%**********************************************************************
%	SIG_GEN is a MATLAB function that generates the SOI and SNOIs
%   signals (with and without Gaussian noise)
%
%	Input Parameters Description
%	----------------------------
%	- N          number of elements
%	- d          inter-element spacing (in wavelength)
%                default value is 0.5 (Nyquist rate)
%   - nn         number of samples
%   - NN         number of samples per cycle or symbol
%   - type       option: 'sinusoid' or 'bpsk'
%   - sig        signals amplitudes and directions of SOI and SNOI
%   - noise      amplitude and variance values
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%
%	Output Parameters Description
%	-----------------------------
%	- S_SOI      SOI reference signal (column vector)
%   - s          matrix containing all signals including noise
%   - m          frequency multiplier for sinusoidal signals
%*************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%--------------------------------------------------------------------------
%---
%%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
k0 = 2*pi;                     % k in free space
n_sig = size(sig,1);           % Total number of signals
m = [];

x = linspace(0,N-1,N);

y = 0;
sig(:,1) = sig(:,1) .* E_pattern(length(E_pattern) - ((length(E_pattern) - 1) / 2 - sig(:,2))).';
sigr = [sig(:,1) pi/180*sig(:,2) zeros(n_sig,1)];       % degrees to radians conversion

[X,Y] = meshgrid(x,y);

err_1 = sprintf('\nSignal type not supported...');
%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%% Generating SOI %%%%%%%%%%%%%%%%%%%%
PSI_SOI = -k0*d*sin(sigr(1,2))*cos(sigr(1,3))*X - k0*d*sin(sigr(1,2))*sin(sigr(1,3))*Y;

switch type
   case 'sinusoid'
      t = linspace(0,nn-1/NN,NN*nn);
      [T1,PSI_SOI] = meshgrid(t,PSI_SOI);
      s = sigr(1,1)*exp(i*(2*pi*T1+PSI_SOI)); 
      S_SOI  = real(s(1,:));
   case 'bpsk'
      STATE1 = sum(100*clock);
      rand('state',STATE1);
      s_rand = round(rand(nn,1));
      s_rand_index = find(s_rand == 0);
      s_rand(s_rand_index) = -1;
      s_NN = ones(1,NN);
      s_rand = (s_rand*s_NN)';
      s_rand = reshape(s_rand,size(s_rand,1)*size(s_rand,2),1);
      [s_rand,PSI_SOI] = meshgrid(s_rand,PSI_SOI);
      s = sigr(1,1)*(s_rand.*exp(i*PSI_SOI));
      S_SOI  = s(1,:);
   otherwise
      error(err_1);
end;

%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%% Generating SNOI %%%%%%%%%%%%%%%%%%%
for k = 2 : n_sig,
   PSI_SNOI = -k0*d*sin(sigr(k,2))*cos(sigr(k,3))*X - k0*d*sin(sigr(k,2))*sin(sigr(k,3))*Y;
   switch type
      case 'sinusoid'
         t = linspace(0,nn-1/NN,NN*nn);
         [T2,PSI_SNOI] = meshgrid(t,PSI_SNOI);
         m = [1; 100*rand(n_sig-1,1)];                        % frequency multipliers
         SNOI = sigr(k,1)*exp(i*(2*pi*m(k)*T2+PSI_SNOI));
         s = s + SNOI;
      case 'bpsk'
         STATE2 = sum(k*100*clock);
         rand('state',STATE2);
         s_rand = round(rand(nn,1));
         s_rand_index = find(s_rand == 0);
         s_rand(s_rand_index) = -1;
         s_NN = ones(1,NN);
         s_rand = (s_rand*s_NN)';
         s_rand = reshape(s_rand,size(s_rand,1)*size(s_rand,2),1);
         [s_rand,PSI_SNOI] = meshgrid(s_rand,PSI_SNOI);
         SNOI = sigr(k,1)*(s_rand.*exp(i*PSI_SNOI));
         s = s + SNOI;
      otherwise
         error(err_1);
   end;
end;
%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%% Thermal Noise %%%%%%%%%%%%%%%%%%%%
if ~isempty(noise)
    for k = 1 : N
        STATE3 = sum(rand(1)*100*clock);
        randn('state',STATE3);
        noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        STATE4 = sum(rand(1)*100*clock);
        randn('state',STATE4);
        noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        noise_data = complex(noise_data_real,noise_data_imag);
        s(k,:) = s(k,:) + noise_data;
    end;
end;
%-------------------------------------------------------%
end


function [w, error] = LMS(w,Mu,x,d)
%***********************************************************************
%	[w, error] = LMS(w,Mu,x,d)
%***********************************************************************
%	LMS is a MATLAB function that computes the weight coefficients
%   using the Least Mean Square algorithm
%
%	Input Parameters Description
%	----------------------------
%	- w          initial weight coefficients (row or column vector)
%   - Mu         convergence factor
%	- x          input data (column vector)
%	- d          sample of desired or reference signal
%
%	Output Parameters Description
%	-----------------------------
%	- w          updated weight coefficients (row or column vector)
%	- error      d - w' * x
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%--------------------------------------------------------------------------
%---
error = d - w' * x;
w = w + 2 * Mu * x * conj(error);
end


function [N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = planar_data_entry
%************************************************************************
%	[N,d,sig,noise,type,nn,NN,AF_thresh,Mu,E_pattern] = planar_data_entry
%************************************************************************
%	DATA_ENTRY is a MATLAB function that collects the
%   user's input data in question form
%
%	Output Parameters Description
%	-----------------------------
%	- N          number of elements in x axis and y axis
%	- d          inter-element spacing (in wavelength) in x and y axis
%                default value is 0.5 (Nyquist rate)
%   - sig        signal parameters (amplitudes and directions)
%                for SOI and SNOI
%   - noise      noise parameter (mean and variance)
%   - type       sinusoidal or BPSK signal
%   - nn         number of samples
%   - NN         number of samples per cycle or symbol
%   - AF_thresh  nulls' depth
%   - Mu         convergence factor for LMS algorithm
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%--------------------------------------------------------------------------
%---
%%%%%%%%%%%%%%%%%%% Default Values %%%%%%%%%%%%%%%%%%
def_N                   = [8 8];
def_d                   = [0.5 0.5];
def_SOI                 = 1;
def_q                   = 1;
def_SNOI                = 1;
def_noise_mean          = 0;
def_noise_var           = 0.1;
def_type_n              = 2;
def_AF_thresh           = -40;
def_Mu                  = 0.001;
def_nn                  = 500;
def_E_pattern_file      = 'planar_isotropic.e';
%---------------------------------------------------%
   
%%%%%%%%%%%%% Strings initialization %%%%%%%%%%%%%%%%
N_string_x              = sprintf('Enter number of elements in the x-axis [%d]: ',def_N(1));
N_string_y              = sprintf('Enter number of elements in the y-axis [%d]: ',def_N(2));
d_string_x              = sprintf('Enter the spacing in the x-axis, dx (lambda) [%2.1f]: ',def_d(1));
d_string_y              = sprintf('Enter the spacing in the y-axis, dy (lambda) [%2.1f]: ',def_d(2));
SOI_string              = sprintf('Enter the Pilot signal (SOI) amplitude [%d]: ',def_SOI);
SOId_string_theta       = sprintf('Enter the Pilot signal (SOI) in theta direction (degrees between 0 and 180): ');
SOId_string_phi         = sprintf('Enter the Pilot singal (SOI) in phi direction (degrees between 0 and 360): ');
q_string                = sprintf('Enter number of interfering signals (SNOI) [%d]: ',def_q);
noise_mean_string       = sprintf('Enter the mean of the noise [%d]: ',def_noise_mean);
noise_var_string        = sprintf('Enter the variance of the noise [%2.1f]: ',def_noise_var);
type                    = sprintf('Type of signal:\n\t[1] sinusoid\n\t[2] BPSK\nEnter number [%d]: ',def_type_n);
nn_string               = sprintf('Enter the number of data samples [%d]: ',def_nn); 
AF_thresh_string        = sprintf('Enter AF threshold (dB) [%d]: ',def_AF_thresh);
Mu_string               = sprintf('Enter a value of Mu in LMS algorithm (0 < Mu < 1) [%4.3f]: ',def_Mu);
E_pattern_string        = sprintf('Enter element pattern filename (*.e) [%s]: ',def_E_pattern_file);
%---------------------------------------------------%

%%%%%%%%%%%%%%%%% Error Messages %%%%%%%%%%%%%%%%%%%%
err_1                   = sprintf('\nSignal type is not supported...');
err_2                   = sprintf('\nAlgorithm type is not supported...');
err_3                   = sprintf('\nMethod is not supported...');
err_4                   = sprintf('\nThis particular weights initialization Method is not supported...');
%---------------------------------------------------%
    
%%%%%%%%%%%%%%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%
% ------------ Number of elements? ------------- %
N_x = input(N_string_x);
if isempty(N_x)
    N_x = def_N(1);
end;
N_y = input(N_string_y);
if isempty(N_y)
    N_y = def_N(2);
end;
N = [N_x N_y];
% ---------------------------------------------- %

% ------------ Inter-element spacing? ---------- %
d_x = input(d_string_x);
if isempty(d_x)
    d_x = def_d(1);
end;
d_y = input(d_string_y);
if isempty(d_y)
    d_y = def_d(2);
end;
d = [d_x d_y];
% ---------------------------------------------- %
    
% -------------- SOI amplitude? ---------------- %
SOI = input(SOI_string);
if isempty(SOI)
    SOI = def_SOI;
end;
% ---------------------------------------------- %
    
% -------------- SOI direction? ---------------- %
SOId_theta = []; SOId_phi = [];
while isempty(SOId_theta)
    SOId_theta = input(SOId_string_theta);
end;
while isempty(SOId_phi)
    SOId_phi = input(SOId_string_phi);
end;
% ---------------------------------------------- %
    
% ------ SNOIs amplitudes and directions? ------ %
q = input(q_string);
if isempty(q)
    q = def_q;
end;
    
sig = [SOI SOId_theta SOId_phi];
    
for k = 1 : q,
    SNOI_k_string           = sprintf('Enter the amplitude of No. %d Interference singal (SNOI_%d)  [%d]: ',k,k, def_SNOI);
    SNOId_k_string_theta    = sprintf('Enter the theta direction of No. %d Interference singal (SNOI_%d) (degrees between 0 and 180): ',k,k);
    SNOId_k_string_phi      = sprintf('Enter the phi direction of No. %d Interference singal (SNOI_%d) (degrees between 0 and 360): ',k,k);
    SNOI_k = input(SNOI_k_string);
    if isempty(SNOI_k)
        SNOI_k = def_SNOI;
    end;
    SNOId_k_theta = []; SNOId_k_phi = [];
    while isempty(SNOId_k_theta)
        SNOId_k_theta = input(SNOId_k_string_theta);
    end;
    while isempty(SNOId_k_phi)
        SNOId_k_phi = input(SNOId_k_string_phi);
    end;
    sig = [sig; SNOI_k SNOId_k_theta SNOId_k_phi];
end;
% ---------------------------------------------- %
    
% ----------------- Noise data? ---------------- %
noise_string = input('Insert noise? ([y]/n): ','s');
if (isempty(noise_string) | noise_string == 'y')
    noise_mean = input(noise_mean_string);
    if isempty(noise_mean)
        noise_mean = def_noise_mean;
    end;
    noise_var = input(noise_var_string);
    if isempty(noise_var)
        noise_var = def_noise_var;
    end;
    noise = [noise_mean noise_var];
else
    noise = [];
end;
% ---------------------------------------------- %
    
% ----------------- Signal type? --------------- %
type_n = [];
if isempty(type_n)
    type_n = def_type_n;
end;
switch type_n
case 1
    type = 'sinusoid';
    def_NN = 100;
case 2
    type = 'bpsk';
    def_NN = 1;
otherwise
    error(err_1);
end;
% ---------------------------------------------- %
    
% ------------ Number of data samples? --------- %
nn = input(nn_string);
if isempty(nn)
    nn = def_nn;
end;
% ---------------------------------------------- %

% ------- Number of samples per symbol? -------- %
NN = [];
if isempty(NN)
   NN = def_NN;
end;
% ---------------------------------------------- %

% ---------------- Mu for LMS? ----------------- %
Mu = input(Mu_string);
if isempty(Mu)
    Mu = def_Mu;
end;
% ---------------------------------------------- %

% ---------------- Nulls depth? ---------------- %
if q==0
    AF_thresh = def_AF_thresh;
else
    AF_thresh = input(AF_thresh_string);
    if isempty(AF_thresh)
        AF_thresh = def_AF_thresh;
    end;
end;
% ---------------------------------------------- %

% -------------- Element pattern? -------------- %
E_pattern_file = [];
if isempty(E_pattern_file)
    E_pattern_file = def_E_pattern_file;
end;
E_pattern = load(E_pattern_file)';
% ---------------------------------------------- %

warning off;

end

function AF=planar_plot_pattern(N,d,w,E_pattern,AF_limit,METHOD,PHI_SCAN)
%************************************************************************
% planar_plot_pattern(N,d,w,E_pattern,AF_limit,METHOD,PHI_SCAN)
%************************************************************************
%	PLOT_PATTERN is a MATLAB function that plots 3-D patterns in
%	rectangular coordinates x axis represents phi and y reresents theta
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%--------------------------------------------------------------------------
%---
H_FIG = 1;

k0 = 2*pi;

switch METHOD
case 'linear'
    %%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
    m     = linspace(0,N-1,N);
    theta = linspace(-pi/2,pi/2,181);
    %-------------------------------------------------------%
    
    %%%%%%%%%%%%%% Total Pattern [Element * AF] %%%%%%%%%%%%%
    AF = E_pattern .* sum(diag(w)*exp(i*(k0*d*m'*sin(theta))));
    AF = 20*log10(abs(AF)./max(abs(AF)));
    %-------------------------------------------------------%
    
    figure(H_FIG);
    semipolar(theta,AF,AF_limit,0,4,'-'); grid on;
    
    AF(find(AF <= AF_limit)) = AF_limit;

    H_FIG = H_FIG + 1; figure(H_FIG);
    plot(theta*180/pi,AF); axis([-90 90 AF_limit 0]); grid on;
    xlabel('\theta (degrees)','fontsize',14);
    ylabel('Magnitude (dB)','fontsize',14);
    fs  = get(gca, 'FontSize');
    set(gca, 'FontSize', 14);
    
case 'planar'
    %%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
    theta = linspace(0,pi/2,361);
    phi   = linspace(0,2*pi,361);
    AF    = zeros(length(theta),length(phi));
    k     = 0;
    %-------------------------------------------------------%
    
    [THETA,PHI] = meshgrid(theta,phi);
    
    PSI_x = k0*d(1)*sin(THETA).*cos(PHI);
    PSI_y = k0*d(2)*sin(THETA).*sin(PHI);
    
    for m = 1:N(1),
        for n = 1:N(2),
            k = k + 1;
            PSI = (m-1)*PSI_x + (n-1)*PSI_y;
            AF  = AF + w(k) * exp(i*PSI);
        end;
    end;
    
    %%%%%%%%%%%%%% Total Pattern [Element * AF] %%%%%%%%%%%%%
    AF = E_pattern .* AF;
    
    AF = 20*log10(abs(AF)./max(max(abs(AF)))); AF = AF';
    %-------------------------------------------------------%

    H_FIG = H_FIG + 1; figure(H_FIG);                % H_FIG = figure(2)
    contourf(180/pi*phi,180/pi*theta,AF,50); hold on; shading('faceted');
    ylabel('Elevation');
    xlabel('Azimuth');
    caxis([AF_limit 0])
    colorbar;
    grid on;
    fs  = get(gca, 'FontSize');
    set(gca, 'FontSize', 14);
 
    if nargin == 7
        H_FIG = H_FIG + 1;
        plot_cut(AF,AF_limit,PHI_SCAN,H_FIG);
    end;

end;
end


function [S_SOI, s, m] = planar_sig_gen(p,d,nn,NN,type,sig,noise,E_pattern)
%************************************************************************
%	[S_SOI, s, m] = planar_sig_gen(p,d,nn,NN,type,sig,noise,E_pattern)
%************************************************************************
%	SIG_GEN is a MATLAB function that generates the SOI and SNOIs
%   signals (with and without Gaussian noise)
%
%	Input Parameters Description
%	----------------------------
%	- p          number of elements in x axis and y axis
%	- d          inter-element spacing (in wavelength)
%                default value is (0.5,0.5) (Nyquist rate)
%   - nn         number of samples
%   - NN         number of samples per cycle or symbol
%   - type       option: 'sinusoid' or 'bpsk'
%   - sig        signals amplitudes and directions of SOI and SNOI
%   - noise      amplitude and variance values
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%
%	Output Parameters Description
%	-----------------------------
%	- S_SOI      SOI reference signal (column vector)
%   - s          matrix containing all signals including noise
%   - m          frequency multiplier for sinusoidal signals
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%
%--------------------------------------------------------------------------
%---
%%%%%%%%%%%%%%% Parameters initialization %%%%%%%%%%%%%%%
k0 = 2*pi;                                                 % k in free space
n_sig = size(sig,1);                                       % Total number of signals
m = [];

x = linspace(0,p(1)-1,p(1));
y = linspace(0,p(2)-1,p(2));

sigr = [sig(:,1) pi/180*sig(:,2:3)];                    % degrees to radians conversion

[X,Y] = meshgrid(x,y);

err_1 = sprintf('\nSignal type not supported...');
%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%% Generating SOI %%%%%%%%%%%%%%%%%%%%
PSI_SOI = -k0*d(1,1)*sin(sigr(1,2))*cos(sigr(1,3))*X - k0*d(1,2)*sin(sigr(1,2))*sin(sigr(1,3))*Y;

PSI_SOI = reshape(PSI_SOI,1,p(1)*p(2));

switch type
   case 'bpsk'
      STATE1 = sum(100*clock);
      rand('state',STATE1);
      s_rand = round(rand(nn,1));
      s_rand_index = find(s_rand == 0);
      s_rand(s_rand_index) = -1;
      s_NN = ones(1,NN);
      s_rand = (s_rand*s_NN)';
      s_rand = reshape(s_rand,size(s_rand,1)*size(s_rand,2),1);
      [s_rand,PSI_SOI] = meshgrid(s_rand,PSI_SOI);
      s = sigr(1,1)*(s_rand.*exp(i*PSI_SOI));
      S_SOI  = s(1,:);
   otherwise
      error(err_1);
end;

%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%% Generating SNOI %%%%%%%%%%%%%%%%%%%
for k = 2 : n_sig,
   PSI_SNOI = -k0*d(1,1)*sin(sigr(k,2))*cos(sigr(k,3))*X - k0*d(1,2)*sin(sigr(k,2))*sin(sigr(k,3))*Y;

   PSI_SNOI = reshape(PSI_SNOI,1,p(1)*p(2));

   switch type
      case 'bpsk'
         STATE2 = sum(k*100*clock);
         rand('state',STATE2);
         s_rand = round(rand(nn,1));
         s_rand_index = find(s_rand == 0);
         s_rand(s_rand_index) = -1;
         s_NN = ones(1,NN);
         s_rand = (s_rand*s_NN)';
         s_rand = reshape(s_rand,size(s_rand,1)*size(s_rand,2),1);
         [s_rand,PSI_SNOI] = meshgrid(s_rand,PSI_SNOI);
         SNOI = sigr(k,1)*(s_rand.*exp(i*PSI_SNOI));
         s = s + SNOI;
      otherwise
         error(err_1);
   end;
end;
%-------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%% Noise %%%%%%%%%%%%%%%%%%%%
if ~isempty(noise)
    n_ant = p(1) * p(2);

    for k = 1 : n_ant
        STATE3 = sum(rand(1)*100*clock);
        randn('state',STATE3);
        noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        STATE4 = sum(rand(1)*100*clock);
        randn('state',STATE4);
        noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        noise_data = complex(noise_data_real,noise_data_imag);
        s(k,:) = s(k,:) + noise_data;
    end;
end;
%-------------------------------------------------------%
end



function polar_plot(p,d,w,AF_limit,E_plane)
%**************************************************************************
% polar_plot.m
%*************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%--------------------------------------------------------------------------
%---
% polar_plot(p,d,w)

k = 0;
k0 = 2*pi;

theta = linspace(0,pi/2,361);
phi   = linspace(0,2*pi,361);

[THETA,PHI] = meshgrid(theta,phi);

PSI_x = -k0*d(1)*sin(THETA).*cos(PHI);
PSI_y = k0*d(2)*sin(THETA).*sin(PHI);

AF = zeros(length(theta),length(phi));
for m = 1:p(1),
    for n = 1:p(2),
	k = k + 1;
        PSI = (m-1)*PSI_x + (n-1)*PSI_y;
        AF = AF + w(k) * exp(i*PSI);
    end;
end;

AF = E_plane .* AF;

AF = 20*log10(abs(AF)./max(max(abs(AF))));

AF(find(AF<=AF_limit)) = AF_limit;

figure(1);
[Y,X]=pol2cart(PHI,THETA);
surf(X,Y,AF(1:1:361,1:1:361));
shading('interp'); axis off; %colorbar;
view(-45,75);

% Draw x, y, and z axes
OX = X(181,181); OY = Y(271,271); s = 0.5;
set(line([OX;OX],[OY;min(min(Y))-s],[min(min(AF));min(min(AF))]),'Color','k'); %PHI=0
set(line([OX;max(max(X))+s],[OY;OY],[min(min(AF));min(min(AF))]),'Color','k'); %PHI=90
set(line([OX;OX],[OY;max(max(Y))+s],[min(min(AF));min(min(AF))]),'Color','k'); %PHI=180
set(line([OX;min(min(X))-s],[OY;OY],[min(min(AF));min(min(AF))]),'Color','k'); %PHI=270

% Label x, y, and z axes
text(min(min(OX)),min(min(Y))-s,min(min(AF)),'\phi=0\circ'  ,'FontSize',14,'Color','k'); %PHI=0
text(max(max(X))+s,min(min(OY)),min(min(AF)),'\phi=90\circ' ,'FontSize',14,'Color','k'); %PHI=90 
text(min(min(OX)),max(max(Y))+s,min(min(AF)),'\phi=180\circ','FontSize',14,'Color','k'); %PHI=180
text(min(min(X))-s,min(min(OY)),min(min(AF)),'\phi=270\circ','FontSize',14,'Color','k'); %PHI=270 

axis normal
end


