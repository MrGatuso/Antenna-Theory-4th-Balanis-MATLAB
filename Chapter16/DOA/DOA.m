function DOA
%*********************************************************************
%    DOA
%*********************************************************************
%    THIS IS A MATLAB BASED PROGRAM THAT ESTIMATES THE DIRECTION OF
%    ARRIVAL USING:
%
%     I.   LINEAR ARRAY
%     II.  RECTANGULAR PLANAR ARRAY
%
%     THE LINEAR ARRAY HAS M ELEMENTS PLACED EQUIDISTANTLY ALONG
%     THE X-AXIS.
%
%     THE RECTANGULAR PLANAR ARRAY HAS M x N ELEMENTS PLACED EQUIDISTANTLY
%     ALONG THE X AND Y AXES. 
%     
%     OPTION I.  LINEAR ARRAY
%
%      ** INPUT PARAMETERS:
%
%         1. NUMBER OF ELEMENTS
%         2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%         3. AMPLITUDE AND DIRECTION (THETA)OF SIGNAL-OF-INTEREST (SNOI)
%         4. AMPLITUDE AND DIRECTION (THETA)OF SIGNALS-NOT-OF-INTEREST
%            (SNOIs)
%         5. NOISE MEAN AND VARIANCE
%
%      ** PROGRAM OUTPUT:
%         1. DOA ESTIMATION  
%
%  OPTION II.  RECTANGULAR PLANAR ARRAY
%
%      ** INPUT PARAMETERS:
%
%         1. NUMBER OF ARRAY ELEMENTS IN X-DIRECTION
%         2. NUMBER OF ARRAY ELEMENTS IN Y-DIRECTION
%         3. SPACING BETWEEN THE ELEMENTS IN X-DIRECTION (IN WAVELENGTHS)
%         4. SPACING BETWEEN THE ELEMENTS IN Y-DIRECTION (IN WAVELENGTHS)
%         5. AMPLITUDE AND DIRECTION (THETA AND PHI) OF SIGNAL-OF-INTEREST
%           (SOI)
%         6. AMPLITUDE AND DIRECTION (THETA AND PHI) OF SIGNALS-NOT-OF-
%            -INTEREST(SNOIs)
%         7. NOISE MEAN AND VARIANCE
%
%      ** PROGRAM OUTPUT:
%         1. DOA ESTIMATION  
%
%     NOTES: ALL THE INPUT PARAMETERS ARE IN TERMS OF THE WAVELENGTH.
%***********************************************************************
%     Written by: Jeff Foutz, Arizona State University, 2/5/2003
%     Revised by: Zhiyong Huang, Arizona State University,  11/14/2004
%***********************************************************************

% Clean-up variables, figures
clear all;
close all;
clc
fprintf('\n *** NOTICE: PLEASE CONTACT DR. C.A. BALANIS IF YOU HAVE ANY QUESTION OR COMMENT. THANKS.\n\n');

%---Choice of array, linear or planar---------------------------------------------
fprintf('Antenna structure option: \n\tOption (1): Linear\n\tOption (2): Planar \n');
ERR = 1;
while(ERR ~= 0)
   array = str2num(input('\n Array structure = ','s'));
   if(array == 1)
       fprintf('\n Linear Array is chosen.\n');
       ERR = 0;
       linear;
   elseif(array == 2)
       fprintf('\n Planar Array is chosen.\n');
       ERR = 0;
       planar;
   else
      fprintf('\n Antenna structure number should be either 1 or 2\n');
   end
end
end

function [theta,phi]=est_doa(X,M,N,dx,dy,num_sig)

%*************************************************************************
%  This function computes the directions of arrival of all signals
%  given the array data, the number of elements, element spacing and
%  number of signals. The function works for a uniform planar or uniform linear array.  
%  If the array is planar, the 2-D Unitary Estimation of Signal Parameters via
%  Rotational Invariance Techniques(ESPRIT) algorithm is used.  If the array 
%  is linear, the total least squares(TLS) ESPRIT is used.

%  Input Variables
%  X        -  MNxL data matrix where L is the number of samples per element.
%              Each row of X contains samples from 1 element.
%  M        -  Number of elememts in x direction.
%  N        -  Number of elements in y direction(for linear array either M=1 or N=1.)
%  dx       -  Spacing between elements in wavelengths in x direction.
%              For linear array, dx is used as interelement spacing and dy is ignored.
%  dy       -  Spacing between elements in wavelengths in y direction.
%  num_sig  -  Number of signals present.

%  Output variables
%  theta    -  Vector containing azimuth direcions of arrival for all signals in radians.
%  phi      -  Vector containing elevation directions of arrival for all signals in radians.
%              If array is linear then phi will be an empty variable.
%*************************************************************************


if N~=1 & M~=1  %if N=1 or M=1 then the array is linear, if not, then use Unitary ESPRIT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations for Unitary ESPRIT for planar array
% See Zoltowski, Haardt and Mathews'  paper(Feb. 1996 IEEE 
% Trans. on Sig. Proc.) for a derivation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set up constants used in Unitary ESPRIT
if rem(N,2)==0  %check if N is even
QN=[eye(N/2) j*eye(N/2) ; fliplr(eye(N/2)) -j*fliplr(eye(N/2))];
QNM1=[eye(N/2-1) zeros(N/2-1,1) j*eye(N/2-1) ; zeros(N/2-1,1)' sqrt(2) zeros(N/2-1,1)' ; fliplr(eye(N/2-1)) zeros(N/2-1,1) -j*fliplr(eye(N/2-1)) ];
else
  QN=[eye((N-1)/2) zeros((N-1)/2,1) j*eye((N-1)/2) ; zeros((N-1)/2,1)' sqrt(2) zeros((N-1)/2,1)' ; fliplr(eye((N-1)/2)) zeros((N-1)/2,1) -j*fliplr(eye((N-1)/2)) ];
  QNM1=[eye((N-1)/2) j*eye((N-1)/2) ; fliplr(eye((N-1)/2)) -j*fliplr(eye((N-1)/2))];
end

if rem(M,2)==0  %check if M is even
QM=[eye(M/2) j*eye(M/2) ; fliplr(eye(M/2)) -j*fliplr(eye(M/2))];
QMM1=[eye(M/2-1) zeros(M/2-1,1) j*eye(M/2-1) ; zeros(M/2-1,1)' sqrt(2) zeros(M/2-1,1)' ; fliplr(eye(M/2-1)) zeros(M/2-1,1) -j*fliplr(eye(M/2-1)) ];
else
QM=[eye((M-1)/2) zeros((M-1)/2,1) j*eye((M-1)/2) ; zeros((M-1)/2,1)' sqrt(2) zeros((M-1)/2,1)' ; fliplr(eye((M-1)/2)) zeros((M-1)/2,1) -j*fliplr(eye((M-1)/2)) ];
  QMM1=[eye((M-1)/2) j*eye((M-1)/2) ; fliplr(eye((M-1)/2)) -j*fliplr(eye((M-1)/2))];
end

%Transform array data
Y=kron(QM',QN')*X;

%Compute a basis for the signal subspace
[U,D,S]=svd([real(Y) imag(Y)]);
Es=U(:,1:num_sig);

%Set up more constants
J2M=eye(M);
J2M=J2M(2:M,:);
J4N=eye(N);
J4N=J4N(2:N,:);

K1=real(QMM1'*J2M*QM);
K2=imag(QMM1'*J2M*QM);
Ku1=kron(K1,eye(N));
Ku2=kron(K2,eye(N));
K3=real(QNM1'*J4N*QN);
K4=imag(QNM1'*J4N*QN);
Kv1=kron(eye(M),K3);
Kv2=kron(eye(M),K4);

PSIu=Ku1*Es\Ku2*Es;
PSIv=Kv1*Es\Kv2*Es;

lambda=eig(PSIu+j*PSIv);
u=2*atan(-real(lambda));
v=2*atan(-imag(lambda));

%Compute azimuth and elevation angles of arrival
phi=atan2(v/dy,u/dx);
theta=asin(u./(2*pi*dx*cos(phi)));

%adjust phi if necessary so that 0<=phi<2*pi
for n=1:num_sig
  if phi(n)<0 
   phi(n)=phi(n)+2*pi;
 end
 if phi(n)>=2*pi 
   phi(n)=phi(n)-2*pi;
 end
end

elseif M==1 & N==1
 disp('Array must have more than 1 element')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Equations for Total Least Squares(TLS) ESPRIT
%   See Roy and Kailath's paper(July 1989 IEEE Trans. on
%   Acoustics, SPeech and Sig. Proc.) for a derivation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else  %IF the array is a uniform linear array, then use the TLS ESPRIT

  R=X*X';                 %Compute an estimate of the spatial covariance matrix
  [U,D,V]=svd(R);         %Do eigendecompositon of R

  S=U(:,1:num_sig);       %Compute a basis for the signal subspace

  %Compute the angles of arrival
  theta = S(1:max(N,M)-1,:)\S(2:max(N,M),:);
  w=angle(eig(theta));
  theta=asin(-w/(dx*2*pi));
  phi=[];  %elevation angle is not used in linear array

end
end

function linear

%*************************************************************************
%	Linear.m
%************************************************************************
%	It is a MATLAB function that estimate the DOA of linear array. 
%************************************************************************
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
[N,d,sig,noise,type,nn,NN,E_pattern] = linear_data_entry;

% Generate signals and estimate DOA.
DOAs = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern);

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
disp('=== DOA Estimations ===');
string = sprintf('        Theta(degrees)'); disp(string);
fprintf('SOI       %3.3f\n', 180/pi*DOAs(1));
for i=2:rowofsig
    fprintf('SNOI%d     %3.3f\n',i-1,180/pi*DOAs(i));
end

disp('*************************************************************************');
% Stop recording

if (save_in_file=='y')
    diary off;
    fprintf('Data saved in file %s.txt\n\n', file_string);
end
end



function [N,d,sig,noise,type,nn,NN,E_pattern] = linear_data_entry
%*************************************************************************
%	[N,d,sig,noise,type,nn,NN,E_pattern] = linear_data_entry
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
     
% -------------- Element pattern? -------------- %
E_pattern_file = [];
if isempty(E_pattern_file)
   E_pattern_file = def_E_pattern_file;
end;
E_pattern = load(E_pattern_file)';
% ---------------------------------------------- %

warning off;
end


function DOAs = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern)
%**********************************************************************
%	DOAs = linear_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern)
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
%   -DOAs        Estimations of DOAs
%************************************************************************
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

if ~isempty(noise)
    s_noise = s;

    for k = 1 : N
        STATE3 = sum(rand(1)*100*clock);
        randn('state',STATE3);
        noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        STATE4 = sum(rand(1)*100*clock);
        randn('state',STATE4);
        noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        noise_data = complex(noise_data_real,noise_data_imag);
        s_noise(k,:) = s_noise(k,:) + noise_data;
    end;
    DOA_theta =est_doa(s_noise,N,1,d,1,1);
else
    DOA_theta =est_doa(s,N,1,d,1,1);
end;    
DOAs = DOA_theta;
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
         if ~isempty(noise)
             s_noise = SNOI;

             for k = 1 : N
                 STATE3 = sum(rand(1)*100*clock);
                 randn('state',STATE3);
                 noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
                 STATE4 = sum(rand(1)*100*clock);
                 randn('state',STATE4);
                 noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
                 noise_data = complex(noise_data_real,noise_data_imag);
                 s_noise(k,:) = s_noise(k,:) + noise_data;
             end;
             DOA_theta = est_doa(s_noise,N,1,d,1,1);
         else
             DOA_theta = est_doa(SNOI,N,1,d,1,1);
         end;    
         DOAs = [DOAs; DOA_theta];
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

function planar


%************************************************************************
%	Planar.m
%************************************************************************
%   It is the Matlab function that estimate DOAs of planar array.
%************************************************************************
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
[N,d,sig,noise,type,nn,NN,E_pattern] = planar_data_entry;

% Generate signals and estimate the DOAs
DOAs = planar_sig_gen(N,d,nn,NN,type,sig,noise,E_pattern);

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
disp('=== DOA Estimations ===');
disp('        Theta(degrees)  Phi(degrees)');
fprintf('SOI       %3.3f            %3.3f\n', 180/pi*DOAs(1,1),180/pi*DOAs(1,2));
for i=2:rowofsig
    fprintf('SNOI%d     %3.3f            %3.3f\n',i-1,180/pi*DOAs(i,1),180/pi*DOAs(i,2));
end
disp('*************************************************************************');

% Stop recording
if (save_in_file=='y')
    diary off;
    fprintf('Data saved in file %s.txt\n\n', file_string);
end
end

function [N,d,sig,noise,type,nn,NN,E_pattern] = planar_data_entry
%************************************************************************
%	[N,d,sig,noise,type,nn,NN,E_pattern] = planar_data_entry
%***********************************************************************
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
%	- E_pattern  samples of element pattern (column vector)
%                if interested in array factor only, enter 1
%                default value is 1
%************************************************************************
%	Credits:
%		S. Bellofiore
%       Zhiyong Huang
%-----------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%% Default Values %%%%%%%%%%%%%%%%%%
def_N                   = [8 8];
def_d                   = [0.5 0.5];
def_SOI                 = 1;
def_q                   = 1;
def_SNOI                = 1;
def_noise_mean          = 0;
def_noise_var           = 0.1;
def_type_n              = 2;
def_nn                  = 500;
def_E_pattern_file      = 'planar_isotropic.e';
%---------------------------------------------------%
   
%%%%%%%%%%%%% Strings initialization %%%%%%%%%%%%%%%%
N_string_x              = sprintf('Enter number of elements in the x-axis [%d]: ',def_N(1));
N_string_y              = sprintf('Enter number of elements in the y-axis [%d]: ',def_N(2));
d_string_x              = sprintf('Enter the spacing in the x-axis, dx (lambda) [%2.1f]: ',def_d(1));
d_string_y              = sprintf('Enter the spacing in the y-axis, dy (lambda) [%2.1f]: ',def_d(2));
SOI_string              = sprintf('Enter the Pilot signal (SOI) amplitude [%d]: ',def_SOI);
SOId_string_theta       = sprintf('Enter the Pilot signal (SOI) in theta direction (degrees between 0 and 90): ');
SOId_string_phi         = sprintf('Enter the Pilot singal (SOI) in phi direction (degrees between 0 and 360): ');
q_string                = sprintf('Enter number of interfering signals (SNOI) [%d]: ',def_q);
noise_mean_string       = sprintf('Enter the mean of the noise [%d]: ',def_noise_mean);
noise_var_string        = sprintf('Enter the variance of the noise [%2.1f]: ',def_noise_var);
type                    = sprintf('Type of signal:\n\t[1] sinusoid\n\t[2] BPSK\nEnter number [%d]: ',def_type_n);
nn_string               = sprintf('Enter the number of data samples [%d]: ',def_nn); 
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
    SNOId_k_string_theta    = sprintf('Enter the theta direction of No. %d Interference singal (SNOI_%d) (degrees between 0 and 90): ',k,k);
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

% -------------- Element pattern? -------------- %
E_pattern_file = [];
if isempty(E_pattern_file)
    E_pattern_file = def_E_pattern_file;
end;
E_pattern = load(E_pattern_file)';
% ---------------------------------------------- %

warning off;

end


function DOAs = planar_sig_gen(p,d,nn,NN,type,sig,noise,E_pattern)
%************************************************************************
%	DOAs = planar_sig_gen(p,d,nn,NN,type,sig,noise,E_pattern)
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
%	-DOAs        Estimations of DOAs
%*************************************************************************
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

if ~isempty(noise)
    s_noise = s;
    n_ant = p(1) * p(2);

    for k = 1 : n_ant
        STATE3 = sum(rand(1)*100*clock);
        randn('state',STATE3);
        noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        STATE4 = sum(rand(1)*100*clock);
        randn('state',STATE4);
        noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
        noise_data = complex(noise_data_real,noise_data_imag);
        s_noise(k,:) = s_noise(k,:) + noise_data;
    end;
    [DOA_theta DOA_phi] = est_doa(s_noise,p(1),p(2),d(1),d(2),1);
else
    [DOA_theta DOA_phi] = est_doa(s,p(1),p(2),d(1),d(2),1);
end;    
DOAs = [DOA_theta DOA_phi];
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
         if ~isempty(noise)
             s_noise = SNOI;
             n_ant = p(1) * p(2);

             for k = 1 : n_ant
                 STATE3 = sum(rand(1)*100*clock);
                 randn('state',STATE3);
                 noise_data_real = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
                 STATE4 = sum(rand(1)*100*clock);
                 randn('state',STATE4);
                 noise_data_imag = noise(1,1) + sqrt(noise(1,2)/2)*randn(1,size(s,2));
                 noise_data = complex(noise_data_real,noise_data_imag);
                 s_noise(k,:) = s_noise(k,:) + noise_data;
             end;
             [DOA_theta DOA_phi] = est_doa(s_noise,p(1),p(2),d(1),d(2),1);
         else
             [DOA_theta DOA_phi] = est_doa(SNOI,p(1),p(2),d(1),d(2),1);
         end;    
         DOAs = [DOAs; [DOA_theta DOA_phi]];
         s = s + SNOI;
      otherwise
         error(err_1);
   end;
end;
%-------------------------------------------------------%
end





