function BEVERAGE
%**********************************************************************************
% This program is about a long wire antenna. Antenna types considered are:
%           SINGLE wire in FREE SPACE;
%           SINGLE wire  above PEC ground: BEVERAGE antenna.
% There are two cases of current distribution considered: 
% a. Traveling wave antenna
% b. Resonant wave antenna
%--------------------------------------------------------------------------------------
% The formulation of the antenna is as follows:
%--------------------------------------------------------------------------------------
%
% CURRENT DISTRIBUTION:
% a. Traveling Wave: Iz(Z')= az Io exp(-jkzZ'), [-L/2,L/2]         Equation (10-1a)
% b. Standing Wave (Resonant):
%                       Io sin[kz(L/2-z')]      [0,+L/2]
%               Iz(Z')=                                             Equation (4-56)
%                       Io sin[kz(L/2+Z')]      [-L/2,0]
% GEOMETRY:
% 
% The coordinate system is follows:
%
% SINGLE WIRE:
% a. The z-axis is along the axis of the wire
% b. The x-y plane is perpendicular to the wire
% 
% BEVERAGE ANTENNA: Wire height above PEC is h.
% a. The z-axis is along the axis of the wire
% b. The x-y plane is perpendicular to the axis of the wire
% c. The x-axis is perpendicular to the interface formed by air and PEC
% d. The y-axis is parallel to the interface formed by air and PEC
%--------------------------------------------------------------------------------------
% FORMULATION:
% 
% SINGLE WIRE
% The formulation of this is given in Chapter 10, Section 10.2.1
%
% BEVERAGE:
% The formulation of this is that given in Chapter 10, Section 10.2.1 along with
% an Array factor based on Image Theory: 
% AF = 2jsin[kL sin(theta) cos(phi)]     
% E(total)=E(single wire) x Array Factor= Equation (10-2b) x 2jsin[kl sin(theta) cos(phi)]
%
% RADIATED POWER Prad:
% a. Single Wire:
%    Integrate radiation intensity U(theta,phi) over full space 
% b. Beverage(traveling, resonant): Integrate radiation intensity U(theta,phi) over half-space above PEC
%
% DIRECTIVITY:
% Do = 4 pi Umax(theta,phi)/Prad
%
% RADIATION RESISTANCE:
% Rr = 2 Prad/|Io|^2
%
% CHARACTERISTIC IMPEDANCE Zo:
% Zo = 60 cosh^-1(h/a); h = height of wire above PEC, a = radius of wire
%
% ----------------------------------------------------------------------------------
% THE PROGRAM CALCULATES THE FOLLOWING:
% a. Characteristic Impedance (only for Beverage)  
% b. Maximum Directivity (dimensionless and in dB)
% c. Radiation Resistance
% d. Amplitude pattern maxima in principal planes for:
%    1. THETA angles (for PHI=0, XZ-plane) 
%    2. PHI angles (for THETA=90,XY-plane)
% e. Amplitude pattern minima in principal planes for:
%    1. THETA angles (for PHI=0, XZ-plane)
%    2. PHI angle (for THETA=90,XY-plane)
% f. Plot normalized amplitude radiation patterns in principal planes for:
%    1. Elevation vs. THETA (for PHI=0,    XZ-plane)
%    2. Azimuth   vs. PHI   (for THETA=90, XY-plane)
%    NOTE: Normalization is done relatively to maximum in the two planes.
% --------------------------------------------------------------------------------- 
% LIMITS OF ANGLES THETA and PHI:
% For SINGLE WIRE:
% a. THETA (in degrees) [0,180]
%    Note: +90 degrees is along XY-plane (perpendicular to the wire)
% b. PHI (in degrees) [0, 360] 
% For BEVERAGE:
% a. THETA (in degrees) [0,180]
%    Note: +90 degrees is along XY-plane (perpendicular to the wire)
% b. PHI (in degrees) [-90, 90]
%    Note: 1. PHI = +90 degrees is along +Y-axis (parallel to PEC/air-PEC interface)
%          2. PHI = 0 degrees   is along  X-axis (perpendicular to PEC/air-PEC interface)
%          3. PHI = -90 degrees is along -Y-axis (parallel to PEC/air-PEC interface)
% ----------------------------------------------------------------------------------   
% INPUT PARAMETERS:
%  Length L of the wire (in wavelengths)
%  Radius a of the wire (in wavelengths) 
%  Height h of the wire above the PEC ground (in wavelengths)
%  Phase velocity ratio K = kz/k: Equation (10-3)
% OPTIONS:
% a. Antenna type:
%    1. Single wire
%    2. Beverage antenna
% b. Current distribution:
%    1. Traveling wave
%    2. Resonant
%---------------------------------------------------------------------------------------
% There are five additional m-files, m-functions, used by the program:
%   URF.m ~ radiation intensity for resonant case in FREE SPACE
%   UTF.m ~ radiation intensity for traveling wave in FREE SPACE
%   UR.m  ~ radiation intensity for resonant wave in the PRESENCE of PEC GROUND
%   UT.m  ~ radiation intensity for traveling wave in the PRESENCE of PEC GROUND
%   polar_dB.m ~ plotting subroutine
%**********************************************************************************

close all;
clear all;
clc

global h l K

%-------------------------------Choice of output-----------------------------------
fprintf('Output device option \n\tOption (1): Screen\n\tOption (2): File \n');
ERR = 1;
while(ERR ~= 0)
   DEVICE = str2num(input('\nOutput device = ','s'));
   if(DEVICE == 1)
      ERR = 0;
   elseif(DEVICE == 2)
      FILNAM = input('Input the desired output filename: ','s');
      ERR = 0;
   else
      fprintf('\nOutputting device number should be either 1 or 2\n');
   end
end

%---------------------Choice of geometry: free space or above PEC ground ----------------
fprintf('\nin Free space or above PEC ground \n\tOption (1): Free space\n\tOption (2): PEC ground \n');
ERR = 1;
while(ERR ~= 0)
   choice2 = str2num(input('\nOption = ','s'));
   if( choice2 == 1)
      ERR = 0;
   elseif(choice2 == 2)
      ERR = 0;
   else
      fprintf('\nchoose 1 or 2\n');
   end
end

%-------------------------Choice of current distribution----------------------------------
fprintf('\nResonant(1) or Traveling wave(2) \n');
ERR = 1;
while(ERR ~= 0)
   choice = str2num(input('\nOption = ','s'));
   if( choice == 1)
      ERR = 0;
   elseif(choice == 2)
      ERR = 0;
   else
      fprintf('\nchoose 1 or 2\n');
   end
end

%----------------------------Read input parameters ---------------------------------------
K = str2num(input('\nPhase velocity ratio kz/k = ','s'));
l = str2num(input('Length of the wire (in wavelengths) = ','s'));
if choice2 == 2,
   a = str2num(input('Radius of the wire (in wavelengths) = ','s'));
   h = str2num(input('Height above the PEC ground (in wavelengths) = ','s'));
   fprintf('Computations may take few minutes. Wait a litlle, please...\n\n');
end;

%********************************* Main part  *******************************'

TH0 = 90; 
PH0 = 0;

%--------------------------------- BEVERAGE  --------------------------------
if choice2 == 2, 
  N  = 377;
  t = zeros(N-1,1);
  p = zeros(N-1,1);
  U = zeros(N-1,N-1);
  dt = pi/(N-1);
  S = 0;
  Umax = 0;
  for n = 1:(N-1),
     t(n) = dt*(n-1);
     if ( (cos(t(n)) == K) | (cos(t(n)) == -K) ),
        t(n) = t(n) + 0.01*dt;
     end;
     for m = 1:(N-1),
         p(m) = 0.5*pi-dt*(m-1);
         if ( choice == 1 ), U(n,m) = UR(t(n),p(m)); end;
         if ( choice == 2 ), U(n,m) = UT(t(n),p(m)); end;        
         if Umax < U(n,m), Umax = U(n,m); end;
         S = S + U(n,m)*sin(t(n));
     end;
 end;
 S = S*dt*dt;
 if ( choice == 1 ), Rr = (120/pi*K^2)*S; end;
 if ( choice == 2 ), Rr = (120/pi    )*S; end;    
 D0 = 4*pi*Umax/S;
 D0dB =10*log10(D0);
 Z0 = 60*acosh(h/a);

 %***************************** THETA angles of maxima when  PHI=0 **********
 M = round( (0.5*pi-pi*PH0/180)/dt ) + 1;
 j=0;
 for n=2:(N-2),
     if ( U(n,M) > U(n-1,M)) & ( U(n,M) > U(n+1,M) ), 
        j=j+1;   
        thm(j) = t(n);
     end;
 end;
 %***************************** THETA angles of minima when PHI = 0  ************
 j=1;
 for n=2:(N-2),
     if ( U(n,M) < U(n-1,M)) & ( U(n,M) < U(n+1,M) ), 
        j=j+1;   
        thn(j) = t(n);
     end;
 end;
 thn(1)=0;
 thn(j+1)=pi;
 thm = 180*thm/pi;
 thn = 180*thn/pi;

 %******************** PHI angles of maxima when THETA=90 (Degree) **********
 MM = round( (pi*TH0/180)/dt )+1;
 j=0;
 for n=2:(N-2),
     if ( U(MM,n) > U(MM,n-1)) & ( U(MM,n) > U(MM,n+1) ), 
         j=j+1;   
         pm(j) = p(n);
      end;
 end;
 while  j == 0, 
     MM = MM + 1; 
     for n=2:(N-2),
         if ( U(MM,n) > U(MM,n-1)) & ( U(MM,n) > U(MM,n+1) ), 
            j=j+1;   
            pm(j) = p(n);
         end;
     end;
 end;
 %******************** PHI angles of minima when THETA=90 (Degree) ************
 j=1;
 for n=2:(N-2),
     if ( U(MM,n) < U(MM,n-1)) & ( U(MM,n) < U(MM,n+1) ), 
        j=j+1;   
        pn(j) = p(n);
     end;
 end;
 pn(1)=pi/2;
 pn(j+1)=-pi/2;
 pm = 180*pm/pi;
 pn = 180*pn/pi;

%***************** maxima in principle planes to normalize the plots **************
 U1max = max(U(:,M));
 U2max = max(U(MM,:));
 if U1max >= U2max, Umax = U1max; else, Umax = U2max; end;

%************************** plotting the pattern, PHI=0 ***************
 UU = U(:,M);
 UUmax = 10*log10( max(UU)/Umax );
 rmin = UUmax-60;
 rmax = UUmax;
 rticks = 4;
 line_style = '-';
 theta = 90-180*t/pi;
 figure(1)
 rho=UU;
 for i=1:length(UU),
      if UU(i) ~= 0,  rho(i)=10*log10(UU(i)/Umax); 
      else, rho(i)= UUmax-60;
      end;
 end;
 polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
 title(strvcat('\it{Elevation amplitude pattern in dB (XZ-plane: PHI=0^0)}','\it{normalized to the maximum in the two principal planes}'),'FontSize',12)
 clear UU;

%********************* plotting the pattern, THETA=90 (Degree) ***************
 UU = U(MM,:)';
 UUmax = 10*log10( max(UU)/Umax );
 rmin = UUmax-60;
 rmax = UUmax;
 rticks = 4;
 line_style = '-';
 theta = 180*p/pi;
 figure(2)
 rho=UU;
 for i=1:length(UU),
      if UU(i) ~= 0,  rho(i)=10*log10(UU(i)/Umax); 
      else
          rho(i)= UUmax-60;
      end;
 end;
 polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
 title(strvcat('\it{Azimuth amplitude pattern in dB (XY-plane: THETA=90^0)}','\it{normalized to the maximum in the two principal planes}'),'FontSize',12)
%---Create output------------
if(DEVICE == 2)
   fid = fopen(FILNAM,'wt');
else
   fid = DEVICE;
   clc
  
end

%---Echo input parameters and output computed parameters---
if( choice == 1 ),
fprintf(fid,'\nBEVEREGE ANTENNA above PEC ground, RESONANT case\n');
end;
if( choice == 2 ),
fprintf(fid,'\nBEVEREGE ANTENNA above PEC ground, TRAVELING WAVE case\n');
end;
fprintf(fid,'\nInput parameters:\n-----------------');
fprintf(fid,'\nLength of the wire (in wavelengths)      = %6.4f',l); 
fprintf(fid,'\nRadius of the wire (in wavelengths)      = %6.4f',a);  
fprintf(fid,'\nHeight above PEC ground (in wavelengths) = %6.4f',h);  
fprintf(fid,'\nPhase velocity ratio                     = %6.4f',K);  
fprintf(fid,'\n\nOutput parameters:\n------------------');
fprintf(fid,'\nCharacteristic impedance (ohms)          = %6.4f',Z0);
fprintf(fid,'\nRadiation resistance (ohms)              = %6.4f',Rr);
fprintf(fid,'\nDirectivity (dimensionless)              = %6.4f',D0);
fprintf(fid,'\nDirectivity (dB)                         = %6.4f',D0dB);
fprintf(fid,'\n');
fprintf(fid,'\nPattern maxima (THETA angles; PHI=0; XZ-plane)');
fprintf(fid,'\nMaxima (degrees)');
for j=1:length(thm), fprintf(fid,'\n = %6.4f',thm(j) ); end;
fprintf(fid,'\n');
fprintf(fid,'\nPattern minima (THETA angles; PHI=0; XZ-plane)');
fprintf(fid,'\nMinima (degrees)');
for j=1:length(thn), fprintf(fid,'\n = %6.4f',thn(j) ); end;
fprintf(fid,'\n');
fprintf(fid,'\nPattern maxima (PHI angles; THETA=90(degrees); XY-plane)');
fprintf(fid,'\nMaxima (degrees)');
for j=1:length(pm), fprintf(fid,'\n = %6.4f', pm(j) ); end;
fprintf(fid,'\n');
fprintf(fid,'\nPattern minima (PHI angles; THETA=90(degrees); XY-plane)');
fprintf(fid,'\nMinima (degrees)');
for j=1:length(pn), fprintf(fid,'\n = %6.4f', pn(j) ); end;
fprintf(fid,'\n\n');
fprintf(fid,'\n\n');
if(DEVICE == 2)
   fclose(fid);
end;
end;

%------------------------------- SINGLE WIRE  --------------------------------
if choice2 == 1, 
  N  = 3*377;
  t = zeros(N-1,1);
  U = zeros(N-1,1);
  dt = pi/(N-1);
  S = 0;
  Umax = 0;
  for n = 1:(N-2),
     t(n) = dt*(n-1);
     if ( (cos(t(n)) == K) | (cos(t(n)) == -K) ),
        t(n) = t(n) + 0.01*dt;
     end;
     if ( choice == 1 ), U(n) = URF(t(n)); end;
     if ( choice == 2 ), U(n) = UTF(t(n)); end;        
     if Umax < U(n), Umax = U(n); end; 
     S = S + U(n)*sin(t(n));
 end;
 S = S*dt;
 if ( choice == 1 ), Rr = (60*K^2)*S; end;
 if ( choice == 2 ), Rr = (60    )*S; end;    
 D0 = 4*pi*Umax/S/2/pi;
 D0dB =10*log10(D0);

 %***************************** THETA angles of maxima **********
 j=0;
 for n=2:(N-2),
     if ( U(n) > U(n-1)) & ( U(n) > U(n+1) ), 
        j=j+1;   
        thm(j) = t(n);
     end;
 end;
 %***************************** THETA angles of minima  ************
 j=1;
 for n=2:(N-2),
     if ( U(n) < U(n-1)) & ( U(n) < U(n+1) ), 
        j=j+1;   
        thn(j) = t(n);
     end;
 end;
 thn(1)=0;
 thn(j+1)=pi;
 thm = 180*thm/pi;
 thn = 180*thn/pi;

%************************** plotting the pattern, PHI=0 ***************

 UU = zeros(2*(N-1),1);
 for n = 1:(N-1),   UU(n) = U(n); end;
 for n = N:(2*N-2), UU(n) = U(2*N-1-n); end;
 UUmax = 10*log10( max(UU)/Umax );
 rmin = -60;
 rmax = 0;
 rticks = 4;
 line_style = '-';
 theta = zeros(2*(N-1),1);
 for n = 1:(N-1),   theta(n) = 90-180*t(n)/pi; end;
 for n = N:(2*N-2), theta(n) = 90-180*(2*pi-t(2*N-1-n ))/pi ; end;
 rho=UU;
 for i=1:length(UU),
      if UU(i) ~= 0,  rho(i)=10*log10(UU(i)/Umax); 
      else, rho(i)= -60;
      end;
 end;
 polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
 title('\it{Normalized elevation amplitude pattern in dB}','FontSize',14)

%---Create output------------
if(DEVICE == 2)
   fid = fopen(FILNAM,'wt');
else
   fid = DEVICE;
   clc
end

%---Echo input parameters and output computed parameters---
if( choice == 1 ),
fprintf(fid,'\nSINGLE WIRE in Free Space, RESONANT case\n');
end;
if( choice == 2 ),
fprintf(fid,'\nSINGLE WIRE in Free Space, TRAVELING WAVE case\n');
end;
fprintf(fid,'\nInput parameters:\n-----------------');
fprintf(fid,'\nLength of the wire (in wavelengths) = %6.4f',l); 
fprintf(fid,'\nPhase velocity ratio                = %6.4f',K);  
fprintf(fid,'\n\nOutput parameters:\n------------------');
fprintf(fid,'\nRadiation resistance (ohms)         = %6.4f',Rr);
fprintf(fid,'\nDirectivity (dimensionless)         = %6.4f',D0);
fprintf(fid,'\nDirectivity (dB)                    = %6.4f',D0dB);
fprintf(fid,'\n');
fprintf(fid,'\nPatern maxima (THETA angles)');
fprintf(fid,'\nMaxima(degrees)');
for j=1:length(thm), fprintf(fid,'\n = %6.4f', thm(j)); end;
fprintf(fid,'\n');
fprintf(fid,'\nPatern minima (THETA angles)');
fprintf(fid,'\nMinima(degrees)');
for j=1:length(thn), fprintf(fid,'\n = %6.4f',thn(j)); end;
fprintf(fid,'\n\n');
if(DEVICE == 2)
   fclose(fid);
end;
end;
end

function y = UR(TH,PHI)

global h l K

% y =  ( sin( 2*pi*h.*sin(TH).*cos(PHI) ).^2 ).* ...
%      (sin(TH).^2).*( cos(pi*l.*cos(TH)) - cos(pi*l*K )  ).^2./(K^2-cos(TH).^2).^2 ;
 
y=( sin( 2*pi*h.*sin(TH).*cos(PHI) ).^2 ).* ...
      (sin(TH).^2).*( cos(pi*l.*cos(TH)) - cos(pi*l*K )  ).^2./(K^2-cos(TH).^2).^2+...
  ( sin( 2*pi*h.*sin(pi-TH).*cos(PHI) ).^2 ).* ... 
     (sin(pi-TH).^2).*( sin( pi*l.*(cos(pi-TH)-K)) ).^2./(cos(pi-TH)-K).^2;
end
function y = URF(TH)

global l K

%y = (sin(TH).^2).*( cos(pi*l.*cos(TH)) - cos(pi*l*K )  ).^2./(K^2-cos(TH).^2).^2 ;

y = (sin(TH).^2).*( sin( pi*l.*(cos(TH)-K)) ).^2./(cos(TH)-K).^2+...
    (sin(pi-TH).^2).*( sin( pi*l.*(cos(pi-TH)-K)) ).^2./(cos(pi-TH)-K).^2;

end
function y = UT(TH,PHI)

global h l K

y =  ( sin( 2*pi*h.*sin(TH).*cos(PHI) ).^2 ).* ... 
     (sin(TH).^2).*( sin( pi*l.*(cos(TH)-K)) ).^2./(cos(TH)-K).^2;
end
function y = UTF(TH)

global l K

y = (sin(TH).^2).*( sin( pi*l.*(cos(TH)-K)) ).^2./(cos(TH)-K).^2;
end

function hpol = polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%************************************************************************
%       polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%***********************************************************************
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
%	Credits:
%		S. Bellofiore
%		S. Georgakopoulos
%		A. C. Polycarpou
%		C. Wangsvick
%		C. Bishop
%
%	Tabulate your data accordingly, and call polar_dB to provide the
%	2-D polar plot
%
%	Note:  This function is different from the polar.m (provided by
%	       MATLAB) because RHO is given in dB, and it can be negative
%-----------------------------------------------------------------------------

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
