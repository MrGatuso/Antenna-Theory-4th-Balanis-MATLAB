%***********************************************************************
% LOOP.m
%**********************************************************************
% This is a MATLAB based program computes the:
%	I.   Maximum directivity (dimensionless and in dB)
%	II.	 Radiation resistance
%   III. Directivity pattern (in dB) in polar form
%	IV.  Normalized elevation far-field pattern (E-phi, in dB)
%
%   and also plots the:
% 	I.   Directivity polar pattern (in dB)
%	II.  Normalized far-field polar amplitude pattern (E-phi, in dB) 
%
% for a circular loop (with constant current).  The loop is radiating
% into free space.
%
% The directivity and radiation resistance are calculated using
% the trailing edge method in increments of 1 degree in theta.
%
%	**Input parameters:
%	1.	A: Loop radius (in wavelengths)
%
%	**Note:
%	The far-zone electric field component E-phi exists for
%	0 < theta < 180 and 0 < phi < 360.
%-----------------------------------------------------------------
% Converted from Fortran to Matlab 3/2002 by Kelly O'Dell.
% Modified by Marios Gkatzianas
%-----------------------------------------------------------------
%
function LOOP
clear all;
close all;
clc

format long;

warning off;

%---Choice of output---

fprintf('Output device option \n\tOption (1): Screen\n\tOption (2): File \n');
ERR = 1;
while(ERR ~= 0)
   DEVICE = input('\nOutput device = ','s');
   DEVICE = str2num(DEVICE);
   if(DEVICE == 1)
      ERR = 0;
   elseif(DEVICE == 2)
      FILNAM = input('Input the desired output filename: ','s');
      ERR = 0;
   else
      error('Outputting device number should be either 1 or 2\n');
   end
end


%---Definition of constants and initialization---

PI = 4.0*atan(1.0);
E = 120.0*PI;
THETA = PI/180.0;
UMAX = 0.0;
PRAD = 0.0;
TOL = 1.0E-5;

%---Input the radius of the loop---

A = input('\nRadius of loop in wavelengths = ','s');
A = str2num(A);
%***Insert input data error loop*****

%---Main program------------------------------------

I = 1;
while(I <= 180)
   XI = I*PI/180.0;
   X = 2.0*PI*A*sin(XI);
   if(abs(X) < TOL)
      F = 0.0;
   else
      F = besselj(1,X);
   end
   U = A^2*(2.0*PI)^2/8.0*E*F^2;
   if(U > UMAX)
      UMAX = U;
   end
   UA = U*sin(XI)*THETA*2*PI;
   PRAD = PRAD+UA;
   I = I+1;
end
D = (4.0*PI*UMAX)/PRAD;
DDB = 10.0*log10(D);
RR = 2.0*PRAD;

%---Calculation of elevation far-field patterns in 1 degree increments---
fid = fopen('ElevPat.dat','w');
fprintf(fid,'\tLoop\n\n\tTheta\t\tH (dB)\n');
fprintf(fid,'\t-----\t\t------\n');
T = zeros(180,1);
HT = zeros(180,1);
HdB = zeros(180,1);

x = 1;
while(x<=180)
   T(x) = x-0.99;
   Y = 2*PI*A*sin(T(x)*THETA);
   HT(x) = besselj(1,Y);
   x = x+1;
end
HT = abs(HT);
HTmax = max(HT);
HdB = 20*log10(abs(HT)/HTmax);

x = 1;
while(x<=180)
   fprintf(fid,'\n %5.0f %12.4f',T(x),HdB(x));
   x = x+1;
end

fclose(fid);


%---Create output------------
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
   clc
end

x=linspace(0,4*pi*A,500);
dx=x(2)-x(1);
den=sum(besselj(2,x)*dx);
theta=linspace(0,pi,300);

Dth=4*pi*A/den*(besselj(1,2*pi*A*sin(theta))).^2;
Dth(1)=0;
Dth_db=10*log10(Dth);


%---Echo input parameters and output computed parameters---
fprintf(fid,'\nLOOP:\n-----');
fprintf(fid,'\n\nInput parameters:\n-----------------');
fprintf(fid,'\nRadius of loop in wavelengths = %6.4f',A);
fprintf(fid,'\n\nOutput parameters:\n------------------');
fprintf(fid,'\nDirectivity (dimensionless) = %6.4f',D);
fprintf(fid,'\nDirectivity (dB) \t= %6.4f',DDB);
fprintf(fid,'\nRadiation resistance (Ohms) = %9.4f',RR);
fprintf(fid,'\n\n***NOTE:\nThe normalized elevation pattern is stored\n');
fprintf(fid,'in an output file called ..........ElevPat.dat\n\n');

if(DEVICE == 2)
   fclose(fid);
end

%---Plot elevation far field pattern------
% plot(T,-abs(HdB),'b');
% axis([0 180 -60 0]);
% grid on;
% xlabel('Theta (degrees)');
% ylabel('Amplitude (dB)');
% legend(['a = ',num2str(A),' \lambda'],0);
% title('Loop Far-Field Elevation Pattern');

T=T'; HdB=HdB';

T=[T T+180];
Hdb=[HdB fliplr(HdB)];

% Figure 1
% ********
polar_dB(T,-abs(Hdb),-60,max(Hdb),4);
title('Elevation plane normalized amplitude pattern (dB)','fontsize',16);

% Figure 2
% ********
figure(2);
polar_dB([theta theta+pi]*180/pi,[Dth_db fliplr(Dth_db)],-60,max(Dth_db),4);
title('Elevation plane directivity pattern (dB)','fontsize',16);
end

%---End of program-----------------------------------------
%***********************************************************************
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

function hpol = polar_dB(theta,rho,rmin,rmax,rticks,line_style) 

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
