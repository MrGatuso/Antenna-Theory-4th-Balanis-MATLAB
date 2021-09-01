
% ****************************************************************************
% DIPOLE.m
%*************************************************************************
% This is a MATLAB based program that computes the: 
%	I.	 Maximum directivity (dimensionless and in dB)
%	II.	 Radiation resistance (Rr)
%   III. Input resistance (Rin)
%	IV.  Reactance relative to current maximum (Xm)
%   V.   Input reactance (Xin)
%	VI.  Normalized current distribution
%   VII. Directivity pattern (in dB) in polar form
%   VIII.Normalized far-field amplitude pattern (E-theta, in dB) in polar form
% for a symmetrical dipole of finite length.  The dipole is radiating
% in free space.
%
% The directivity, resistances and resistances are calculated using the trailing 
% edge method in increments of 1 degree in theta.
%
% **Input parameters
% 1.	L:	Dipole length (in wavelengths)
% 2.    a:  Dipole radius (in wavelengths)
%
% **Note:
% The far zone electrif field component, E-theta, exists for
% 0 < theta < 180 and 0 < phi < 360.
%------------------------------------------------------------------------
% Converted from Fortran to Matlab by Kelly O'Dell 3/2002
% Modified by Marios Gkatzianas
%------------------------------------------------------------------------
%
function []=dipole;

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
TOL = 1.0E-6;

%---Input the length of the dipole---

L = input('\nLength of dipole in wavelengths = ','s');
L = str2num(L);
%***Insert input data error loop***

r=input('Radius of dipole in wavelengths = ');


%---Main program----------------------

A = L*PI;
I = 1;
while(I <= 180)
   XI = I*PI/180.0;
   if(XI ~= PI)
      U = ((cos(A*cos(XI))-cos(A))/sin(XI))^2*(E/(8.0*PI^2));
      if(U > UMAX)
         UMAX = U;
      end
   end
   UA = U*sin(XI)*THETA*2.0*PI;
   PRAD = PRAD+UA;
   I = I+1;
end
D = (4.0*PI*UMAX)/PRAD;
DDB = 10.0*log10(D);
RR = 2.0*PRAD;
if(A ~= PI)
   RIN = RR/(sin(A))^2;
end


%---Calculation of elevation far-field patterns in 1 degree increments---
fid = fopen('ElevPat.dat','w');
fprintf(fid,'\tDipole\n\n\tTheta\t\tE (dB)\n');
fprintf(fid,'\t----\t\t------');
T = zeros(180,1);
ET = zeros(180,1);
EdB = zeros(180,1);

x = 1;
while(x<=180)
   T(x) = x-0.99;
   ET(x) = (cos(PI*L*cos(T(x)*THETA))-cos(PI*L))/sin(T(x)*THETA);
   x = x+1;
end
ET = abs(ET);
ETmax = max(abs(ET));
EdB = 20*log10(abs(ET)/ETmax);

x = 1;
while(x<=180)
   fprintf(fid,'\n %5.4f %12.4f',T(x),EdB(x));
   x = x+1;
end

fclose(fid);


n=120*pi;
k=2*pi;

if exist('cosint')~=2,
   disp(' ');
   disp('Symbolic toolbox is not installed. Switching to numerical computation of sine and cosine integrals.');
   Xm=30*(2*si(k*L)+cos(k*L)*(2*si(k*L)-si(2*k*L))- ...
          sin(k*L)*(2*ci(k*L)-ci(2*k*L)-ci(2*k*r^2/L)));
   Xin=Xm/(sin(k*L/2))^2;
      
elseif exist('cosint')==2,
   Xm=30*(2*sinint(k*L)+cos(k*L)*(2*sinint(k*L)-sinint(2*k*L))- ...
          sin(k*L)*(2*cosint(k*L)-cosint(2*k*L)-cosint(2*k*r^2/L)));
   Xin=Xm/(sin(k*L/2))^2;
end;   



%---Create output------------
if(DEVICE == 2)
   fid = fopen(FILNAM,'w');
else
   fid = DEVICE;
end
clc
%---Echo input parameters and output computed parameters---
fprintf(fid,'\nDIPOLE:\n-------');
fprintf(fid,'\n\nInput parameters:\n-----------------');
fprintf(fid,'\nLength of dipole in wavelengths = %6.4f',L);
fprintf(fid,'\nRadius of dipole in wavelengths = %6.7f',r);
fprintf(fid,'\n\nOutput parameters:\n------------------');
fprintf(fid,'\nDirectivity (dimensionless) = %6.4f',D);
fprintf(fid,'\nDirectivity (dB) \t= %6.4f\n',DDB);
fprintf(fid,'\nRadiation resistance based on current maximum (Ohms) = %10.4f',RR);
fprintf(fid,'\nReactance based on current maximum (Ohms) = %10.4f\n',Xm);

if(abs(sin(A)) < TOL)
   fprintf(fid,'\nInput resistance = INFINITY');
   fprintf(fid,'\nInput reactance = INFINITY\n\n');   
else
   fprintf(fid,'\nInput resistance (Ohms) = %10.4f',RIN);
 % fprintf(fid,'\nInput reactance based on current maximum (Ohms) = %10.4f',Xm);
   fprintf(fid,'\nInput reactance (Ohms) = %10.4f',Xin);
   fprintf(fid,'\n\n***NOTE:\nThe normalized elevation pattern is stored\n');
   fprintf(fid,'in an output file called ..........ElevPat.dat\n\n');
end


if(DEVICE == 2)
   fclose(fid);
end

%---Plot elevation far field pattern------
% plot(T,EdB,'b');
% axis([0 180 -60 0]);
% grid on;
% xlabel('Theta (degrees)');
% ylabel('Amplitude (dB)');
% legend(['L = ',num2str(L),' \lambda'],0);
% title('Dipole Far-Field Elevation Pattern');


% Figure 1
% ********
z=linspace(-L/2,L/2,500);
k=2*pi;
I=sin(k*(L/2-abs(z)));
plot(z,abs(I));
xlabel('z^{\prime}/\lambda','fontsize',12);
ylabel('Normalized current distribution','fontsize',12);


% Figure 2
% ********
figure(2);
T=T'; EdB=EdB';
EdB=[EdB fliplr(EdB)];
T=[T T+180];

polar_dB(T,EdB,-60,0,4);
title('Elevation plane normalized amplitude pattern (dB)','fontsize',16);

% Figure 3
% ********
figure(3);
theta=linspace(0,2*pi,300);

Eth=(cos(k*L/2*cos(theta))-cos(k*L/2))./sin(theta);
Dth=4*pi*120*pi/(8*pi^2)*Eth.^2/PRAD;
Dth_db=10*log10(Dth);
Dth_db(Dth_db<=-60)=-60;

polar_dB(theta*180/pi,Dth_db,-60,max(Dth_db),4);
title('Elevation plane directivity pattern (dB)','fontsize',16);

%---End program----------------------------------------------
end

function [y]=si(x);

v=linspace(0,x/pi,500);
dv=v(2)-v(1);

y=pi*sum(sinc(v)*dv);
end

function [y]=ci(x);

v=linspace(0,x/(2*pi),500);
dv=v(2)-v(1);

y1=2*pi*sum(sinc(v).*sin(pi*v)*dv);

y=.5772+log(x)-y1;
end
%**************************************************************************
%       polar_dB(theta,rho,rmin,rmax,rticks,line_style) 
%**************************************************************************
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

function y=sinc(x)
i=find(x==0);                                                              
x(i)= 1;                             
y = sin(pi*x)./(pi*x);                                                     
y(i) = 1;   
end
