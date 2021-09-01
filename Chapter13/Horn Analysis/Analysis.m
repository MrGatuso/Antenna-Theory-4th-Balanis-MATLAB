%*******************************************************************************
%      PYRAMIDAL HORN: ANALYSIS
%**************************************************************************
%
%     THIS PROGRAM COMPUTES FOR A PYRAMIDAL HORN THE:
%
%     I.    FAR-FIELD E- AND H-PLANE AMPLITUDE PATTERNS* BASED ON THE
%           THEORY OF SECTION 13.4 EQUATIONS (13-46) - (13-48c)
%     II.   DIRECTIVITY (IN dB) BASED ON EQUATION (13-52)
%     III.  DIRECTIVITY (IN dB) OF THE CORRESPONDING E-PLANE SECTORAL
%           HORN BASED ON EQUATION (13-19)
%     IV.   DIRECTIVITY (IN dB) OF THE CORRESPONDING H-PLANE SECTORAL
%           HORN BASED ON EQUATION (13-41)
%
%        ** INPUT PARAMETERS:
%
%       1.  RHO1 (IN WAVELENGTHS)
%       2.  RHO2 (IN WAVELENGTHS)
%       3.  WAVEGUIDE DIMENSIONS a & b (IN WAVELENGTHS)
%       4.  HORN APERTURE DIMENSIONS a1 & b1 (IN WAVELENGTHS)
%
%
%        ** NOTE: REFER TO FIGURE 13.18 FOR THE GEOMETRY.
%                 THE E- AND H-PLANE AMPLITUDE PATTERNS ARE STORED IN
%                 TWO DATA FILES NAMELY E-Plane Horn.dat AND H-Plane Horn.dat,
%                 RESPECTIVELY.
%**************************************************************************
%     Written by: Keith B. Washington, Arizona State University
%
%**************************************************************************

function []=horn

close all;
clear all
clc
% Input parameters

disp(' ')
disp('E-Plane and H-Plane Horn Specifications');
disp('--------------------------------------------')

R1=[]; R2=[];
while isempty(R1)|~isa(R1,'double')|R1<=0,
   R1 = input('rho1(in wavelengths) = ');
end;

while isempty(R2)|~isa(R2,'double')|R2<=0,
   R2 = input('rho2(in wavelengths) = ');
end;

disp(' ');
disp('Waveguide Dimensions');
disp('--------------------------------------------');

a=[]; b=[];
while isempty(a)|~isa(a,'double')|a<=0,
   a = input('a(in wavelengths) = ');
end;

while isempty(b)|~isa(b,'double')|b<=0,
   b = input('b(in wavelengths) = ');
end;  

disp(' ');
disp('Horn Aperture Dimensions')
disp('--------------------------------------------')

a1=[]; b1=[];
while isempty(a1)|~isa(a1,'double')|a1<=0,
   a1 = input('a1(in wavelengths) = ');
end;

while isempty(b1)|~isa(b1,'double')|b1<=0,
   b1 = input('b1(in wavelengths) = ');
end;

disp(' ')

opt=[];
while isempty(opt)|(opt~=1&opt~=2),
   opt=input(['Select output option\n','--------------------------------------------\n', ...
              '1. Screen\n','2. Output file\n','->']);
end;

if opt==2,
   fn = input('Input the Desired Filename:\n','s');
end;

disp(' ')
disp('please wait')


% compute internal functions

u = (1/sqrt(2))*((sqrt(R2)/a1)+(a1/sqrt(R2)));

v = (1/sqrt(2))*((sqrt(R2)/a1)-(a1/sqrt(R2)));

u = Fresnel(u);

v = Fresnel(v);

DH = 4*pi*b*R2/a1*((real(u)-real(v))^2 + (imag(u)-imag(v))^2);

w = Fresnel(b1/sqrt(2*R1));

DE = 64*a*R1/(pi*b1)*((real(w))^2 + (imag(w))^2);

DP = pi/(32*a*b)*DE*DH;

k = 2*pi;

Emax = 0;

Hmax = 0;


% E and H plane Outputs


% E-Plane Amplitude

fid1 = fopen('E-Plane Horn.dat','wt');

fprintf(fid1,'                         E-Plane Horn\n');
fprintf(fid1,'Theta     E-Theta\n');


for(theta = 0:0.5:360);

    I = theta*2 + 1;

    theta = theta*pi/180;

    phi = pi/2;

    ky = k*sin(theta);

    kxp = pi/a1;

    kxdp = -pi/a1;

    t1 = sqrt(1/(pi*k*R1))*(-k*b1/2-ky*R1);

    t2 = sqrt(1/(pi*k*R1))*(k*b1/2-ky*R1);

    t1p = sqrt(1/(pi*k*R2))*(-k*a1/2-pi/a1*R2);

    t2p = sqrt(1/(pi*k*R2))*(k*a1/2-pi/a1*R2);

    t1dp = -t2p;

    t2dp = -t1p;

    I1 = .5*sqrt(pi*R2/k)*(exp(j*R2/(2*k)*kxp^2)*(Fresnel(t2p) - Fresnel(t1p)) + exp(j*R2/(2*k)*kxdp^2)*(Fresnel(t2dp) - Fresnel(t1dp)));

    I2 = sqrt(pi*R1/k) * exp(j*R1/(2*k)*ky^2) * (Fresnel(t2) - Fresnel(t1));

    y(I) = (1 + cos(theta))*I1*I2;

    y(I) = abs(y(I));

end

for(I = 1:721)
    if(y(I) > Emax)
       Emax = y(I);
    end
end

for(I = 1:721)
    if(y(I) <= 0)
       Edb = -100;
    else
       Edb = 20*log10(abs(y(I))/Emax);
    end

     theta = (I-1)/2;
     x(I)=theta;
     q1(I)=Edb;
    fprintf(fid1,'%6.1f %12.6f\n',theta,Edb);

end

fclose(fid1);

%subplot(2,1,1)
%plot(x,q);
%xlabel('Theta (degrees)');
%ylabel('Field Pattern (dB)');
%title('E-Plane');



% H-Plane Amplitude

fid2 = fopen('H-Plane Horn.dat','wt');

fprintf(fid2,'                         H-Plane Horn\n');
fprintf(fid2,'Theta     E-Phi\n');

for(theta = 0:0.5:360);

    I = theta*2 + 1;

    theta = theta*pi/180;

    phi = 0;

    kxp = k*sin(theta) + pi/a1;

    kxdp = k*sin(theta) - pi/a1;

    t1 = sqrt(1/(pi*k*R1))*(-k*b1/2);

    t2 = sqrt(1/(pi*k*R1))*(k*b1/2);

    t1p = sqrt(1/(pi*k*R2))*(-k*a1/2-kxp*R2);

    t2p = sqrt(1/(pi*k*R2))*(k*a1/2-kxp*R2);

    t1dp = sqrt(1/(pi*k*R2))*(-k*a1/2-kxdp*R2);

    t2dp = sqrt(1/(pi*k*R2))*(k*a1/2-kxdp*R2);

    I1 = .5*sqrt(pi*R2/k)*(exp(j*R2/(2*k)*kxp^2)*(Fresnel(t2p) - Fresnel(t1p)) + exp(j*R2/(2*k)*kxdp^2)*(Fresnel(t2dp) - Fresnel(t1dp)));

    I2 = sqrt(pi*R1/k) * exp(j*R1/(2*k)*ky^2) * (Fresnel(t2) - Fresnel(t1));

    y(I) = (1 + cos(theta))*I1*I2;

    y(I) = abs(y(I));

end

for(I = 1:721)
    if(y(I) > Hmax)
       Hmax = y(I);
    end
end

for(I = 1:721)
    if(y(I) <= 0)
        Hdb = -100;
    else
        Hdb = 20*log10(abs(y(I))/Hmax);
     end

     theta = (I-1)/2;
     x(I)=theta;
     q2(I)=Hdb;
    fprintf(fid2,'%6.1f %12.6f\n',theta,Hdb);

end

fclose(fid2);

% Figure 1
% ********
ha=plot(x,q1); set(ha,'linestyle','-','linewidth',2);
hold on; hb=plot(x,q2,'r--'); set(hb,'linewidth',2);
xlabel('Theta (degrees)');
ylabel('Field Pattern (dB)');
title('Horn Analysis');
legend('E-Plane','H-Plane');
grid on;
axis([0 360 -60 0]);

% Figure 2
% *********
figure(2);
ht1=elevation(x*pi/180,q1,-60,0,4,'b-');
hold on;
ht2=elevation(x*pi/180,q2,-60,0,4,'r--');
set([ht1 ht2],'linewidth',2);

legend([ht1 ht2],{'E-plane','H-plane'});
title('Field patterns');


% Directivity Output

if exist('fn'),
   fid3 = [1 fopen(fn,'wt')];
else
   fid3=1;
end;

for fin=1:length(fid3),
   fprintf(fid3(fin),'*******************************************************************************\n');
   fprintf(fid3(fin),'PROGRAM OUTPUT\n');
   fprintf(fid3(fin),'*******************************************************************************\n\n');

   fprintf(fid3(fin),'\n');
   fprintf(fid3(fin),'Pyramidal Horn\n');
   fprintf(fid3(fin),'--------------\n');
   fprintf(fid3(fin),'Directivity = %.2f dB\n',10*log10(DP));
   fprintf(fid3(fin),'Directivity = %.2f dimensionless\n\n',DP);
   fprintf(fid3(fin),'\n');
   fprintf(fid3(fin),'E-Plane Sectoral Horn\n');
   fprintf(fid3(fin),'---------------------\n');
   fprintf(fid3(fin),['Directivity = ',num2str(10*log10(DE)),' dB\n']);
   fprintf(fid3(fin),['Directivity = ',num2str(DE),' dimensionless\n']);
   fprintf(fid3(fin),'\n');
   fprintf(fid3(fin),'H-Plane Sectoral Horn\n');
   fprintf(fid3(fin),'---------------------\n');
   fprintf(fid3(fin),['Directivity = ',num2str(10*log10(DH)),' dB\n']);
   fprintf(fid3(fin),['Directivity = ',num2str(DH),' dimensionless\n']);
   fprintf(fid3(fin),'\n');
   fprintf(fid3(fin),'*** Note:\n');
   fprintf(fid3(fin),'    The E-Plane Amplitude Pattern is stored in E-Plane Horn.dat\n');
   fprintf(fid3(fin),'    The H-Plane Amplitude Pattern is stored in H-Plane Horn.dat\n');
end;

fclose('all');

end
% Fresnel Subfunction

function[y] = Fresnel(x);

A(1) = 1.595769140;
A(2) = -0.000001702;
A(3) = -6.808508854;
A(4) = -0.000576361;
A(5) = 6.920691902;
A(6) = -0.016898657;
A(7) = -3.050485660;
A(8) = -0.075752419;
A(9) = 0.850663781;
A(10) = -0.025639041;
A(11) = -0.150230960;
A(12) = 0.034404779;

B(1) = -0.000000033;
B(2) = 4.255387524;
B(3) = -0.000092810;
B(4) = -7.780020400;
B(5) = -0.009520895;
B(6) = 5.075161298;
B(7) = -0.138341947;
B(8) = -1.363729124;
B(9) = -0.403349276;
B(10) = 0.702222016;
B(11) = -0.216195929;
B(12) = 0.019547031;

CC(1) = 0;
CC(2) = -0.024933975;
CC(3) = 0.000003936;
CC(4) = 0.005770956;
CC(5) = 0.000689892;
CC(6) = -0.009497136;
CC(7) = 0.011948809;
CC(8) = -0.006748873;
CC(9) = 0.000246420;
CC(10) = 0.002102967;
CC(11) = -0.001217930;
CC(12) = 0.000233939;

D(1) = 0.199471140;
D(2) = 0.000000023;
D(3) = -0.009351341;
D(4) = 0.000023006;
D(5) = 0.004851466;
D(6) = 0.001903218;
D(7) = -0.017122914;
D(8) = 0.029064067;
D(9) = -0.027928955;
D(10) = 0.016497308;
D(11) = -0.005598515;
D(12) = 0.000838386;


if(x==0)
   y=0;
   return
elseif(x<0)
   x=abs(x);
   x=(pi/2)*x^2;
   F=0;
   if(x<4)
      for(k=1:12)
         F=F+(A(k)+j*B(k))*(x/4)^(k-1);
      end
      y = F*sqrt(x/4)*exp(-j*x);
      y = -y;
      return
   else
      for(k=1:12)
         F=F+(CC(k)+j*D(k))*(4/x)^(k-1);
      end
      y = F*sqrt(4/x)*exp(-j*x)+(1-j)/2;
      y =-y;
      return
   end
else
   x=(pi/2)*x^2;
   F=0;
   if(x<4)
      for(k=1:12)
         F=F+(A(k)+j*B(k))*(x/4)^(k-1);
      end
      y = F*sqrt(x/4)*exp(-j*x);
      return
   else
      for(k=1:12)
         F=F+(CC(k)+j*D(k))*(4/x)^(k-1);
      end
      y = F*sqrt(4/x)*exp(-j*x)+(1-j)/2;
      return
   end
end

end

%***********************************************************************
%       elevation(theta,rho,rmin,rmax,rticks,line_style) 
%**********************************************************************
%       GAINPLOT makes an antenna gain plot using polar coordinates of
%       the angle THETA, in radians, versus the radius RHO, where RHO may 
%       be negative.
%
%     - rmin sets the value of the center of the plot. (for ex. -40 dB)
%     - rmax sets the value of the outer ring of the plot (for ex. 0 dB)
%     - rticks is the # of radial ticks (or circles) you want. 
%              it is NOT THE SPACING BETWEEN THEM. Also, if rticks is an even
%              number > 5, then it will be divided by 2, or else if rticks > 5 
%              and divisible by 3, rticks will be divided by 3.  
%     - linestyle is solid or dashed etc. (default is a solid yellow line) 
%
%       POLAR(THETA,RHO,S) uses the linestyle specified in string S.
%       See PLOT for a description of legal linestyles.
%       See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.
%********************************************************************
%       Credits:
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

function hpol = elevation(theta,rho,rmin,rmax,rticks,line_style) 

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
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
	'DefaultTextFontName',   get(cax, 'FontName'), ...
	'DefaultTextFontSize',   get(cax, 'FontSize'), ...
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
	'DefaultTextFontName',   fName , ...
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
