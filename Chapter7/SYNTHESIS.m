%*************************************************************************************************
% SYNTHESIS.m
%*************************************************************************
% This is a MATLAB based program that:
%
% Implements the various synthesis methods presented in Chapter 7 of Antenna Theory and Design. 
% 
% Specifically, the:
%
% A. Schelkunoff
% B. Fourier
% C. Woodward-Lawson
% D. Taylor (both Tschebyscheff-Error and One-Parameter) 
% 
% methods are used to synthesize line sources and/or linear arrays (where applicable). 
%
% The user is guided in entering the correct input data through interactive questions 
% with built-in error checking. 
%
% The output data:
% A. Are plotted
% B. Can be saved both in MATLAB's binary .mat format as well as user-specified text files. 
% .........................................................................
% Written by Marios Gkatzianas
% Arizona State University, September 2002
% .........................................................................
% Revised by Bo Yang
% Arizona State University, September 2004
% ........................................................................
function SYNTHESIS
close all;
clear all;
clc

figure(1)
%%% disp('Please read the following carefully:');
%%% disp('The MATLAB prompt is in pause mode right now and an empty MATLAB figure is shown.');
%%% disp('Resize the size of the figure with the mouse (DON''T PRESS ANY KEY) until you are');
%%% disp('satisfied and THEN PRESS any key. The figure size you have selected will be used');
%%% disp('throughout this MATLAB program for all subsequent figures.');
%%% pause;

pos=get(gcf,'position');
close all;

disp(' ');
disp('Figure position acquired.'); disp(' ');
load bal.mat;

form=5;
while (form~=1)&(form~=2)&(form~=3)&(form~=4),
    form=input(['Select format for output data\n','*****************************\n', ...
         '1. short (5 decimal digits)\n','2. long (15 decimal digits)\n','3. short e\n','4. long e\n', ...
         'Select desired format (ENTER for default=short):']);
    form(isempty(form))=1; 
    
    if (form~=1)&(form~=2)&(form~=3)&(form~=4),
       hform=msgbox('Specified number must be either 1,2,3 or 4','Invalid choice','custom',x,map,'modal'); 
    end;                
end;

switch form
   case 1,
      format short;
   case 2
      format long;
   case 3
      format short e;
   otherwise
      format long e;
end;


method=[];
while isempty(method)|((method~=1)&(method~=2)&(method~=3)&(method~=4)),
   method=input(['Choose one of the following synthesis methods\n', ...
          '*********************************************\n' ...
          '1. Schelkunoff method\n', '2. Fourier transform method\n', ...
          '3. Woodward-Lawson method\n', '4. Taylor method\n','->']);
        
    if isempty(method)|((method~=1)&(method~=2)&(method~=3)&(method~=4)),
        hm=msgbox('Specified number must be either 1,2,3 or 4','Invalid choice','custom',x,map,'modal'); 
    end;   
end;

switch method
   case 1,    % Schelkunoff

        output_mode=[];
        while isempty(output_mode)|((output_mode~=1)&(output_mode~=2)),
            output_mode=input(['Select output method\n', ...
                     '********************\n', '1. Screen\n', ...
                     '2. Output file\n', '->']);
            if isempty(output_mode)|(output_mode~=1&output_mode~=2);
                hf=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal');  
            end;      
        end;		
	
        outfile=[]; outpath=[];
   
%%%	     if (output_mode==2),
%%%	         while isempty(outfile)|(outfile==0), 

% don't comment out the next two lines
%	      outfile=input(['Give name of output file to save Array Factor (don''t', ...
%                     ' use quotes)\n'],'s');

%%%                 if exist('hout'),
%%%                    waitfor(hout);
%%%                 end;  
%%%                 [outfile,outpath]=uiputfile('*.txt','Save As');
 
%%%                 if isa(outfile,'double')|isa(outpath,'double'),     
%%%                     delete(gco);                    
%%%                     hout=msgbox('Incorrect output file specification','Invalid choice', ...
%%%                                 'custom',x,map,'modal');                 
%%%                     set(hout,'interruptible','off','busyaction','queue');
%%%                 end;

%%%             end;
             
%%% 	     end;

        if output_mode==2,
             outfile=input(['Give name of output file to save Array Factor (don''t use quotes)\n'],'s');
             while isempty(outfile),
                outfile=input(['Filename must be non-empty. Re-enter name of output file ' ...
                              '(don''t use quotes).\n'],'s');
             end;     
        end;
   
        outfile=[pwd '\' outfile];
        schel(output_mode,outfile,pos);

   case 2     % Fourier
      
        fourier;

   case 3,    % Woodward

	wood_main(pos);

   case 4,    % Taylor

	taylor_main(pos);

end;  % switch
end



function [st]=chap(SFAF,str_mod);
st=0;
if norm(SFAF-fliplr(SFAF))>=1e-3,
    disp(' ');
    disp(['Warning: Specified ' str_mod ' is not symmetrical with respect to 90 deg.' ...
          ' Phase information is missing and']); 
    disp(['cannnot be uniquely reconstructed. Program will continue but another synthesis method ' ...
         'is strongly suggested.']);
    st=1;    
end;   

if abs(SFAF(1))>=1e-3|abs(SFAF(end))>=1e-3,
   disp(' ');
   disp(['Warning: Specified ' str_mod ' does not have compact support. Invisible region cannot be ' ...
         'uniquely defined.']);
   disp('Program will continue but another synthesis method is strongly suggested.');
   st=1;
end;   
end
function [stnice]=cosm(x);

stnice=num2str(x,'%8.4f');
stnice(stnice=='i')='j';
stnice=cellstr(stnice);                      
           
for ik=1:length(stnice),
   pp=findstr(stnice{ik},'+');  
   if ~isempty(pp),
      stnice{ik}=[stnice{ik}(1:pp-1) ' + ' stnice{ik}(pp+1:end)];
   end;
   pm=findstr(stnice{ik},'-');
   pm(pm==1)=[];     % ignore minus on the real part 
   if ~isempty(pm),
      stnice{ik}=[stnice{ik}(1:pm-1) ' - ' stnice{ik}(pm+1:end)];
   end;   
end;   
end

%************************************************************************
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
function [z_nul,psi_nul,z_invis,z_visi,nul_final]=file_format(z_nul,beta,d, ...
                    filename,output_mode,fid,poly_opt)


z_nul=z_nul(abs(abs(z_nul)-1)<1e-4);  % use tolerance to get all
						 % nulls on unit circle
psi_nul=angle(z_nul);

disp('Phase of roots on the unit circle (in degrees):')
if ~isempty(psi_nul),
   disp(strjust(num2str(sort(psi_nul*180/pi),'%6.2f'),'right'));
else 
   disp('none');
end;   
disp(' '); 
        	
disp(['Visible region extends from ', num2str(beta*180/pi-360*d) ...
      ' degrees to ', num2str(beta*180/pi+360*d) ' degrees']);
disp(' ');

disp('Roots in invisible region:');
z_invis=z_nul(psi_nul<(beta-2*pi*d)|psi_nul>(beta+2*pi*d));

if isempty(z_invis),
   disp('none'); disp(' ');
else
   disp(char(cosm(z_invis))); disp(' ');
end; 

disp('Roots in visible region:');            
z_visi=z_nul(psi_nul>(beta-2*pi*d)&psi_nul<(beta+2*pi*d));
           
if isempty(z_visi),              
   disp('none'); disp(' ');
else            
  disp(char(cosm(z_visi))); disp(' ');           
end;

if (output_mode==2),             
    fprintf(fid,'\n');              
    fprintf(fid,'%% Phase of roots on the unit circle (in degrees):\n');
    if ~isempty(psi_nul),
       fprintf(fid,'%% %-+6.2f\n',psi_nul*180/pi);              
    else
       fprintf(fid,'%% none\n');
    end;   
    fprintf(fid,'\n');              
    fprintf(fid,['%% Visible region extends from %6.2f degrees to ' ...
                 '%6.2f degrees\n'], beta*180/pi-360*d,beta*180/pi+360*d);
    fprintf(fid,'\n');
    fprintf(fid,'%% Roots in invisible region:\n'); 
              
    if isempty(z_invis),         
        fprintf(fid,'%% none\n');                  
    else
        fprintf(fid,'%% %8.3f%+8.3f j\n',[real(z_invis.'); imag(z_invis.')]);
    end;
                           
    fprintf(fid,'\n%% Roots in visible region:\n');

    if isempty(z_visi),	         
        fprintf(fid,'%% none\n');
    else
        fprintf(fid,'%% %8.3f%+8.3f j\n',[real(z_visi.'); imag(z_visi.')]);   
    end;

end;

psi_nul_vis=psi_nul((abs(psi_nul-beta)<=(2*pi*d)+1e-3));   
                       % psi's in visible region (again use tolerance)
	   
if (isempty(psi_nul_vis)),
    disp('No nulls exist');

    if (output_mode==2),
        fprintf(fid,'\n%% No nulls exist \n');
    end;
	   
else

% get only the real part of nulls 	       
theta_nul=real(acos((psi_nul_vis-beta)/(2*pi*d)));

    for mu=1:ceil(d*(1-min(cos(theta_nul)))), % extra code to cover 
                                              % the case d>0.5 lambda  
                 
        theta_nul1(mu,:)=acos(cos(theta_nul')+mu/d);
        theta_nul2(mu,:)=acos(cos(theta_nul')-mu/d);
    end;

    theta_nul=[theta_nul; theta_nul1(:); theta_nul2(:)];
    theta_nul(abs(imag(theta_nul))>=1e-4)=[];
    nul_tot=[theta_nul*180/pi];

	       
    for w=1:length(nul_tot),   % extra code to remove multiple null
		           % occurences
       nul_final(w)=nul_tot(1);
       nul_tot(abs(nul_tot-nul_final(w))<=1e-2)=[];
	          
       if (isempty(nul_tot)),
          break;
       end;
    end;	

	       
    if (isempty(nul_final)),
       disp('No Array Factor nulls exist');
    else
       ats=nul_final(find(nul_final(nul_final<0|nul_final>180)));
       disp(ats);
    end;  
       
    if ~isempty(ats),
        disp('Null found in invisible region');
    end;

    disp('Array Factor nulls inside visible region were found at:');
    disp(strcat(num2str(sort(nul_final'),'%6.2f'),  ...
         repmat(' deg.',length(nul_final),1)));
    
    if exist('fid')&(output_mode==2),  
       fprintf(fid,'\n%% Array Factor nulls inside visible region were found at:\n');  
       fprintf(fid,'%% %6.2f deg.\n',sort(nul_final'));
    end;
             
end;

end
function [prec,HPBW_app,tm]=find_hpbw(theta,U)

% U is radiation intensity vector
% *******************************
U=U(theta<=pi);
theta=theta(theta<=pi);

[Um,im]=max(U);
th_3db=[];
prec=1e-3;

while isempty(th_3db),
   th_3db=theta(abs(U-0.5*Um)<=prec);
%   th_low=max(th_3db(th_3db<=theta(im)));
%   th_high=min(th_3db(th_3db>=theta(im)));
   if isempty(th_3db)|length(th_3db)<2,
      prec=prec*10;
   end;     
end;   

if theta(im)<min(th_3db)|theta(im)>max(th_3db),
   HPBW_app=(min(th_3db)+abs(pi-max(th_3db)))*180/pi;      
else
   th_low=max(th_3db(th_3db<theta(im)));
   th_high=min(th_3db(th_3db>theta(im)));
   HPBW_app=abs(th_low-th_high)*180/pi;   
end;

tm=theta(im);

end

% this program synthesizes a linear array using the Schelkunoff
% polynomial method

% Written by Marios Gkatzianas	
% Arizona State University, September 2002

function []=schel(output_mode,filename,pos);

% load bal.mat;
[x,map]=imread('bal.tiff');

close all;
warning off;

mode_schel=3;
while (mode_schel~=1)&(mode_schel~=2),
    mode_schel=input(['Schelkunoff synthesis options\n', ...
           '*****************************\n',... 
           '1. Compute polynomial for specific nulls\n', ...
           '2. Compute nulls for a specified polynomial\n', ...
           'Choose one of the above:\n','->']);
     
    if isempty(mode_schel)|((mode_schel~=1)&(mode_schel~=2)),
        hso=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal'); 
    end;   
end;

d=0;
while ~isreal(d)|(isempty(d)|(d<=0)),
   d=input('Specify spacing between elements (in wavelengths)\n');
   d=d+1e-12;
   if ~isreal(d)|isempty(d)|d<=0,
      hss=msgbox('Spacing must be positive','Invalid choice','custom',x,map,'modal');
   end;   
end;

if (mode_schel==1),
   disp('It is assumed that beta=0'); beta=0;
else

   beta=[];
   while ~isreal(beta)|isempty(beta),	
      beta=input('Specify progressive phase shift beta (in degrees)\n');
      if ~isreal(beta)|isempty(beta),
         create.WindowStyle='modal';
         create.Interpreter='tex';
         hsp=msgbox(strvcat('Progressive phase \beta must be non-empty','and real'), ...
                           'Invalid choice','custom',x,map,create);
      end;   
   end;

end;
beta=beta*pi/180;
poly_opt=4;
switch mode_schel
    case 1,  % nulls -> poly

        Nnul=0;
        while ~isreal(Nnul)|isempty(Nnul)|Nnul<=0|mod(real(Nnul),1)~=0,
           Nnul=input('Specify number of desired nulls\n');
           if ~isreal(Nnul)|isempty(Nnul)|isempty(Nnul)|Nnul<=0|mod(real(Nnul),1)~=0,
              hnnul=msgbox('Number of nulls must be positive integer','Invalid choice', ...
                           'custom',x,map,'modal');
           end;   
        end;

	Nel=Nnul+1;
	nul=zeros(size(Nnul));

	for n=1:Nnul,
	    aux=[];

	    while ~isreal(aux)|isempty(aux)|aux<0|aux>180,
     	     aux=input(['Specify position of null #',num2str(n), ...
                   ' (in degrees)\n']);
           if ~isreal(aux)|isempty(aux)|aux<0|aux>180,
              hposn=msgbox('Position of null must be between 0 and 180 degrees','Invalid choice', ...
                            'custom',x,map,'modal'); 
           end;     
       end;

	    nul(n)=aux;
	end;

	nul=nul*pi/180;   % convert to rads

    case 2, % poly -> nulls   
     

	while ((poly_opt~=1)&(poly_opt~=2)&(poly_opt~=3)),
	    poly_opt=input([ 'Polynomial definition options\n', ...         
                 '*****************************\n', ...
                 '1. Enter vector of polynomial coefficients in descending form\n', ...
                 '---------------------------------------------------\n', ...
                 'Example: z^3-z^2+z-1 should be input as [1 -1 1 -1]\n', ...
                 '---------------------------------------------------\n\n', ...
                 '2. Specify polynomial as a function of z\n', ... 
                 '--------------------\n', ...
                 'Example: z^3-z^2+z-1\n', ...
                 '--------------------\n\n', ...
                 '3. Specify roots of the desired polynomial\n' ...
                 '--------------------\n', ...
                 'Example(4 roots): +j, +1, -j, 0.707 + 0.707j\n',... 
                 '--------------------\n', ...
                 'Choose one of the above\n','->' ]);
           
       if isempty(poly_opt)|((poly_opt~=1)&(poly_opt~=2)&(poly_opt~=3)),
          hpoly=msgbox('Specified number must be either 1,2 or 3','Invalid choice', ...
                        'custom',x,map,'modal');
       end;   
   end;

	if (poly_opt==1),
	   a=[];
	   while isempty(a),
	      a=input(strcat('Enter row vector of polynomial coefficients', ...
                        ' (in descending powers of z)\n','Example: z^3-z^2+z-1 should be input as [1 -1 1 -1]\n'));	
         if isempty(a)|~isa(a,'double'),
            hrow=msgbox(strvcat('Polynomial coefficients must be inside square brackets'), ...
                        'Invalid choice', 'custom',x,map,'modal');
         end;   
	   end;

	   a=fliplr(a);   % flip a in ascending powers of z   
       
   elseif (poly_opt==2),
       pstring=[];
       while isempty(pstring)|~isempty(find(pstring(isletter(pstring))~='z'&pstring(isletter(pstring))~='i')),
          
           pstring=input(strcat('Specify polynomial as a function of z. ', ...
              ' Use ^ for all powers and i for the imaginary unit and don''t put quotes\n', ...
              'Example: z^3-z^2+z-1\n'),'s');
           if isempty(pstring)|~isempty(find(pstring(isletter(pstring))~='z'&pstring(isletter(pstring))~='i')),                                                      
              hpstr=msgbox(strvcat('String must be non empty and must only contain','a z variable (apart from i)'), ...
                           'Invalid choice','custom',x,map,'modal');             
           end;   
       end;
       pz=inline(vectorize(pstring));

   else  % define roots

       Nroot=[];   
       while (isempty(Nroot)|Nroot<=0|~isreal(Nroot)|mod(Nroot,1)~=0),    
          Nroot=input('Specify number of roots for the polynomial\n');
          if isempty(Nroot)|((Nroot<=0)|~isreal(Nroot)|mod(Nroot,1)~=0),
              hnumroot=msgbox('Number or roots must be positive integer','Invalid choice', ...
                           'custom',x,map,'modal');   
          end;
       end;

       root=[];

       for n=1:Nroot,
           temp_root=[];               
           while isempty(temp_root),
              temp_root=input(['Specify root #', num2str(n),' of the polynomial\n']);
              
              if isempty(temp_root),
                 hroot=msgbox(['Root #' num2str(n) ' must be non-empty'],'Invalid choice', ...
                               'custom',x,map,'modal');
              end;   
           end;

           root=[root temp_root];
           
% the following portion of the code includes the conjugate roots to make the polynomial real. 
% To be used optionally           
           
%           if imag(temp_root)~=0,
%               disp(['The root ', num2str(conj(temp_root)) ' is also included to' ...
%                     ' make the polynomial have real coefficients']);
%               root=[root conj(temp_root)];
%           end;
        end;
                
        polycoef=poly(root);
        a=polycoef;
        
        spoly=[];
        for n=1:length(a)-1,
           stemp=['(' num2str(real(a(n)),'%+6.3f') num2str(imag(a(n)),'%+6.3f') 'j)' '*z^' num2str(length(a)-n) ' '];
           spoly=[spoly stemp];
        end;
        
        spoly=[spoly '(' num2str(real(a(end)),'%+6.3f') num2str(imag(a(end)),'%+6.3f') 'j)'];
                
    end;

end;  % end of switch 

Ntheta=0;

while ~isreal(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0,
   Ntheta=input(['Specify number of sample points for the angle theta', ...
         ' for the pattern plots\n', '[ENTER for default=360]: ']); 
   Ntheta(isempty(Ntheta))=360;
   
   if ~isreal(Ntheta)|Ntheta<=0|mod(Ntheta,1)~=0,
      create.WindowStyle='modal';
      create.Interpreter='tex';
      hNthe=msgbox('Number of sample points for \theta must be integer positive','Invalid choice', ...
         'custom',x,map,create);      
   end;   
end;

thres=0;
while ~isreal(thres)|thres>=0,
   thres=input(['Specify lowest threshold of Array Factor for plotting' ...
         ' purposes (in -dB)\n','[ENTER for default=-120]:']);
   thres(isempty(thres))=-120;
   
   if ~isreal(thres)|thres>=0,
      hthres=msgbox('Threshold must be negative','Invalid choice','custom',x,map,'modal');
   end;
   
end;
thres(isempty(thres))=-120;

theta=linspace(0,2*pi,2*Ntheta);
psi=2*pi*d*cos(theta)+beta;


switch mode_schel,
    case 1   % nulls -> polynomial	

	     psi_nul=2*pi*d*cos(nul)+beta;
	     a=poly(exp(j*psi_nul));  % coefficients of z polynomial
        
        spoly=[];
        for n=1:length(a)-1,
           stemp=['(' num2str(real(a(n)),'%+6.3f') num2str(imag(a(n)),'%+6.3f') 'j)' '*z^' num2str(length(a)-n) ' '];
           spoly=[spoly stemp];
        end;
        
        spoly=[spoly '(' num2str(real(a(end)),'%+6.3f') num2str(imag(a(end)),'%+6.3f') 'j)'];       
        
	     z=exp(j*psi);

	     for n=1:Nel,
	        temp(n,:)=z.^(n-1);
	     end;

	     AF=fliplr(a)*temp;  % final array factor
        	
        for mu=1:ceil(d*(1-min(cos(nul)))),    % extra code to cover the case
					       % d>0.5 lambda  
           theta_nul1(mu,:)=acos(cos(nul)+mu/d)*180/pi;
           theta_nul2(mu,:)=acos(cos(nul)-mu/d)*180/pi;
        end;

	     theta_nul=[theta_nul1(:) theta_nul2(:)];
        theta_nul(abs(imag(theta_nul))>=1e-4)=[];
	     nul_tot=[nul*180/pi theta_nul];

	     for w=1:length(nul_tot),   % extra code to remove multiple null
				   % occurences
	         nul_final(w)=nul_tot(1);
            nul_tot(abs(nul_tot-nul_final(w))<=1e-2)=[];
	         if (isempty(nul_tot)),
                break;
	         end;
	     end;
        
        disp(' '); disp('Schelkunoff polynomial found as:');
        disp(spoly); disp(' ');
	     disp('Array Factor nulls found at:');
		  disp(strcat(num2str(sort(nul_final)','%6.2f'), ...
             repmat(' deg.',length(nul_final),1)))                  
        fprintf(2,'\nExcitation coefficients:\n');   % save output on screen too
        fprintf(2,'   Real part \t    Imag part\n');
        fprintf(2,'%10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);
  
          
        if output_mode==2,
           if ~exist('fid'),
              fid=fopen(filename,'wt');
           end;
           fprintf(fid,'%% Schelkunoff synthesis method\n\n');
           fprintf(fid,['%% Spacing d=' num2str(d) ' lambda\n']);
           fprintf(fid,['%% Progressive phase beta=' num2str(beta*180/pi) ' degrees\n']);             
           fprintf(fid,'\n%% User specified nulls:\n');
           fprintf(fid,'%% %5.2f\n',nul'*180/pi);
           fprintf(fid,'\n%% Schelkunoff polynomial:\n');
           fprintf(fid,'%% %s\n',spoly);  
           fprintf(fid,'\n%% Array Factor nulls found at:\n');
           fprintf(fid,'%% %5.2f\n',sort(nul_final'));
           fprintf(fid,'\n%% Excitation coefficients:\n');
           fprintf(fid,'%% Real part \t Imag part\n');
           fprintf(fid,'%% %10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);                                  
        end;

    case 2  % polynomial -> nulls

        if (poly_opt==3),
            psi=2*pi*d*cos(theta)+beta;
            z=exp(j*psi);
            AF=polyval(polycoef,z);
            strpoly=[];
            z_nul=root.';
            Nel=length(z_nul)+1;            

            for n=1:length(polycoef)-1,
               if polycoef(n)~=0
                  tempstr=[num2str(polycoef(n),'%+5.2f') '*z^' ...
                           num2str(length(polycoef)-n)];
               end; 
  	            strpoly=[strpoly tempstr];                
            end;
            
            if polycoef(end)~=0, 
               strpoly=[strpoly num2str(polycoef(end),'%+5.2f')];
            end;  

            disp('Roots of the polynomial:');
            disp(char(cosm(root.')));
            disp('The Schelkunoff polynomial is:');           
            disp(spoly);

	         if (output_mode==2),
                fid=fopen(filename,'wt');
                fprintf(fid,'%% Schelkunoff synthesis method\n');
                fprintf(fid,['%% Spacing d=' num2str(d) ' lambda\n']);
                fprintf(fid,['%% Progressive phase beta=' num2str(beta*180/pi) ' degrees\n\n']);
                fprintf(fid,'%% Roots of the polynomial:\n');
                fprintf(fid,'%% %8.3f%+8.3f j\n',[real(z_nul.'); imag(z_nul.')]);         
                fprintf(fid,'\n%% The Schelkunoff polynomial is:\n');
                fprintf(fid,'%% %s\n',spoly);
                fprintf(fid,'\n%% Excitation coefficients:\n');
                fprintf(fid,'%% Real part \t Imag part\n');
                fprintf(fid,'%10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);
            end;      

            if ~exist('fid'),
	             fid=[];          
            end;

            [z_nul,psi_nul,z_invis,z_visi,nul_final]=file_format(z_nul,beta,d, ...
                     filename,output_mode,fid,poly_opt);

        end;  % poly_opt=3

	     if (poly_opt==1),  % row vector
	        z=exp(j*psi);
           Nel=length(a);

	        for n=1:Nel,
	           temp(n,:)=z.^(n-1);
           end;                      
           
           AF=fliplr(a)*temp;
 
	        z_nul=roots(a);
           a=fliplr(a);
           
           spoly=[];
           for n=1:length(a)-1,
              stemp=['(' num2str(real(a(n)),'%+6.3f') num2str(imag(a(n)),'%+6.3f') 'j)' ...
                     '*z^' num2str(length(a)-n) ' '];
               spoly=[spoly stemp];
           end;
        
           spoly=[spoly '(' num2str(real(a(end)),'%+6.3f') num2str(imag(a(end)),'%+6.3f') 'j)'];
           
           disp(' ');
           disp('Schelkunoff polynomial is:');
           disp(spoly); disp(' ');                                                  
           
           disp('Roots of the polynomial:');
           disp(char(cosm(z_nul)));                      
                                          
           if (output_mode==2),
              fid=fopen(filename,'wt');
              fprintf(fid,'%% Schelkunoff synthesis method\n\n');
              fprintf(fid,['%% Spacing d= ' num2str(d) ' lambda\n']);
              fprintf(fid,['%% Progressive phase beta= ' num2str(beta*180/pi) ' degrees\n\n']);                      
           
              strpoly=[];
              for n=length(a):-1:2,
                 if a(n)~=0
                    tempstr=[num2str(a(n),'%+5.2f') '*z^' ...
                             num2str(n-1)];
                 end; 
  	              strpoly=[strpoly tempstr];                
              end;
           
              if a(1)~=0, 
                 strpoly=[strpoly num2str(a(1),'+%6.2f')];
              end; 
                     
              fprintf(fid,'%% The Schelkunoff polynomial is:\n%% %s\n\n',spoly);
              fprintf(fid,'%% Roots of the polynomial:\n');
              fprintf(fid,'%% %8.3f%+8.3f j\n',[real(z_nul.'); imag(z_nul.')]);  
              fprintf(fid,'\n%% Excitation coefficients:\n');
              fprintf(fid,'%% Real part \t Imag part\n');
			     fprintf(fid,'%% %10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);        
           end;                     

           if ~exist('fid'),
	          fid=[];          
           end;

           [z_nul,psi_nul,z_invis,z_visi,nul_final]=file_format(z_nul,beta,d,filename,output_mode,fid,poly_opt);                                 
                 
        elseif (poly_opt==2),   % inline function
                                % attention:  add more code for output file       
                                                                        
	        z=exp(j*psi);
           AF=pz(z);
           
           strpz=char(pz);
           strpz(strpz=='.')=[]; 
           
           if exist('sym')~=2,
              disp(['Symbolic toolbox not installed. Skipping symbolic null computations.' ...
                    ' You have to determine the nulls visually.']);
              disp(' '); disp('Schelkunoff polynomial is:');
              disp(strpz); disp(' ');
              
           else   
              zs=sym('zs');                                          
              pzs=strrep(char(pz),'z','zs');
              pzs(pzs=='.')=[];
              zpoly=sym2poly(sym(pzs));
              z_nul=roots(zpoly);           
              a=zpoly;                                         
                                         
              disp(' '); disp('Schelkunoff polynomial is:');
              disp(strpz); disp(' ');            
              disp('Excitation coefficients:');
              disp(char(cosm(num2str(a.','%8.4f')))); disp(' ');
                            
              if ~exist('fid')&(output_mode==2),
                 fid=fopen(filename,'wt');
                 fprintf(fid,'%% Schelkunoff synthesis method\n\n');
                 fprintf(fid,['%% Spacing d= ' num2str(d) ' lambda\n']);
                 fprintf(fid,['%% Progressive phase beta= ' num2str(beta*180/pi) ' degrees\n\n']);                 
                 fprintf(fid,'%% The Schelkunoff polynomial is:\n%% %s\n',strpz);
                 fprintf(fid,'\n%% Excitation coefficients:\n');
                 fprintf(fid,'%% Real part \t Imag part\n');
                 fprintf(fid,'%% %10.4f \t %10.4f\n',[real(fliplr(a)); imag(fliplr(a))]);
              end;
              
              if ~exist('fid'),
                 fid=[];
              end;
              
              [z_nul,psi_nul,z_invis,z_visi,nul_final]=file_format(z_nul,beta,d, ...
                     filename,output_mode,fid,poly_opt);                                 
              Nel=length(a);
           end;                     
           
	     end;   % if 

end;

	AF_db=20*log10(abs(AF)/max(abs(AF)));
	AF_db(AF_db<=thres)=thres;
   
% directivity computation
% ***********************   

dth=theta(2)-theta(1);   % theta in radians
U=(abs(AF(theta<pi))).^2;
Prad=2*pi*sum(U.*sin(theta(theta<pi))*dth);
D=4*pi*max(U)/Prad;
D_db=10*log10(D);
   
%if (mode_schel==1),
if exist('a'),

      % Figure 1
      % ********
      figure(1);
      set(gcf,'units','pixels','position',pos);
      subplot(2,1,1); stem(fliplr(abs(a))); grid on; 
      set(gca,'fontname','timesnewroman','fontsize',14,'xtick',[1:Nel]);	
      ylabel('Amplitude excitation','fontname','timesnewroman','fontsize',18);
      title(['Amplitude and phase excitation (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'],'fontsize',18);        
      set(findobj(findobj(gca,'type','line'),'color',[0 0 1]),'linewidth',2);

      subplot(2,1,2); stem(fliplr(angle(a))*180/pi); grid on;
      set(gca,'fontname','timesnewroman','fontsize',14);   
      set(gca,'xtick',[1:Nel]);	set(gca,'ylim',[-180 180]);
%      set(gca,'ytick',[-180 -90 0 90 180]);
      xlabel('element number','fontname','timesnewroman','fontsize',18);
      ylabel('Phase excitation (in degrees)','fontname','timesnewroman','fontsize',18);  
      set(findobj(findobj(gca,'type','line'),'color',[0 0 1]),'linewidth',2);

end;

% Figure 2
% ********
figure(2);
set(gcf,'units','pixels','position',pos);
set(gca,'fontname','timesnewroman','fontsize',14); 
hl=plot(theta*180/pi,AF_db); grid; zoom; thres=get(gca,'ylim'); grid on;
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[0 180]);
set(hl,'linewidth',2); 
xlabel('\theta (in degrees)','fontname','timesnewroman','fontsize',18, ...
       'verticalalignment','top');
ylabel('Normalized Array Factor (in dB)','fontname','timesnewroman', ...
       'fontsize',18,'verticalalignment','bottom');
title(['Synthesized Array Factor using Schelkunoff polynomial(N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
      'fontname','timesnewroman','fontsize',18);


% Figure 3
% ********
figure(3);

if (mod(thres(1),4)~=0),  
   thres(1)=ceil(thres(1)/4)*4;   
end;

set(gcf,'units','pixels','position',pos);

elevation(theta,AF_db,thres(1),0,4); 
set(gca,'fontname','timesnewroman','fontsize',14);
set(findobj(gca,'type','text'),'fontname','timesnewroman','fontsize',14);
set(findobj(findobj(gca,'type','line'),'color',[0 0 1]),'linewidth',2);

title(['Synthesized Array Factor using Schelkunoff polynomial(N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
      'fontname','timesnewroman','fontsize',18);
   
disp(' ');
disp(['Directivity (dimensionless) as computed from synthesized Array Factor = ' num2str(D,'%6.2f')]);
disp(['Directivity (dB) as computed from synthesized Array Factor = ' num2str(D_db,'%6.2f\n')]);
           
if (output_mode==2),

   if exist('poly_opt')&poly_opt==1&~isempty(nul_final),
       fprintf(fid,'\n%% Array Factor nulls found at:\n');
       fprintf(fid,'%% %6.2f deg. \n',sort(nul_final'));
   end;

   if ~exist('fidm'),
      fidm=fopen('AF_array.txt','wt');
   end;
   
%% WARNING: add more code to write relevant information in second file   
  

   fprintf(fidm,'%% Schelkunoff polynomial:\n');
   if poly_opt==2,
      fprintf(fidm,['%% ' strpz '\n']);
   elseif mode_schel==1|poly_opt==1|poly_opt==3,
      fprintf(fidm,['%% ' spoly '\n']);
   end;      
   fprintf(fidm,'\n%s\n','% Theta (deg.)    AF (dB)');    
   fprintf(fidm,'   %6.2f \t %6.2f\n', [(theta*180/pi); AF_db]);
   fprintf(fidm,'\n');   
   fclose(fidm);
   fprintf(fid,'\n%% Directivity (dimensionless) as computed from Array Factor = %6.2f\n',D);
   fprintf(fid,'%% Directivity (dB) as computed from Array Factor = %6.2f\n',D_db);  
   fclose(fid);
   
   disp(['Output saved in file ''',filename,'''.']);
   disp(['Array Factor saved in file ''' pwd '\AF_array.txt''']);   
end;  

disp(['Workspace saved in file ''' pwd '\synth.mat''']);
save synth.mat;

warning on;
end

function []=taylor_cheby;

% this program synthesizes a line source using the Taylor method based on the 
% Chebyshev error

% Written by Marios Gkatzianas
% Arizona State University, September 2002

close all;
clear all;

nbar=0;
while (nbar<=1),
   if (nbar==1),
      disp('nbar must be greater than 1');
   end;   
   nbar=input('Specify nbar\n');
end;

side_lobe=0;
while (side_lobe<=0),
   side_lobe=input([ 'Specify desired side lobe level for the first ',... 
                    num2str(nbar),' lobes (in dB)\n']);
end;   

R0=10^(side_lobe/20);  % voltage ratio
A=1/pi*acosh(R0);
sig=nbar/sqrt(A^2+(nbar-0.5)^2);

len=0;
while (len<=0),
   len=input('Specify the length of the line source (in wavelengths)\n');
end;   

N=0;
while (N<=0),
   N=input('Select number of sample points for line source\n');   
end;
zs=linspace(-len/2,len/2,N);

n=1:(nbar-1);
u_con=pi*sig*sqrt(A^2+(n-0.5).^2);

if (len>=nbar),
   u_decr=(nbar:floor(len))*pi;
end;   

u=[u_con u_decr];
theta_null=acos(u/(pi*len))*180/pi; 
theta_null_sort=sort([theta_null 180-theta_null]);

p=1:(nbar-1);
c1=(fact2(nbar-1))^2./(fact2(nbar-1+p).*fact2(nbar-1-p));
temp=ones(size(p));

for m=1:(nbar-1),
   temp=temp.*(1-(p*pi/u(m)).^2);
end;   
SF_coef=c1.*temp;

Is=ones(size(zs));

Ntheta=0;
while (Ntheta<=0),
   Ntheta=input('Give number of sample points for the angle domain [0 180]\n');
end;   
theta=linspace(0,pi,Ntheta);

u2=pi*len*cos(theta);
SF=sinc(u2/pi);

for q=1:(nbar-1),
   Is=Is+2*SF_coef(q)*cos(2*pi*q*zs/len);
   SF=SF.* (1-(u2/u(q)).^2)./(1-(u2/(q*pi)).^2);
end;
Is=Is/len;

thres=0;
while (thres>=0),   
   thres=input('Give threshold for SF (in -dB)\n');
end;   

SF_db=20*log10(abs(SF)/max(abs(SF)));
SF_db(SF_db<=thres)=thres;

% Figure 1
% ********
stem(p,abs(SF_coef));
set(gca,'xtick',[1:(nbar-1)]);

% Figure 2
% ********
figure(2);
plot(zs,abs(Is)/max(abs(Is))); grid;
title('Normalized current distribution');
set(gca,'xlim',[-len/2 len/2]); %,'xtick',linspace(-len/2,len/2,6));
xlabel('z/\lambda','fontsize',12);

% Figure 3
% ********
figure(3);
plot(theta*180/pi,SF_db); grid; hold on;
xlabel('\theta (in \circ)'); ylabel('SF (dB)');
title('Normalized space factor');
h=line([0 180],[-side_lobe -side_lobe]); set(h,'linestyle','--','color','r');

disp('The SF nulls (in degrees) are:');
disp(num2str(theta_null_sort','%6.2f'));
end

function [y]=fact(x);

for n=1:length(x),
   y(n)=prod(1:x(n));
end;   
end
% this program synthesizes a line source using the Taylor method 
% variations of Chebyshev error and one parameter

% Written by Marios Gkatzianas  
% Arizona State University, September 2002

function []=taylor_main(pos);

% load bal.mat; 
[x,map]=imread('bal.tiff');

close all;

taylor_mode=[];
while ~isreal(taylor_mode)|isempty(taylor_mode)|((taylor_mode~=1)&(taylor_mode~=2)),
   taylor_mode=input( ['Select method of Taylor line source synthesis\n' ...
                     '*********************************************\n' ...
                     '1. Tschebyscheff error\n','2. One-Parameter\n','->'] );
   if ~isreal(taylor_mode)|isempty(taylor_mode)|(taylor_mode~=1&taylor_mode~=2),           
      htaymod=msgbox('Specified number must be either 1 or 2' ,'Invalid choice','custom',x,map,'modal');     
   end;   
end;

output_mode=[];
while ~isreal(output_mode)|isempty(output_mode)|((output_mode~=1)&(output_mode~=2)),
    output_mode=input(['Select output method\n','********************\n', ...
          '1. Screen\n','2. Output file\n','->']);
    if ~isreal(output_mode)|isempty(output_mode)|(output_mode~=1&output_mode~=2),
       houtmod=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal');
    end;   
end;


outfile=[]; outpath=[];

%%% if (output_mode==2),
%%%	 while isempty(outfile)|(outfile==0), 

% don't comment in the next two lines for GUI style

%	      outfile=input(['Give name of output file to save Array Factor (don''t', ...
%                     ' use quotes)\n'],'s');

%%%       if exist('hout'),
%%%           waitfor(hout);
%%%       end;  
%%%       [outfile,outpath]=uiputfile('*.txt','Save As');
 
%%%       if isa(outfile,'double')|isa(outpath,'double'),     
%%%           delete(gco);                    
%%%           hout=msgbox('Incorrect output file specification','Invalid choice', ...
%%%                       'custom',x,map,'modal');                 
%%%           set(hout,'interruptible','off','busyaction','queue');
%%%        end;

%%%     end;
             
%%% end;
         
%%% outfile=[outpath outfile];

if output_mode==2,
   outfile=input('Give name of output file to save Space Factor (don''t use quotes)\n','s');
   while isempty(outfile),
      outfile=input('Filename must be non-empty. Re-enter name of output file (don''t use quotes).\n','s');
   end;      
   outpath=pwd;
   outfile=[outpath '\' outfile];   
end;

len=[];
while ~isreal(len)|isempty(len)|len<=0,
   len=input('Specify length of line source (in wavelengths)\n');
   if ~isreal(len)|isempty(len)|len<=0,
      hlen=msgbox('Length must be non-empty and real positive','Invalid choice','custom',x,map,'modal');
   end;   
end;

side_lobe=[];
while ~isreal(side_lobe)|isempty(side_lobe)|side_lobe>=0,
   if taylor_mode==1,
      side_lobe=input('Specify desired constant side lobe level (in -dB)\n');
   else
      side_lobe=input('Specify desired maximum side lobe level (in -dB)\n');
   end;   
   if ~isreal(side_lobe)|isempty(side_lobe)|side_lobe>=0,
      hslobe=msgbox('Side lobe leve must be real negative','Invalid choice','custom',x,map,'modal');
   end;   
end;

if (taylor_mode==1),
    nbar=[];
    while ~isreal(nbar)|isempty(nbar)|nbar<=0|mod(nbar,1)~=0|real(nbar)>len, 
       nbar=input('Specify number of side lobes of the same constant level (nbar)\n');
       if ~isreal(nbar)|isempty(nbar)|nbar<=0|mod(real(nbar),1)~=0|real(nbar)>len,          
          hnbar=msgbox(strvcat('nbar must be a positive integer less than the' , ... 
                    'length of the line source'),'Invalid choice','custom',x,map,'modal');
       end;   
    end;
end;

Ntheta=[];
while ~isreal(Ntheta)|isempty(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0, 
   Ntheta=input(['Specify number of samples for the angle theta for the pattern plots\n' ...
                 '[ENTER for default=360]:']);
   Ntheta(isempty(Ntheta))=360;
   if ~isreal(Ntheta)|isempty(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0,        
      hNthe=msgbox('Specified number must be positive integer','Invalid choice','custom',x,map,'modal');
   end;   
end;

theta=linspace(0,pi,Ntheta);
u=pi*len*cos(theta);

Nz=[];
while ~isreal(Nz)|isempty(Nz)|Nz<=0|mod(real(Nz),1)~=0,
   Nz=input(['Specify number of sample points for the line source for plotting purposes\n' ...
            '[ENTER for default=200]:']);
   Nz(isempty(Nz))=200;
   if ~isreal(Nz)|isempty(Nz)|Nz<=0|mod(real(Nz),1)~=0,   
      hNz=msgbox('Specified number must be positive integer','Invalid choice','custom',x,map,'modal');
   end;   
end;

thres=[];
while ~isreal(thres)|isempty(thres)|thres>=0,
   thres=input(['Specify the lowest threshold for dB scale (in -dB) for plotting purposes ' ...
               '[ENTER for default=-120]:']);
   thres(isempty(thres))=-120;
   if ~isreal(thres)|isempty(thres)|thres>=0,
      hthre=msgbox('Specified number must be real negative','Invalid choice','custom',x,map,'modal');
   end;   
end;

R0=10^(-side_lobe/20);

switch taylor_mode 
    case 1   % Chebyshev error

	A=1/pi*acosh(R0);
	sig=nbar/sqrt(A^2+(nbar-0.5)^2);

	n=1:nbar-1;
	costh1=sig/len*sqrt(A^2+(n-0.5).^2);
	un=pi*len*costh1(abs(costh1)<=1);
	costh1=[-costh1 costh1];

	ng=nbar:floor(len);
	costh2=ng/len;
	costh2=[-costh2 costh2];

	costh=[costh1 costh2];
	costh(abs(costh)>1)=[];

	theta_nul=sort(acos(costh))*180/pi;  % side lobes with constant level

	SF_coef=ones(size(un));

	for p=1:nbar-1,
	    SF_coef(p)=(fact2(nbar-1))^2/(fact2(nbar-1+p)*fact2(nbar-1-p))* ...
                        prod(1-(p*pi./un).^2);
	end

	zs=linspace(-len/2,len/2,Nz);
	Is=ones(size(zs));

	for p=1:nbar-1,
	   Is=Is+2*SF_coef(p)*cos(2*pi*p*zs/len);
	end;
	Is=Is/len;

SF_rec=sinc(u/pi);   

   
	for p=1:nbar-1,
	   SF_rec=SF_rec.*(1-(u/un(p)).^2)./(1-(u/(p*pi)).^2);
	end;

     case 2   % one parameter

	     lobe_look=[-10:-5:-40];  % lookup table
	     B=[j*0.4597 0.3558 0.7386 1.0229 1.2761 1.5136 1.7415];
        
        B=B(find(lobe_look==side_lobe));

        if isempty(B),   % find B numerically           
	        f=inline(['sinh(pi*x)/(pi*x)-', num2str(R0), '/4.603']); 
	        x=fzero(f,1,1e-3);
	        B=x;
        end;
  
	     zs=linspace(-len/2,len/2,Nz);
	     Is=besselj(0,j*pi*B*sqrt(1-(2*zs/len).^2));

	     SF_rec=zeros(size(theta));

 	     SF_rec(abs(u)<pi*B)=len*sinh(sqrt((pi*B)^2-u(abs(u)<pi*B).^2))./ ...
                            sqrt((pi*B)^2-u(abs(u)<pi*B).^2);

        SF_rec(abs(u)>pi*B)=len*sin(sqrt(u(abs(u)>pi*B).^2-(pi*B)^2))./ ...
			    sqrt(u(abs(u)>pi*B).^2-(pi*B)^2);

end;



SF_recdb=20*log10(abs(SF_rec)/max(abs(SF_rec)));
SF_recdb(SF_recdb<=thres)=thres;

if (taylor_mode==1),
   disp(' ');
   disp(['Scaling factor sigma = ' num2str(sig)]); disp(' ');
   disp('Space Factor nulls located at:');
   disp([num2str(theta_nul','%6.2f') repmat(' deg.',length(theta_nul),1)]);
   HPBW_app=2*asin(sig/(pi*len)*sqrt((acosh(R0)^2-acosh(R0/sqrt(2))^2)))*180/pi;
   
   [y,ind]=min(abs(SF_recdb+3));
   HPBW_graph=2*abs(theta(ind)*180/pi-90);
   disp(' '); disp(['HPBW based on approximate formula = ', num2str(HPBW_app,'%6.2f') ' deg.']);
   disp(['HPBW as computed from the Space Factor graph = ' num2str(HPBW_graph,'%6.2f') ' deg.']);
   
   if output_mode==2,
      if ~exist('fid'),
         fid=fopen(outfile,'wt');
      end;
      fprintf(fid,'%% Taylor synthesis method based on Tschebysheff error\n\n');
      fprintf(fid,['%% Length of line source = ' num2str(len) ' lambda\n']);
      fprintf(fid,['%% Desired constant side lobe level = ' num2str(side_lobe) ' dB\n']);
      fprintf(fid,['%% Desired number of equilevel lobes nbar = ' num2str(nbar) '\n\n']);
      fprintf(fid,['%% Scaling factor sigma = ' num2str(sig) '\n\n']);
      fprintf(fid,'%% Space Factor nulls located at:\n');   
      fprintf(fid,'%% %6.2f deg.\n',theta_nul);
      fprintf(fid,'\n\n%% HPBW based on approximate formula = %6.2f deg.\n',HPBW_app);            
      fprintf(fid,'%% HPBW as computed from the Space Factor graph = %6.2f deg.\n\n',HPBW_graph);
   end;
   
end;

if taylor_mode==2&output_mode==2,
   if ~exist('fid'),
      fid=fopen(outfile,'wt');
      fprintf(fid,'%% Taylor synthesis method based on One-Parameter\n\n');
      fprintf(fid,['%% Length of line source = ' num2str(len) ' lambda\n']);      
      fprintf(fid,['%% Maximum side lobe level = ' num2str(side_lobe) ' dB\n\n']);
   end;   
end;   


% Figure 1
% *********
set(gcf,'units','pixels','position',pos)
plot(zs,abs(Is/max(abs(Is))),'linewidth',2); grid;
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[-len/2, len/2]); 
xlabel('z^{\prime}/\lambda','fontname','timesnewroman','fontsize',18,'verticalalign','top');
ylabel('Current amplitude','fontname','timesnewroman','fontsize',18, ...
       'verticalalign','bottom');

if (taylor_mode==1),
   title(['Current distribution for the Taylor method (Tschebyscheff error, L = ',num2str(len),'\lambda)',], ...
         'fontname','timesnewroman','fontsize',18);
else
   title(['Current distribution for the Taylor method (One-Parameter, L = ',num2str(len),'\lambda)',], ...
         'fontname','timesnewroman','fontsize',18);
end;

% Figure 2
% ********
figure(2);
set(gcf,'units','pixels','position',pos);
plot(theta*180/pi,SF_recdb,'linewidth',2); grid;
set(gca,'xlim',[0 180],'fontname','timesnewroman','fontsize',14); 

xlabel('\theta (in degrees)','fontname','timesnewroman','fontsize',18,'verticalalign','top');
ylabel('Normalized Space Factor (dB)','fontname','timesnewroman','fontsize',18, ...
       'verticalalign','bottom');

h=line([0 180],[-3 -3]); set(h,'color','r','linestyle','--','linewidth',2);
hth=line([0 180],[side_lobe side_lobe],'color',[81 164 109]/255,'linestyle','--','linewidth',2);

set(gca,'units','normalized');
text(0.85,0.95,'-3 dB level','units','normalized','fontname','timesnewroman','fontsize',14);

if (taylor_mode==1),
    title(['Synthesized SF: Taylor method based on Tschebyscheff error (L = ',num2str(len),'\lambda)',], ...
          'fontname','timesnewroman','fontsize',18);
else
   title(['Synthesized SF: Taylor method based on One-Parameter (L = ',num2str(len),'\lambda)',], ...
         'fontname','timesnewroman','fontsize',18);
end;

% Figure 3
% ********
figure(3);
set(gcf,'units','pixels','position',pos);
elevation([theta theta+pi],[SF_recdb fliplr(SF_recdb)],thres(1),0,4);
set(gca,'fontname','timesnewroman','fontsize',14);
set(findobj(gca,'type','text'),'fontname','timesnewroman','fontsize',14);
set(findobj(findobj(gca,'type','line'),'color',[0 0 1]),'linewidth',2);

if (taylor_mode==1),
   title(['Synthesized SF: Taylor method based on Tschebyscheff error (L = ',num2str(len),'\lambda)',], ...
          'fontname','timesnewroman','fontsize',16)
else   
   title(['Synthesized SF: Taylor method based on One-Parameter (L = ',num2str(len),'\lambda)',], ...
         'fontname','timesnewroman','fontsize',16);

end;

if ~exist('HPBW_app'),  % one-parameter 
   [y,in]=min(abs(SF_recdb+3));
   HPBW_app=2*abs(90-theta(in)*180/pi);
end;

U=(abs(SF_rec(theta<pi))).^2;
Prad=sum(U.*sin(theta(theta<pi)))*(theta(2)-theta(1));
D=2*max(U)/Prad;
D_mcd=101/(HPBW_app-0.0027*HPBW_app^2);
D_db=10*log10(D);
D_poz=-172.4+191*sqrt(0.818+1/HPBW_app);

if (output_mode==2),
                 
    if taylor_mode==1,
       fidsf=fopen('SF_Taylor.txt','wt');
       fprintf(fidsf,'%% Space factor for the Taylor method based on the Tschebyscheff error\n\n');  
       fprintf(fidsf,'  %% Theta (deg.)      SF (-dB) \n');    
       fprintf(fidsf,'    %6.2f \t      %6.2f\n',[theta*180/pi; SF_recdb]);               
    else
       fidsf=fopen('SF_Taylor.txt','wt');       
       fprintf(fidsf,'%% Space factor for the Taylor method based on One-Parameter\n \n');       
       fprintf(fidsf,'  %% Theta (deg.)      SF (-dB) \n');    
       fprintf(fidsf,'    %6.2f \t      %6.2f\n',[theta*180/pi; SF_recdb]); 
       fprintf(fid,'\n\n%% HPBW based on approximate formula = %6.2f deg.\n\n',HPBW_app);
    end;    
    
    fprintf(fid,'%% Directivity (dimensionless) as computed from synthesized Space Factor = %6.2f\n',D);
    fprintf(fid,'%% Directivity (dB) as computed from synthesized Space Factor = %6.2f\n\n',D_db);       
    fprintf(fid,'%% Directivity (dimensionless) from McDonald''s formula = %6.2f\n',D_mcd);
    fprintf(fid,'%% Directivity (dB) from McDonald''s formula = %6.2f\n\n',10*log10(D_mcd));
    fprintf(fid,'%% Directivity (dimensionless) from Pozar''s formula = %6.2f\n',D_poz);      
    fprintf(fid,'%% Directivity (dB) from Pozar''s formula = %6.2f\n',10*log10(D_poz));
        
end;

disp(' ');
if taylor_mode==2,
   disp(['HPBW as computed from the Space Factor graph = ' num2str(HPBW_app,'%6.2f') ' deg.']);
   disp(' ');
end;   
disp(['Directivity (dimensionless) as computed from synthesized Space Factor = ' num2str(D,'%6.2f')]);
disp(['Directivity (dB) as computed from synthesized Space Factor = ' num2str(D_db,'%6.2f')]); 
disp(' ');
disp(['Directivity (dimensionless) from McDonald''s formula = ' num2str(D_mcd,'%6.2f')]);
disp(['Directivity (dB) from McDonald''s formula = ' num2str(10*log10(D_mcd),'%6.2f')]); 
disp(' ');
disp(['Directivity (dimensionless) from Pozar''s formula = ' num2str(D_poz,'%6.2f')]);
disp(['Directivity (dB) from Pozar''s formula = ' num2str(10*log10(D_poz),'%6.2f')]);
disp(' ');  
if output_mode==2,
   disp(['Output data saved in file ''' outfile '''']);
   disp(['Space Factor saved in file ''' pwd '\SF_Taylor.txt''']);
end;   

disp(['Workspace saved in file ''' pwd '\synth.mat''']);
save synth.mat;
fclose('all');

end

function [y]=fact2(x);

for i=1:length(x),
   if (x(i)==0),
      y(i)=1;
   else
      y(i)=prod(1:x(i));
   end;
end;
end
function [theta,theta_samp,b,varargout]=wood(str,d,len,Nel,sfaf_mode,filename,pstring,Nsamp,Ntheta,theta,SFAF);

close all;

if ~isempty(pstring),
   px=inline(vectorize(pstring));
else
   px=[];
end;

if (mod(Nsamp,2)==0),   % even number    
    m=1:(Nsamp/2);
    theta_samp=acos([(2*m-1)/(2*len) (1-2*m)/(2*len)])*180/pi;
    costhm=[(2*m-1)/(2*len) (1-2*m)/(2*len)];

%    costhm(find(abs(imag(theta_samp))>=1e-2))=[];
%    costhm(abs(costhm)>1)=[];

    theta_samp(abs(imag(theta_samp))>=1e-2)=[];
    theta_samp=real(theta_samp);

    costhm=cos(theta_samp*pi/180);

    switch sfaf_mode,
       case 1  % file 
          b=interp1(theta,SFAF,theta_samp,'linear');  % linear interpolation
       case 2  % function 
          b=px(theta_samp*pi/180);
    end;

else     % odd number
    m=-(Nsamp-1)/2:(Nsamp-1)/2;
    theta_samp=acos(m/len)*180/pi;
%     theta_samp(abs(imag(theta_samp))>=1e-2)=[];
%     costhm=m/len;
% %    costhm(find(abs(imag(theta_samp))>=1e-2))=[];
% 
     theta_samp=real(theta_samp);
    costhm=cos(theta_samp*pi/180);

    switch sfaf_mode,
       case 1  % file 
          b=interp1(theta,SFAF,theta_samp,'linear');  % linear interpolation
       case 2  % function 
          b=px(theta_samp*pi/180);
    end;
  
end;  

SF_rec=zeros(size(theta));
AF_rec=zeros(size(theta));


for m=1:length(costhm),
   if str=='Space Factor',
      SF_rec=SF_rec+b(m)*sin(pi*len*(cos(theta*pi/180)-costhm(m)))./ ...
             (pi*len*(cos(theta*pi/180)-costhm(m)));
   else
      AF_rec=AF_rec+b(m)*sin(Nel*pi*d*(cos(theta*pi/180)-costhm(m)))./ ...
             (Nel*sin(pi*d*(cos(theta*pi/180)-costhm(m))));
   end;
end;

if str=='Space Factor',
   zs=linspace(-len/2,len/2,length(theta));
   Is=zeros(size(zs));

   for m=1:length(theta_samp),
      Is=Is+1/len*b(m)*exp(-j*2*pi*zs*costhm(m));
   end;

   varargout={zs,Is,SF_rec};

else

   if mod(Nel,2)==0,   % even samples
      za=((1:Nel/2)-0.5)*d;
      za=[-fliplr(za) za];
   else   % odd samples
      za=(-(Nel-1)/2:(Nel-1)/2)*d;
   end;

   Ia=zeros(size(za));
   
   for m=1:Nel,
      Ia(m)=Ia(m)+1/Nel*sum(b.*exp(-j*2*pi*za(m)*costhm));
   end;

   varargout={za,Ia,AF_rec};

end;

end
% this program synthesizes a line source or a linear array using the 
% Woodward-Lawson method

% Written by Marios Gkatzianas	
% Arizona State University, September 2002

function []=wood_main(pos);

% load bal.mat;
[x,map]=imread('bal.tiff');

close all;
warning off;

wood_mode=[];
while ~isreal(wood_mode)|isempty(wood_mode)|((wood_mode~=1)&(wood_mode~=2)),
    wood_mode=input( ['Select method of Woodward-Lawson synthesis\n' ...
                      '******************************************\n' ...
                      '1. Line source\n','2. Linear array\n','->'] );
    if ~isreal(wood_mode)|isempty(wood_mode)|(wood_mode~=1&wood_mode~=2),
       hwomod=msgbox('Specifed number must be either 1 or 2','Invalid choice','custom',x,map,'modal');            
    end;               
end;

output_mode=[];
while ~isreal(output_mode)|isempty(output_mode)|((output_mode~=1)&(output_mode~=2)),
    output_mode=input([ 'Select output method\n','********************\n', ...
          '1. Screen\n', '2. Output file\n', '->']);
    if ~isreal(output_mode)|isempty(output_mode)|(output_mode~=1&output_mode~=2),
       houtmod=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal');
    end;   
end;

if (wood_mode==1),
   str_mod='Space Factor';
else
   str_mod='Array Factor';
end;

outfile=[]; outpath=[];

%%% if (output_mode==2),
%%%	 while isempty(outfile)|(outfile==0), 

%%% don't comment out the next two lines

%	      outfile=input(['Give name of output file to save Array Factor (don''t', ...
%                     ' use quotes)\n'],'s');

%%%       if exist('hout'),
%%%           waitfor(hout);
%%%       end;  
%%%       [outfile,outpath]=uiputfile('*.txt','Select output file');
 
%%%       if isa(outfile,'double')|isa(outpath,'double'),     
%%%           delete(gco);                    
%%%           hout=msgbox('Incorrect output file specification','Invalid choice', ...
%%%                       'custom',x,map,'modal');                 
%%%           set(hout,'interruptible','off','busyaction','queue');
%%%        end;

%%%     end;             
%%% end;

if output_mode==2,
   outfile=input(['Give name of output file to save ' str_mod ' (don''t use quotes)\n'],'s');
   while isempty(outfile),
      outfile=input('Filename must be non-empty. Re-enter name of output file (don''t use quotes).\n','s');
   end;      
   outpath=pwd;
   outfile=[outpath '\' outfile];   
end;


switch wood_mode
    case 1       
	   len=[];
	   while ~isreal(len)|isempty(len)|len<=0,
         len=input('Give length of line source (in wavelengths)\n');
         if ~isreal(len)|isempty(len)|len<=0,
            hlen=msgbox('Length must be real positive','Invalid choice','custom',x,map,'modal');
         end;  
	   end;

	   M=round(len);
      d=[];
      Nel=[];

    case 2
        Nel=[];
        while ~isreal(Nel)|isempty(Nel)|(Nel<=0)|(mod(real(Nel),1)~=0),
           Nel=input('Specify number of elements for linear array\n');            
           if ~isreal(Nel)|isempty(Nel)|Nel<=0|mod(real(Nel),1)~=0,
              hNel=msgbox('Specified number must be positive integer','Invalid choice', ...
                          'custom',x,map,'modal');
           end;   
        end;
        
	     d=[];
        while ~isreal(d)|isempty(d)|(d<=0),
           d=input('Specify spacing between the elements (in wavelengths)\n');
           if ~isreal(d)|isempty(d)|d<=0,
              hd=msgbox('Specified number must be real positive','Invalid choice','custom',x,map,'modal');
           end;
        end;
        if d>0.5,
           disp(['Element spacing is greater that lambda/2. Insufficient spectral resolution ' ...
                    '(excessive aliasing expexted).']); disp(' ');
        end;
	     len=Nel*d;
        M=round(len);
end;  % switch

sfaf_mode=[];
while ~isreal(sfaf_mode)|isempty(sfaf_mode)|((sfaf_mode~=1)&(sfaf_mode~=2)),
   sfaf_mode=input([str_mod,' definition options\n','********************************\n', ...
         '1. Read ', str_mod,' from a file (angles must be in degrees and ' str_mod ...
         ' in linear scale)\n',...
         '(Note: SF function cannot be explicitly expressed as a function of theta. For example: A rectangular pulse.)\n','\n' ...
         '2. Define ', str_mod,' as a function\n', ...
         '(Note: SF function can be explicitly expressed as a function of theta.)\n', ...
         'Choose one of the above\n->']);
   if ~isreal(sfaf_mode)|isempty(sfaf_mode)|(sfaf_mode~=1&sfaf_mode~=2),
      hsfafmod=msgbox('Specified number must be either 1 or 2','Invalid choice','custom',x,map,'modal');
   end;   
end;

switch sfaf_mode,
    case 1  % read SF/AF from file. Angles must be in degrees, abs(SF/AF) in 
	    % absolute units   
%        filename=[];
%        while isempty(filename),	 	
%           filename=input(['Specify the name of the file containing ', ...
%                            str_mod,' (don''t use quotes)\n'],'s');
%        end;

        inpath=[]; infile=[];
               
%%%        while isempty(infile)|isa(infile,'double'),
%%%           if exist('hout'),
%%%              waitfor(hout);
%%%           end;   
%%%           [infile,inpath]=uigetfile('*.txt',['Select text file containing ' str_mod]);
%%%           if isa(inpath,'double')|isa(infile,'double'),  
%%%              hout=msgbox('Incorrect file specification','Invalid choice','custom',x,map,'modal');
%%%           end; 
%%%         end;

        infile=input(['Give name of file to read ' str_mod ' (don''t use quotes)\n', ...
              '(Note: The file should have two lines. The first line should contain the values of theta (in degrees),\n', ...
              'and the second line should have the values of SF.\n', ...
              'Read the file ''sfinput.m'' as an example how to create the SF file.\n' ...
              '(''sf1.m'' for rectangular ',str_mod, ','' sf2.m'' for triangular ', str_mod,'))\n'],'s');
        while isempty(infile),
           infile=input('Filename must be non-empty. Re-enter name of output file (don''t use quotes).\n','s');
        end;      
        inpath=pwd;
        filename=[inpath '\' infile];         
        
        ftid=fopen(filename,'rt');
        if ftid==-1,
           disp('File does not exist or cannot be opened. Program terminated.');
           return;
        else   
           fclose(ftid);
        end;   
           
        SFAF=load(filename);		
%         SFAF=SFAF';
               
        theta=SFAF(1,:);  % degrees
        SFAF=SFAF(2,:);	
	
	     Ntheta=length(theta);	    
	     pstring=[];	          
        
    case 2  

        pstring=[];
        while isempty(pstring),
           pstring=input(['Give ',str_mod,' as a function of x (x=theta).' ...
              ' Use ^ for all powers and don''t use quotes.\n' ...
             '-------------------------------------\n' ... 
             'For example: sin(x); x must be in radians\n' ...
             '-------------------------------------\n'],'s');
           if isempty(pstring),
              hpstr=msgbox('String must be non-empty and contain only a x variable','Invalid choice', ...
                           'custom',x,map,'modal');
           end;              
        end;

        px=inline(vectorize(pstring));

        Ntheta=[];
        while ~isreal(Ntheta)|isempty(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0,
           Ntheta=input(['Specify number of sample points for the angle theta' ...
                         ' for the pattern plots\n', ...
                       '(NOT the sample points for the pattern) [ENTER for default=360]:']);
           Ntheta(isempty(Ntheta))=360;	
           if ~isreal(Ntheta)|isempty(Ntheta)|Ntheta<=0|mod(real(Ntheta),1)~=0,
              hNthe=msgbox('Specified number must be positive integer','Invalid choice', ...
                           'custom',x,map,'modal');
           end;   
        end;

	     theta=linspace(0,360,2*Ntheta);
        SFAF=px(theta*pi/180);

	filename=[];

end;   % end of switch

eo_switch=[];
while ~isreal(eo_switch)|isempty(eo_switch)|((eo_switch~=1)&(eo_switch~=2)&(eo_switch~=3)),
    eo_switch=input( ['Select even and/or odd samples\n' ...
                      '******************************************\n' ...
                      '1. Even samples\n','2. Odd samples\n','3. Even & Odd samples\n','->'] );
    if ~isreal(eo_switch)|isempty(eo_switch)|(eo_switch~=1&eo_switch~=2&eo_switch~=3),
       hwomod=msgbox('Specifed number must be either 1,2 or 3','Invalid choice','custom',x,map,'modal');            
    end;               
end;

[theta_even,theta_samp1,b1,z1,I1,SFAF_rec_even]=wood(str_mod,d,len,Nel,sfaf_mode, ...
     filename,pstring,2*M,Ntheta,theta,SFAF);

[theta_odd,theta_samp2,b2,z2,I2,SFAF_rec_odd]=wood(str_mod,d,len,Nel,sfaf_mode, ...
     filename,pstring,2*M+1,Ntheta,theta,SFAF);

U=(abs(SFAF(theta<=180))).^2;  U(isnan(U))=0;
U_even=abs(SFAF_rec_even(theta<=180)).^2; U_even(isnan(U_even))=0;
U_odd=abs(SFAF_rec_odd(theta<=180)).^2; U_odd(isnan(U_odd))=0;

Prad=sum(U.*sin(theta(theta<=180)*pi/180))*(theta(2)-theta(1))*pi/180;
Prad_even=sum(U_even.*sin(theta(theta<=180)*pi/180))*(theta(2)-theta(1))*pi/180;
Prad_odd=sum(U_odd.*sin(theta(theta<=180)*pi/180))*(theta(2)-theta(1))*pi/180;

D=2*max(U)/Prad;
D_even=2*max(U_even)/Prad_even;
D_odd=2*max(U_odd)/Prad_odd;

       Ltheta=[];
       while isempty(Ltheta)|Ltheta<=-180|Ltheta>=180|mod(Ntheta,1)~=0,
          Ltheta=input(['Specify the angle limits (0 <= theta <= 180 degrees, default angle limits is [0 180]): ']);
          if isempty(Ltheta)
          Ltheta=[0,180];
          break;
          end
       end;    


% Figure 1
% ********
figure(1);

switch eo_switch,
    case 1 % Even samples
switch wood_mode
    case 1 
        set(gcf,'units','pixels','position',pos);
 	     plot(z1,abs(I1),'k','linewidth',2); grid; hold on; 
        set(gca,'fontname','timesnewroman','fontsize',14);
 	     legend(['even (' num2str(2*M),') samples']);
        title(['Current amplitude using the W-L method for the line source',' (L = ',num2str(len),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);
   grid on;
        
    case 2
        set(gcf,'units','pixels','position',pos);
 	     hs1=plot(z1,abs(I1),'k--','marker','s'); grid; hold on; 
        set(gca,'fontname','timesnewroman','fontsize',14);
        set([hs1],'linewidth',2,'markersize',6);
        set(hs1,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'k' 'k' 8});
 	     hl=legend([hs1(1)],['even (' num2str(2*M),') samples']);
        title(['Amplitude distribution using the W-L method for the array (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);        
end;  % switch        
        
    case 2 % Odd samples
switch wood_mode
    case 1 
        set(gcf,'units','pixels','position',pos);
	     plot(z2,abs(I2),'linewidth',2,'color','r'); 
        set(gca,'fontname','timesnewroman','fontsize',14);
 	     legend(['odd (' num2str(2*M+1),') samples']);
        title(['Current amplitude using the W-L method for the line source',' (L = ',num2str(len),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);
   grid on;
        
    case 2
        set(gcf,'units','pixels','position',pos);
	     hs2=plot(z2,abs(I2),'r--','marker','d'); grid;
        set(gca,'fontname','timesnewroman','fontsize',14);
        set([hs2],'linewidth',2,'markersize',6);
        set(hs2,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'r' 'r' 8});
 	     hl=legend([hs2(1)],['odd (' num2str(2*M+1),') samples']);
        title(['Amplitude distribution using the W-L method for the array (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);        
end;  % switch        
    case 3 % Both
switch wood_mode
    case 1 
        set(gcf,'units','pixels','position',pos);
 	     plot(z1,abs(I1),'k','linewidth',2); grid; hold on; 
	     plot(z2,abs(I2),'linewidth',2,'color','r'); 
        set(gca,'fontname','timesnewroman','fontsize',14);
 	     legend(['even (' num2str(2*M),') samples'],['odd (' num2str(2*M+1),') samples']);
        title(['Current amplitude using the W-L method for the line source',' (L = ',num2str(len),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);
   grid on;
        
    case 2
        set(gcf,'units','pixels','position',pos);
 	     hs1=plot(z1,abs(I1),'k--','marker','s'); grid; hold on; 
	     hs2=plot(z2,abs(I2),'r--','marker','d'); grid;
        set(gca,'fontname','timesnewroman','fontsize',14);
        set([hs1 hs2],'linewidth',2,'markersize',6);
        set(hs1,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'k' 'k' 8});
        set(hs2,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'r' 'r' 8});
 	     hl=legend([hs1(1) hs2(1)],['even (' num2str(2*M),') samples'], ...
                   ['odd (' num2str(2*M+1),') samples']);
        title(['Amplitude distribution using the W-L method for the array (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'],'fontname','timesnewroman', ...
       'fontsize',18);        
end;  % switch        
end
grid on;
xlabel('z^{\prime}/\lambda (normalized)','fontname','timesnewroman','fontsize',18, ... 
       'verticalalign','top'); 
ylabel('|I(z^{\prime})|','fontname', 'timesnewroman','fontsize',18,'verticalalign','bottom');

% Figure 2
% ********
figure(2);

switch eo_switch,
    case 1 % Even samples
set(gcf,'units','pixels','position',pos);
plot(theta,abs(SFAF),'linewidth',2); hold on; grid; 
plot(theta_even,abs(SFAF_rec_even),'k','linewidth',2); 
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[Ltheta(1) Ltheta(2)]);
hl=legend('Desired',['even (' num2str(2*M),') samples']);
switch wood_mode
   case 1
       hl=legend('Desired',['Line-source, even (',num2str(2*M),') samples']);
   case 2
       hl=legend('Desired',['Linear array, even (',num2str(2*M),') samples']);
end
list1=findall(hl,'type','line','color','k');
set(list1,'marker','s','markersize',8,'markerfacecolor','k','markeredgecolor','k');
h1=plot(theta_samp1,abs(b1),'marker','s','linestyle','none');
set(h1,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'k' 'k' 8});
switch wood_mode,
    case 1
      title(['Synthesized ' str_mod ' using the W-L method (L = ',num2str(len),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
    case 2
       title(['Synthesized ' str_mod ' using the W-L method (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
end
    case 2 % Odd samples
set(gcf,'units','pixels','position',pos);
plot(theta,abs(SFAF),'linewidth',2); hold on; grid; 
plot(theta_odd,abs(SFAF_rec_odd),'linewidth',2,'color','r');
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[Ltheta(1) Ltheta(2)]);
hl=legend('Desired',['odd (' num2str(2*M+1),') samples']);
switch wood_mode
   case 1
       hl=legend('Desired',['Line-source, odd (',num2str(2*M+1),') samples']);
   case 2
       hl=legend('Desired',['Linear array, odd (',num2str(2*M+1),') samples']);
end
list2=findall(hl,'type','line','color','r');
set(list2,'marker','d','markersize',8,'markerfacecolor','r','markeredgecolor','r');
h2=plot(theta_samp2,abs(b2),'marker','d','linestyle','none');
set(h2,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'r' 'r' 8});
switch wood_mode,
    case 1
      title(['Synthesized ' str_mod ' using the W-L method (L = ',num2str(len),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
    case 2
       title(['Synthesized ' str_mod ' using the W-L method (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
end
    case 3 % Both
set(gcf,'units','pixels','position',pos);
plot(theta,abs(SFAF),'linewidth',2); hold on; grid; 
plot(theta_even,abs(SFAF_rec_even),'k','linewidth',2); 
plot(theta_odd,abs(SFAF_rec_odd),'linewidth',2,'color','r');
set(gca,'fontname','timesnewroman','fontsize',14,'xlim',[Ltheta(1) Ltheta(2)]);
hl=legend('Desired',['even (' num2str(2*M),') samples'], ...
           ['odd (' num2str(2*M+1),') samples']);
switch wood_mode
   case 1
       hl=legend('Desired',['Line-source, even (',num2str(2*M),') samples'],['Line-source, odd (',num2str(2*M+1),') samples']);
   case 2
       hl=legend('Desired',['Linear array, even (',num2str(2*M),') samples'],['Linear array, odd (',num2str(2*M+1),') samples']);
end
list1=findall(hl,'type','line','color','k');
list2=findall(hl,'type','line','color','r');
set(list1,'marker','s','markersize',8,'markerfacecolor','k','markeredgecolor','k');
set(list2,'marker','d','markersize',8,'markerfacecolor','r','markeredgecolor','r');
h1=plot(theta_samp1,abs(b1),'marker','s','linestyle','none');
h2=plot(theta_samp2,abs(b2),'marker','d','linestyle','none');
set(h1,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'k' 'k' 8});
set(h2,{'markerfacecolor','markeredgecolor','markersize'}, ...
             {'r' 'r' 8});
switch wood_mode,
    case 1
      title(['Synthesized ' str_mod ' using the W-L method (L = ',num2str(len),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
    case 2
       title(['Synthesized ' str_mod ' using the W-L method (N = ',num2str(Nel),', d = ',num2str(d),'\lambda)'], ...
       'fontname','timesnewroman','fontsize',18);
end
      
end % end switch

xlabel('\theta (in degrees)','fontname','timesnewroman','fontsize',18, ...
       'verticalalign','top'); 
ylabel(['|' str_mod '|'],'fontname','timesnewroman','fontsize',18,'verticalalign','bottom');
grid on;


if (output_mode==2),

    fid=fopen(outfile,'wt');
     
    if (wood_mode==1),
       
       fidsaf=fopen('SF_wood.txt','wt');
       fprintf(fid,'%% Woodward-Lawson synthesis method for line source\n\n');
       fprintf(fid,['%% Length of line source = ',num2str(len) ' lambda\n\n']); 
       fprintf(fid,'%% Desired Directivity (dimensionless) = %6.2f\n',D);
       fprintf(fid,'%% Directivity (dimensionless) as computed from synthesized Space Factor = %6.2f (even samples)\n',D_even);     
       fprintf(fid,['%%' blanks(34) '>>' blanks(35) '= %6.2f (odd samples)\n\n'],D_odd);         
       fprintf(fid,'%% Desired Directivity (dB) = %6.2f\n',10*log10(D));
       fprintf(fid,'%% Directicity (dB) as computed from synthesized Space Factor = %6.2f (even samples)\n',10*log10(D_even));
       fprintf(fid,['%%' blanks(29) '>>' blanks(29) '= %6.2f (odd samples)\n'],10*log10(D_odd));
                     
       fprintf(fidsaf,'%% Woodward-Lawson synthesis method for line source\n\n');                           
       fprintf(fidsaf,['%% Length of line source = ',num2str(len) ' lambda\n\n']);        
       fprintf(fidsaf,['%%   Theta (deg.)      SF_des      SF_rec (even)' ...
              '      SF_rec (odd)\n']);
       fprintf(fidsaf,'    %8.4f \t     %8.4f \t     %8.4f \t     %8.4f \n',[theta; ...
               SFAF; SFAF_rec_even; SFAF_rec_odd]);
       fclose(fid);	            
    else
       fidsaf=fopen('AF_wood.txt','wt');
                     
       fprintf(fid,'%% Woodward-Lawson synthesis method for linear array\n\n');
       fprintf(fid,['%% Number of array elements = ' num2str(Nel) '\n']);
       fprintf(fid,['%% Spacing = ' num2str(d) ' lambda\n\n']);
       fprintf(fid,'%% Excitation coefficients:\n');              
       fprintf(fid,['%% Real part (' num2str(2*M+1) ' pts) \t Imag part (' num2str(2*M+1) ' pts) \t' ...
                    ' Real part (' num2str(2*M) ' pts) \t Imag part (' num2str(2*M) ' pts)\n']);
       fprintf(fid,'%%   %8.4f \t\t    %8.4f \t\t    %8.4f \t\t   %8.4f\n', ...
               [real(I2); imag(I2); real(I1); imag(I1)]);                    
       fprintf(fid,'\n%% Desired Directivity (dimensionless) = %6.2f\n',D);
       fprintf(fid,'%% Directivity (dimensionless) as computed from synthesized Array Factor = %6.2f (even samples)\n',D_even);     
       fprintf(fid,['%%' blanks(34) '>>' blanks(35) '= %6.2f (odd samples)\n\n'],D_odd);         
       fprintf(fid,'%% Desired Directivity (dB) = %6.2f\n',10*log10(D));
       fprintf(fid,'%% Directicity (dB) as computed from synthesized Array Factor = %6.2f (even samples)\n',10*log10(D_even));
       fprintf(fid,['%%' blanks(29) '>>' blanks(29) '= %6.2f (odd samples)\n'],10*log10(D_odd));
              
       fprintf(fidsaf,'%% Woodward-Lawson synthesis method for linear array\n\n');
       fprintf(fidsaf,['%% Number of array elements = ' num2str(Nel) '\n']);
       fprintf(fidsaf,['%% Spacing = ' num2str(d) ' lambda\n\n']);       
       fprintf(fidsaf,['\n%%   Theta (deg.)      AF_des      AF_rec (even)' ...
              '      AF_rec (odd)\n']);
       fprintf(fidsaf,'    %8.4f \t     %8.4f \t     %8.4f \t     %8.4f \n',[theta; ...
               SFAF; SFAF_rec_even; SFAF_rec_odd]);
       fclose(fid);     
    end;
    
end;

disp(' ');
disp(['Desired Directivity (dimensionless) = ' num2str(D,'%6.2f')]);   
disp(['Directivity (dimensionless) as computed from synthesized ' str_mod ' = ' num2str(D_even,'%6.2f') ...
      ' (even samples)']);
disp([blanks(34) '>>' blanks(34) '= ' num2str(D_odd,'%6.2f') ' (odd samples)']);
disp(' ');
disp(['Desired Directivity (dB) = ' num2str(10*log10(D),'%6.2f')]);   
disp(['Directivity (dB) as computed from synthesized ' str_mod ' = ' num2str(10*log10(D_even),'%6.2f') ...
      ' (even samples)']);
disp([blanks(29) '>>' blanks(28) '= ' num2str(10*log10(D_odd),'%6.2f') ' (odd samples)']);
      

if output_mode==2,
   disp(' ');
   disp(['Output data saved in file ' '''' outfile '''']);
   if wood_mode==1,
      disp([str_mod ' saved in file ''' pwd '\SF_wood.txt''']);
   else
      disp([str_mod ' saved in file ''' pwd '\AF_wood.txt''']);
   end;   
end;   

disp(' ');
disp(['Workspace saved in file ''' pwd '\synth.mat''']);
save synth.mat;

fclose('all');



end

function y=sinc(x)
hh=find(x==0);
x(hh)= 1;
y = sin(pi*x)./(pi*x);
y(hh) = 1;
end



