%--------------------------------------------------------------------------
% Program written by SIVASEETHARAMAN PANDI
%--------------------------------------------------------------------------
% This program animates the polarized Electric Field of a wave propagating  
% in +z direction as shown in Fig. 2.23(a). It also animates the polarization 
% ellipse shown in Fig. 2.23(b).

% The user specifies the normalized magnitude and phase of Electric Field
% components. The normalization of the field magnitude is made relative to 
% the amplitude of the largest of the two Electric Field components.

% Apart from the animation, the program also computes the Axial Ratio(AR), 
% sense of rotation and tilt angle relative to Y-axis.

%  Notations used in the program:
% 1) Ex --> x component of the Electric field.
% 2) Ey --> y component of the Electric field.
% 3) Ex(t) --> Instantaneous x component of the Electric field.
% 4) Ey(t) --> Instantaneous y component of the Electric field.
% 5) Ex + Ey --> Total Electric field.
% 6) Ex(t) + Ey(t) --> Instantaneous total Electric field.
% 7) PE --> Polarization Ellipse.
% Input:
% 1. Normalized magnitude of Ex.
% 2. Phase of Ex in degrees.
% 3. Normalized magnitude of Ey.
% 4. Phase of Ey in degrees.
% Output:
% 1. 3D Polarization Animation.
% 2. Polarization Ellipse animation.
% 3. Axial Ratio.
% 4. Sense of rotation.
% 5. Tilt Angle measured relative to Y-axis(in degrees).
% -------------------------------------------------------------------------
function Polarization_Diagram_Ellipse_Animation
clc;
clear all;
close all;
disp('*------------------------- CHAPTER - 4 ------------------------------');
disp('*');
disp('*     The program animates 3D polarization diagram of a rotating E-field vector.');
disp('*');
disp('*     The wave travels in the +z direction.');
disp('*');
disp('*     The program uses exp(jwt) time convention.');
disp('*');
disp('*     The symbols used in the program are as follows:');
disp('*');
disp('*     1) Ex --> x component of the Electric field');
disp('*     2) Ey --> y component of the Electric field');
disp('*     3) Ex(t) --> Instantaneous x component of the Electric field');
disp('*     4) Ey(t) --> Instantaneous y component of the Electric field');
disp('*     5) Ex + Ey --> Total Electric field');
disp('*     6) Ex(t) + Ey(t) --> Instantaneous total Electric field');
disp('*     7) PE --> Polarization Ellipse');
disp('*');
disp('*     Orthogonal field components are represented by: ');
disp('*     Magnitude and phase of:');
disp('*     1) Ex');
disp('*     2) Ey');
disp('*');
disp('*     NOTE: Normalization is made relative to the magnitude of the');
disp('*     largest of the two components');
disp('*');
disp('*     Enter the normalized magnitude of Ex(maximum is 1):');
disp('*');
mag_x=input('|Ex|: ');
disp('*');
disp('*     Enter the normalized magnitude of Ey(maximum is 1): ');
disp('*');
mag_y=input('|Ey|: ');
disp('*');
disp('*     Enter the phase of Ex (in degrees):');
disp('*');
phase_x=input('/_Ex: ');
disp('*');
disp('*     Enter the phase of Ey (in degrees):');
disp('*');
phase_y=input('/_Ey: ');
disp('*');
disp('*     Enter the number of cycles of animation to be run');
disp('*'); 
c=input('--> ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('*');
disp('*     SELECT ONE OF THE OPTIONS BELOW FOR THE DISPLAY OF THE OUTPUT');
disp('*');
disp('*     1) 3D plot of rotating polarization vector(See Fig. 2.23(a))');
disp('*');
disp('*     2) 2D Polarization ellipse(See Fig. 2.23(b))');
disp('*');
disp('*     3) Both');
disp('*');
n=input('*      Your option is: ');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('*');
disp('*     How do you want the output?');
disp('*');
disp('*     1) On the screen');
disp('*');
disp('*     2) To be saved in an output file(Graphics Interchange Format-GIF)');
disp('*');
m=input('*      Your option is: ');
if (n==1 || n==2) 
    disp('*');
    %disp('*     Enter the file name: ');
    disp('*');
    filename=input('Enter the file name --> ','s');
elseif n==3 
    disp('*');
    %disp('*     Enter the filename of 3D rotating polarization vector');
    disp('*');
    filename=input('Enter the filename of 3D rotating polarization vector --> ','s');
    disp('*');
    %disp('*     Enter the filename of polarization ellipse');
    disp('*');
    filename_ellipse=input('Enter the filename of polarization ellipse --> ','s');
end
disp('*');
disp('*     NOTE: The computed data will be displayed once the animation is over');
disp('*');
disp('*     NOTE: Do not disturb the animation while it is being saved as an output GIF file');
disp('*');
if n==1
     [rot,ar,ta]=three_dim(mag_x,mag_y,deg2rad(phase_x),deg2rad(phase_y),m,filename);
elseif n==2
       [rot,ar,ta]=polar(mag_x,mag_y,deg2rad(phase_x),deg2rad(phase_y),m,c,filename);
elseif n==3
        [rot,ar,ta]=three_dim(mag_x,mag_y,deg2rad(phase_x),deg2rad(phase_y),m,filename);
close all;
pause(2);
   [rot,ar,ta]=polar(mag_x,mag_y,deg2rad(phase_x),deg2rad(phase_y),m,c,filename_ellipse);
else
    disp('*');
    disp('*     ENTER THE RIGHT VALUE');
end

function [rot,axialratio,tilt]=three_dim(mag_x,mag_y,phase_x,phase_y,m,filename)
T=5;
z=1e-2:1e-2:600e-2;
wavelen=3;
n=1;
const_z = 2;
z_ind= length(z)*const_z/max(z);
for t=1e-2:5e-2:c*500e-2
    
    E_x = mag_x * cos((2*pi*t/T)-(2*pi*z/wavelen)+phase_x);
    E_y = mag_y * cos((2*pi*t/T)-(2*pi*z/wavelen)+phase_y);
    
    const_x = E_x(const_z);
    const_y = E_y(const_z);
    ar(n)=sqrt((const_x^2)+(const_y^2));
    n=n+1;
    if t<=11e-2
        x(n)=const_x;
        y(n)=const_y;
    end
     
end
for n=1:3
theta(n)=rad2deg(atan(x(n)/y(n)));
if theta(n)<0 && y(n)<0
    theta(n)=theta(n)+180;
elseif theta(n)>0 && y(n)<0
    theta(n)=theta(n)+180;
elseif theta(n)<0 && y(n)>0
    theta(n)=theta(n)+360;
end
end
if theta(3)>theta(2)
    if (theta(3)-theta(2))>180
        rot=-1;
    else
    rot=1;
    end
elseif theta(3)<theta(2)
    if (theta(2)-theta(3))>180
        rot=1;
    else
    rot=-1;
    end
elseif theta(3)==theta(2)
    rot=0;
end
if abs(rad2deg(phase_x-phase_y))==180
    rot=0;
end
axialratio=max(ar)/min(ar);
fname=num2str(axialratio);
n=1;
% tilt=rad2deg(acos(mag_x/max(ar)));
for t=1e-2:8e-2:500e-2
    
    E_x = mag_x * cos((2*pi*t/T)-(2*pi*z/wavelen)+phase_x);
    E_y = mag_y * cos((2*pi*t/T)-(2*pi*z/wavelen)+phase_y);
    
    const_x = E_x(z_ind);
    const_y = E_y(z_ind);
    
    empty=zeros(1,length(E_x));


                fig=plot3(z,E_y,empty,'b-',const_z,const_y,0,'ks',...
                z,empty,E_x,'r-',const_z,0,const_x,'kd',...
                z,E_y,E_x,'k-',const_z,const_y,const_x,'rs',...
                'linewidth',1.5);
         
            delph=phase_x-phase_y;
            tilt=rad2deg((pi/2)-0.5*atan((2*mag_x*mag_y*cos(delph))/(mag_x^2-mag_y^2)));
            h=legend(fig,'Ey','Ey(t)','Ex','Ex(t)','Ex+Ey',...
                'Ex(t)+Ey(t)','Location','northeast','Orientation','vertical');

            for k=1:5:600
                line([z(k) z(k)], [0 E_y(k)],[0 0],'color','b','linewidth',1.25,'HandleVisibility','off');
                line([z(k) z(k)], [0 0],[0 E_x(k)],'color','r','linewidth',1.25,'HandleVisibility','off');
                line([z(k) z(k)], [0 E_y(k)],[0 E_x(k)], 'color','k','linewidth',1.25,'HandleVisibility','off');
            end
        if (rot==1 && mag_x==mag_y)
        title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-8]'];[' Axial ratio: ',fname]; [' Sense of rotation - Left Hand Circularly Polarized']});
        elseif rot==-1 && mag_x == mag_y
        title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-8]'];[' Axial ratio: -',fname] ;[' Sense of rotation - Right Hand Circularly Polarized']});
        elseif rot==1 && mag_x~=mag_y
        title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-8]'];[' Axial ratio: ',fname]; [' Sense of rotation - Left Hand Elliptically Polarized'];[' Tilt Angle: ',num2str(tilt)]});
        elseif rot==-1 && mag_x~=mag_y
        title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-8]'];[' Axial ratio: -',fname]; [' Sense of rotation - Right Hand Elliptically Polarized'];[' Tilt Angle: ',num2str(tilt)]});
        else
            title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-8]'];['Linearly Polarized wave'];[' Tilt Angle: ',num2str(tilt)]});
        end
            axis([0 6 -1 1 -1 1]);
            set(gca,'YDir','reverse');
            grid on;
            xlabel('z --> ','fontsize',8);
            ylabel('y --> ','fontsize',8);
            zlabel('x --> ','fontsize',8);
            drawnow;
            if m==2
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if t == 1e-2
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
            end
            n=n+1;
end
end
function [rot,axialratio,tilt]=polar(mag_x,mag_y,phase_x,phase_y,m,c,filename)
T=5;
z=1e-2:1e-2:600e-2;
wavelen=3;
const_z = 2;
z_ind = length(z)*const_z/max(z);
n=1;

delph=phase_x-phase_y;
tilt=rad2deg((pi/2)-0.5*atan((2*mag_x*mag_y*cos(delph))/(mag_x^2-mag_y^2)));

l1 = line([0 0], [0 0]);
l2 = line([0 0], [0 0]);
l3 = line([0 0], [0 0]);
for t=1e-2:5e-2:c*500e-2
    
    E_x = mag_x * cos((2*pi*t/T)-(2*pi*z/wavelen)+phase_x);
    E_y = mag_y * cos((2*pi*t/T)-(2*pi*z/wavelen)+phase_y);
    
    const_x = E_x(const_z);
    const_y = E_y(const_z);
    ar(n)=sqrt((const_x^2)+(const_y^2));
    n=n+1;
    if t<=11e-2
        x(n)=const_x;
        y(n)=const_y;
    end
     
end
for n=1:3
theta(n)=rad2deg(atan(x(n)/y(n)));
if theta(n)<0 && y(n)<0
    theta(n)=theta(n)+180;
elseif theta(n)>0 && y(n)<0
    theta(n)=theta(n)+180;
elseif theta(n)<0 && y(n)>0
    theta(n)=theta(n)+360;
end
end
if theta(3)>theta(2)
    if (theta(3)-theta(2))>180
        rot=-1;
    else
    rot=1;
    end
elseif theta(3)<theta(2)
    if (theta(2)-theta(3))>180
        rot=1;
    else
    rot=-1;
    end
elseif theta(3)==theta(2)
    rot=0;
end
if abs(rad2deg(phase_x-phase_y))==180
    rot=0;
end
axialratio=max(ar)/min(ar);
fname=num2str(axialratio);
n=1;
% tilt=rad2deg(acos(mag_x/max(ar)));
for t=1e-2:5e-2:c*500e-2
    
    E_x = mag_x * cos((2*pi*t/T)-(2*pi*z/wavelen)+phase_x);
    E_y = mag_y * cos((2*pi*t/T)-(2*pi*z/wavelen)+phase_y);
    
    const_x = E_x(const_z);
    const_y = E_y(const_z);
    delete(l1); delete(l2); delete(l3);
        l1 = line([0 const_y], [0 0],'color','b','linewidth',1.5,'marker','*','HandleVisibility','on');
        l2 = line([0 0], [const_x 0],'color','r','linewidth',1.5,'marker','*','HandleVisibility','on');
        l3 = line([0 const_y], [0 const_x],'color','k','linewidth',1.5,'marker','*','HandleVisibility','on');
        
%         if t==1e-2
        %h=legend('Ey(t)','Ex(t)','Ex(t)+Ey(t)','Polarization Ellipse','Location','NorthEastOutside');
        h=legend('Ey(t)','Ex(t)','Ex(t)+Ey(t)','Location','NorthEastOutside');
        set(h,'FontAngle','italic');
%         end
        
        l4 = line([0 const_y], [0 const_x],'color','g','linewidth',1.35,'marker','.', 'linestyle','none','HandleVisibility','off');
        if (rot==1 && mag_x==mag_y)
        title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-19]'];[' Axial ratio: ',fname]; ['. Sense of rotation - Left Hand Circularly Polarized']});
        elseif rot==-1 && mag_x == mag_y
        title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-19]'];[' Axial ratio: -',fname]; ['. Sense of rotation - Right Hand Circularly Polarized']});
        elseif rot==1 && mag_x~=mag_y
        title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-19]'];[' Axial ratio: ',fname]; ['. Sense of rotation - Left Hand Elliptically Polarized'];[' Tilt Angle: ',num2str(tilt)]});
        elseif rot==-1 && mag_x~=mag_y
        title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-19]'];[' Axial ratio: -',fname]; ['. Sense of rotation - Right Hand Elliptically Polarized'];[' Tilt Angle: ',num2str(tilt)]});
        else
            title({['filename: ',filename,' [Chapter - 4 ',' Figure 4-19]'];['Linearly Polarized wave'];[' Tilt Angle: ',num2str(tilt)]});
        end
        axis equal;
        axis([-1.5 1.5 -1.5 1.5]);
        grid on;
        xlabel('Ey -->');
        ylabel('Ex -->');
        
        % This is where the legend was initially 
        drawnow;
        if m==2
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
            if t == 1e-2;
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
        end
   n=n+1;
   
end
end
s2=('.txt');
filename_txt=strcat(filename,s2);
del_phi=phase_x-phase_y;
major=sqrt(.5*(mag_x^2+mag_y^2+sqrt(mag_x^4+mag_y^4+2*mag_x^2*mag_y^2*cosd(2*del_phi))));
minor=sqrt(.5*(mag_x^2+mag_y^2-sqrt(mag_x^4+mag_y^4+2*mag_x^2*mag_y^2*cosd(2*del_phi))));
if major > minor
    AR=major/minor;
else 
    AR=minor/major;
end
if m==1
    disp('*');
    disp(['*    File Name:',filename]);
    disp('*');
    disp('*--------------- Input Parameters -----------------------------');
    disp('*');
    disp(['*    Magnitude of Ex: ',num2str(mag_x)]);
    disp('*');
    disp(['*    Magnitude of Ey: ',num2str(mag_y)]);
    disp('*');
    disp(['*    Phase of Ex(in degrees): ',num2str(phase_x)]); 
    disp('*');
    disp(['*    Phase of Ey(in degrees): ',num2str(phase_y)]);
    disp('*');
    disp('*    Direction of wave travel: +z');
    disp('*');
    disp('*--------------- Output Parameters ----------------------------');
    disp('*');
    if rot==-1
    disp(['*    Axial ratio calculated from the animation = ',num2str(-ar)]);
    disp('*');
    disp(['*    Axial ratio calculated from the input parameters = ',num2str(-AR)]);
    disp('*');
    else
    disp(['*    Axial ratio calculated from the animation = ',num2str(ar)]);
    disp('*');
    disp(['*    Axial ratio calculated from the input parameters = ',num2str(AR)]);
    disp('*');
    end
    if (rot==1 && mag_x==mag_y)
        disp('*    Sense of rotation: Left Hand Circular Polarization');
        elseif rot==-1 && mag_x == mag_y
        disp('*    Sense of rotation: Right Hand Circular Polarization');
        elseif rot==1 && mag_x~=mag_y
        disp('*    Sense of rotation: Left Hand Elliptical Polarization');
        disp('*');
        disp(['*    Tilt Angle measured relative to Y - axis: ',num2str(ta)]);
        elseif rot==-1 && mag_x~=mag_y
        disp('*    Sense of rotation: Right Hand Elliptical Polarization');
        disp('*');
        disp(['*    Tilt Angle measured relative to Y - axis (in degrees): ',num2str(ta)]);
    else
        disp('*     Linear Polarization');
        disp('*');
        disp(['*    Tilt Angle measured relative to Y - axis (in degrees): ',num2str(ta)]);    
    end
elseif m==2
    fid=fopen(filename_txt,'wt+');
    fprintf(fid,['*********************  ',filename,'  ******************\n']);
    fprintf(fid,'\n');
    fprintf(fid,'--------------- Input Parameters -----------------------------\n');
    fprintf(fid,'\n');
    fprintf(fid,['    Magnitude of Ex: ',num2str(mag_x),' \n']);
    fprintf(fid,'\n');
    fprintf(fid,['    Magnitude of Ey:: ',num2str(mag_y),' \n']);
    fprintf(fid,'\n');
    fprintf(fid,['    Phase of Ex(in degrees): ',num2str(phase_x),' \n']);
    fprintf(fid,'\n');
    fprintf(fid,['    Phase of Ey(in degrees): ',num2str(phase_y),' \n']);
    fprintf(fid,'\n');
    fprintf(fid,'    Direction of wave travel: +z \n');
    fprintf(fid,'\n');
    fprintf(fid,'--------------- Output Parameters ----------------------------\n');
    fprintf(fid,'\n');
    if rot==-1
    fprintf(fid,['    Axial ratio calculated from the animation = ',num2str(-ar),'\n']);
    fprintf(fid,'\n');
    fprintf(fid,['    Axial ratio calculated from the input parameters = ',num2str(-AR),'\n']);
    fprintf(fid,'\n');
    else
    fprintf(fid,['    Axial ratio calculated from the animation = ',num2str(ar),'\n']);
    fprintf(fid,'\n');
    fprintf(fid,['    Axial ratio calculated from the input parameters = ',num2str(AR),'\n']);
    fprintf(fid,'\n');
    end
    if (rot==1 && mag_x==mag_y)
    fprintf(fid,'     Sense of rotation: Left Hand Circular Polarization \n');
    elseif rot==-1 && mag_x == mag_y
    fprintf(fid,'    Sense of rotation: Right Hand Circular Polarization \n');
    elseif rot==1 && mag_x~=mag_y
    fprintf(fid,'    Sense of rotation: Left Hand Elliptical Polarization \n');
    fprintf(fid,'\n');
    fprintf(fid,['    Tilt Angle: ',num2str(ta),'\n']);
    elseif rot==-1 && mag_x~=mag_y
    fprintf(fid,'    Sense of rotation: Right Hand Elliptical Polarization \n');
    fprintf(fid,'\n');
    fprintf(fid,['    Tilt Angle measured relative to Y - axis (in degrees): ',num2str(ta),'\n']);
    else
    fprintf(fid,'    Linear Polarization \n');
    fprintf(fid,'\n');
    fprintf(fid,['    Tilt Angle measured relative to Y - axis (in degrees): ',num2str(ta),'\n']);
    fclose(fid);
    end
end
% pause
end

    
    
        
        


