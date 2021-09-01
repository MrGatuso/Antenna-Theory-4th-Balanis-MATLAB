%***********************    Polarization_Propag    ************************
%
%   This program computes the Poincare sphere angles, and thus the
%   polarization of a plane wave travelling in an infinite homogeneous
%   medium, based on the equations (4-57c), (4-58), and (4-61a).
%
%   Input (assuming a wave travelling in positive z-direction):
%       1. Magnitude and phase of x-component of E-field
%       2. Magnitude and phase of y-component of E-field
%   Output:
%       1. Poincare sphere angles: gamma, delta, epsilon, and tau
%       2. Polarization and sense of rotation
%       3. Axial Ratio
%
%           Written by: Manpreet Saini. June 24, 2010
%**************************************************************************

function Polarization_Propag
clear all
clc

ERR = 1;
save_in_file = input('\nDo you want to save results in file? ([n]/y):','s');
while(ERR ~= 0)
   if (isempty(save_in_file) || save_in_file == 'n')
       fprintf('------------Results will not be saved in file.-----------\n\n');
       ERR=0;
   elseif(save_in_file =='y')
       file_string = input('\nInput the desired output filename: ','s');
       out_file = sprintf('%s.txt',file_string);
       ERR=0;
   else
       save_in_file=input('\nDo you want to save the results in file or not? ([n]/y):','s');
   end;
end

% User input
fprintf('\nAssuming the wave is travelling in positive z direction...');
Ex=input('\nEnter magnitude of the x-component of Electric field:');
Ey=input('\nEnter magnitude of the y-component of Electric field:');
PHIx=input('\nEnter phase angle of the x-component of Electric field (in degrees):');
PHIy=input('\nEnter phase angle of the y-component of Electric field (in degrees):');
if(Ex<0)
    PHIx=PHIx+180;
    Ex=abs(Ex);
end
if(Ey<0)
    PHIy=PHIy+180;
    Ey=abs(Ey);
end
del_PHI=PHIx-PHIy;
%Computations
gamma=atand(Ey/Ex);
delta=PHIy-PHIx;
while(delta < -180)
    delta=delta+360;
end
while(delta > 180)
    delta=delta-360;
end
epsilon=(asind(sind(2*gamma)*sind(delta)))/2;
if(epsilon < -45 || epsilon > 45)
    gamma=atand(Ex/Ey);
    epsilon=(asind(sind(2*gamma)*sind(delta)))/2;
end
tau=(pi-atan(2*Ex*Ey*cosd(del_PHI)/(Ex.^2-Ey.^2 )))/2*180/pi;
if(Ex==0)
    tau=90;
elseif(Ey==0)
    tau=0;
end

AR=cotd(epsilon);
if(AR > 0)
    SOR='Counter-clockwise/Left-hand';
elseif(AR < 0)
    SOR='Clockwise/Right-hand';
end

if(isinf(AR))
    polar='Linear';
    SOR='None';
elseif(abs(AR)-1==-eps)
    polar='Circular';
else
    polar='Elliptical';
end

%Output
if (save_in_file=='y')
    diary(out_file);
end


fprintf('\n\n\t\t\t         Polarization_Propag');
if (save_in_file =='y')
fprintf(['\n\n\t\t\t------------ ' out_file ' ------------']);
end
fprintf('\n\n\t\t\t----------------Input----------------');
fprintf('\nFor a wave travelling in positive z direction...');
fprintf('\nx-component of Electric field: %g at an angle of %g degrees',Ex,PHIx);
fprintf('\ny-component of Electric field: %g at an angle of %g degrees',Ey,PHIy);
fprintf('\n\n\t\t\t----------------Output----------------');
fprintf('\nGamma= %6.4f degrees',gamma);
fprintf('\nDelta= %6.4f degrees',delta);
fprintf('\nEpsilon= %6.4f degrees',epsilon);
fprintf('\nTau= %6.4f degrees',tau);
fprintf('\nAxial Ratio= %6.4f',AR);
fprintf('\nPolarization= %s',polar);
fprintf('\nSense of Rotation= %s\n\n',SOR);


% Stop recording
if (save_in_file=='y')
    diary off;
    fprintf('\nData saved in file %s.txt\n\n', file_string);
end

pause

end