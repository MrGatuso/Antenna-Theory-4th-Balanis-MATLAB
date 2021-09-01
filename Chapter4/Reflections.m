%--------------------------------------------------------------------------
%                     Infinitesimal Electric Dipole
%
%  This MATLAB M-file calculates the reflection coefficient and the                            
%amplitude patterns of the fields radiated when an infinitesimal electric 
%dipole is placed a height h above an infinite and flat reflecting surface. 
%The radiating element can be either a vertical or horizontal 
%infinitesimal electric dipole. 
%
%The inputs of the program are: 
%    a. Type of surface (PEC, Lossless or Lossy)
%    b. Dipole (Vertical, Horizontal)
%    c. Reflecting surface characteristics
%       [relative permittivity/dielectric constant, conductivity (S/m)] 
%    d. Frequency, f (in GHz) 
%    e. Height, h (in wavelengths), of dipole above reflecting surface
%
%The outputs of the program are: 
%   a. Reflection Coefficient Magnitude for phi=0 and phi=90 degree planes
%   b. Reflection Coefficient Phase for phi=0 and phi=90 degree planes
%   c. Amplitude Patterns for phi=0 and phi=90 degree planes
%
%  Author: Alix Rivera-Albino
%  Date: February 2014
%--------------------------------------------------------------------------

clc
clear all
close all
rep=1;
ii=1;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')

fprintf ('           INFINITESIMAL ELECTRIC DIPOLE \n')
fprintf ('-------------------------------------------------------\n')

while rep==1
    fprintf('-> Select one of the following alternatives \n')
    fprintf('     1. PEC Case\n')
    fprintf('     2. Lossless Case\n')
    fprintf('     3. Lossy Case\n')
    loss= input (' Selected Case = ');
    fprintf ('-> Select one of the following alternatives \n')
    fprintf('     1. Vertical Dipole\n')
    fprintf('     2. Horizontal Dipole\n')
    pol= input (' Selected Case = ');
    
    if loss==1
        sigma=Inf;
        er=1;
        if pol==1
            legd{ii}=char('Vertical: PEC');
        else
            legd{ii}=char('Horizontal: PEC');
        end
    elseif loss==2
        sigma=0;
        er= input('-> Enter the relative permitivity, eps_r = ');
        if pol==1
            legd{ii}=[char([ 'Vertical: Lossless, \it{\epsilon}\rm{_r =} ' num2str(er)])];
        else
            legd{ii}=[char([ 'Horizontal: Lossless, \it{\epsilon}\rm{_r =} ' num2str(er)])];
        end
        
    else
        sigma=input('-> Enter the conductivity in S/m, sigma = ');
        er= input('-> Enter the relative permitivity, eps_r = ');
        if pol==1
            legd{ii}=char(['Vertical: Lossy (\it{\sigma}\rm{ = }' num2str(sigma) ' S/m),\it{\epsilon}\rm{_r =}' num2str(er)]);
        else
          legd{ii}=char(['Horizontal: Lossy (\it{\sigma}\rm{ = }' num2str(sigma) ' S/m),\it{\epsilon}\rm{_r =}' num2str(er)]);
        end
    end
    
    freq= input ('-> Enter the frequency in GHz, f = ')*10^9;
    h=input ('-> Enter the height in wavelengths, h =' );
    lambda= 3e8/freq/sqrt(er);
    h=h*lambda;
    mu0=4*pi*10^-7;
    eps0= 8.854187817e-12;
    omega= 2*pi*freq;
    eps1=er*eps0;
    mu1=mu0;
    eta0=sqrt(mu0/eps0);
    eta1=sqrt(1j*omega*mu1/(sigma+1j*omega*eps1));
    if sigma==Inf
        eta1=0;
    end
    
    %---------------------------------------------
    alpha1= omega*sqrt(mu1*eps1)*(1/2*(sqrt(1+(sigma/(omega*eps1))^2)-1))^1/2;
    beta1=  omega*sqrt(mu1*eps1)*(1/2*(sqrt(1+(sigma/(omega*eps1))^2)+1))^1/2;
    k1=2*pi/lambda;
    if pol==1 %%%%vertical dipole
        thetaid=-90:0.001:90;
        thetai=thetaid*pi/180;
        Gamma(ii,:)= (eta0.*cos(thetai)-eta1.*sqrt(1-(sin(thetai./er)).^2))./(eta0.*cos(thetai)+eta1.*sqrt(1-(sin(thetai./er)).^2));
        E=abs(sin(thetai).*(exp(1j.*k1*h.*cos(thetai))+Gamma(ii,:).*exp(-1j.*k1*h.*cos(thetai))));
        E=E/max(E);
        E=20*log10(E);
        E=E+40;
        E(E<0)=0;
        Er(ii,:)= E;
       Gamma2(ii,:)= Gamma(ii,:);
        Er2=Er;  %%%%phi=0 and phi=90 are the same
    else       %%%%%%horizontal dipole
        %Plane phi =90 degrees
        thetaid=-90:0.001:90;
        thetai=thetaid*pi/180;
        Gamma(ii,:)= -(eta0.*cos(thetai)-eta1.*sqrt(1-(sin(thetai./er)).^2))./(eta0.*cos(thetai)+eta1.*sqrt(1-(sin(thetai./er)).^2));
        E=abs(cos(thetai).*(exp(1j.*k1*h.*cos(thetai))+Gamma(ii,:).*exp(-1j.*k1*h.*cos(thetai))));
        E=E/max(E);
        E=20*log10(E);
        E=E+40;
        E(E<0)=0;
        Er(ii,:)= E;
        
        %Plane phi=0 degrees
        thetaid2=-90:0.001:90;
        thetai2=thetaid*pi/180;
        Gamma2(ii,:)= (eta1.*cos(thetai)-eta0.*sqrt(1-(sin(thetai./er)).^2))./(eta1.*cos(thetai)+eta0.*sqrt(1-(sin(thetai./er)).^2));
        E2=abs(1.*(exp(1j.*k1*h.*cos(thetai))+ Gamma2(ii,:).*exp(-1j.*k1*h.*cos(thetai))));
        E2=E2/max(E2);
        E2=20*log10(E2);
        E2=E2+40;
        E2(E2<0)=0;
        Er2(ii,:)= E2;
        
    end
  
    fprintf('-> Do you want to include another case in the analysis? \n')
    fprintf('   (The results will be displayed in the same plots) \n')
    rep=input ('           (1) Yes (2) No : ');
    
    if rep==1
        ii=ii+1;
    else
    end
    
end

%-------------------------------------------------------------------------
%                             PLOTS
%-------------------------------------------------------------------------
%----------------Magnitude of the reflection coefficient------------------
figure (1)
h=plot(thetaid, abs(Gamma), 'LineWidth', 2);
h_legend=legend(legd,'Location', 'Best');
xlabel('\it{\theta_i} \rm{(degrees)}', 'FontSize', 16)
ylabel('Magnitude of R', 'FontSize', 16)
title(['Reflection Coefficient (\it{\phi}\rm{ = 90}' char(176) ')'],'FontSize',16)
axis([0 90 0 1.2])
box on

%----------------Phase of the reflection coefficient-----------------------
figure (2)
h=plot(thetaid, 180/pi*angle(Gamma), 'LineWidth', 2);
h_legend=legend(legd,'Location', 'Best');
xlabel('\it{\theta_i} \rm{(degrees)}', 'FontSize', 16)
ylabel('Phase of R', 'FontSize', 16)
axis([0 90 -200 200])
title(['Reflection Coefficient (\it{\phi}\rm{ = 90}' char(176) ')'],'FontSize',16)
box on


%----------------Magnitude of the reflection coefficient------------------
figure (3)
h=plot(thetaid, abs(Gamma2), 'LineWidth', 2);
h_legend=legend(legd,'Location', 'Best');
xlabel('\it{\theta_i} \rm{(degrees)}', 'FontSize', 16)
ylabel('Magnitude of R', 'FontSize', 16)
title(['Reflection Coefficient (\it{\phi}\rm{ = 0}' char(176) ')'],'FontSize',16)
axis([0 90 0 1.2])
box on

%----------------Phase of the reflection coefficient-----------------------
figure (4)
h=plot(thetaid, 180/pi*angle(Gamma2), 'LineWidth', 2);
h_legend=legend(legd,'Location', 'Best');
xlabel('\it{\theta_i} \rm{(degrees)}', 'FontSize', 16)
ylabel('Phase of R', 'FontSize', 16)
axis([0 90 -200 200])
title(['Reflection Coefficient (\it{\phi}\rm{ = 0}' char(176) ')'],'FontSize',16)
box on


%------------- Radiation Pattern: Phi= 90 degrees -------------------------
figure (5)

h=polar(repmat(thetai,ii,1).',Er.');
view([90 -90])
h_legend=legend(legd,'Location', 'South');
set (h,'linewidth',2);
hold on
%%%%%Set Limits of POLAR PLOT
t = 0 : .01 : pi/2;
P = polar(t, 40 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
t = 0 : .01 :  pi/2;
P = polar(t, 0 * ones(size(t)));

set(P, 'Visible', 'off')
set(0, 'ShowHiddenHandles', 'on')
set(findobj(gca, 'String', '  40'),'String', '  0'  );  
set(findobj(gca, 'String', '  30'),'String', '-10');  
set(findobj(gca, 'String', '  20'),'String', '-20');  
set(findobj(gca, 'String', '  10'),'String', '-30'); 
set(0, 'ShowHiddenHandles', 'off')

set(findall(gcf, 'String', '330') ,'String', '30');
set(findall(gcf, 'String', '300') ,'String', '60');
set(findall(gcf, 'String', '270') ,'String', '90');
set(findall(gcf, 'String', '240') ,'String', '120');
set(findall(gcf, 'String', '210') ,'String', '150');

pos=get(gca,'Position')-[0 0.05 0 0];
set(gca,'Position',pos)
t=title(['Elevation Plane (\it{\phi}\rm{ = 90}' char(176) ')'],'FontSize',14);
post=get(t,'Position')+[5 0 0];
set(t,'Position', post);

%------------- Radiation Pattern: Phi= 90 degrees -------------------------
figure (6)

h=polar(repmat(thetai,ii,1).',Er2.');
view([90 -90])
h_legend=legend(legd,'Location', 'South');
set (h,'linewidth',2);
hold on
%%%%%Set Limits of POLAR PLOT
t = 0 : .01 : pi/2;
P = polar(t, 40 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
t = 0 : .01 :  pi/2;
P = polar(t, 0 * ones(size(t)));

set(P, 'Visible', 'off')
set(0, 'ShowHiddenHandles', 'on')
set(findobj(gca, 'String', '  40'),'String', '  0'  );  
set(findobj(gca, 'String', '  30'),'String', '-10');  
set(findobj(gca, 'String', '  20'),'String', '-20');  
set(findobj(gca, 'String', '  10'),'String', '-30'); 
set(0, 'ShowHiddenHandles', 'off')

set(findall(gcf, 'String', '330') ,'String', '30');
set(findall(gcf, 'String', '300') ,'String', '60');
set(findall(gcf, 'String', '270') ,'String', '90');
set(findall(gcf, 'String', '240') ,'String', '120');
set(findall(gcf, 'String', '210') ,'String', '150');

pos=get(gca,'Position')-[0 0.05 0 0];
set(gca,'Position',pos)
t=title(['Elevation Plane (\it{\phi}\rm{ = 0}' char(176) ')'],'FontSize',14);
post=get(t,'Position')+[5 0 0];
set(t,'Position', post);
