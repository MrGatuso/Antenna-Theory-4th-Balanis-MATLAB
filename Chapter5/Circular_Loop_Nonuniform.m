%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Circular Loop Antenna with Non-Uniform Current                    %
%                                                                         %
% This program computes the characteristics of a circular loop antenna    %
% with non-uniform current.                                               %
%  The inputs are:                                                        %
%       * The radius of the loop a (in wavelengths) and the radius of the %
%         wire b (in wavelengths), or Omega and the circumference (in     %
%         wavelengths)                                                    %
%       * For a conical radiation pattern: the angle theta (in degrees)   %
%         or for an elevation (great circle) radiation pattern: the angle %
%         phi (in degrees)                                                %
%                                                                         %
%   The outputs are:                                                      %
%       * Current distribution                                            %
%       * Input impedance                                                 %
%       * 2-D conical pattern for a given theta angle                     %
%         or 2-D elevation (great circle) pattern for a given elevation   %
%         angle                                                           %
%       * Directivity pattern (in dB)                                     %
%       * Maximum directivity (dimensionless and in dB)                   %
%  By: Alix Rivera-Albino                                                 %
%     May 2014                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
warning('off','all');
set(0,'DefaultAxesFontName', 'Times New Roman');

fprintf('   CIRCULAR LOOP ANTENNA WITH NON-UNIFORM CURRENT\n')
fprintf('-----------------------------------------------------------\n\n')
fprintf(' Inputs: \n')
fprintf('----------\n')

fprintf('Select one of the following options:\n')
fprintf(' (1)Enter the radius of the loop a (in wavelengths), and the radius of the wire b (in wavelengths)\n')
fprintf(' (2)Enter Omega (Omega= 2*ln(2*pi*a/b)) and the circumference in wavelengths (C/lambda = ka)\n');
opt= input ('      Selected option: ');

if opt==1
    a=input('Enter the radius of the loop a (in wavelengths)= ');
    b=input('Enter the radius of the wire b (in wavelengths)= ');
    omega= 2*log(2*pi*a/b);
    C=2*pi*a;
else
    omega=input('Enter the value of Omega = 2*ln(2*pi*a/b))= ');
    C=input('Enter the loop circumference (in wavelengths), C/lambda = ka  = ');
    a = C/(2*pi);
    b= (2*pi*a/exp(omega/2));
end

fprintf('Select one of the following options: \n')
fprintf(' (1) Conical Radiation Pattern \n (2) Elevation (great circle) Radiation Pattern \n');

pat= input ('      Selected option: ');
if pat==1
    thetade= input('Enter the angle theta in degrees = ');
else
    phide= input ('Enter the angle phi in degrees = ');
end

fprintf('Besides the amplitude and directivity patterns, do you want also a plot of the impedance for circumference sizes from 0 to 2.5 wavelengths?\n')
imp= input('(1) Yes (2) No : ' );
t1=tic;
%-------------------------------------------------------------------------%
%                           Variables                                     %
%-------------------------------------------------------------------------%
gammae=0.577215664901532;
tol=1e-6;
eta=120*pi;
n0=2*a/b*exp(-gammae);
JO =log(n0/4.5)/(n0/4.5)*real(-expint(-log(n0/4.5)));
modemax=5;
modemaxc=21;
K0=besselk(0,(1:modemaxc+1).*b./a);
I0=besseli(0,(1:modemaxc+1).*b./a);
ac=zeros(1,modemax);
N=zeros(modemax+1,1);

%-------------------------------------------------------------------------%
%            Impedance for different loop circumferences                  %
%-------------------------------------------------------------------------%

if imp==1
    kav=0:0.01:2.5;
    Z=zeros(1, length(kav));
    psis=2.*pi./log(n0./4.5).*kav./4.5.*(JO+1/3.*(kav./4.5).^2);
    for kk=1:length(kav)
        ka=kav(kk);
        
        N0=(1./pi.*log(8.*a./b)-1./2.*1/pi...
            .*dblquad(@(phi,x) sin(x.*sin(phi)),0, pi, 0,2*ka,tol) ...
            -1./2.*quad(@(x)1j.*besselj(0,x), 0, 2.*ka,tol));
        
        for m=1:6
            Cm=log(4.*m)+gammae-2.*sum((2.*(0:m-1)+1).^-1);
            N(m,1)= (1./pi.*(K0(m).*I0(m)+Cm)-1/2.*1/pi.*...
                dblquad(@(phi,x) sin(x.*sin(phi)-2*m.*phi),0, pi, 0,2*ka,tol) ...
                -1/2.*quad(@(x)1j.*besselj(2.*m, x), 0, 2.*ka,tol));
        end
        a0=ka*N(1,1);
        ac(1)=ka.*(N(2,1)+N0)/2-1/ka.*N(1,1);
        for m=2:5
            ac(m)=ka.*(N(m+1,1)+N(m-1,1))/2-m^2/ka.*N(m,1);
        end
        Z(kk)=1j*pi*eta./ (1./a0+2*sum(1./ac)-psis(kk));
    end
    %IMPEDANCE
    figure
    plot (kav, real(Z), kav,imag(Z),'r--','LineWidth', 2.2)
    xlabel( '\it C/\lambda = ka \rm( circumference in wavelengths)', 'FontSize',12)
    ylabel('Input Resistance, Input Reactance (Ohms)', 'FontSize',12)
    legend('Input Resistance','Input Reactance')
    s=char(['Circular loop Impedance, \Omega = ', num2str(omega)]);
    title (s, 'FontSize',12)
    
end
%-------------------------------------------------------------------------%
%            Impedance for the desired loop circumference                 %
%-------------------------------------------------------------------------%
ka=C;
psis=2.*pi./log(n0./4.5).*ka./4.5.*(JO+1/3.*(ka./4.5).^2);
N0=(1./pi.*log(8.*a./b)-1./2.*1/pi...
    .*dblquad(@(phi,x) sin(x.*sin(phi)),0, pi, 0,2*ka,tol) ...
    -1./2.*quad(@(x)1j.*besselj(0,x), 0, 2.*ka,tol));

for m=1:6
    Cm=log(4.*m)+gammae-2.*sum((2.*(0:m-1)+1).^-1);
    N(m,1)= (1./pi.*(K0(m).*I0(m)+Cm)-1/2.*1/pi.*...
        dblquad(@(phi,x) sin(x.*sin(phi)-2*m.*phi),0, pi, 0,2*ka,tol) ...
        -1/2.*quad(@(x)1j.*besselj(2.*m, x), 0, 2.*ka,tol));
end
a0=ka*N(1,1);
ac(1)=ka.*(N(2,1)+N0)/2-1/ka.*N(1,1);
for m=2:5
    ac(m)=ka.*(N(m+1,1)+N(m-1,1))/2-m^2/ka.*N(m,1);
end
Zc=1j*pi*eta./ (1./a0+2*sum(1./ac)-psis);


%-------------------------------------------------------------------------%
%              Current for the desired loop circumference                 %
%-------------------------------------------------------------------------%
phid=0:360;
I=zeros(1,length(phid));

for m=7:modemaxc+1
    Cm=log(4.*m)+gammae-2.*sum((2.*(0:m-1)+1).^-1);
    N(m,1)= (1./pi.*(K0(m).*I0(m)+Cm)-1/2.*1/pi.*dblquad(@(phi,x) sin(x.*sin(phi)-2*m.*phi),0, pi, 0,2*ka,tol) ...
        -1/2.*quad(@(x)1j.*besselj(2.*m, x), 0, 2.*ka,tol));
end

for m=6:modemaxc
    ac(m)=ka.*(N(m+1,1)+N(m-1,1))/2-m^2/ka.*N(m,1);
end
for kk=1:length(phid)
    I(kk)=1/(1j*pi*eta)*(1/a0+2*sum(cosd((1:modemaxc)*phid(kk))./ac));
end


%CURRENT
figure
plot(phid,abs(I)/10^-3, 'LineWidth', 2)
xlabel( '\it{\phi\prime} \rm (degrees)', 'FontSize',12)
ylabel('Current amplitude (mA)', 'FontSize',12)
s=char(['Circular Loop Current (\Omega = ', num2str(omega), ',\it C/\lambda = ka =\rm ', num2str(C), ' )' ]);
title (s, 'FontSize',12)
xlim([0 360])
set(gca,'XTick',[0:45:360])
ylim([floor(min(abs(I)/10^-3) ) ceil(max(abs(I)/10^-3))]);


figure
plot(phid,(angle(I))*180/pi, 'LineWidth', 2)
xlabel( '\it{\phi\prime} \rm (degrees)', 'FontSize',12)
ylabel('Current phase (degrees)', 'FontSize',12)
s=char(['Circular Loop Current (\Omega = ', num2str(omega), ',\it C/\lambda = ka =\rm ', num2str(C), ' )' ]);
title (s, 'FontSize',12)
axis([0 360 -180 180])
set(gca,'XTick',[0:45:360]);
set(gca,'YTick',[-180:45:180]);
%-------------------------------------------------------------------------%
%                            Electric Field                               %
%-------------------------------------------------------------------------%
In=1j./(pi.*eta).*[1./a0,2./ac];
n=0:modemaxc;

if pat==1
    phid=0:360;
    theta=thetade*pi/180;
    for jj=1:length (phid)
        phi=phid(jj).*pi/180;
        w=C.*sin(theta);
        ETheta(jj)=-eta/2*cot(theta).*sum(n.*(1j).^n.*In.*sin(n.*phi).*besselj(n,w));
        EPhi(jj)=-eta/2*C.*sum((1j).^n.*In.*cos(n.*phi).* 0.5.*(besselj(n-1,w) - besselj(n+1,w)));
    end
    ang=phid*pi/180;
else
    thetad=-180:180;
    phi=phide*pi/180;
    for jj=1:length (thetad)
        theta=thetad(jj).*pi/180;
        w=C.*sin(theta);
        ETheta(jj)=-eta/2*cot(theta).*sum(n.*(1j).^n.*In.*sin(n.*phi).*besselj(n,w));
        EPhi(jj)=-eta/2*C.*sum((1j).^n.*In.*cos(n.*phi).* 0.5.*(besselj(n-1,w) - besselj(n+1,w)));
        ang=thetad*pi/180;
    end
end

E=ETheta +EPhi;
Emag=abs(E);
Emag=Emag/max(max(Emag));
Emag=20*log10(Emag);
Emag=Emag+40;
Emag(Emag<0)=0;
ange=ang;
ange(isnan(Emag))=[];
Emag(isnan(Emag))=[];
%Electric Field
figure
h=polar(ange, Emag);

view([90 -90])
set (h,'linewidth',2);
hold on
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
if pat ==2
    set(0, 'ShowHiddenHandles', 'off')
    set(findall(gcf, 'String', '330') ,'String', '30');
    set(findall(gcf, 'String', '300') ,'String', '60');
    set(findall(gcf, 'String', '270') ,'String', '90');
    set(findall(gcf, 'String', '240') ,'String', '120');
    set(findall(gcf, 'String', '210') ,'String', '150');
end
pos=get(gca,'Position')-[0 0.05 0 0];
set(gca,'Position',pos)
if pat==1
    t=title(['Conical Radiation Pattern (\it{\theta}\rm{ = }' num2str(thetade) char(176), ', \it C/\lambda = ka =\rm ', num2str(C),', \Omega =', num2str(omega),  ')'],'FontSize',12);
else
    t=title(['Elevation (great circle) Radiation Pattern (\it{\phi}\rm{ = }' num2str(phide) char(176) , ',\it C/\lambda = ka =\rm ', num2str(C),', \Omega =', num2str(omega),  ')'],'FontSize',12);
end

post=get(t,'Position')+[5 0 0];
set(t,'Position', post);

%-------------------------------------------------------------------------%
%                 Directivity at desired plane
%-------------------------------------------------------------------------%

U=1./(2.*eta).*(abs(ETheta).^2+abs(EPhi).^2);

int1 = dblquad(@(phi,th) abs(cot(th)).^2.*...
    abs(1.*(1j).^1.*In(2).*sin(1.*phi).*besselj(1,C.*sin(th))+...
    2.*(1j).^2.*In(3).*sin(2.*phi).*besselj(2,C.*sin(th))+...
    3.*(1j).^3.*In(4).*sin(3.*phi).*besselj(3,C.*sin(th))+...
    4.*(1j).^4.*In(5).*sin(4.*phi).*besselj(4,C.*sin(th))+...
    5.*(1j).^5.*In(6).*sin(5.*phi).*besselj(5,C.*sin(th))+...
    6.*(1j).^6.*In(7).*sin(6.*phi).*besselj(6,C.*sin(th))+...
    7.*(1j).^7.*In(8).*sin(7.*phi).*besselj(7,C.*sin(th))+...
    8.*(1j).^8.*In(9).*sin(8.*phi).*besselj(8,C.*sin(th))+...
    9.*(1j).^9.*In(10).*sin(9.*phi).*besselj(9,C.*sin(th))+...
    10.*(1j).^10.*In(11).*sin(10.*phi).*besselj(10,C.*sin(th))+...
    11.*(1j).^11.*In(12).*sin(11.*phi).*besselj(11,C.*sin(th))+...
    12.*(1j).^12.*In(13).*sin(12.*phi).*besselj(12,C.*sin(th))+...
    13.*(1j).^13.*In(14).*sin(13.*phi).*besselj(13,C.*sin(th))+...
    14.*(1j).^14.*In(15).*sin(14.*phi).*besselj(14,C.*sin(th))+...
    15.*(1j).^15.*In(16).*sin(15.*phi).*besselj(15,C.*sin(th))+...
    16.*(1j).^16.*In(17).*sin(16.*phi).*besselj(16,C.*sin(th))+...
    17.*(1j).^17.*In(18).*sin(17.*phi).*besselj(17,C.*sin(th))+...
    18.*(1j).^18.*In(19).*sin(18.*phi).*besselj(18,C.*sin(th))+...
    19.*(1j).^19.*In(20).*sin(19.*phi).*besselj(19,C.*sin(th))+...
    20.*(1j).^20.*In(21).*sin(20.*phi).*besselj(20,C.*sin(th))+...
    21.*(1j).^21.*In(22).*sin(21.*phi).*besselj(21,C.*sin(th))).^2.*...
    sin(th), 0,2.*pi,0 ,pi);


int2 = dblquad(@(phi,th) ...
    abs((1j).^0.*In(1).*cos(0.*phi).*0.5.*(besselj(-1,C.*sin(th))- besselj(1,C.*sin(th)))+...
    (1j).^1.*In(2).*cos(1.*phi).*0.5.*(besselj(0,C.*sin(th)) - besselj(2,C.*sin(th)))+...
    (1j).^2.*In(3).*cos(2.*phi).*0.5.*(besselj(1,C.*sin(th)) - besselj(3,C.*sin(th)))+...
    (1j).^3.*In(4).*cos(3.*phi).*0.5.*(besselj(2,C.*sin(th)) - besselj(4,C.*sin(th)))+...
    (1j).^4.*In(5).*cos(4.*phi).*0.5.*(besselj(3,C.*sin(th)) - besselj(5,C.*sin(th)))+...
    (1j).^5.*In(6).*cos(5.*phi).*0.5.*(besselj(4,C.*sin(th)) - besselj(6,C.*sin(th)))+...
    (1j).^6.*In(7).*cos(6.*phi).*0.5.*(besselj(5,C.*sin(th))- besselj(7,C.*sin(th)))+...
    (1j).^7.*In(8).*cos(7.*phi).*0.5.*(besselj(6,C.*sin(th)) - besselj(8,C.*sin(th)))+...
    (1j).^8.*In(9).*cos(8.*phi).*0.5.*(besselj(7,C.*sin(th)) - besselj(9,C.*sin(th)))+...
    (1j).^9.*In(10).*cos(9.*phi).*0.5.*(besselj(8,C.*sin(th)) - besselj(10,C.*sin(th)))+...
    (1j).^10.*In(11).*cos(10.*phi).*0.5.*(besselj(9,C.*sin(th)) - besselj(11,C.*sin(th)))+...
    (1j).^11.*In(12).*cos(11.*phi).*0.5.*(besselj(10,C.*sin(th)) - besselj(12,C.*sin(th)))+...
    (1j).^12.*In(13).*cos(12.*phi).*0.5.*(besselj(11,C.*sin(th))- besselj(13,C.*sin(th)))+...
    (1j).^13.*In(14).*cos(13.*phi).*0.5.*(besselj(12,C.*sin(th)) - besselj(14,C.*sin(th)))+...
    (1j).^14.*In(15).*cos(14.*phi).*0.5.*(besselj(13,C.*sin(th)) - besselj(15,C.*sin(th)))+...
    (1j).^15.*In(16).*cos(15.*phi).*0.5.*(besselj(14,C.*sin(th)) - besselj(16,C.*sin(th)))+...
    (1j).^16.*In(17).*cos(16.*phi).*0.5.*(besselj(15,C.*sin(th)) - besselj(17,C.*sin(th)))+...
    (1j).^17.*In(18).*cos(17.*phi).*0.5.*(besselj(16,C.*sin(th)) - besselj(18,C.*sin(th)))+...
    (1j).^18.*In(19).*cos(18.*phi).*0.5.*(besselj(17,C.*sin(th))- besselj(19,C.*sin(th)))+...
    (1j).^19.*In(20).*cos(19.*phi).*0.5.*(besselj(18,C.*sin(th)) - besselj(20,C.*sin(th)))+...
    (1j).^20.*In(21).*cos(20.*phi).*0.5.*(besselj(19,C.*sin(th)) - besselj(21,C.*sin(th)))+...
    (1j).^21.*In(22).*cos(21.*phi).*0.5.*(besselj(20,C.*sin(th)) - besselj(22,C.*sin(th)))).^2.*...
    sin(th), 0,2.*pi,0 ,pi);

Prad=eta/8*int1+eta*(C)^2/8*int2;

De= 4*pi*U./Prad;
Dedb= 10*log10(De);
ang(isnan(Dedb))=[];
Dedb(isnan(Dedb))=[];


figure
plot(ang*180/pi, Dedb, 'LineWidth',2);
if pat ==1
    xlabel( '\it \phi \rm(degrees)', 'FontSize',12)
else
    xlabel( '\it \theta \rm(degrees)', 'FontSize',12)
end

ylabel('Directivity (dB)', 'FontSize',12)
if pat==1
    s= char(['Circular Loop Directivity (\it{\theta}\rm{ = }' num2str(thetade) char(176), ', \it C/\lambda = ka =\rm ', num2str(C),', \Omega =', num2str(omega),  ')']);
else
    s=char(['Circular Loop Directivity (\it{\phi}\rm{ = }' num2str(phide) char(176) , ',\it C/\lambda = ka =\rm ', num2str(C),', \Omega =', num2str(omega),  ')']);
end
title (s, 'FontSize',12)
xlim([min(ang*180/pi) max(ang*180/pi)])
if pat==1
    set(gca,'XTick',[0:45:360])
else
    set(gca,'XTick',[-180:45:180])
end
%-------------------------------------------------------------------------%
%                       Maximum Directivity
%-------------------------------------------------------------------------%
thetadl= 0:180;
phidl=0:360;

for kk= 1: length (thetadl)
    theta=thetadl(kk)*pi/180;
    for jj=1:length (phidl)
        phi=phidl(jj).*pi/180;
        w=C.*sin(theta);
        ETheta(jj,kk)=-eta/2*cot(theta).*sum(n.*(1j).^n.*In.*sin(n.*phi).*besselj(n,w));
        EPhi(jj,kk)=-eta/2*C.*sum((1j).^n.*In.*cos(n.*phi).* 0.5.*(besselj(n-1,w) - besselj(n+1,w)));
    end
    ang=phid*pi/180;
end

U=1./(2.*eta).*(abs(ETheta).^2+abs(EPhi).^2);

Dm= 4*pi*U./Prad;
Dmdb= 10*log10(Dm);

%-------------------------------------------------------------------------%
%                       Print Outs
%-------------------------------------------------------------------------%
fprintf('\n\n Outputs: \n')
fprintf('----------\n')
if imag(Zc)>0
    fprintf('The input impedance is: Zin = %.2f + j %.2f\n', real(Zc), imag(Zc))
else
    fprintf('The input impedance is: Zin = %.2f - j %.2f \n', real(Zc), abs(imag(Zc)))
end
if pat==1
    fprintf('The maximum directivity [in the plane theta = %g (degrees)] is: D =  %.2f (dimensionless)=  %.2f dB\n', thetade,max (De), max(Dedb))
else
    fprintf('The maximum directivity [in the plane phi =  %g (degrees)] is: D =  %.2f (dimensionless)=  %.2f dB\n', phide, max(De), max(Dedb))
end

fprintf('The maximum directivity overall is: D = %.2f (dimensionless)=  %.2f dB \n', max(max(Dm)), max(max((Dmdb))))

[pdm, tdm]= find(Dmdb(1:360, 1:181)==max(max((Dmdb))));

for gg=1:length(pdm)
    fprintf('  The maximum directivity overall occurs at: (phi, theta) = ( %g degrees, %g degrees)\n', phidl(pdm(gg)), thetadl(tdm(gg)));
end

t2= toc(t1);
