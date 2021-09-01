% A Matlab computer program, designated as DRA_Analysis_Design, has been 
% developed which performs:
%   •	Analysis  for the four DRAs of Figures 14.57–14.60
%   •	Design for the cylindrical DRA of Figure 14.58
% _______________________________________________________________________
%
%   •  The Analysis part of the program, once the radius a, height h and 
% dielectric constant er are specified, computes the:
% 1. Dominant mode resonant frequencies, of any of the four DRA 
% geometries of Figures 14.57-14.60, based on the modal solution expressions 
% of (14-108a)-(14-114b).
% 2. The resonant frequencies and Qs of the cylindrical resonator of
% Figure 14.58 for all three modes (TE01d, TM01d or HE11d), 
% based on (14-115a)-(14-117b) 
%    •	The Design part of the program is for a cylindrical DRA of 
% Figure 14.58 which, once the:
%   1. Mode (TE01d, TM01d or HE11d)
%   2. Fractional bandwidth (BW, in %)
%   3. VSWR 
%   4. Resonant frequency (fr, in GHz)
% are specified, the program performs the following: 
%   1. Computes the Q of the cavity [based on (14-88a)]. 
%   2. Prompts the user to select a desired dielectric constant er from 
%      the following:
%    •	From a range of values for the TE01d mode.
%    •	Greater than some minimum value for the TM01d mode.
%    •	From a range of values for the HE11d  mode.
% The above ranges or values for the dielectric constant are determined, 
% once the Q is computed using (14-88a), by solving using a 
% nonlinear procedure for:
%    •	TE01d  mode:  (14-115b) for 0.5 < a/h < 5.  
%    •	TM01d  mode (14-116b) for 0.125 < a/h < 5.  
%    •	HE11d  mode (14-117b) for 0.5 < a/h < 5.                                                                         
% Otherwise the stated design specification cannot be met.
% 3.	Finally, once the Q, resonant frequency fr and dielectric constant 
% er have been decided in the previous steps, to complete the design based 
% on the stated specifications, the computer program, using a nonlinear 
% procedure, then determines the dimensions of the cylindrical DRA 
% (radius a and height h, both in cm) by solving simultaneously:
%    •	TE01d  mode:  (14-115a) and (14-115b).
%    •	TM01d  mode: (14-116a) and (14-116b).
%    •  HE11d  mode: (14-117a) and  (14-117b).                                                                       

%%DRA_Analysis_Design
close all
clc
clear all
 
c0=3e8;
    fprintf('--------------------------------------------------------\n')
    fprintf('                      Input                               \n')
    fprintf('--------------------------------------------------------\n')
fprintf('Select the desired geometry: \n')
fprintf('   1. Cubic resonator\n')
fprintf('   2. Cylindrical\n')
fprintf('   3. Hemicylindrical\n')
fprintf('   4. Hemispherical\n')
geo=input ('     Selected DRA: ');
 
if geo==1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           Cubic Resonator                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    a= input('Enter the length a (in cm): ')/100;
    b= input('Enter the width b (in cm): ')/100;
    c= input('Enter the height c (in cm): ')/100;
    er= input('Enter the relative permittivity, er: ');
    
    %TE Modes
    [A,B,C] = ndgrid(1:5,1:5,0:4);
    X=  [A(:),B(:),C(:)];
    m=X(:,1);
    n=X(:,2);
    p=X(:,3);
    
    fTE=c0./(2.*pi.*sqrt(er)).*sqrt((m.*pi./a).^2+(n.*pi./b).^2+((2.*p-1).*(pi./(2.*c))).^2);
    
    %TM modes
    [A,B,C] = ndgrid(0:4,0:4,1:5);
    Y=  [A(:),B(:),C(:)];
    [w,y,z]= find(Y(:,1)==0 & Y(:, 2)==0);
    Y(w,:,:)=[];
    m=Y(:,1);
    n=Y(:,2);
    p=Y(:,3);
    
    fTM=c0./(2.*pi.*sqrt(er)).*sqrt((m.*pi./a).^2+(n.*pi./b).^2+(((2*p-1).*pi./(2.*c))).^2);
    
    %Find the first 5 modes
    
    R= [ones(size(fTE)), X, fTE; zeros(size(fTM)),Y, fTM];
    R= sortrows(R, 5);
    
    modes= R(1:5, :,:);
    
    %Find 5 first TE modes
    R_TE = [ones(size(fTE)), X, fTE];
    R_TE = sortrows(R_TE,5);
    TE_modes = R_TE(1:5, :, :);

    %Find 5 first TM modes
    R_TM = [ones(size(fTM)), Y, fTM];
    R_TM = sortrows(R_TM,5);
    TM_modes = R_TM(1:5, :, :);
    
    %Display 5 first dominant modes
    
    fprintf('\n\n--------------------------------------------------------\n')
    fprintf('                      Output                               \n')
    fprintf('--------------------------------------------------------\n')
        
    
    fprintf('The first five modes for a cubic resonator are: \n');
    
    for jj=1:5
        if modes(jj, 1)==0
            fprintf(' fTE(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz\n',TE_modes(jj,2),TE_modes(jj,3),TE_modes(jj,4),TE_modes(jj,5)/10^9,TM_modes(jj,2),TM_modes(jj,3),TM_modes(jj,4),TM_modes(jj,5)/10^9,modes(jj,2),modes(jj,3),modes(jj,4),modes(jj,5)/10^9) ;
            
        else
            fprintf(' fTE(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz   fTE(%g,%g,%g) = %.2f GHz\n',TE_modes(jj,2),TE_modes(jj,3),TE_modes(jj,4),TE_modes(jj,5)/10^9,TM_modes(jj,2),TM_modes(jj,3),TM_modes(jj,4),TM_modes(jj,5)/10^9,modes(jj,2), modes(jj,3),modes(jj,4),modes(jj,5)/10^9);
        end
        
    end
    
elseif geo==2
    fprintf('Select one of the following:\n');
    fprintf('   1. Analysis\n');
    fprintf('   2. Design\n');
    DoA = input('     Selected one:');
    if DoA == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       Cylindrical Resonator                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    a= input('Enter the radius a (in cm): ')/100;
    h= input('Enter the height h (in cm): ')/100;
    er= input('Enter the relative permittivity, er: ');
    
    %Zeros Bessel
    for mm= 0:5
        zeroXmn(1,mm+1) = fzero(@(x) besselj(mm,x),mm+2);
        y=mm+2;
        for nn= 2:5
            zeroXmn(nn,mm+1)=zeroXmn(nn-1,mm+1) ;
            while zeroXmn(nn,mm+1)<=zeroXmn(nn-1,mm+1)+0.1
                y=y+1;
                zeroXmn(nn,mm+1) = fzero(@(x) besselj(mm,x),y);
            end
        end
    end
    
    %Zeros derivative bessel
    for mm= 0:11
        if mm==0
            zeroXmnp(1,mm+1) = fzero(@(x) 0.5*(besselj(mm-1,x)-besselj(mm+1,x)),3);
            y=7;
        else
            zeroXmnp(1,mm+1) = fzero(@(x) 0.5*(besselj(mm-1,x)-besselj(mm+1,x)),mm+1);
            y=mm+2;
        end
        for nn= 2:5
            zeroXmnp(nn,mm+1)=zeroXmnp(nn-1,mm+1) ;
            while zeroXmnp(nn,mm+1)<=zeroXmnp(nn-1,mm+1)+0.1
                y=y+1;
                zeroXmnp(nn,mm+1) = fzero(@(x) 0.5*(besselj(mm-1,x)-besselj(mm+1,x)),y);
            end
        end
    end
    
    %TE Modes
    [A,B,C] = ndgrid(0:4,1:5,0:4);
    X=  [A(:),B(:),C(:)];
    m=X(:,1);
    n=X(:,2);
    p=X(:,3);
    
    idx = sub2ind(size(zeroXmn), [n], [m+1]);
    Xmn=zeroXmn(idx);
    fTE=c0./(2.*pi.*sqrt(er)).*sqrt((Xmn/a).^2+((2.*p+1).*(pi./(2.*h))).^2);
    
    %TM modes
    [A,B,C] = ndgrid(0:4,1:5,0:4);
    Y=  [A(:),B(:),C(:)];
    m=Y(:,1);
    n=Y(:,2);
    p=Y(:,3);
    
    idx = sub2ind(size(zeroXmnp), [n], [m+1]);
    Xmnp=zeroXmnp(idx);
    fTM=c0./(2.*pi.*sqrt(er)).*sqrt((Xmnp/a).^2+((2.*p+1).*(pi./(2.*h))).^2);
    
    %Find the first 5 modes
    
    R= [ones(size(fTE)), X, fTE; zeros(size(fTM)),Y, fTM];
    R= sortrows(R, 5);
    modes= R(1:5, :,:);
     %Find 5 first TE modes
    R_TE = [ones(size(fTE)), X, fTE];
    R_TE = sortrows(R_TE,5);
    TE_modes = R_TE(1:5, :, :);

    %Find 5 first TM modes
    R_TM = [ones(size(fTM)), Y, fTM];
    R_TM = sortrows(R_TM,5);
    TM_modes = R_TM(1:5, :, :);
    
       
    %Find the resonant frequency and Q
    
    fr_TE10d= c0/(2*pi*a)*(2.327/sqrt(er+1))*(1+0.2123*(a/h)-0.00898*(a/h)^2);
    
    Q_TE10d=0.078192*(er)^1.27*(1+17.31*(h/a)-21.57*(h/a)^2+10.86*(h/a)^3- 1.98*(h/a)^4);
    
    fr_TM01d= c0/(2*pi*a)*(1/sqrt(er+2))*sqrt(3.83^2+(pi*a/(2*h))^2);    
    Q_TM01d=0.008721*(er)^0.888413*exp(0.0397475*er)*(1-(0.3-0.2*(a/h))*((38-er)/28))*(9.498186*(a/h) + 2058.33*(a/h)^(4.322261)*exp(-3.50099*(a/h)));
    
    fr_HE11d= (c0/(2*pi*a))*(6.324/sqrt(er+2))*(0.27 + 0.18*(a/h) + 0.005*((a/h)^2));
    Qrad_HE11d= 0.01007*(er^1.3)*(a/h)*(1 + 100 * exp(-2.05*(0.5*(a/h) - 0.0125 * (a/h)^2)));
    
    
   %Display 5 first dominant modes
    
    fprintf('\n\n--------------------------------------------------------\n')
    fprintf('                      Output                               \n')
    fprintf('--------------------------------------------------------\n')
        
    
    fprintf('The first five modes for a cylindrical resonator are: \n');
    
    for jj=1:5
        if modes(jj, 1)==0
            fprintf(' fTE(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz\n',TE_modes(jj,2),TE_modes(jj,3),TE_modes(jj,4),TE_modes(jj,5)/10^9,TM_modes(jj,2),TM_modes(jj,3),TM_modes(jj,4),TM_modes(jj,5)/10^9,modes(jj,2),modes(jj,3),modes(jj,4),modes(jj,5)/10^9) ;
            
        else
            fprintf(' fTE(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz   fTE(%g,%g,%g) = %.2f GHz\n',TE_modes(jj,2),TE_modes(jj,3),TE_modes(jj,4),TE_modes(jj,5)/10^9,TM_modes(jj,2),TM_modes(jj,3),TM_modes(jj,4),TM_modes(jj,5)/10^9,modes(jj,2), modes(jj,3),modes(jj,4),modes(jj,5)/10^9);
        end
        
    end
    
    fprintf('\n\n-----------------------------------------------------------\n')
    fprintf( '   Mode       Resonant frequency, fr            Qrad \n')
    fprintf( '-----------------------------------------------------------\n')
    fprintf( '  TE_{01d}          %07.4f GHz                %07.4f \n', fr_TE10d/10^9,Q_TE10d)
    fprintf( '  TM_{01d}          %07.4f GHz                %07.4f \n', fr_TM01d/10^9,Q_TM01d)
    fprintf( '  HE_{11d}          %07.4f GHz                %07.4f \n', fr_HE11d/10^9,Qrad_HE11d)
    
    elseif DoA == 2

%%%%%Design of the cylindrical DRA%%%%%
 fprintf('Select one of the following modes:\n')
 fprintf('   1. Transverse Electric (TE_01d)\n');
 fprintf('   2. Transverse Magnetic (TM_01d)\n');
 fprintf('   3. Hybrid (HE_11d)\n');
 Mode = input('   Selected mode:');
 %%%%%%%%%%% TE Mode %%%%%%%%%%%%%%%%%%
 if Mode == 1
 
 prompt = 'Enter the fractional bandwidth (in percent):';
 BW = input(prompt)/100;
 
 prompt = 'Enter VSWR:';
 VSWR = input(prompt);


prompt = 'Enter the resonant frequency (in GHz):';
f_r = input(prompt)*10^9;

 %%%%%Finding the ratio of h/a%%%%%
c = 3*10^8;
z = 0.5:0.01:5;
Q_spec = (VSWR-1)/(BW*sqrt(VSWR));
fprintf('\n--------------------------------------------------------\n')
fprintf('(a)Q(specified)= %.4f \n',Q_spec)
fprintf('--------------------------------------------------------\n')
er_range = ((Q_spec/0.0782)^(1/1.27))*(1+ (17.31*z.^-1) - 21.57*(z.^-2) + 10.86*(z.^-3)-1.98*(z.^-4)).^(-1/1.27);
er_max = max(er_range);
er_min = min(er_range);

fprintf('(b)The dielectric constant should be within the range of: ');
fprintf('%.4f < er < %.4f\n',er_min,er_max)
fprintf('--------------------------------------------------------')
prompt = '\n(c)Enter your dielectric constant:';
er = input(prompt);
Q_rad = 0.078192*(er^1.27)*(1+ (17.31*z.^-1) - 21.57*(z.^-2) + 10.86*(z.^-3)-1.98*(z.^-4)); 
Diff =  Q_rad - Q_spec;
j = 1;
for i = 1:1:450
    if (abs(Diff(i))<0.05 && j<3)
        r = 100*(c/(2*pi*f_r))*(2.327/((er+1)^0.5))*(1+ 0.2123*(z(i))-0.00898*((z(i))^2));
        a = r;
        h = a*(z(i)^-1);
        j = j+1;     
    end
end

    if j== 1
    fprintf('It is impossible to design a DRA with the specified parameters.\n')
    else
    fprintf('\n--------------------------------------------------------\n')
    fprintf('                      Output                               \n')
    fprintf('--------------------------------------------------------\n')
    fprintf( '(d)  a(cm)= %.4f                         h(cm)= %.4f\n',a,h)
    fprintf('--------------------------------------------------------\n')
    end
%%%%%%%%%%%%%%%%% TM Mode %%%%%%%%%%%%%%%%    
elseif Mode == 2
prompt = 'Enter the fractional bandwidth (in percent):';
BW = input(prompt)/100;
 
prompt = 'Enter VSWR:';
VSWR = input(prompt);

prompt = 'Enter the resonant frequency (in GHz):';
f_r = input(prompt)*10^9;

%%%%%Finding the ratio of h/a%%%%%
c = 3*10^8;
z = 0.5:0.01:5;
Q_spec = (VSWR-1)/(BW*sqrt(VSWR));
fprintf('\n--------------------------------------------------------\n')
fprintf('(a)Q(specified)= %.4f \n',Q_spec)
fprintf('--------------------------------------------------------\n')
er_range = 1:0.01:100;

for m = 1:9900
    Q_range = 0.00872*(er_range(m)^0.888)*(exp(0.03975*er_range(m)))*(1- (0.3-0.2*z)*((38-er_range(m))/28)).*(9.498*z + 2058.33*(z.^4.3226).*(exp(-3.501*z)));
        if (max(Q_range) > Q_spec)
            fprintf('(b)Your dielectirc constant should be greater than %.4f\n', er_range(m));
            fprintf('--------------------------------------------------------')
            break;
        end

end


prompt = '\n(c)Enter your dielectric constant:';
er = input(prompt);
Q_rad = 0.008721*(er^0.888)*(exp(0.03975*er))*(1- (0.3-0.2*z)*((38-er)/28)).*(9.498*z + 2058.33*(z.^4.3226).*(exp(-3.501*z))); 
Diff =  Q_rad - Q_spec;
j = 1;
for i = 1:1:450
    if (abs(Diff(i))<1 && j<6)
        r = 100*(c/(2*pi*f_r))*(1/sqrt(er + 2))*sqrt(3.83^2 + (pi*z(i)/2)^2);
        a = r;
        h = a*(z(i)^-1);
        j = j+1;     
    end
end


    fprintf('\n--------------------------------------------------------\n')
    fprintf('                      Output                               \n')
    fprintf('--------------------------------------------------------\n')
    fprintf( '(d)  a(cm)= %.4f                      h(cm)= %.4f\n',a,h)
    fprintf('--------------------------------------------------------\n')
%%%%%%%%%%%%%%%%%%% HE Mode %%%%%%%%%%%%%%%%%% 
elseif Mode == 3
prompt = 'Enter the fractional bandwidth (in percent):';
BW = input(prompt)/100;
 
prompt = 'Enter VSWR:';
VSWR = input(prompt);

prompt = 'Enter the resonant frequency (in GHz):';
f_r = input(prompt)*10^9;

c = 3*10^8;
z = 0.5:0.01:5;
Q_spec = (VSWR-1)/(BW*sqrt(VSWR));
fprintf('\n--------------------------------------------------------\n')
fprintf('(a)Q(specified)= %.4f \n',Q_spec)
fprintf('--------------------------------------------------------\n')

er_range = ((0.01007/Q_spec)*z.*(1 + 100 * exp(-2.05 * (0.5*z - 0.0125*(z.^2))))).^(-1/1.3);
er_max = max(er_range);
er_min = min(er_range);

fprintf('(b)The dielectric constant should be within the range of: ');
fprintf('%.4f < er < %.4f\n',er_min,er_max)
fprintf('--------------------------------------------------------\n')
prompt = '(c)Enter your dielectric constant:';
er = input(prompt);
Q_rad = (0.01007*(er^1.3)* z.*(1 + 100 * exp(-2.05 * (0.5*z - 0.0125*(z.^2))))); 
Diff =  Q_rad - Q_spec;
j = 1;
for i = 1:1:450
    if (abs(Diff(i))<0.1 && j<3)
        r = 100*(c/(2*pi*f_r))*6.324/sqrt(er + 2)*(0.27 + 0.18*z(i) + 0.005 * z(i)^2);
        a = r;
        h = a*(z(i)^-1);
        j = j+1;     
    end
end

    if j== 1
    fprintf('It is impossible to design a DRA with the specified parameters.\n')
    else
    fprintf('\n--------------------------------------------------------\n')
    fprintf('                           Output                               \n')
    fprintf('--------------------------------------------------------\n')
    fprintf( '(d)  a(cm)= %.4f                         h(cm)= %.4f\n',a,h)
    fprintf('--------------------------------------------------------\n')
    end
end

end
elseif geo==3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   Hemicylindrical Resonator                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a= input('Enter the radius a (in cm): ')/100;
    h= input('Enter the length h (in cm): ')/100;
    er= input('Enter the relative permittivity, er: ');
    
    %Zeros Bessel
    for mm= 0:5
        zeroXmn(1,mm+1) = fzero(@(x) besselj(mm,x),mm+2);
        y=mm+2;
        for nn= 2:5
            zeroXmn(nn,mm+1)=zeroXmn(nn-1,mm+1) ;
            while zeroXmn(nn,mm+1)<=zeroXmn(nn-1,mm+1)+0.1
                y=y+1;
                zeroXmn(nn,mm+1) = fzero(@(x) besselj(mm,x),y);
            end
        end
    end
    
    %Zeros derivative bessel
    for mm= 0:11
        if mm==0
            zeroXmnp(1,mm+1) = fzero(@(x) 0.5*(besselj(mm-1,x)-besselj(mm+1,x)),3);
            y=7;
        else
            zeroXmnp(1,mm+1) = fzero(@(x) 0.5*(besselj(mm-1,x)-besselj(mm+1,x)),mm+1);
            y=mm+2;
        end
        for nn= 2:5
            zeroXmnp(nn,mm+1)=zeroXmnp(nn-1,mm+1) ;
            while zeroXmnp(nn,mm+1)<=zeroXmnp(nn-1,mm+1)+0.1
                y=y+1;
                zeroXmnp(nn,mm+1) = fzero(@(x) 0.5*(besselj(mm-1,x)-besselj(mm+1,x)),y);
            end
        end
    end
    
    %TE Modes
    [A,B,C] = ndgrid(0:4,1:5,0:4);
    X=  [A(:),B(:),C(:)];
    m=X(:,1);
    n=X(:,2);
    p=X(:,3);
    
    idx = sub2ind(size(zeroXmn), [n], [m+1]);
    Xmn=zeroXmn(idx);
    fTE=c0./(2.*pi.*sqrt(er)).*sqrt((Xmn/a).^2+(((p)*pi./(2.*h))).^2);
    
    %TM modes
    [A,B,C] = ndgrid(1:5,1:5,0:4);
    Y=  [A(:),B(:),C(:)];
    m=Y(:,1);
    n=Y(:,2);
    p=Y(:,3);
    
    idx = sub2ind(size(zeroXmnp), [n], [m+1]);
    Xmnp=zeroXmnp(idx);
    fTM=c0./(2.*pi.*sqrt(er)).*sqrt((Xmnp/a).^2+(((p)*pi./(2.*h))).^2);
    
    %Find the first 5 modes
    
    R= [ones(size(fTE)), X, fTE; zeros(size(fTM)),Y, fTM];
    R= sortrows(R, 5);
    modes= R(1:5, :,:);
    
    %Find 5 first TE modes
    R_TE = [ones(size(fTE)), X, fTE];
    R_TE = sortrows(R_TE,5);
    TE_modes = R_TE(1:5, :, :);

    %Find 5 first TM modes
    R_TM = [ones(size(fTM)), Y, fTM];
    R_TM = sortrows(R_TM,5);
    TM_modes = R_TM(1:5, :, :);
    
    %Display 5 first dominant modes
    
    fprintf('\n\n--------------------------------------------------------\n')
    fprintf('                      Output                               \n')
    fprintf('--------------------------------------------------------\n')
        
    
    fprintf('The first five modes for a Hemicylindrical resonator are: \n');
    
    for jj=1:5
        if modes(jj, 1)==0
            fprintf(' fTE(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz\n',TE_modes(jj,2),TE_modes(jj,3),TE_modes(jj,4),TE_modes(jj,5)/10^9,TM_modes(jj,2),TM_modes(jj,3),TM_modes(jj,4),TM_modes(jj,5)/10^9,modes(jj,2),modes(jj,3),modes(jj,4),modes(jj,5)/10^9) ;
            
        else
            fprintf(' fTE(%g,%g,%g) = %.2f GHz   fTM(%g,%g,%g) = %.2f GHz   fTE(%g,%g,%g) = %.2f GHz\n',TE_modes(jj,2),TE_modes(jj,3),TE_modes(jj,4),TE_modes(jj,5)/10^9,TM_modes(jj,2),TM_modes(jj,3),TM_modes(jj,4),TM_modes(jj,5)/10^9,modes(jj,2), modes(jj,3),modes(jj,4),modes(jj,5)/10^9);
        end
        
    end
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     Hemispherical Resonator                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a= input('Enter the radius a (in cm): ')/100;
    er= input('Enter the relative permittivity, er: ');
    
    X= [1, 1,1];
    Y= [1, 1,1];
    
    
    %TE Modes
    fTE=2.744*c0/(2*pi*a*sqrt(er));
    
    %TM modes
    
    fTM=4.493*c0/(2*pi*a*sqrt(er));
    
    %Order the modes
    R= [ones(size(fTE)), X, fTE; zeros(size(fTM)),Y, fTM];
    R= sortrows(R, 5);
    modes= R(1, :,:);
    
    %Display the modes
    
    fprintf('\n\n--------------------------------------------------------\n')
    fprintf('                      Output                               \n')
    fprintf('--------------------------------------------------------\n')
    
    fprintf('The dominant mode for a hemispherical resonator is: \n');
    
    if modes(1, 1)==0
        fprintf('   fTM(%g,%g,%g) = %.4f GHz\n', modes(1,2),modes(1,3),modes(1,4),modes(1,5)/10^9) ;
        
    else
        fprintf('   fTE(%g,%g,%g) = %.4f GHz\n', modes(1,2),modes(1,3),modes(1,4),modes(1,5)/10^9)
    end
    
end


