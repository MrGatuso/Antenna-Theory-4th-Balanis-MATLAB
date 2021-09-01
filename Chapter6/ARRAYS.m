%************************************************************************
%     ARRAYS
%************************************************************************
%    This is a MATLAB based program that computes the
%    radiation characteristics of:
%
%     I.   LINEAR ARRAYS(UNIFORM & BROADSIDE NONUNIFORM)
%     II.  PLANAR ARRAY (BROADSIDE UNIFORM)
%     II.  CIRCULAR ARRAY (BROADSIDE UNIFORM)
%
%     THE UNIFORM AND BROADSIDE NONUNIFORM LINEAR ARRAYS HAVE N ELEMENTS
%     PLACED EQUIDISTANTLY ALONG THE Z-AXIS.
%
%     BROADSIDE PLANAR UNIFORM ARRAY HAS M x N ELEMENTS PLACED
%     EQUIDISTANTLY ALONG THE X AND Y AXES
%
%     OPTION I.  LINEAR ARRAYS
%
%       OPTION A.  UNIFORM
%         ** CHOICES: ARRAY TYPE
%
%         1. BROADSIDE (MAXIMUM ALONG THETA = 0 DEGREES)
%         2. ORDINARY END-FIRE
%            (a). MAXIMUM ALONG THETA = 0 DEGREES
%            (b). MAXIMUM ALONG THETA = 180 DEGREES
%         3. HANSEN WOODYARD END-FIRE
%            (a). MAXIMUM ALONG THETA = 0 DEGREES
%            (b). MAXIMUM ALONG THETA = 180 DEGREES
%         4. SCANNING (MAXIMUM ALONG THETA = THETA_MAX)
%
%         ** INPUT PARAMETERS:
%
%            1. NUMBER OF ELEMENTS
%            2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%            3. DIRECTION OF ARRAY MAXIMUM (THETA_MAX IN DEGREES)
%
%         ** PROGRAM OUTPUT:
%
%            1. NORMALIZED ARRAY FACTOR
%            2. DIRECTIVITY (DIMENSIONLESS & IN dB) USING NUMERICAL
%               INTEGRATION OF THE ARRAY FACTOR
%            3. HALF-POWER BEAMWIDTH (IN DEGREES) USING AN ITERATIVE
%               METHOD (FOR ALL MAXIMA IN THE PATTERN)
%
%       OPTION B.  NONUNIFORM BROADSIDE
%         ** CHOICES: ARRAY TYPE
%
%         1. BINOMIAL
%         2. DOLPH-TSCHEBYSCHEFF
%
%         ** BINOMIAL ARRAY INPUT PARAMETERS:
%
%         1. NUMBER OF ELEMENTS
%         2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%
%         ** DOLPH-TSCHEBYSCHEFF INPUT PARAMETERS:
%
%         1. NUMBER OF ELEMENTS
%         2. SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS)
%         3. SIDE LOBE LEVEL (IN POSITIVE dB; i.e., 30 dB)
%
%         ** PROGRAM OUTPUT:
%
%         1. NORMALIZED EXCITATION COEFFICIENTS (An)
%         2. NORMALIZED ARRAY FACTOR
%         3. DIRECTIVITY (IN dB) USING NUMERICAL INTEGRATION OF THE
%            ARRAY FACTOR
%         4. HALF-POWER BEAMWIDTH (IN DEGREES) USING AN ITERATIVE
%            METHOD (FOR ALL MAXIMA THAT OCCUR IN THE PATTERN)
%
%     OPTION II.  PLANAR ARRAY
%
%         ** ARRAY INPUT PARAMETERS:
%
%         1. NUMBER OF ARRAY ELEMENTS IN X-DIRECTION
%         2. SPACING BETWEEN THE ELEMENTS IN X-DIRECTION (IN
%            WAVELENGTHS)
%         3. NUMBER OF ARRAY ELEMENTS IN Y-DIRECTION
%         4. SPACING BETWEEN THE ELEMENTS IN Y-DIRECTION (IN
%            WAVELENGTHS)
%         5. MAXIMUM BEAM DIRECTION - ANGLE THETA (IN DEGREES)
%         6. MAXIMUM BEAM DIRECTION - ANGLE PHI (IN DEGREES)
%         7. THE PHI ANGLE (IN DEGREES) AT WHICH THE 2-D ANTENNA
%            PATTERN NEEDS TO BE EVALUATED (PHIEVAL IN DEG.)
%         8. TYPE OF ARRAY IN X-DIRECTION
%         9. TYPE OF ARRAY IN Y-DIRECYION
%
%       OPTION A.  UNIFORM
%
%       OPTION B.  BINOMIAL
%
%       OPTION C.  DOLPH-TSCHEBYSCHEFF
%
%         ** INPUT PARAMETERS:
%
%         1. SIDE LOBE LEVEL (IN POSITIVE dB; i.e., 30 dB)
%
%
%      NOTE: ONLY THE ELEVATION ANTENNA PATTERN IS EVALUATED.  THIS
%      ----- PATTERN RANGES FROM THETA=0 DEG. TO THETA=180 DEG. WHEREAS
%            PHI REMAINS CONSTANT AT PHIEVAL.  IF THE PATTERN NEEDS TO
%            BE EVALUATED IN THE BACKSIDE REGION OF THE 2-D ARRAY,
%            THEN THE PROGRAM NEEDS TO BE RE-RUN FOR A NEW PHIEVAL
%            WHICH MUST BE EQUAL TO THE PREVIOUS VALUE PLUS 180 DEG.
%
%         ** PROGRAM OUTPUT:
%
%         1. NORMALIZED ARRAY FACTOR EVALUATED AT A GIVEN ANGLE
%            (PHIEVAL)
%         2. UNIFORM PROGRESSIVE PHASE SHIFT IN X AND Y DIRECTIONS
%            (IN DEGREES)
%         3. DIRECTIVITY (IN dB) USING NUMERICAL INTEGRATION OF THE
%            ARRAY FACTOR
%         4. HALF-POWER BEAMWIDTH (IN DEGREES) FOR ALL MAXIMA THAT
%            OCCUR IN THE ELEVATION PLANE OF THE 2-D ANTENNA PATTERN
%         5. DIRECTIVITY (DIMENSIONLES AND IN dB) USING (6-100) AND (6-101)
%         6. HPBW (IN DEGREES) USING (6-95) AND (6-96)
%         7. BEAM SOLID ANGLE
%
%     OPTION III.  CIRCULAR ARRAY
%         ** CHOICES: EQUATION TYPE
%
%         1. EXPONENTIAL FUNCTION FORM
%         2. BESSEL FUNCTION FORM
%
%
%     ALL THE INPUT PARAMETERS ARE IN TERMS OF THE WAVELENGTH.
%     *****************************************************************
%
%     Written by: Seunghwan Yoon, Arizona State University,  08/12/2002
%
%     Revised by: Seunghwan Yoon, Arizona State University,  12/03/2002
%
%     Updated (to include non-uniform planar arrays) by: Manpreet Saini,
%                                 Arizona State University,  10/27/2009
%
%     *****************************************************************

function[]=ARRAYS
close all;
clc;
option_a=0;


while((option_a~=1)&&(option_a~=2)),
    disp(strvcat('OUTPUT DEVICE OPTION FOR THE OUTPUT PARAMETERS',...
        'OPTION (1):SCREEN','OPTION (2):OUTPUT FILE'));
    option_a=input('OUTPUT DEVICE =');
end
if option_a==2, %OUTPUT FILE
    filename=input('INPUT THE DESIRED OUTPUT FILENAME <in single quotes> = ',...
        's');
    fid=fopen(filename,'wt');
end

theta_low=0;


% CHOICE OF LINEAR OR PLANAR ARRAY OR CIRCULAR ARRAY

option_b=0;
while ((option_b~=1)&&(option_b~=2)&&(option_b~=3)),
    disp(strvcat('LINEAR, PLANAR OR CIRCULAR ARRAY','OPTION (1):LINEAR ARRAY','OPTION (2):PLANAR ARRAY','OPTION (3):CIRCULAR ARRAY'));
    option_b=input('OPTION NUMBER =');
end;
rephi=1;
recalc=1;
while rephi==1;
    
    if option_b==1, % LINEAR ARRAY
        
        patt=1;
        
        theta_up=180;
        disc=181;
        M=1800;
        
        
        phi_low=0;
        phi_up=360;
        MM=disc;
        NN=disc;
        
        
        % CHOICE OF UNIFORM OR NONUNIFORM
        option_c=0;
        while ((option_c~=1)&&(option_c~=2)),
            disp(strvcat('UNIFORM OR NONUNIFORM ARRAY','OPTION (1):UNIFORM ARRAY','OPTION (2):NONUNIFORM ARRAY'));
            option_c=input('OPTION NUMBER =');
        end;
        
        %MM=180;
        %NN=180;
        k=2*pi;
        if patt==1
            theta=linspace(0,pi,M+1);
        else
            theta=linspace(0,pi/2,M+1);
        end
        theta3=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
        phi3=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);
        [THETA,PHI]=meshgrid(theta3,phi3);
        if patt==1
            dtheta=pi/M;
        else
            dtheta=pi/M/2;
        end
        
        if option_c==1, % UNIFORM ARRAY
            % CHOICE OF ARRAY TYPE
            % (BROADSIDE,ORDINARY END-FIRE,HANSEN-WOODYARD END-FIRE,SCANNING)
            option_d=0;
            while ((option_d~=1)&&(option_d~=2)&&(option_d~=3)&&(option_d~=4)),
                disp(strvcat('ARRAY NAMES','OPTION (1):BROADSIDE ARRAY(MAXIMUM ALONG THETA = 90 DEGREES)','OPTION (2):ORDINARY END-FIRE ARRAY','OPTION (3):HANSEN-WOODYARD END-FIRE ARRAY','OPTION (4):SCANNING ARRAY'));
                option_d=input('OPTION NUMBER =');
            end;
            
            if option_d==1, % BROADSIDE
                Nelem=0;
                while (Nelem<1),
                    Nelem=floor(input('NUMBER OF ELEMENTS ='));
                end
                d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
                beta=0;
                psi=k.*d.*cos(theta)+beta;
                psi3=k.*d.*cos(THETA)+beta;
                AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);   % I used sinc function.
                AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);  % sinc(x)=sin(pi*x)/(pi*x).  not, sin(x)/x.
                
            elseif option_d==2, % ORDINARY END-FIRE
                Nelem=0;
                while (Nelem<1),
                    Nelem=floor(input('NUMBER OF ELEMENTS ='));
                end
                d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
                thmax=90;
                while((thmax~=0)&&(thmax~=180)),
                    thmax=input('ANGLE WHERE MAXIMUM OCCURS (THETA = 0 OR 180 DEG.)=');
                end
                if abs(thmax)<eps,beta=-k*d;
                elseif abs(thmax-180)<eps,beta=k*d;
                end
                psi=k*d*cos(theta)+beta;
                psi3=k*d*cos(THETA)+beta;
                AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
                AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);
                
            elseif option_d==3, % HANSEN-WOODYARD END-FIRE
                Nelem=0;
                while (Nelem<1),
                    Nelem=floor(input('NUMBER OF ELEMENTS ='));
                end
                d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
                thmax=90;
                while((thmax~=0)&&(thmax~=180)),
                    thmax=input('ANGLE WHERE MAXIMUM OCCURS (THETA = 0 OR 180 DEG.)=');
                end
                if abs(thmax)<eps,beta=-(k*d+pi/(Nelem-0.00001)); %what is the value(0.00001)??
                elseif abs(thmax-180)<eps,beta=k*d+pi/Nelem;
                end
                psi=k*d*cos(theta)+beta;
                psi3=k*d*cos(THETA)+beta;
                AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
                AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);
                
            elseif option_d==4, % SCANNING
                Nelem=0;
                while (Nelem<1),
                    Nelem=floor(input('NUMBER OF ELEMENTS ='));
                end
                d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
                thmax=input('ANGLE WHERE MAXIMUM OCCURS (THETA = 0 OR 180 DEG.)=');
                beta=-k*d*cos(thmax*pi/180);
                psi=k*d*cos(theta)+beta;
                psi3=k*d*cos(THETA)+beta;
                AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
                AF3=sinc((Nelem.*psi3./2)/pi)./sinc((psi3./2)/pi);
            end;
            
        elseif option_c==2, % NONUNIFORM ARRAY
            % CHOICE OF ARRAY TYPE
            % (BINOMIAL,DOLPH-TSCHEBYSCHEFF)
            option_e=0;
            while ((option_e~=1)&&(option_e~=2)&&(option_e~=3)),
                disp(strvcat('ARRAY NAMES','OPTION (1):BINOMIAL - BROADSIDE','OPTION (2):DOLPH-TSCHEBYSCHEFF - BROADSIDE','OPTION (3):DOLPH-TSCHEBYSCHEFF - SCANNING'));
                option_e=input('OPTION NUMBER =');
            end;
            
            if option_e==1, % BINOMIAL BROADSIDE
                Nelem=0;
                while (Nelem<1),
                    Nelem=floor(input('NUMBER OF ELEMENTS ='));
                end
                d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
                beta=0;
                [AF,Ncoef,Coef]=bin(cos(theta),Nelem,d,beta);
                [AF3,Ncoef3,Coef3]=bin(cos(THETA),Nelem,d,beta);
                % for i=1:M+1
                % [AF,Ncoef,Coef]=bin(X,Nelem,d,beta);
                % end
                % for i=1:MM+1
                % [AF3,Ncoef3,Coef3]=bin(X,Nelem,d,beta);
                % end
                
            elseif option_e==2, % DOLPH-TSCHEBYSCHEFF BROADSIDE
                Nelem=0;
                while (Nelem<1),
                    Nelem=floor(input('NUMBER OF ELEMENTS ='));
                end
                d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
                beta=0;
                RdB=input('SIDE LOBE LEVEL (IN dB; NEGATIVE NUMBER) =');
                RdB=-RdB;
                Ro=10^(RdB/20);
                P=Nelem-1;
                Zo=0.5*((Ro+sqrt(Ro^2-1))^(1/P)+(Ro-sqrt(Ro^2-1))^(1/P));
                dmax=(1/pi)*acos(-1/Zo);
                gamma=cosh((1/P)*log(Ro+sqrt((Ro^2)-1)));
                dopt=1-(acos(1/gamma)/pi);
                [AF,Ncoef,Coef,Zo]=tscheby(cos(theta),Nelem,d,RdB,beta);
                [AF3,Ncoef3,Coef3,Zo]=tscheby(cos(THETA),Nelem,d,RdB,beta);
                [AFmax,Ncoefmax,Coefmax,Zomax]=tscheby(cos(theta),Nelem,dmax,RdB,beta);
                [AFopt,Ncoefopt,Coefopt,Zoopt]=tscheby(cos(theta),Nelem,dopt,RdB,beta);
                % for i=1:M+1;
                % [AF,Ncoef,Coef]=tscheby(cos(theta),Nelem,d,RdB,beta);
                % end
                % for i=1:MM+1;
                % [AF3,Ncoef3,Coef3]=tscheby(cos(THETA),Nelem,d,RdB,beta);
                % end
                          R=10^(RdB/20);
                           f=1+0.636*(2/R*cosh(sqrt((acosh(R))^2-pi^2)))^2;
                      Dir=2*R^2/(1+(R^2-1)*f/(d*Nelem));
        
           HP_uni=acos(-0.443/(Nelem*d))-acos(0.443/(Nelem*d));
           HPT=HP_uni*f*180/pi;
           % DIRECTIVITY FOR MAXIMUM SPACING
           Umax=(abs(AFmax)./max(abs(AFmax))).^2;
        Pradmax=2*pi*sum(Umax.*sin(theta).*dtheta);
        
        Dmax=4*pi*Umax/Pradmax;
        DdBmax=10.*log10(Dmax+eps);
        Domax=max(Dmax);
        DodBmax=max(DdBmax);
        
        % DIRECTIVITY FOR OPTIMUM SPACING
        
        Uopt=(abs(AFopt)./max(abs(AFopt))).^2;
        Pradopt=2*pi*sum(Uopt.*sin(theta).*dtheta);
        
        Dopt=4*pi*Uopt/Pradopt;
        DdBopt=10.*log10(Dopt+eps);
        Doopt=max(Dopt);
        DodBopt=max(DdBopt);
                
%             end
%             Coef=Coef(1:Ncoef);
%             Ncoef=Coef(1:Ncoef)/Coef(Ncoef);
%         end
%        
%         U=(abs(AF)./max(abs(AF))).^2;
%         Prad=2*pi*sum(U.*sin(theta).*dtheta);
%         
%         D=4*pi*U/Prad;
%         DdB=10.*log10(D+eps);
%         Do=max(D);
%         DodB=max(DdB);
%        
%         if Nelem==1||max(AF)<2*min(AF)
%             hp=0;
%             thmax=0;
%         else
%             [hp,thmax]=hpbw(U,M);
%              
%         end
%         No_maxima=length(thmax);
        
        elseif option_e==3, % DOLPH-TSCHEBYSCHEFF SCANNING
                Nelem=0;
                while (Nelem<1),
                    Nelem=floor(input('NUMBER OF ELEMENTS ='));
                end
                d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
                beta=0;
                RdB=input('SIDE LOBE LEVEL (IN dB; NEGATIVE NUMBER) =');
                theta0=input('ENTER THE SCAN ANGLE IN DEGREES =');
                RdB=-RdB;
                Ro=10^(RdB/20);
                P=Nelem-1;
                Zo=0.5*((Ro+sqrt(Ro^2-1))^(1/P)+(Ro-sqrt(Ro^2-1))^(1/P));
                dmax=(1/(pi*(1+abs(cosd(theta0)))))*acos(-1/Zo);
                gamma=cosh((1/P)*log(Ro+sqrt((Ro^2)-1)));
                dopt=dmax;
                [AF,Ncoef,Coef,Zo]=tscheby(cos(theta)-cosd(theta0),Nelem,d,RdB,beta);
                [AF3,Ncoef3,Coef3,Zo]=tscheby(cos(THETA)-cosd(theta0),Nelem,d,RdB,beta);
                [AFmax,Ncoefmax,Coefmax,Zomax]=tscheby(cos(theta)-cosd(theta0),Nelem,dmax,RdB,beta);
                [AFopt,Ncoefopt,Coefopt,Zoopt]=tscheby(cos(theta)-cosd(theta0),Nelem,dopt,RdB,beta);
                % for i=1:M+1;
                % [AF,Ncoef,Coef]=tscheby(cos(theta),Nelem,d,RdB,beta);
                % end
                % for i=1:MM+1;
                % [AF3,Ncoef3,Coef3]=tscheby(cos(THETA),Nelem,d,RdB,beta);
                % end
                          R=10^(RdB/20);
                           f=1+0.636*(2/R*cosh(sqrt((acosh(R))^2-pi^2)))^2;
                      Dir=2*R^2/(1+(R^2-1)*f/(d*Nelem));
        
           HP_uni=acos(-0.443/(Nelem*d))-acos(0.443/(Nelem*d));
           HPT=HP_uni*f*180/pi;
           % DIRECTIVITY FOR MAXIMUM SPACING
           Umax=(abs(AFmax)./max(abs(AFmax))).^2;
        Pradmax=2*pi*sum(Umax.*sin(theta).*dtheta);
        
        Dmax=4*pi*Umax/Pradmax;
        DdBmax=10.*log10(Dmax+eps);
        Domax=max(Dmax);
        DodBmax=max(DdBmax);
        
        % DIRECTIVITY FOR OPTIMUM SPACING
        
        Uopt=(abs(AFopt)./max(abs(AFopt))).^2;
        Pradopt=2*pi*sum(Uopt.*sin(theta).*dtheta);
        
        Dopt=4*pi*Uopt/Pradopt;
        DdBopt=10.*log10(Dopt+eps);
        Doopt=max(Dopt);
        DodBopt=max(DdBopt);
                
            end
            Coef=Coef(1:Ncoef);
            Ncoef=Coef(1:Ncoef)/Coef(Ncoef);
        end
       
        U=(abs(AF)./max(abs(AF))).^2;
        Prad=2*pi*sum(U.*sin(theta).*dtheta);
        
        D=4*pi*U/Prad;
        DdB=10.*log10(D+eps);
        Do=max(D);
        DodB=max(DdB);
       
        if Nelem==1||max(AF)<2*min(AF)
            hp=0;
            thmax=0;
        else
            [hp,thmax]=hpbw(U,M);
             
        end
        No_maxima=length(thmax);

    elseif option_b==2, % PLANAR ARRAY
        
       
        disp(strvcat('PLOTS OPTIONS',...
            'OPTION (1): ASSUMES RADIATION ABOVE AND BELOW XY PLANE ',...
            'OPTION (2): ASSUMES RADIATION ONLY ABOVE XY PLANE'));
        patt=input('PLOTS OPTIONS =');
        
        
        if patt==1
            theta_up=180;
            disc=181;
            M=1800;
        else
            theta_up=90;
            disc=91;
            M=900;
        end
        
        phi_low=0;
        phi_up=360;
        MM=disc;
        NN=disc;
        %M=180;
        %N=180;
        k=2*pi;
        dphi=2*pi/NN;
        
        if patt==1;
            dtheta=pi/MM;
        else
            dtheta=pi/MM/2;
        end
        if recalc==1
            Mx=0;
            
            while (Mx<2),
                Mx=floor(input('NUMBER OF ELEMENTS IN THE X-DIRECTION ='));
            end
            Ny=0;
            while (Ny<2),
                Ny=floor(input('NUMBER OF ELEMENTS IN THE Y-DIRECTION ='));
            end
            
            dx=input('SPACING dx BETWEEN THE ELEMENTS IN X-DIRECTION (IN WAVELENGTHS) =');
            dy=input('SPACING dy BETWEEN THE ELEMENTS IN Y-DIRECTION (IN WAVELENGTHS) =');
            
            % CHOICE OF ARRAY TYPE
            % (UNIFORM,BINOMIAL,DOLPH-TSCHEBYSCHEFF)
            option_fx=0;
            while ((option_fx~=1)&&(option_fx~=2)&&(option_fx~=3)),
                disp(strvcat('ARRAY TYPE IN X-DIRECTION','OPTION (1):UNIFORM','OPTION (2):BINOMIAL','OPTION (3):DOLPH-TSCHEBYSCHEFF'));
                option_fx=input('OPTION NUMBER =');
            end;
            option_fy=0;
            while ((option_fy~=1)&&(option_fy~=2)&&(option_fy~=3)),
                disp(strvcat('ARRAY TYPE IN Y-DIRECTION','OPTION (1):UNIFORM','OPTION (2):BINOMIAL','OPTION (3):DOLPH-TSCHEBYSCHEFF'));
                option_fy=input('OPTION NUMBER =');
            end;
            
            thmax2=input('MAXIMUM BEAM DIRECTION - ANGLE THETA (IN DEGREES - INTEGER #)=');
            phimax2=input('MAXIMUM BEAM DIRECTION - ANGLE PHI (IN DEGREES - INTEGER #)=');
        end
        phieval=input('THE PATTERN IS EVALUATED AT AN ANGLE PHI (IN DEGREES - INTEGER #)=');
        dtor=pi/180;
        betax=-k*dx*sin(dtor*thmax2)*cos(dtor*phimax2);
        betay=-k*dy*sin(dtor*thmax2)*sin(dtor*phimax2);
        
        theta3=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
        phi3=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);
        [THETA,PHI]=meshgrid(theta3,phi3);
        if patt==1
            theta=linspace(0,pi,10*MM+1);
        else
            theta=linspace(0,pi/2,10*MM+1);
        end
        
        phi=phieval*dtor;
        phi180=phieval*dtor-pi;
        
        HP_uni_x=acos(-0.443/(Mx*dx))-acos(0.443/(Mx*dx));          %HPBW for Uniform broadside array using (6-22)
        HP_uni_y=acos(-0.443/(Ny*dy))-acos(0.443/(Ny*dy));
        
        if option_fx==1                                             %Uniform
            AF3x=af10x(THETA,PHI,Mx,dx,betax);
            AFx=af10x(theta,phi,Mx,dx,betax);
            AF180x=af10x(theta,phi180,Mx,dx,betax);
            Coefx=ones(1,Mx);
            th_x0=HP_uni_x;
            D0_x=2*Mx*dx;
        elseif option_fx==2                                         %Binomial
            X=sin(THETA).*cos(PHI);
            [AF3x,Ncoef3x,Coef3x]=bin(X,Mx,dx,betax);
            X=sin(theta).*cos(phi);
            [AFx,Ncoefx,Coefx]=bin(X,Mx,dx,betax);
            X=sin(theta).*cos(phi180);
            [AF180x,Ncoef180x,Coef180x]=bin(X,Mx,dx,betax);
            
            if dx==0.5                      %Cannot calculate HPBW for spacing other than lambda/2
                th_x0=1.06/sqrt(Mx-1);
                D0_x=1.77*sqrt(Mx);
            else
                th_x0=0;
                D0_x=0;
            end
        elseif option_fx==3                                         %Dolph-Tschebyscheff
            if recalc==1
                RxdB=input('SIDE LOBE LEVEL FOR TSCHEBYSCHEFF ARRAY IN X-DIRECTION(IN dB; NEGATIVE NUMBER) =');
                RxdB=-RxdB;
            end
            
            [AF3x,Ncoef3x,Coef3x,Zox]=tscheby(sin(THETA).*cos(PHI),Mx,dx,RxdB,betax);
            [AFx,Ncoefx,Coefx,Zo]=tscheby(sin(theta).*cos(phi),Mx,dx,RxdB,betax);
            [AF180x,Ncoef180x,Coef180x,Zo]=tscheby(sin(theta).*cos(phi180),Mx,dx,RxdB,betax);
            Rx=10^(RxdB/20);
            Px=Mx-1;
            Zx=0.5*((Rx+sqrt(Rx^2-1))^(1/Px)+(Rx-sqrt(Rx^2-1))^(1/Px));
            dmaxx=(1/pi)*acos(-1/Zx);
            HP_uni_xmax=acos(-0.443/(Mx*dmaxx))-acos(0.443/(Mx*dmaxx));
            fx=1+0.636*(2/Rx*cosh(sqrt((acosh(Rx))^2-pi^2)))^2;
            th_x0=fx*HP_uni_x;
            th_x0max=fx*HP_uni_xmax;
            D0_x=2*Rx^2/(1+(Rx^2-1)*fx/(Mx*dx));
            D0_xmax=2*Rx^2/(1+(Rx^2-1)*fx/(Mx*dmaxx));
        end
        
        if option_fy==1                                 %Uniform
            AF3y=af10y(THETA,PHI,Ny,dy,betay);
            AFy=af10y(theta,phi,Ny,dy,betay);
            AF180y=af10y(theta,phi180,Ny,dy,betay);
            Coefy=ones(1,Ny);
            th_y0=HP_uni_y;
            D0_y=2*Ny*dy;
        elseif option_fy==2                             %Binomial
            [AF3y,Ncoef3y,Coef3y]=bin(sin(THETA).*sin(PHI),Ny,dy,betay);
            [AFy,Ncoefy,Coefy]=bin(sin(theta).*sin(phi),Ny,dy,betay);
            [AF180y,Ncoef180y,Coef180y]=bin(sin(theta).*sin(phi180),Ny,dy,betay);
            if dy==0.5                      %Cannot calculate HPBW for spacing other than lambda/2
                th_y0=1.06/sqrt(Ny-1);
                D0_y=1.77*sqrt(Ny);
            else
                th_y0=0;
                D0_y=0;
            end
        elseif option_fy==3                             %Dolph-Tschebyscheff
            if recalc==1
                RydB=input('SIDE LOBE LEVEL FOR TSCHEBYSCHEFF ARRAY IN Y-DIRECTION(IN dB; NEGATIVE NUMBER) =');
                RydB=-RydB;
            end
            [AF3y,Ncoef3y,Coef3y,Zoy]=tscheby(sin(THETA).*sin(PHI),Ny,dy,RydB,betay);
            [AFy,Ncoefy,Coefy,Zo]=tscheby(sin(theta).*sin(phi),Ny,dy,RydB,betay);
            [AF180y,Ncoef180y,Coef180y,Zo]=tscheby(sin(theta).*sin(phi180),Ny,dy,RydB,betay);
            Ry=10^(RydB/20);
            Py=Ny-1;
            Zy=0.5*((Ry+sqrt(Ry^2-1))^(1/Py)+(Ry-sqrt(Ry^2-1))^(1/Py));
            dmaxy=(1/pi)*acos(-1/Zy);
            fy=1+0.636*(2/Ry*cosh(sqrt((acosh(Ry))^2-pi^2)))^2;
            th_y0=fy*HP_uni_y;
            D0_y=2*Ry^2/(1+(Ry^2-1)*fy/(Ny*dy));
            D0_ymax=2*Ry^2/(1+(Ry^2-1)*fy/(Ny*dmaxy));
        end
        
        Ncoefx=Coefx./Coefx(end);             %Normalize w.r.t the edge element
        Ncoefy=Coefy./Coefy(end);
        
        if option_fx==option_fy && Mx==Ny              %Calculate HPBW using (6-95) and (6-96)
            theta_h=th_x0*sec(thmax2*pi/180);
            psi_h=th_x0;
        else
            theta_h=1/sqrt(((cos(thmax2*pi/180))^2)*((th_x0^(-2))*(cos(phimax2*pi/180))^2+...
                (th_y0^(-2))*(sin(phimax2*pi/180))^2));
            psi_h=1/sqrt((th_x0^(-2))*(sin(phimax2*pi/180))^2+(th_y0^(-2))*(cos(phimax2*pi/180))^2);
        end
        theta_h=theta_h*180/pi;
        psi_h=psi_h*180/pi;
        
        omega_A=theta_h*psi_h;                  %Calculate Beam Solid Angle using (6-97)
        
        D0_2=pi*cos(thmax2*pi/180)*D0_x*D0_y;       %Directivity using (6-100)
        if option_fx==3 && option_fy==3
        D0_2max=pi*cos(thmax2*pi/180)*D0_xmax*D0_ymax;
        end
        if option_fx==3 && option_fy~=3
        D0_2max=pi*cos(thmax2*pi/180)*D0_xmax*D0_y;
        end
        if option_fx~=3 && option_fy==3
        D0_2max=pi*cos(thmax2*pi/180)*D0_x*D0_ymax;
        end
        D0_3=32400/omega_A;                         %Directivity using Kraus' formula(6-101)
        
        AF3=AF3x.*AF3y;
        Prad=sum(sum((abs(AF3).^2).*sin(THETA)*dtheta*dphi));
        D1=4*pi*abs(AF3).^2/Prad;
        D1dB=10.*log10(D1);
        Do=4*pi*max(max(abs(AF3).^2))/Prad;
        DodB=10.*log10(Do);
        
        AF=AFx.*AFy;
        D=4*pi*abs(AF).^2/Prad;
        DdB=10.*log10(D);
        U=(abs(AF)./max(abs(AF))).^2;
        
        AF180=AF180x.*AF180y;
        D180=4*pi*abs(AF180).^2/Prad;
        D180dB=10.*log10(D180);
        U180=(abs(AF180)./max(abs(AF180))).^2;
        
        
        
    elseif option_b==3, % CIRCULAR ARRAY
        disp(strvcat('PLOTS OPTIONS',...
            'OPTION (1): ASSUMES RADIATION ABOVE AND BELOW XY PLANE ',...
            'OPTION (2): ASSUMES RADIATION ONLY ABOVE XY PLANE'));
        patt=input('PLOTS OPTIONS =');
        
        
        if patt==1
            theta_up=180;
            disc=181;
            M=1800;
        else
            theta_up=90;
            disc=91;
            M=900;
        end
        
        phi_low=0;
        phi_up=360;
        MM=disc;
        NN=disc;
        option_f=0;
        if recalc==1
            while ((option_f~=1)&&(option_f~=2)),
                disp(strvcat('FUNCTIONS','OPTION (1):EXPONENTIAL FUNCTION FORM','OPTION (2):BESSEL FUNCTION FORM'));
                option_f=input('OPTION NUMBER =');
            end;
            
            %M=180;
            %N=180;
            k=2*pi;
            if patt==1
                dtheta=pi/MM;
            else
                dtheta=pi/MM/2;
            end
            dphi=2*pi/NN;
            
            Nelem=0;
            while (Nelem<2),
                Nelem=floor(input('NUMBER OF ELEMENTS ='));
            end
            rad=input('RADIUS (IN WAVELENGTHS) =');
            if option_f==2, % BESSEL
                m=-1;
                while (m<0),
                    m=floor(input('Residuals[>2]='));
                end
            end
            
            theta0=input('MAXIMUM BEAM DIRECTION - ANGLE THETA (IN DEGREES - INTEGER #)=');
            phi0=input('MAXIMUM BEAM DIRECTION - ANGLE PHI (IN DEGREES - INTEGER #)=');
        end
        
        phieval=input('THE PATTERN IS EVALUATED AT AN ANGLE PHI (IN DEGREES - INTEGER #)=');
        
        if option_f==1, % EXPONENTIAL
            dtor=pi/180;
            theta=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
            phi=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);
            [THETA,PHI]=meshgrid(theta,phi);
            AF3=afc(THETA,PHI,theta0,phi0,Nelem,rad);
            Prad=sum(sum(abs(AF3).^2.*sin(THETA)*dtheta*dphi));
            D1=4*pi*abs(AF3).^2/Prad;
            D1dB=10.*log10(D1);
            Do=4*pi*max(max(abs(AF3).^2))/Prad;
            DodB=10.*log10(Do);
            if patt==1
                theta=linspace(0,pi,10*MM+1);
            else
                theta=linspace(0,pi/2,10*MM+1);
            end
            phi=phieval*dtor;
            AF=afc(theta,phi,theta0,phi0,Nelem,rad);
            D=4*pi*abs(AF).^2/Prad;
            DdB=10.*log10(D);
            U=(abs(AF)./max(abs(AF))).^2;
            AFdB=10.*log10(U);
            
            phi180=phieval*dtor-pi;
            AF180=afc(theta,phi180,theta0,phi0,Nelem,rad);
            D180=4*pi*abs(AF180).^2/Prad;
            D180dB=10.*log10(D180);
            U180=(abs(AF180)./max(abs(AF180))).^2;
            
            
        end
        
        if option_f==2, % BESSEL
            dtor=pi/180;
            theta=linspace(theta_low*pi/180,theta_up*pi/180,MM+1);
            phi=linspace(phi_low*pi/180,phi_up*pi/180,NN+1);
            [THETA,PHI]=meshgrid(theta,phi);
            AF3=afc3(THETA,PHI,theta0,phi0,Nelem,rad,m);
            Prad=sum(sum(abs(AF3).^2.*sin(THETA)*dtheta*dphi));
            D1=4*pi*abs(AF3).^2/Prad;
            D1dB=10.*log10(D1);
            Do=4*pi*max(max(abs(AF3).^2))/Prad;
            DodB=10.*log10(Do);
            
            AF30=afc3(THETA,PHI,theta0,phi0,Nelem,rad,0);
            Prad0=sum(sum(abs(AF30).^2.*sin(THETA)*dtheta*dphi));
            D10=4*pi*abs(AF30).^2/Prad0;
            D1dB0=10.*log10(D10);
            Do0=4*pi*max(max(abs(AF30).^2))/Prad0;
            DodB0=10.*log10(Do0);
            
            if patt==1
                theta=linspace(0,pi,10*MM+1);
            else
                theta=linspace(0,pi/2,10*MM+1);
            end
            phi=phieval*dtor;
            AF=afc3(theta,phi,theta0,phi0,Nelem,rad,m);
            D=4*pi*abs(AF).^2/Prad;
            DdB=10.*log10(D);
            U=(abs(AF)./max(abs(AF))).^2;
            AFdB=10.*log10(U);
            
            AF0=afc3(theta,phi,theta0,phi0,Nelem,rad,0);
            D0=4*pi*abs(AF0).^2/Prad0;
            DdB0=10.*log10(D0);
            U0=(abs(AF0)./max(abs(AF0))).^2;
            AFdB0=10.*log10(U0);
            
            phi180=phieval*dtor-pi;
            AF180=afc3(theta,phi180,theta0,phi0,Nelem,rad,m);
            D180=4*pi*abs(AF180).^2/Prad;
            D180dB=10.*log10(D180);
            U180=(abs(AF180)./max(abs(AF180))).^2;
            
            
        end
    end
    if option_a==2
        diary(filename);
    end
    
    if recalc==1
        scale=0;
        while ((scale~=1)&&(scale~=2)),
            disp(strvcat('DIMENSIONLESS OR dB SCALE IN 3D DIRECTIVITY PLOT','OPTION (1):DIMENSIONLESS SCALE','OPTION (2):dB SCALE'));
            scale=input('OPTION NUMBER =');
        end;
    end
    
    clc
    
    
    
    % Let's go output!!!
    disp(strvcat('********************************************************'));
    disp(strvcat('PROGRAM OUTPUT'));
    disp(strvcat('********************************************************'));
    disp(strvcat('INPUT SPECIFICATION'));
    disp(strvcat('--------------------------------------------------------'));
    
    if option_b==1
        if option_c==1&&option_d==1
            disp(strvcat('UNIFORM BROADSIDE LINEAR ARRAY'));
        end
        if option_c==1&&option_d==2
            disp(strvcat('UNIFORM ORDINARY END-FIRE LINEAR ARRAY'));
        end
        if option_c==1&&option_d==3
            disp(strvcat('UNIFORM HANSEN-WOODYARD END-FIRE LINEAR ARRAY'));
        end
        if option_c==1&&option_d==4
            disp(strvcat('UNIFORM SCANNING LINEAR ARRAY'));
        end
        if option_c==2&&option_e==1
            disp(strvcat('NONUNIFORM BINOMIAL(BROADSIDE) LINEAR ARRAY'));
        end
        if option_c==2&&option_e==2
            disp(strvcat('NONUNIFORM DOLPH-TSCHEBYSCHEFF(BROADSIDE) LINEAR ARRAY'));
        end
        if option_c==2&&option_e==3
            disp(strvcat('NONUNIFORM DOLPH-TSCHEBYSCHEFF(SCANNING) LINEAR ARRAY'));
        end
        disp(['NUMBER OF ARRAY ELEMENTS = ',num2str(Nelem)]);
        disp(['SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS) = ',num2str(d)]);
        if option_c==1&&option_d~=1
            disp(['MAXIMUM NEEDS TO OCCUR AT = ',num2str(thmax)]);
        end
        if option_c==2&&option_e==2
            disp(['SIDE LOBE LEVEL (IN dB) = ',num2str(-RdB)]);
        end
        if option_c==2&&option_e==3
            disp(['SIDE LOBE LEVEL (IN dB) = ',num2str(-RdB)]);
            disp(['SCAN ANGLE (IN DEGREES) = ',num2str(theta0)]);
        end
        disp('OUTPUT CHARACTERISTICS OF THE ARRAY');
        disp('--------------------------------------------------------');
        
       
        disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
        disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
        disp(['NUMBER OF MAXIMA BETWEEN 0 AND 180 DEGREES = ',num2str(No_maxima)]);
        for i=1:No_maxima;
            disp(['HPBW using 0.5*AF^2 FOR MAXIMUM # =',num2str(i),'    ',num2str(hp(i,:)),' degrees     THMAX = ',num2str(thmax(i)),' degrees']);
        end
        if option_b==1&&option_c==2
            
            disp('TOTAL EXCITATION COEFFICIENTS FOR THE ARRAY DESIGN');
            disp(Coef);
            disp('NORMALIZED TOTAL EXCITATION COEFFICIENTS (RELATIVE TO EDGE)');
            Ncoef=Coef/Coef(round(Nelem/2));
            disp(Ncoef);
            disp('NORMALIZED TOTAL EXCITATION COEFFICIENTS (RELATIVE TO CENTER)');
            Ncoef=Ncoef./Ncoef(1);
            disp(Ncoef);
        end
        
        if option_c==2  
            if option_e==2
        fprintf('\nBEAM BROADENING FACTOR (BBF), f = %6.4f \n', f);
        fprintf('DIRECTIVITY using (6-79) = %6.4f dB \n', 10*log10(Dir));
        fprintf('DIRECTIVITY using (6-79) = %6.4f dimensionless\n', Dir);
        fprintf('HALF POWER BEAMWIDTH for uniform array, using (6-22a)= %6.4f degrees\n', HP_uni*180/pi);
        fprintf('HALF POWER BEAMWIDTH for Dolph-Tschebyscheff array computed using\n')
        fprintf('      HPBW(uniform)*f(BBF) = %6.4f degrees\n', HPT);
     disp('--------------------------------------------------------');
     disp(['OPTIMUM SPACING (dopt IN WAVELENGTHS) = ',num2str(dmax)]);
     disp(['MAXIMUM SPACING (dmax IN WAVELENGTHS) = ',num2str(dmax)]);
     disp('MAXIMUM SPACING dmax = OPTIMUM SPACING dopt');
     %disp(['DIRECTIVITY (d=dopt) = ',num2str(DodBopt),' dB']);
        %disp(['DIRECTIVITY (d=dopt)= ',num2str(Doopt),' dimensionless']);
     disp(['DIRECTIVITY (d=dmax=dopt) = ',num2str(DodBmax),' dB']);
        disp(['DIRECTIVITY (d=dmax=dopt)= ',num2str(Domax),' dimensionless']);
        disp('--------------------------------------------------------');
        end
        end
        
    if option_c==2  
            if option_e==3
        fprintf('\nBEAM BROADENING FACTOR (BBF), f = %6.4f \n', f);
        fprintf('DIRECTIVITY using (6-79) = %6.4f dB \n', 10*log10(Dir));
        fprintf('DIRECTIVITY using (6-79) = %6.4f dimensionless\n', Dir);
        fprintf('HALF POWER BEAMWIDTH for uniform array, using (6-22a)= %6.4f degrees\n', HP_uni*180/pi);
        fprintf('HALF POWER BEAMWIDTH for Dolph-Tschebyscheff array computed using\n')
        fprintf('      HPBW(uniform)*f(BBF) = %6.4f degrees\n', HPT);
     disp('--------------------------------------------------------');
     disp(['OPTIMUM SPACING (dopt IN WAVELENGTHS) = ',num2str(dmax)]);
     disp(['MAXIMUM SPACING (dmax IN WAVELENGTHS) = ',num2str(dmax)]);
     disp('MAXIMUM SPACING dmax = OPTIMUM SPACING dopt');
     %disp(['DIRECTIVITY (d=dopt) = ',num2str(DodBopt),' dB']);
        %disp(['DIRECTIVITY (d=dopt)= ',num2str(Doopt),' dimensionless']);
     disp(['DIRECTIVITY (d=dmax=dopt) = ',num2str(DodBmax),' dB']);
        disp(['DIRECTIVITY (d=dmax=dopt)= ',num2str(Domax),' dimensionless']);
        disp('--------------------------------------------------------');
        end
        end    
        
    end
   
    if option_b==2
        disp('PLANAR ARRAY WITH');
        if option_fx==1
        disp(strvcat('UNIFORM DISTRIBUTION IN X-DIRECTION'));
        elseif option_fx==2
        disp(strvcat('BINOMIAL DISTRIBUTION IN X-DIRECTION'));  
        elseif option_fx==3
        disp(strvcat('DOLPH-TSCHEBYSCHEFF DISTRIBUTION IN X-DIRECTION'));
        end
        if option_fy==1
        disp(strvcat('UNIFORM DISTRIBUTION IN Y-DIRECTION'));
        elseif option_fy==2
        disp(strvcat('BINOMIAL DISTRIBUTION IN Y-DIRECTION'));  
        elseif option_fy==3
        disp(strvcat('DOLPH-TSCHEBYSCHEFF DISTRIBUTION IN Y-DIRECTION'));
        end
        disp(['NUMBER OF ELEMENTS IN X-DIRECTION = ',num2str(Mx)]);
        disp(['SPACING BETWEEN THE ELEMENTS IN X-DIRECTION (IN WAVELENGTHS) = ',num2str(dx)]);
        disp(['NUMBER OF ELEMENTS IN Y-DIRECTION = ',num2str(Ny)]);
        disp(['SPACING BETWEEN THE ELEMENTS IN Y-DIRECTION (IN WAVELENGTHS) = ',num2str(dy)]);
        disp(['MAXIMUM BEAM DIRECTION - THETA (IN DEGREES) = ',num2str(thmax2)]);
        disp(['MAXIMUM BEAM DIRECTION - PHI (IN DEGREES) = ',num2str(phimax2)]);
        disp(['THE 2D ANTENNA PATTERN IS EVALUATED AT AN ANGLE PHI (IN DEGREES) = ',num2str(phieval)]);
        if option_fx==3
            disp(['SIDE LOBE LEVEL IN THE X-DIRECTION (IN dB) =', num2str(-RxdB)]);
        end
        if option_fy==3
            disp(['SIDE LOBE LEVEL IN THE Y-DIRECTION (IN dB) =', num2str(-RydB)]);
        end
        disp(strvcat('OUTPUT CHARACTERISTICS OF THE ARRAY'));
        disp(strvcat('--------------------------------------------------------'));
        if option_fx==3
            disp(['z_0 IN THE X-DIRECTION =', num2str(Zox)]);
            disp(['BEAM BROADENING FACTOR IN THE X-DIRECTION =', num2str(fx)]);
        end
        if option_fy==3
            
            disp(['z_0 IN THE Y-DIRECTION =', num2str(Zoy)]);
            disp(['BEAM BROADENING FACTOR IN THE Y-DIRECTION =', num2str(fy)]);
            %disp(['Maximum Directivity in the Y-direction =',num2str(D0_ymax)]);
        end
        disp(['PROGRESSIVE PHASE SHIFT IN X-DIRECTION = ',num2str(betax/dtor),' degrees']);
        disp(['PROGRESSIVE PHASE SHIFT IN Y-DIRECTION = ',num2str(betay/dtor),' degrees']);
        disp('DIRECTIVITY BASED ONLY ON THE FIELDS ABOVE THE XY-PLANE')
        if patt==1
            disp(['DIRECTIVITY = ',num2str(DodB+10*log10(2)),' dB']);
            disp(['DIRECTIVITY = ',num2str(Do*2),' dimensionless']);
        else
            disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
            disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
        end
        
        %     if max(AF3)<2*min(AF3)
        %         hp=0;
        %         thmax=0;
        %     else
        [hp,thmax]=hpbw(U,10*MM);
        %     end
        No_maxima=length(thmax);
        disp('DIRECTIVITY BASED ON THE FIELDS ABOVE AND BELOW THE XY-PLANE')
        
        if patt==1
            disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
            disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
        else
            disp(['DIRECTIVITY = ',num2str(DodB-10*log10(2)),' dB']);
            disp(['DIRECTIVITY = ',num2str(Do/2),' dimensionless']);
        end
        
        disp(['EVALUATION PLANE: NUMBER OF MAXIMA BETWEEN 0 AND 180 DEGREES = ',num2str(No_maxima)]);
        for i=1:No_maxima;
            disp(['HPBW FOR MAXIMUM #',num2str(i),'    ',num2str(hp(i,:)),' degrees     THMAX = ',num2str(thmax(i)),' degrees']);
        end
        fprintf('\nDIRECTIVITY using (6-100) (dimensionless):  %6.4f \n', D0_2);
        fprintf('DIRECTIVITY using (6-100) (dB):  %6.4f dB \n', 10*log10(D0_2));
        fprintf('DIRECTIVITY using (6-101) (dimensionless):  %6.4f \n', D0_3);
        fprintf('DIRECTIVITY using (6-101) (dB):  %6.4f dB \n', 10*log10(D0_3));
        fprintf('HALF POWER BEAMWIDTH, using (6-95) : %6.4f degrees\n', theta_h);
        fprintf('HALF POWER BEAMWIDTH, using (6-96): %6.4f degrees\n', psi_h);
        fprintf('BEAM SOLID ANGLE: %6.4f square degrees\n',omega_A);
        disp('--------------------------------------------------------');
        if (option_fx==3 && option_fy==3)
            disp(['MAXIMUM SPACING IN THE X-DIRECTION =',num2str(dmaxx)]);
            disp(['MAXIMUM SPACING IN THE Y-DIRECTION =',num2str(dmaxy)]);
            disp(['DIRECTIVITY (dx = dmax) =',num2str(D0_xmax),' dimensionless']);
            disp(['DIRECTIVITY (dy = dmax) =',num2str(D0_ymax),' dimensionless']);
            D0_xmaxdb=10*log10(D0_xmax);
            D0_ymaxdb=10*log10(D0_ymax);
            D0_2maxdb=10*log10(D0_2max);
            disp(['DIRECTIVITY IN THE X-DIRECTION (dx = dmax) =',num2str(D0_xmaxdb),' dB']);
            disp(['DIRECTIVITY IN THE Y-DIRECTION (dy = dmax) =',num2str(D0_ymaxdb),' dB']);
            disp(['OVERALL DIRECTIVITY WITH MAXIMUM SPACING IN BOTH X & Y DIRECTIONS =',num2str(D0_2maxdb),' dB']);
        elseif (option_fx==3 && option_fy~=3)
            disp(['MAXIMUM SPACING IN THE X-DIRECTION =',num2str(dmaxx)]);
            disp(['DIRECTIVITY (dx = dmax) =',num2str(D0_xmax),' dimensionless']);
            D0_xmaxdb=10*log10(D0_xmax);
            D0_2maxdb=10*log10(D0_2max);
            disp(['DIRECTIVITY IN THE X-DIRECTION (dx = dmax) =',num2str(D0_xmaxdb),' dB']);
            disp(['OVERALL DIRECTIVITY WITH MAXIMUM SPACING IN X-DIRECTION =',num2str(D0_2maxdb),' dB']);
        elseif (option_fx~=3 && option_fy==3)
            disp(['MAXIMUM SPACING IN THE Y-DIRECTION =',num2str(dmaxy)]);
            disp(['DIRECTIVITY (dy = dmax) =',num2str(D0_ymax),' dimensionless']);
            D0_ymaxdb=10*log10(D0_ymax);
            D0_2maxdb=10*log10(D0_2max);
            disp(['DIRECTIVITY IN THE Y_DIRECTION (dy = dmax) =',num2str(D0_ymaxdb),' dB']);
            disp(['OVERALL DIRECTIVITY WITH MAXIMUM SPACING IN Y-DIRECTION =',num2str(D0_2maxdb),' dB']);
        end
        disp('--------------------------------------------------------');
        if (option_fx==2 && dx~=0.5) || (option_fy==2 && dy~=0.5)
            fprintf('HPBW COULD NOT BE COMPUTED FOR BINOMIAL DISTRIBUTION FOR SPACING OTHER THAN HALF-WAVELENGTH');
        end
    end
    
    if option_b==3
        disp(strvcat('UNIFORM CIRCULAR ARRAY'));
        disp(['NUMBER OF ELEMENTS = ',num2str(Nelem)]);
        disp(['RADIUS (IN WAVELENGTHS) = ',num2str(rad)]);
        disp(['MAXIMUM BEAM DIRECTION - THETA (IN DEGREES) = ',num2str(theta0)]);
        disp(['MAXIMUM BEAM DIRECTION - PHI (IN DEGREES) = ',num2str(phi0)]);
        disp(['THE 2D ANTENNA PATTERN IS EVALUATED AT AN ANGLE PHI (IN DEGREES) = ',num2str(phieval)]);
        disp(strvcat('OUTPUT CHARACTERISTICS OF THE ARRAY'));
        disp(strvcat('--------------------------------------------------------'));
        disp('DIRECTIVITY BASED ONLY ON THE FIELDS ABOVE THE XY-PLANE')
        if patt==1
            disp(['DIRECTIVITY = ',num2str(DodB+10*log10(2)),' dB']);
            disp(['DIRECTIVITY = ',num2str(Do*2),' dimensionless']);
        else
            disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
            disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
        end
        if max(AF3)<2*min(AF3)
            hp=0;
            thmax=0;
        else
            [hp,theta0]=hpbw(U,10*MM);
        end
        No_maxima=length(theta0);
        disp('DIRECTIVITY BASED ON THE FIELDS ABOVE AND BELOW THE XY-PLANE')
        if patt==1
            disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
            disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
        else
            disp(['DIRECTIVITY = ',num2str(DodB-10*log10(2)),' dB']);
            disp(['DIRECTIVITY = ',num2str(Do/2),' dimensionless']);
        end
        disp(['EVALUATION PLANE: NUMBER OF MAXIMA BETWEEN 0 AND 180 DEGREES = ',num2str(No_maxima)]);
        for i=1:No_maxima;
            disp(['HPBW FOR MAXIMUM #',num2str(i),'    ',num2str(hp(i,:)),' degrees     THMAX = ',num2str(theta0(i)),' degrees']);
        end
    end
    
    
    
    
    
    
    
    
    
    
    disp(' *** NOTE:');
    disp(' THE NORMALIZED ARRAY FACTOR (in dB) IS STORED IN');
    disp(' AN OUTPUT FILE CALLED  ............   ArrFac.dat');
    disp(' ================================================');
    
    diary off;
    
    AFdB=10.*log10(U);
    if option_b==1,
        for i=1:901
            thetarec(i)=theta(i*2-1);
            AFdBrec(i)=AFdB(i*2-1);
        end
        for i=902:1801
            thetarec(i)=2*pi-theta(3601-i*2+3);
            AFdBrec(i)=AFdB(3601-i*2+3);
            
        end
    elseif option_b==2||option_b==3,
        for i=1:MM+1
            thetarec(i)=theta(i*10-9);
            AFdBrec(i)=AFdB(i*10-9);
        end
        for i=MM+2:2*MM+1
            thetarec(i)=2*pi-thetarec(2*MM+2-i);
            AFdBrec(i)=AFdBrec(2*MM+2-i);
        end
        
    end
    
    fidaf=fopen('ArrFac.dat','wt');
    fprintf(fidaf,'%7.3f        %9.5f\n',[thetarec.*180/pi; AFdBrec]);
    fclose(fidaf);
    
    
    
    
    
    
    
    
    
    
    
    
    
    % PLOT THE GRAPHS
    % ARRAY FACTOR
    
    clf;
    
    plot(theta*180/pi,AFdB,'m','linewidth',2);
    if patt==1
        axis([0 180 max(min(AFdB)-1,-60) 1]);
    else
        axis([0 90 max(min(AFdB)-1,-60) 1]);
    end
    xlabel(['\theta',' (degrees)']),ylabel('ARRAY FACTOR(dB)')
    grid on;
    
    t1=text(1,1,['HPBW = ',num2str(max(hp)),' (degrees)']);
    set(t1,'units','normalized','position',[1 1.05],'horizontalalign','right');
  
    if option_b==1&&option_c==1&&option_d==1
        s1=title('LINEAR UNIFORM BROADSIDE','Fontsize',15);
        set(gca,'units','normalized');
        set(s1,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==2
        s2=title('LINEAR UNIFORM ORDINARY END-FIRE','Fontsize',15);
        set(gca,'units','normalized');
        set(s2,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==3
        s3=title('LINEAR UNIFORM HANSEN-WOODYARD END-FIRE','Fontsize',15);
        set(gca,'units','normalized');
        set(s3,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==4
        s4=title('LINEAR UNIFORM SCANNING','Fontsize',15);
        set(gca,'units','normalized');
        set(s4,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==2&&option_e==1
        s5=title('LINEAR NONUNIFORM BINOMIAL','Fontsize',15);
        set(gca,'units','normalized');
        set(s5,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==2&&option_e==2
        s6=title('LINEAR NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
        set(gca,'units','normalized');
        set(s6,'position',[0 1],'horizontalalign','left');
    end
    if option_b==2
        switch option_fx
            case 1
                strX='UNIFORM';
            case 2
                strX='BINOMIAL';
            case 3
                strX='DOLPH-TSCHEBYSCHEFF';
        end
        switch option_fy
            case 1
                strY='UNIFORM';
            case 2
                strY='BINOMIAL';
            case 3
                strY='DOLPH-TSCHEBYSCHEFF';
        end
        s7=title(['PLANAR, ',strX,' in X- and ',strY,' in Y- Direction'],'Fontsize',15);
        set(gca,'units','normalized');
        set(s7,'position',[0 1],'horizontalalign','left');
    end
    if option_b==3&&option_f==1
        s8=title('CIRCULAR UNIFORM (EXPONENTIAL)','Fontsize',15);
        set(gca,'units','normalized');
        set(s8,'position',[0 1],'horizontalalign','left');
    end
    if option_b==3&&option_f==2
        hold on;
        plot(theta*180/pi,AFdB0,'b--','linewidth',2);
        legend('principal+residuals','principal')
        s9=title('CIRCULAR UNIFORM (BESSEL)','Fontsize',15);
        set(gca,'units','normalized');
        set(s9,'position',[0 1],'horizontalalign','left');
    end
    
    
    %Array Factor linear
    figure;
       plot(theta*180/pi,abs(AF)/max(abs(AF)),'m','linewidth',2);
    if patt==1
        axis([0 180 0  1]);
    else
        axis([0 90 0 1]);
    end
    xlabel(['\theta',' (degrees)']),ylabel('NORMALIZED ARRAY FACTOR')
    grid on;
    
    t1=text(1,1,['HPBW = ',num2str(max(hp)),' (degrees)']);
    set(t1,'units','normalized','position',[1 1.05],'horizontalalign','right');
  
    if option_b==1&&option_c==1&&option_d==1
        s1=title('LINEAR UNIFORM BROADSIDE','Fontsize',15);
        set(gca,'units','normalized');
        set(s1,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==2
        s2=title('LINEAR UNIFORM ORDINARY END-FIRE','Fontsize',15);
        set(gca,'units','normalized');
        set(s2,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==3
        s3=title('LINEAR UNIFORM HANSEN-WOODYARD END-FIRE','Fontsize',15);
        set(gca,'units','normalized');
        set(s3,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==4
        s4=title('LINEAR UNIFORM SCANNING','Fontsize',15);
        set(gca,'units','normalized');
        set(s4,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==2&&option_e==1
        s5=title('LINEAR NONUNIFORM BINOMIAL','Fontsize',15);
        set(gca,'units','normalized');
        set(s5,'position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==2&&option_e==2
        s6=title('LINEAR NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
        set(gca,'units','normalized');
        set(s6,'position',[0 1],'horizontalalign','left');
    end
    if option_b==2
        switch option_fx
            case 1
                strX='UNIFORM';
            case 2
                strX='BINOMIAL';
            case 3
                strX='DOLPH-TSCHEBYSCHEFF';
        end
        switch option_fy
            case 1
                strY='UNIFORM';
            case 2
                strY='BINOMIAL';
            case 3
                strY='DOLPH-TSCHEBYSCHEFF';
        end
        s7=title(['PLANAR, ',strX,' in X- and ',strY,' in Y- Direction'],'Fontsize',15);
        set(gca,'units','normalized');
        set(s7,'position',[0 1],'horizontalalign','left');
    end
    if option_b==3&&option_f==1
        s8=title('CIRCULAR UNIFORM (EXPONENTIAL)','Fontsize',15);
        set(gca,'units','normalized');
        set(s8,'position',[0 1],'horizontalalign','left');
    end
    if option_b==3&&option_f==2
        hold on;
        plot(theta*180/pi,AFdB0,'b--','linewidth',2);
        legend('principal+residuals','principal')
        s9=title('CIRCULAR UNIFORM (BESSEL)','Fontsize',15);
        set(gca,'units','normalized');
        set(s9,'position',[0 1],'horizontalalign','left');
    end
    
    
    % DIRECTIVITY
    figure;
    
    diff=Do-min(D);
    subplot(2,1,1)
    
    
    plot(theta*180/pi,D,'r','linewidth',2);
    
    xlabel(['\theta',' (degrees)']),ylabel('DIRECTIVITY(dimensionless)')
    grid on;
    if patt==1
        axis([0 180 floor(min(D)-0.1*diff-.1) ceil(Do+0.1*diff+.1)]);
        
    else
        axis([0 90 floor(min(D)-0.1*diff-.1) ceil(Do+0.1*diff+.1)]);
        
    end
    t2=text(1,1,['D_0 = ',num2str(Do),' (dimensionless)']);
    set(t2,'units','normalized','position',[1 1.05],'horizontalalign','right');
    
    if option_b==1&&option_c==1&&option_d==1
        s1=title('LINEAR UNIFORM BROADSIDE','Fontsize',15);
        set(gca,'units','normalized');
        set(s1,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==2
        s2=title('LINEAR UNIFORM ORDINARY END-FIRE','Fontsize',15);
        set(gca,'units','normalized');
        set(s2,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==3
        s3=title('LINEAR UNIFORM HANSEN-WOODYARD END-FIRE','Fontsize',15);
        set(gca,'units','normalized');
        set(s3,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==1&&option_d==4
        s4=title('LINEAR UNIFORM SCANNING','Fontsize',15);
        set(gca,'units','normalized');
        set(s4,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==2&&option_e==1
        s5=title('LINEAR NONUNIFORM BINOMIAL','Fontsize',15);
        set(gca,'units','normalized');
        set(s5,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    if option_b==1&&option_c==2&&option_e==2
        s6=title('LINEAR NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
        set(gca,'units','normalized');
        set(s6,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    if option_b==2
        switch option_fx
            case 1
                strX='UNIFORM';
            case 2
                strX='BINOMIAL';
            case 3
                strX='DOLPH-TSCHEBYSCHEFF';
        end
        switch option_fy
            case 1
                strY='UNIFORM';
            case 2
                strY='BINOMIAL';
            case 3
                strY='DOLPH-TSCHEBYSCHEFF';
        end
        s7=title(['PLANAR, ',strX,' in X- and ',strY,' in Y- Direction'],'Fontsize',15);
        set(gca,'units','normalized');
        set(s7,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    if option_b==3&&option_f==1
        s8=title('UNIFORM CIRCULAR(EXPONENTIAL)','Fontsize',15);
        set(gca,'units','normalized');
        set(s8,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    if option_b==3&&option_f==2
        %         hold on;
        %         plot(theta*180/pi,D0,'b--','linewidth',2);
        %         legend('principal+residuals','principal')
        s9=title('UNIFORM CIRCULAR(BESSEL)','Fontsize',15);
        set(gca,'units','normalized');
        set(s9,'units','normalized','position',[0 1],'horizontalalign','left');
    end
    
    
    
    
    diffdB=DodB-min(DdB);
    subplot(2,1,2)
    plot(theta*180/pi,DdB,'b','linewidth',2);
    
    t3=text(1,1,['D_0 = ',num2str(DodB),' (dB)']);
    
    
    %if option_b==3&option_f==2
    %         hold on;
    %         plot(theta*180/pi,D0,'r--','linewidth',2);
    %         legend('principal+residuals','principal');
    %         end
    
    set(t3,'units','normalized','position',[1 1.05],'horizontalalign','right');
    xlabel(['\theta',' (degrees)']),ylabel('DIRECTIVITY(dB)')
    grid on;
    
    if patt==1
        axis([0 180 max(-50,10*floor(min(DdB)/10)) 10*ceil(DodB/10)]);
    else
        axis([0 90 max(-50,10*floor(min(DdB)/10)) 10*ceil((DodB)/10)]);
    end
    
    
    % AMPLITUDE AND PHASE OF EXCITATION COEFFICIENTS
    
    if option_b==1&&option_c~=2;
        for i=1:Nelem;
            Ncoef(i)=1;
        end
    end
    if option_b==1;
        figure;
        x=(1-Nelem)*d/2:d:(Nelem-1)*d/2;
        y=(1-Nelem)/2:1:(Nelem-1)/2;
        [AX,H1,H2]=plotyy(x,Ncoef(ceil(abs(y)+0.1)),x,beta.*y.*180./pi);
        set(get(AX(1),'Ylabel'),'String','AMPLITUDE','color','r');
        set(get(AX(2),'Ylabel'),'String','PHASE (degrees)','color','b');
        set(AX(1),'ycolor','r');
        set(AX(2),'ycolor','b');
        set(H1,'Linestyle','-','color','r','linewidth',2,'marker','s');
        set(H2,'Linestyle',':','color','b','linewidth',2,'marker','o');
        
        xlabel(['ARRAY LENGTH',' (\lambda)']);
        %axis([],'color','r');
        %axis([(1-Nelem)*d/2 (Nelem-1)*d/2 0 max(Ncoef)+0.1],'color','r');
        grid on;
        [hle,l3]=legend('AMPLITUDE','PHASE','Location','NorthEast'); set(hle,'color',[1 1 1]);
        
        if option_b==1&&option_c==1&&option_d==1
            s1=title('LINEAR UNIFORM BROADSIDE','Fontsize',15);
            % set(gca,'units','normalized');
            set(s1,'units','normalized','position',[0 1],'horizontalalign','left');
        end
        if option_b==1&&option_c==1&&option_d==2
            s2=title('LINEAR UNIFORM ORDINARY END-FIRE','Fontsize',15);
            set(gca,'units','normalized');
            set(s2,'units','normalized','position',[0 1],'horizontalalign','left');
        end
        if option_b==1&&option_c==1&&option_d==3
            s3=title('LINEAR UNIFORM HANSEN-WOODYARD END-FIRE','Fontsize',15);
            set(gca,'units','normalized');
            set(s3,'units','normalized','position',[0 1],'horizontalalign','left');
        end
        if option_b==1&&option_c==1&&option_d==4
            s4=title('LINEAR UNIFORM SCANNING','Fontsize',15);
            set(gca,'units','normalized');
            set(s4,'units','normalized','position',[0 1],'horizontalalign','left');
        end
        if option_b==1&&option_c==2&&option_e==1
            s5=title('LINEAR NONUNIFORM BINOMIAL','Fontsize',15);
            set(gca,'units','normalized');
            set(s5,'units','normalized','position',[0 1],'horizontalalign','left');
        end
        if option_b==1&&option_c==2&&option_e==2
            s6=title('LINEAR NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
            set(gca,'units','normalized');
            set(s6,'units','normalized','position',[0 1],'horizontalalign','left');
        end
    end
    
    if option_b==2
        
        figure;
        x=(1-Mx)*dx/2:dx:(Mx-1)*dx/2;
        y=(1-Mx)/2:1:(Mx-1)/2;
        subplot(2,1,1)
        [AX,H1,H2]=plotyy(x,Ncoefx(ceil(abs(y)+0.1)),x,betax.*y.*180./pi);
        
        
        set(get(AX(1),'Ylabel'),'String','AMPLITUDE(X)','color','r');
        set(get(AX(2),'Ylabel'),'String','PHASE (degrees)','color','b');
        set(AX(1),'ycolor','r');
        set(AX(2),'ycolor','b');
        set(H1,'Linestyle','-','color','r','linewidth',2,'marker','s');
        set(H2,'Linestyle',':','color','b','linewidth',2,'marker','o');
        xlabel(['ARRAY LENGTH ',' (\lambda)']);
        %axis([(1-Mx)*dx/2 (Mx-1)*dx/2 0 max(Ncoef)+0.1]);
        grid on;
        [hle,l3]=legend('AMPLITUDE','PHASE','Location','NorthEast'); set(hle,'color',[1 1 1]);
        
        switch option_fx
            case 1
                strX='UNIFORM';
            case 2
                strX='BINOMIAL';
            case 3
                strX='DOLPH-TSCHEBYSCHEFF';
        end
        switch option_fy
            case 1
                strY='UNIFORM';
            case 2
                strY='BINOMIAL';
            case 3
                strY='DOLPH-TSCHEBYSCHEFF';
        end
        s7=title(['PLANAR, ',strX,' in X- and ',strY,' in Y- Direction'],'Fontsize',15);
        set(gca,'units','normalized');
        set(s7,'units','normalized','position',[0 1],'horizontalalign','left');
        
        x1=(1-Ny)*dy/2:dy:(Ny-1)*dy/2;
        y1=(1-Ny)/2:1:(Ny-1)/2;
        subplot(2,1,2)
        [AX,H1,H2]=plotyy(x1,Ncoefy(ceil(abs(y1)+0.1)),x1,betay.*y1.*180./pi);
        set(get(AX(1),'Ylabel'),'String','AMPLITUDE(Y)','color','r');
        set(get(AX(2),'Ylabel'),'String','PHASE (degrees)','color','b');
        set(AX(1),'ycolor','r');
        set(AX(2),'ycolor','b');
        set(H1,'Linestyle','-','color','r','linewidth',2,'marker','s');
        set(H2,'Linestyle',':','color','b','linewidth',2,'marker','o');
        xlabel(['ARRAY LENGTH ',' (\lambda)']);
        %axis([(1-Ny)*dy/2 (Ny-1)*dy/2 0 max(Ncoef)+0.1]);
        grid on;
        [hle,l3]=legend('AMPLITUDE','PHASE','Location','NorthEast'); set(hle,'color',[1 1 1]);
        
    end
    % Directivity, absolute
    figure;
    
    if option_b==1,
        
        polar_dB(theta*180/pi,DdB,-60+max(DdB),max(DdB),12,'-')
        hold on
        polar_dB(-theta*180/pi,DdB,-60+max(DdB),max(DdB),12,'-')
        
        title('Polar plot of Directivity (0< \phi <360 degrees)','Fontsize',15)
    else
        polar_dB(theta*180/pi,DdB,-60+max(DdB),max(DdB),12,'-')
        
        hold on
        
        polar_dB(-theta*180/pi,D180dB,-60+max(DdB),max(DdB),12,'-')
        
        title(['Polar plot of Directivity (\phi= ',num2str(phieval),' degrees)'],'Fontsize',15)
    end
    
    
    
    %Directivity, normalized
    figure;
    if option_b==1,
        polar_dB(theta*180/pi,DdB-max(DdB),max(-60,6*floor(min(DdB)/6)),0,12,'-')
        hold on
        polar_dB(-theta*180/pi,DdB-max(DdB),max(-60,6*floor(min(DdB)/6)),0,12,'-')
        title('Polar plot of Relative Directivity (0< \phi <360 degrees)','Fontsize',15)
    else
        %manu=min(DdB)-max(max(DdB),max(D180dB))
        polar_dB(theta*180/pi,DdB-max(max(DdB),max(D180dB)),max(-60),0,12,'-')
        
        hold on
        polar_dB(-theta*180/pi,D180dB-max(max(DdB),max(D180dB)),max(-60),0,12,'-')
        
        %polar_dB(-theta*180/pi,D180dB-max(max(DdB),max(D180dB)),max(-60,6*floor(min(DdB)/6)),0,12,'-')
        title(['Polar plot of Relative Directivity (\phi= ',num2str(phieval),' degrees)'],'Fontsize',15)
    end
    
    %Spherical Plot3D
    D3=4*pi*(abs(AF3).^2)/Prad;
    
    D3(D3==0)=1e-20;
    D3dB=10.*log10(D3);
    D3dB1=D3dB-min(min(D3dB));
    
    %figure;
    %[x,y,z]=sph2cart(PHI,pi/2-THETA,D3dB);
    %surf(x,y,z);
    %title('Directivity(Spherical Plot3D)','Fontsize',15);
    
    
    if scale==1
        DD=D3;
    end
    if scale==2
        DD=D3dB1;
    end
    disc=size(DD,1);
    spherical_plot(DD,THETA,PHI,disc)
    if scale==1
        ss=title('3D Spherical plot of Directivity (Dimensionless)','Fontsize',15);
    else
        ss=title('3D Spherical plot of Directivity (dB)','Fontsize',15);
    end
    %title('3D Spherical plot of Directivity','Fontsize',15)
    
    
    %Amplitude Pattern, normalized
    
    figure;
    if option_b==1,
        polar_dB(theta*180/pi,20*log10(abs(AF/max( abs(AF)))),-60,0,12,'-')
        hold on
        polar_dB(-theta*180/pi, 20*log10(abs(AF/max(abs(AF)))),-60,0,12,'-')
        
        title('Polar Plot of Normalized Amplitude Pattern (0< \phi <360 degrees)','Fontsize',15)
    else
        
        polar_dB(theta*180/pi,20*log10(abs(AF/max( abs(AF)))),-60,0,12,'-')
        hold on
        polar_dB(-theta*180/pi, 20*log10(abs(AF180/max(abs(AF)))),-60,0,12,'-')
        
        title(['Polar Plot of Normalized Amplitude Pattern (\phi= ',num2str(phieval),' degrees)'],'Fontsize',15)
    end
    
    
    
    
    if option_b==1
        rephi=2;
    else
        disp(strvcat('PLOT PATTERN IN A DIFFERENT PLANE?','OPTION (1):YES','OPTION (2):NO'));
        rephi=input('OPTION NUMBER =');
        if rephi==2
        else
            recalc=2;
            close all
            clc
        end
    end
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HPBWCALC
function[hp,thmax]=hpbw(U,M)
tol=0.001;
imax=0;
j=0;


for i=1:M+1;
    
    if abs(U(i)-1)<tol && floor((j+2)/2)==imax+1,
        imax=imax+1;
        thmax(imax)=(i-1)/10;
    end
    if i>1 && abs(U(i)-1)<tol && U(i)>U(i-1) && j~=0,
        thmax(imax)=(i-1)/10;
    end
    if i>1,
        y(1)=U(i)-0.5;
        y(2)=U(i-1)-0.5;
        x(1)=(i-1)/10;
        x(2)=(i-2)/10;
        sign=y(1)*y(2);
        if sign<0,
            j=j+1;
            root(j)=x(2)-y(2)*(x(2)-x(1))/(y(2)-y(1));
            if j>=2 && y(2)>y(1),
                hp(imax,1)=root(j)-root(j-1);
            elseif j==1 && y(2)>y(1),
                hp(imax,1)=2.*root(j);
            end
        end
    end
end
if j==0
    hp(imax,:)='N/A';
elseif thmax(imax)>root(j),
    hp(imax,1)=2.*(180-root(j));
else
end

% BINOMIAL
function[AF,Ncoef,Coef]=bin(X,Nelem,d,beta)
if 2*floor(Nelem/2)==Nelem,Ncoef=Nelem/2;
else Ncoef=(Nelem+1)/2;
end
Coef=zeros(1,Ncoef);
for i=1:Ncoef;
    Coef(i)=1;
    for j=1:Ncoef-i;
        Coef(i)=Coef(i).*(Nelem-j)./j;
    end
end
if 2*floor(Nelem/2)~=Nelem,Coef(1)=Coef(1)/2;
end
u=2*pi*d*X+beta;
if 2*floor(Nelem/2)==Nelem,
    AF=0;
    for i=1:Ncoef;
        AF=AF+Coef(i).*cos((2.*i-1)/2.*u);
    end
else AF=0;
    for i=1:Ncoef;
        AF=AF+Coef(i).*cos((i-1).*u);
    end
end
if 2*round(Nelem/2)~=Nelem
    Coef(1)=2*Coef(1);
end

% TSCHEBY(THETA,NELEM,D,RDB)
function[AF,Ncoef,Coef,Zo]=tscheby(X,Nelem,d,RdB,beta)
Ro=10^(RdB/20);
P=Nelem-1;
Zo=0.5*((Ro+sqrt(Ro^2-1))^(1/P)+(Ro-sqrt(Ro^2-1))^(1/P));
if 2*floor(Nelem/2)==Nelem,
    Ncoef=Nelem/2;
    M=Ncoef;
    Coef=zeros(1,M);
    for i=1:M;
        Coef(i)=0;
        for j=i:M;
            Coef(i)=Coef(i)+(-1)^(M-j)*Zo^(2*j-1)*fact(j+M-2)*(2*M-1)/(fact(j-i)*fact(j+i-1)*fact(M-j));
        end
    end
elseif 2*floor((Nelem+1)/2)==Nelem+1,
    Ncoef=(Nelem+1)/2;
    M=Ncoef-1;
    Coef=zeros(1,M);
    for i=1:M+1;
        Coef(i)=0;
        for j=i:M+1;
            if i==1,EN=2;
            else EN=1;
            end
            Coef(i)=Coef(i)+(-1)^(M-j+1)*Zo^(2*(j-1))*fact(j+M-2)*2*M/(EN*fact(j-i)*fact(j+i-2)*fact(M-j+1));
        end
    end
    
end
u=2*pi*d*X+beta;
if 2*floor(Nelem/2)==Nelem,
    AF=0;
    for i=1:Ncoef;
        AF=AF+Coef(i)*cos((2*i-1)/2*u);
    end
elseif 2*floor((Nelem+1)/2)==Nelem+1,
    AF=0;
    for i=1:Ncoef;
        AF=AF+Coef(i)*cos((i-1)*u);
    end
end
if 2*round(Nelem/2)~=Nelem
    Coef(1)=2*Coef(1);
end

% FACT(IARG)
function[f7]=fact(iarg)
f7=1;
for j=1:iarg;
    f7=j*f7;
end

% PLANAR UNIFORM(THETA,PHI,M,D,BETA) in X-Direction
function[AFx]=af10x(theta,phi,Mx,dx,betax)
k=2*pi;
psix=k.*dx.*sin(theta).*cos(phi)+betax;
AFx=sinc((Mx.*psix./2)./pi)./sinc((psix./2)./pi);

% PLANAR UNIFORM(THETA,PHI,M,D,BETA) in Y-Direction
function[AFy]=af10y(theta,phi,Ny,dy,betay)
k=2*pi;
psiy=k.*dy.*sin(theta).*sin(phi)+betay;
AFy=sinc((Ny.*psiy./2)./pi)./sinc((psiy./2)./pi);

% CIRCULAR(THETA,PHI,THETA0,PHI0,NELEM,RAD)
function[fc]=afc(theta,phi,theta0,phi0,Nelem,rad)
k=2*pi;
dtor=pi/180;
%n=linspace(1,Nelem,Nelem);
%AF=sum(exp(i.*(k.*rad.*sin(theta).*cos(phi-phin(n))+alpha(n))));
AF=0;
for n=1:Nelem
    phin(n)=2*pi*n/Nelem;
    alpha(n)=-k.*rad.*sin(dtor.*theta0).*cos(dtor.*phi0-phin(n));
    AF=AF+exp(1i.*(k.*rad.*sin(theta).*cos(phi-phin(n))+alpha(n)));
end
AFabs=abs(AF);
fc=AFabs;


% bessel with residuals
function[fc3]=afc3(theta,phi,theta0,phi0,Nelem,rad,m)
if theta0==0||phi0==0,
    theta0=theta0+0.000001;
    phi0=phi0-0.000001;
end
k=2*pi;
dtor=pi/180;
rho0=rad.*sqrt((sin(theta).*cos(phi)-sin(theta0).*cos(phi0)).^2+(sin(theta).*sin(phi)-sin(theta0).*sin(phi0)).^2);
zeta=atan((sin(theta).*sin(phi)-sin(theta0).*sin(phi0))./(sin(theta).*cos(phi)-sin(theta0).*cos(phi0)));
AFpri=besselj(0,k.*rho0);
AFres=0;
for n=1:m
    AFres=AFres+(1i^(n.*Nelem))*besselj(n.*Nelem,k.*rho0).*cos(n.*Nelem.*zeta);
end
AFres=2*AFres;
AF=AFpri+AFres;
AFabs=abs(AF);
fc3=AFabs;


% spherical plot
function spherical_plot(r,THETA,PHI,disc)
%theta = linspace(theta_low,theta_up,disc);
%phi   = linspace(phi_low,phi_up,disc);

%[THETA,PHI] = meshgrid(theta,phi);

% spherical to rectangular conversion
x = abs(r).*sin(THETA).*cos(PHI);
y = abs(r).*sin(THETA).*sin(PHI);
z = abs(r).*cos(THETA);

% do the plot
figure; surf(x,y,z); view(135,20);
C = [.8 .8 .8]; colormap(C); axis off equal;

% Draw x, y, and z axes
set(line([1e-8;max(max(x))+3],[1e-8;1e-8],[1e-8;1e-8]),'Color','r');
set(line([1e-8;1e-8],[1e-8;max(max(y))+3],[1e-8;1e-8]),'Color','r');
set(line([1e-8;1e-8],[1e-8;1e-8],[1e-8;max(max(z))+3]),'Color','r');

% Label x, y, and z axes
text(max(max(x))+4,0,0,'x','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');
text(0,max(max(y))+4,0,'y','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');
text(0,0,max(max(z))+4,'z','FontSize',14,'FontName','Times','FontAngle','italic','Color','r');

% Fill surface using patches
patch_1 = zeros(3,disc+1);  patch_2 = zeros(3,disc+1);
patch_1(1,1:disc) = x(1,:); patch_2(1,1:disc) = x(disc,:);
patch_1(2,1:disc) = y(1,:); patch_2(2,1:disc) = y(disc,:);
patch_1(3,1:disc) = z(1,:); patch_2(3,1:disc) = z(disc,:);
patch(patch_1(1,:),patch_1(2,:),patch_1(3,:),C);
patch(patch_2(1,:),patch_2(2,:),patch_2(3,:),C);



%-----------------------------------------------------------------------------
%       polar_dB(theta,rho,rmin,rmax,rticks,line_style)
%
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
if isstr(theta) || isstr(rho)
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
        if i==ceil(i)
            text(0,is+rinc/20,['  ' num2str(i)],'verticalalignment','bottom' );
        else
            text(0,is+rinc/20,['  ' num2str(i,'%0.2f')],'verticalalignment','bottom' );
        end
        
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
        if theta(i)*180/pi >=0 && theta(i)*180/pi <=90
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

function y=sinc(x)
hh=find(x==0);
x(hh)= 1;
y = sin(pi*x)./(pi*x);
y(hh) = 1;

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
