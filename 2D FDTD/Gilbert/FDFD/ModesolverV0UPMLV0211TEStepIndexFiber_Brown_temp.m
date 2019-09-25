% This is Finite-Difference Frequency-Domain (FDFD) code
% This program will calculate the eigen mode with propagation constant
% (beta)
% And plot the E^2 intensity pattern, H^2 pattern
% This is a TE code, and is still in progress

clear all;
jmax=30;
imax=30;
nmax=jmax*imax;
     orderbc=2;
 	 iebc=10;
     jebc=10;
     ie=imax-1;
     ib=imax;
     je=jmax-1;
     jb=jmax;
     jbbc=iebc+1 ;
     ibbc=iebc+1;
   
     
cc=3.0e+8;
muz=4.0*pi*1.0D-7;
epsz0=1.0/(muz*cc*cc);
dx=1;
dt=cc/2/dx;
wavelength=10;
omega=2*pi*cc/wavelength;
k0=2.*pi/wavelength;
% calculate cut off freq
i=0;
llowmode=zeros(10,10);
XX=1:imax;
YY=1:jmax;
[Xp,Yp]=meshgrid(XX,YY);


TEon=1;
for mm=0:2;
for nn=0:2;
LL=(imax-1)*dx;
lamdaC=1/sqrt((mm/2/LL)^2+(nn/2/LL)^2);
if (wavelength <= lamdaC)
    i=i+1;
    %allowmode(mm,nn)=1;
    if (TEon ==1) 
%    betaz(i)=2*pi*sqrt(1/wavelength^2-1/lamdaC^2);
    figure(10000+100*mm+nn);
    subplot(2,3,1);pcolor(cos(mm*pi/LL.*Xp).*sin(nn*pi/LL.*Yp));shading flat;colorbar;
    subplot(2,3,2);pcolor(sin(mm*pi/LL.*Xp).*cos(nn*pi/LL.*Yp));shading flat;colorbar;
    subplot(2,3,4);pcolor(sin(mm*pi/LL.*Xp).*cos(nn*pi/LL.*Yp));shading flat;colorbar;
    subplot(2,3,5);pcolor(cos(mm*pi/LL.*Xp).*sin(nn*pi/LL.*Yp));shading flat;colorbar;
    subplot(2,3,6);pcolor(cos(mm*pi/LL.*Xp).*cos(nn*pi/LL.*Yp));shading flat;colorbar;
    else
    figure(i);
    subplot(2,3,1);pcolor(cos(mm*pi/LL.*Xp).*sin(nn*pi/LL.*Yp));shading flat;colorbar;
    subplot(2,3,2);pcolor(sin(mm*pi/LL.*Xp).*cos(nn*pi/LL.*Yp));shading flat;colorbar;
    subplot(2,3,3);pcolor(sin(mm*pi/LL.*Xp).*sin(nn*pi/LL.*Yp));shading flat;colorbar; 
    subplot(2,3,4);pcolor(sin(mm*pi/LL.*Xp).*cos(nn*pi/LL.*Yp));shading flat;colorbar;
    subplot(2,3,5);pcolor(cos(mm*pi/LL.*Xp).*sin(nn*pi/LL.*Yp));shading flat;colorbar;
    end
       
end
end
end


%stop

%betaz/k0;



% create the eps  square sparse matrix
epsr=1.00^2;
epsm=1.00^2;

material=2;  % choose material =1 (Au) and  =2 (Ag)
if (material ==1)
%~~~~~~~~~~~Fit for Au from 400~750nm~~~~~~~~~~~~
% ----- Drude pole -----------------------------
     %~~~~~~~~~ fit from 450~800nm~~~~~~~~~~~
	 % epsInf=6.147013866
     % Wev_Dp=4.443163516
	 % DampEV_Dp=1.03E-01
     %~~~~~~~~~ fit from 400~750nm~~~~~~~~~~~
      epsInf=5.8640944;
      Wev_Dp=8.9103231;
	  DampEV_Dp=0.0500110128;

      W_Dp=Wev_Dp/27.2214D0/(2.41888D-17);
	  gama_Dp=DampEV_Dp/27.2214D0/(2.41888D-17);

    

%------- Lorentzian 1st pole -------------------
	  delta_eps_Lp1=0.1038229*6.748518709;
	  Wev_Lp1=2.742915731;
	  DampEV_Lp1=0.30074702;
      W_Lp1=Wev_Lp1/27.2214D0/(2.41888D-17);
      gama_Lp1=DampEV_Lp1/27.2214D0/(2.41888D-17);
%------- Lorentzian 2nd pole -------------------      
	  delta_eps_Lp2=0.6639161*2.024737979;
	  Wev_Lp2=3.297578769;
	  DampEV_Lp2=0.4281952;
      W_Lp2=Wev_Lp2/27.2214D0/(2.41888D-17);
      gama_Lp2=DampEV_Lp2/27.2214D0/(2.41888D-17);

  elseif (material ==2)     
 %   Ag Drude Lorentz fit 300_700nm
% epsinf	3.32597307
% omp	8.796904782
% gamma	7.83E-02
% DeltaEps1	2.26E-01
% GL1	1.983587109
% WL1	4.291700676
% RL1	1.62E-01
% DeltaEps2	5.621772015
% GL2	4.25E-02
% WL2	2.814031556
% RL2	9.70E-01
%~~~~~~~~~ fit from Ag 300~700nm~~~~~~~~~~~
      epsInf=3.32597307;
      Wev_Dp=8.796904782;
      DampEV_Dp=7.83E-02;
	 W_Dp=Wev_Dp/27.2214D0/(2.41888D-17);
	 gama_Dp=DampEV_Dp/27.2214D0/(2.41888D-17);
%------- Lorentzian 1st pole -------------------
	  delta_eps_Lp1=2.26E-01*1.983587109;
	  Wev_Lp1=4.291700676;
	  DampEV_Lp1=1.62E-01;
      W_Lp1=Wev_Lp1/27.2214D0/(2.41888D-17);
      gama_Lp1=DampEV_Lp1/27.2214D0/(2.41888D-17);
%------- Lorentzian 2nd pole -------------------      
	  delta_eps_Lp2=5.621772015*4.25E-02;
	  Wev_Lp2=4.25E-02;
	  DampEV_Lp2=9.70E-01;
      W_Lp2=Wev_Lp2/27.2214D0/(2.41888D-17);
      gama_Lp2=DampEV_Lp2/27.2214D0/(2.41888D-17);
%===================================================
      
end   
	  

% check dielectric constant
ri=complex(0,1);
%   for iff=1:iffmax
       omega2=omega^2;
  epslon=  epsInf- W_Dp^2/(omega2+ ri* gama_Dp*omega) ...
           -delta_eps_Lp1*(W_Lp1^2)/(omega2-W_Lp1^2+2*ri*gama_Lp1*omega) ...
           -delta_eps_Lp2*(W_Lp2^2)/(omega2-W_Lp2^2+2*ri*gama_Lp2*omega);
  
 
   
 


epsrij=ones(imax,jmax);
Z=ones(imax,jmax);
X=ones(imax,jmax);
Y=ones(imax,jmax);
radius=60;
icenter= (imax)/2;
jcenter= (jmax)/2;

for j=1:jmax;
    for i=1:imax;
  
       check=radius^2-(i-icenter)^2-(j-jcenter)^2;
     
        if (check >= 0  )
            %epsrij(i,j)=epslon*epsz0;
            epsrij(i,j)=epsr;   % *epsz0;
        else
            epsrij(i,j)=epsm;  % *epsz0;
        end
    end
end
% for symmetric assignment of dielectric constant

Z(2:imax,2:jmax)=0.25*(epsrij(2:imax,2:jmax)+epsrij(1:imax-1,2:jmax)+epsrij(2:imax,1:jmax-1)+epsrij(1:imax-1,1:jmax-1));
% Z(1,1:jmax)=Z(2,1:jmax);
% Z(1:imax,1)=Z(1:imax,2);
% Z(1,1)=Z(2,2);

X(1:imax,2:jmax)=0.5*(epsrij(1:imax,2:jmax)+epsrij(1:imax,1:jmax-1));
% X(1:imax,1)=X(1:imax,2);

Y(2:imax,1:jmax)=0.5*(epsrij(2:imax,1:jmax)+epsrij(1:imax-1,1:jmax));
% Y(1,1:jmax)=Y(2,1:jmax);


%figure(1);
%pcolor(abs(Z));shading flat;daspect([1 1 1]);colorbar;

% create the sparse matrix of the material eps matrix
epsrz=sparse(Z);
clear Z;
epsrx=sparse(X);
clear X;
epsry=sparse(Y);
clear Y;
% reshape the material eps matrix
epsrz1=reshape(epsrz,nmax,1);
clear epsrz;
epsrx1=reshape(epsrx,nmax,1);
clear epsrx;
epsry1=reshape(epsry,nmax,1);
clear epsry;
% create the sparse diagonal matrix of eps
epsz=sparse(1:nmax,1:nmax,epsrz1,nmax, nmax);
clear epsrz1;
epsx=sparse(1:nmax,1:nmax,epsrx1,nmax, nmax);
clear epsrx1;
epsy=sparse(1:nmax,1:nmax,epsry1,nmax, nmax);
clear epsry1;

% UPML grid
%.......Fill grid with ether......
      
      axe(1:imax)=1.0;
	 
	  aye(1:jmax)=1.0;
	
      axh(1:imax)=1.0;

      ayh(1:jmax)=1.0;

      upml=0;
      if (upml ==1 )
%....................Fill PML Grid.........................
%!----------------------------
    
     Kmax=1.0;
     rmax=1.e-6; 
      delbc=(iebc)*dx   ;

	  eta=sqrt(muz/epsz0);
      sigmam=-log(rmax)*(orderbc+1.0)/(2.0*delbc*eta)   ;
% by changing the coefficient of bacfactor to 0.*, this code will use PEC
% as boundary condition. If 1.*, this code use UPML boundary condition
      bcfactor=1.*sigmam/(dx*(delbc^orderbc)*(orderbc+1.0)) ;
	  kfactor=(Kmax-1.0)/(dx*(delbc^orderbc)*(orderbc+1.0)) ;

%      .........

        for i=1:iebc;
          x1=dx*i   ;
          x2=dx*(i-1.0)   ;
          sigma=bcfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0))  ;
		  kapa=1.0+kfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0))   ;

		  axhbc(i)=1./complex(1,-(sigma/omega/epsz0/epsm));
         
      end
          x1=dx*0.50   ;
          sigma=bcfactor*x1^(orderbc+1.0)   ;
		  kapa=1.+kfactor*x1^(orderbc+1.0)  ;
		 
          axebc(1)= 1./complex(1,-(sigma/omega/epsz0/epsm));
                 
        for i=2:iebc    ;
          x1=dx*(i-0.50)  ;
          x2=dx*(i-1.50)   ;
          sigma=bcfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0))    ;
		  kapa=1.0+kfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0))    ;
	
          axebc(i)=1./complex(1,-(sigma/omega/epsz0/epsm));
          
      end
 
% ....FRONT REGION.....

% .......Fill EZ and HX with PML Media in Front Region Everywhere.......
	  for  i=1:iebc   ;
		 i0=iebc+1-i  ; 
	   	 i1=ibbc+1-i   ;
		 i2=ie-iebc+i    ;
		 i3=ie-iebc+i    ;
    % remove i0 i1 for symmetric case
    %    axh(i0) =axhbc(i)  ;
		 axh(i2) =axhbc(i)   ;
		
	%	 axe(i1) =axebc(i)   ;
		 axe(i3) =axebc(i)   ;
		
     end
% axh(1)=1.0 ;

	    for j=1:jebc  ;
	   	 j0=jebc+1-j    ;
		 j1=jbbc+1-j   ;
		 j2=je-jebc+j   ;
		 j3=je-jebc+j   ;
	% remove j0 j1 for symmetric case 
    %  	 ayh(j0) =axhbc(j)   ;
		 ayh(j2) =axhbc(j)   ;
		
	%	 aye(j1) =axebc(j)   ;
		 aye(j3) =axebc(j)   ;
	
     end
     
 end     % end if upml=1
     
     
 %    ayh(1)=1. ;
    for m=1:jmax;
        for i=1:imax;
            ii=(m-1)*jmax+i ;
     Axee(ii)=axe(i);
     Axhh(ii)=axh(i);
 end
end


        for j=1:jmax;
             for m=1:imax;
            jj=(j-1)*jmax+m ;
     Ayee(jj)=aye(j);
     Ayhh(jj)=ayh(j);
 end
end



     Sxe=sparse(1:nmax,1:nmax,Axee,nmax,nmax);
     Sye=sparse(1:nmax,1:nmax,Ayee,nmax,nmax);
     Sxh=sparse(1:nmax,1:nmax,Axhh,nmax,nmax);
     Syh=sparse(1:nmax,1:nmax,Ayhh,nmax,nmax);
    
%=======================================
% Create the Ux, Uy, Vx,  Vy square sparse matrix
A=sparse(1:nmax,1:nmax,-1,nmax,nmax);
A_Ux=A;
A_Ux(1:imax,1:imax)=0; %PEC BC for Hy
B=sparse(1:nmax-1,2:nmax,1,nmax,nmax);
B_Ux=B;
B_Ux(1:imax,2:imax+1)=0; %PEC BC for Hy
B_Ux(imax:imax:nmax-imax,imax+1:imax:nmax+1-imax)=0; % correction for ...
% Hy last column at i=imax, no Ez (:,imax+1) out of bound

A1=sparse(imax:imax:nmax,imax:imax:nmax,1,nmax,nmax);
B1=sparse(imax:imax:nmax-imax,imax+1:imax:nmax-imax+1,-1,nmax,nmax);
UUx=A_Ux+B_Ux; %+A1+B1;
% add A1 for Ux every imax be zero --> Hy Right edge Ez=0
% add B1 for Ux every (imax, imax+1) be zero --> Hy Right edge =0
Ux=(1/dx)*Sxe*UUx;
%Ux=(1/dx)*(A+B);
%Bx=Sxh*UUx;
A_Uy=A;
A_Uy(1:imax:nmax-imax,1:imax:nmax-imax)=0;
C=sparse(1:nmax-imax,imax+1:nmax,1,nmax,nmax);
C_Uy=C;
C_Uy(1:imax:nmax-imax,imax+1:imax:nmax-imax+1)=0;
%Uy=(1/dx)*(A+C);
A2=sparse(nmax-imax+1:nmax,nmax-imax+1:nmax,1,nmax,nmax);
UUy=A_Uy+C_Uy; %+A2;
% add A2 for Uy every upper edge be zero
Uy=(1/dx)*Sye*UUy;
%By=Sxh*Uy;
%Uy(nmax,nmax)=1;
%Vx=-Ux';
%Vy=-Uy';
%VVx=-UUx';
%VVy=-UUy';
D=sparse(2:nmax,1:nmax-1,-1,nmax,nmax);
D_Vx=D;
D_Vx(imax+1:imax:nmax+1-imax,jmax:jmax:nmax-jmax)=0; % correct Ey first colum at i=1...
%         no Hz(0:j) out of bound

D1=sparse(imax+1:imax:nmax-imax+1,imax:imax:nmax-imax,1,nmax,nmax);
% add D1 for Vx Hy left edge be zero
A3=sparse(1:imax:nmax-imax+1,1:imax:nmax-imax+1,1,nmax,nmax);
% add A3 (+1) for Vy Ez even mode along y-axis
% add A3 (-1) for Vy Ez odd  mode along y-axis
A_Vx=-A;
A_Vx(1:imax:nmax-imax,1:imax:nmax-imax)=0;  % correction for Ey=0 at i=1 PEC BC
VVx=A_Vx+D_Vx; % +D1;  %+A3;




A_Vy=-A;
A_Vy(1:imax,1:imax)=0;    %PEC BC for Ey at lower edge(start) to correct Hz(j,l)
A_Vy(nmax-imax+1:nmax,nmax-imax+1:nmax)=0;%PEC BC for Ey at upper edge(over) to correct Hz(j,l)
E=sparse(imax+1:nmax,1:nmax-imax,-1,nmax,nmax);
E_Vy=E;
E_Vy(nmax-imax+1:nmax,nmax-2*imax+1:nmax-imax)=0;%PEC BC for Ey at upper edge(over) to correct Hz(j-1,l)
%A4=sparse(1:imax,1:imax,1,nmax,nmax);
% add A4 (+1) for Ez even mode along x-axis 
% add A4 (-1) for Ez odd mode along x-axis
VVy=A_Vy+E_Vy; % -A4;
Vx=(1/dx)*Sxh*VVx;
Vy=(1/dx)*Syh*VVy;

     
%=======================================
% create eigenvalue matrix  -- eqn (12) from  Opt Exp v10 p857 (2002)
%Qxx=(-1/k0^2)*Vx*Uy*Ux*(epsz\Vy)+(epsy+(1/k0^2)*Vx*Ux)*((k0^2)*(-A)+Uy*(epsz\Vy));
%Qyy=(-1/k0^2)*Vy*Ux*Uy*(epsz\Vx)+(epsx+(1/k0^2)*Vy*Uy)*((k0^2)*(-A)+Ux*(epsz\Vx));
%Qxy=-(epsy+(1/k0^2)*Vx*Ux)*Uy*(epsz\Vx)+(1/k0^2)*Vx*Uy*(k0^2*(-A)+Ux*(epsz\Vx));
%Qyx=-(epsx+(1/k0^2)*Vy*Uy)*Ux*(epsz\Vy)+(1/k0^2)*Vy*Ux*(k0^2*(-A)+Uy*(epsz\Vy));
%=======================================
% solve for the eigenvalue and Eigen vector
%[H,beta]=eigs([Qxx Qxy; Qyx Qyy]);

%Pxx=(1/k0^2)*((-1/k0^2)*Ux*(epsz\Vy)*Vx*Uy+((k0^2)*(-A)+Ux*(epsz\Vx))*(epsx+(1/k0^2)*Vy*Uy));
%Pyy=(1/k0^2)*((-1/k0^2)*Uy*(epsz\Vx)*Vy*Ux+((k0^2)*(-A)+Uy*(epsz\Vy))*(epsy+(1/k0^2)*Vx*Ux));
%Pxy=(1/k0^2)*(Ux*(epsz\Vy)*(epsy+(1/k0^2)*Vx*Ux)-(1/k0^2)*(k0^2*(-A)+Ux*(epsz\Vx))*Vy*Ux);
%Pyx=(1/k0^2)*(Uy*(epsz\Vx)*(epsx+(1/k0^2)*Vy*Uy)-(1/k0^2)*(k0^2*(-A)+Uy*(epsz\Vy))*Vx*Uy);
Pxx=((-1/k0^2)*Ux*(epsz\Vy)*Vx*Uy+((k0^2)*(-A)+Ux*(epsz\Vx))*(epsx+(1/k0^2)*Vy*Uy));
Pyy=((-1/k0^2)*Uy*(epsz\Vx)*Vy*Ux+((k0^2)*(-A)+Uy*(epsz\Vy))*(epsy+(1/k0^2)*Vx*Ux));
Pxy=(Ux*(epsz\Vy)*(epsy+(1/k0^2)*Vx*Ux)-(1/k0^2)*(k0^2*(-A)+Ux*(epsz\Vx))*Vy*Ux);
Pyx=(Uy*(epsz\Vx)*(epsx+(1/k0^2)*Vy*Uy)-(1/k0^2)*(k0^2*(-A)+Uy*(epsz\Vy))*Vx*Uy);
%=======================================
% solve for the eigenvalue and Eigen vector
%[H,beta]=eigs([Qxx Qxy; Qyx Qyy]);
%[E,beta]=eigs([Pxx Pxy; Pyx Pyy],6,'lm');
Z=Pxx;
Pxx=eye(size(Z))*Z;
Z=Pxy;
Pxy=eye(size(Z))*Z;
Z=Pyx;
Pyx=eye(size(Z))*Z;
Z=Pyy;
Pyy=eye(size(Z))*Z;
[E,beta]=eig([Pxx Pxy; Pyx Pyy]);
M=[Pxx Pxy; Pyx Pyy];

n=60;
for i=1:length(M)/n    
AAAA(i,:)=eig(M(1+n*(i-1):n*i,1+n*(i-1):n*i));
end


for i=1:length(M)/n;
for j=1:n;
    if AAAA(i,j)>0
        AAAA(i,j)=AAAA(i,j);
        
    elseif AAAA(i,j)<=0
        AAAA(i,j)=0;
    end

end
end


%=======================================
% plot the E field
istart=1;  %iebc+1;
jstart=1;  %jebc+1;
iend=imax;  %ie-iebc;
jend=jmax;  %je-jebc;

for i=1:100
    beta00(i)=beta(i,i);
end

nmode=1;
figure(61);
Ex=reshape(E(1:nmax,1),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,1),imax,jmax);
H=-Uy*E(1:nmax,1)+Ux*E(nmax+1:2*nmax,1);
Hz=reshape(H,imax,jmax);
Hx=-i/sqrt(beta(1))*(Vx*H-i*k0*epsy*E(nmax+1:2*nmax,nmode));
Hy=-i/sqrt(beta(1))*(Vy*H+i*k0*epsx*E(1:nmax,nmode)); 
Ez=-i/k0*epsz\(Vy*Hx-Vx*Hy);


figure(62); 
Ex=reshape(E(1:nmax,2),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,2),imax,jmax);
H=-Uy*E(1:nmax,2)+Ux*E(nmax+1:2*nmax,2);
Hz=reshape(H,imax,jmax);
Hx=-i/sqrt(beta(2))*(Vx*H-i*k0*epsy*E(nmax+1:2*nmax,2));
Hy=-i/sqrt(beta(2))*(Vy*H+i*k0*epsx*E(1:nmax,2)); 
Ez=-i/k0*epsz\(Vy*Hx-Vx*Hy);





subplot(2,4,1);pcolor((abs(Ex')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,3,2);pcolor((abs(Ey')));colorbar;shading flat;daspect([1 1 1]);title('Ey');
subplot(2,4,3);pcolor((abs(reshape(Ez,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,4);pcolor((abs(Ex').^2+abs(Ey').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(2,3,4);pcolor((abs(reshape(Hx,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);title('Hx');
subplot(2,3,5);pcolor((abs(reshape(Hy,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);title('Hy');
subplot(2,3,6);pcolor((abs(Hz')));colorbar;shading flat;daspect([1 1 1]);title('H^2');
subplot(2,4,8);pcolor((abs(reshape(Hx,imax,jmax)').^2+(abs(reshape(Hy,imax,jmax)')).^2));colorbar;shading flat;daspect([1 1 1]);

figure(62);
Ex=reshape(E(1:nmax,2),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,2),imax,jmax);
H=-Uy*E(1:nmax,2)+Ux*E(nmax+1:2*nmax,2);
Hz=reshape(H,imax,jmax);
Hx=-i/sqrt(beta(2))*(Vx*H-i*k0*epsy*E(nmax+1:2*nmax,2));
Hy=-i/sqrt(beta(2))*(Vy*H+i*k0*epsx*E(1:nmax,2)); 
Ez=-i/k0*epsz\(Vy*Hx-Vx*Hy);

subplot(2,4,1);pcolor((abs(Ex')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,2);pcolor((abs(Ey')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,3);pcolor((abs(reshape(Ez,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,4);pcolor((abs(Ex').^2+abs(Ey').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,5);pcolor((abs(reshape(Hx,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,6);pcolor((abs(reshape(Hy,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,7);pcolor((abs(Hz')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,8);pcolor((abs(reshape(Hx,imax,jmax)').^2+(abs(reshape(Hy,imax,jmax)')).^2));colorbar;shading flat;daspect([1 1 1]);

figure(63);
Ex=reshape(E(1:nmax,3),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,3),imax,jmax);
H=-Uy*E(1:nmax,3)+Ux*E(nmax+1:2*nmax,3);
Hz=reshape(H,imax,jmax);
Hx=-i/sqrt(beta(3))*(Vx*H-i*k0*epsy*E(nmax+1:2*nmax,3));
Hy=-i/sqrt(beta(3))*(Vy*H+i*k0*epsx*E(1:nmax,3)); 
Ez=-i/k0*epsz\(Vy*Hx-Vx*Hy);

subplot(2,4,1);pcolor((abs(Ex')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,2);pcolor((abs(Ey')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,3);pcolor((abs(reshape(Ez,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,4);pcolor((abs(Ex').^2+abs(Ey').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,5);pcolor((abs(reshape(Hx,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,6);pcolor((abs(reshape(Hy,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,7);pcolor((abs(Hz')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,8);pcolor((abs(reshape(Hx,imax,jmax)').^2+(abs(reshape(Hy,imax,jmax)')).^2));colorbar;shading flat;daspect([1 1 1]);
 
figure(64);
Ex=reshape(E(1:nmax,4),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,4),imax,jmax);
H=-Uy*E(1:nmax,4)+Ux*E(nmax+1:2*nmax,4);
Hz=reshape(H,imax,jmax);
Hx=-i/sqrt(beta(4))*(Vx*H-i*k0*epsy*E(nmax+1:2*nmax,4));
Hy=-i/sqrt(beta(4))*(Vy*H+i*k0*epsx*E(1:nmax,4)); 
Ez=-i/k0*epsz\(Vy*Hx-Vx*Hy);

subplot(2,4,1);pcolor((abs(Ex')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,2);pcolor((abs(Ey')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,3);pcolor((abs(reshape(Ez,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,4);pcolor((abs(Ex').^2+abs(Ey').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,5);pcolor((abs(reshape(Hx,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,6);pcolor((abs(reshape(Hy,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,7);pcolor((abs(Hz')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,8);pcolor((abs(reshape(Hx,imax,jmax)').^2+(abs(reshape(Hy,imax,jmax)')).^2));colorbar;shading flat;daspect([1 1 1]);



figure(65);
Ex=reshape(E(1:nmax,5),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,5),imax,jmax);
H=-Uy*E(1:nmax,5)+Ux*E(nmax+1:2*nmax,5);
Hz=reshape(H,imax,jmax);
Hx=-i/sqrt(beta(5))*(Vx*H-i*k0*epsy*E(nmax+1:2*nmax,5));
Hy=-i/sqrt(beta(5))*(Vy*H+i*k0*epsx*E(1:nmax,5)); 
Ez=-i/k0*epsz\(Vy*Hx-Vx*Hy);

subplot(2,4,1);pcolor((abs(Ex')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,2);pcolor((abs(Ey')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,3);pcolor((abs(reshape(Ez,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,4);pcolor((abs(Ex').^2+abs(Ey').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,5);pcolor((abs(reshape(Hx,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,6);pcolor((abs(reshape(Hy,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,7);pcolor((abs(Hz')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,8);pcolor((abs(reshape(Hx,imax,jmax)').^2+(abs(reshape(Hy,imax,jmax)')).^2));colorbar;shading flat;daspect([1 1 1]);
 

figure(66);
Ex=reshape(E(1:nmax,6),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,6),imax,jmax);
H=-Uy*E(1:nmax,6)+Ux*E(nmax+1:2*nmax,6);
Hz=reshape(H,imax,jmax);
Hx=-i/sqrt(beta(6))*(Vx*H-i*k0*epsy*E(nmax+1:2*nmax,6));
Hy=-i/sqrt(beta(6))*(Vy*H+i*k0*epsx*E(1:nmax,6)); 
Ez=-i/k0*epsz\(Vy*Hx-Vx*Hy);

subplot(2,4,1);pcolor((abs(Ex')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,2);pcolor((abs(Ey')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,3);pcolor((abs(reshape(Ez,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,4);pcolor((abs(Ex').^2+abs(Ey').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,5);pcolor((abs(reshape(Hx,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,6);pcolor((abs(reshape(Hy,imax,jmax)')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,7);pcolor((abs(Hz')));colorbar;shading flat;daspect([1 1 1]);
subplot(2,4,8);pcolor((abs(reshape(Hx,imax,jmax)').^2+(abs(reshape(Hy,imax,jmax)')).^2));colorbar;shading flat;daspect([1 1 1]);





TMPlotOn=0;
TEPlotOn=0;
if (TMPlotOn==1)


% TM mode plot
figure(62);
Ex=reshape(E(1:nmax,1),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,1),imax,jmax);
subplot(2,3,1);pcolor((abs(Ex(istart:iend,jstart:jend)').^2+abs(Ey(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Ex=reshape(E(1:nmax,2),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,2),imax,jmax);
subplot(2,3,2);pcolor((abs(Ex(istart:iend,jstart:jend)').^2+abs(Ey(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Ex=reshape(E(1:nmax,3),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,3),imax,jmax);
subplot(2,3,3);pcolor((abs(Ex(istart:iend,jstart:jend)').^2+abs(Ey(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Ex=reshape(E(1:nmax,4),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,4),imax,jmax);
subplot(2,3,4);pcolor((abs(Ex(istart:iend,jstart:jend)').^2+abs(Ey(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Ex=reshape(E(1:nmax,5),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,5),imax,jmax);
subplot(2,3,5);pcolor((abs(Ex(istart:iend,jstart:jend)').^2+abs(Ey(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Ex=reshape(E(1:nmax,6),imax,jmax);
Ey=reshape(E(nmax+1:2*nmax,6),imax,jmax);
subplot(2,3,6);pcolor((abs(Ex(istart:iend,jstart:jend)').^2+abs(Ey(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);


figure(63);
H=-Uy*E(1:nmax,1)+Ux*E(nmax+1:2*nmax,1);
Hz=reshape(H,imax,jmax);
subplot(2,3,1);pcolor((abs(Hz(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

H=-Uy*E(1:nmax,2)+Ux*E(nmax+1:2*nmax,2);
Hz=reshape(H,imax,jmax);
subplot(2,3,2);pcolor((abs(Hz(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

H=-Uy*E(1:nmax,3)+Ux*E(nmax+1:2*nmax,3);
Hz=reshape(H,imax,jmax);
subplot(2,3,3);pcolor((abs(Hz(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

H=-Uy*E(1:nmax,4)+Ux*E(nmax+1:2*nmax,4);
Hz=reshape(H,imax,jmax);
subplot(2,3,4);pcolor((abs(Hz(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

H=-Uy*E(1:nmax,5)+Ux*E(nmax+1:2*nmax,5);
Hz=reshape(H,imax,jmax);
subplot(2,3,5);pcolor((abs(Hz(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

H=-Uy*E(1:nmax,6)+Ux*E(nmax+1:2*nmax,6);
Hz=reshape(H,imax,jmax);
subplot(2,3,6);pcolor((abs(Hz(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);



end


if (TEPlotOn ==1)
% TE plot
Hx=reshape(H(1:nmax,1),imax,jmax);
Hy=reshape(H(nmax+1:2*nmax,1),imax,jmax);
figure(52);
subplot(2,3,1);pcolor((abs(Hx(istart:iend,jstart:jend)').^2+abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
figure(54);
subplot(3,4,1);pcolor((abs(Hx(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(3,4,2);pcolor((abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);


Hx=reshape(H(1:nmax,2),imax,jmax);
Hy=reshape(H(nmax+1:2*nmax,2),imax,jmax);
figure(52);
subplot(2,3,2);pcolor((abs(Hx(istart:iend,jstart:jend)').^2+abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
figure(54);
subplot(3,4,3);pcolor((abs(Hx(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(3,4,4);pcolor((abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Hx=reshape(H(1:nmax,3),imax,jmax);
Hy=reshape(H(nmax+1:2*nmax,3),imax,jmax);
figure(52);
subplot(2,3,3);pcolor((abs(Hx(istart:iend,jstart:jend)').^2+abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
figure(54);
subplot(3,4,5);pcolor((abs(Hx(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(3,4,6);pcolor((abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Hx=reshape(H(1:nmax,4),imax,jmax);
Hy=reshape(H(nmax+1:2*nmax,4),imax,jmax);
figure(52);
subplot(2,3,4);pcolor((abs(Hx(istart:iend,jstart:jend)').^2+abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
figure(54);
subplot(3,4,7);pcolor((abs(Hx(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(3,4,8);pcolor((abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Hx=reshape(H(1:nmax,5),imax,jmax);
Hy=reshape(H(nmax+1:2*nmax,5),imax,jmax);
figure(52);
subplot(2,3,5);pcolor((abs(Hx(istart:iend,jstart:jend)').^2+abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
figure(54);
subplot(3,4,9);pcolor((abs(Hx(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(3,4,10);pcolor((abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

Hx=reshape(H(1:nmax,6),imax,jmax);
Hy=reshape(H(nmax+1:2*nmax,6),imax,jmax);
figure(52);
subplot(2,3,6);pcolor((abs(Hx(istart:iend,jstart:jend)').^2+abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
figure(54);
subplot(3,4,11);pcolor((abs(Hx(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);
subplot(3,4,12);pcolor((abs(Hy(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);


figure(53);
E=epsz\(-Vy*H(1:nmax,1)+Vx*H(nmax+1:2*nmax,1));  %/epsz;
Ez=reshape(E,imax,jmax);
subplot(2,3,1);pcolor((abs(Ez(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

E=epsz\(-Vy*H(1:nmax,2)+Vx*H(nmax+1:2*nmax,2));  %/epsz;
Ez=reshape(E,imax,jmax);
subplot(2,3,2);pcolor((abs(Ez(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

E=epsz\(-Vy*H(1:nmax,3)+Vx*H(nmax+1:2*nmax,3));   %/epsz;
Ez=reshape(E,imax,jmax);
subplot(2,3,3);pcolor((abs(Ez(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

E=epsz\(-Vy*H(1:nmax,4)+Vx*H(nmax+1:2*nmax,4));   %/epsz;
Ez=reshape(E,imax,jmax);
subplot(2,3,4);pcolor((abs(Ez(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

E=epsz\(-Vy*H(1:nmax,5)+Vx*H(nmax+1:2*nmax,5));    %/epsz;
Ez=reshape(E,imax,jmax);
subplot(2,3,5);pcolor((abs(Ez(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

E=epsz\(-Vy*H(1:nmax,6)+Vx*H(nmax+1:2*nmax,6));    %/epsz;
Ez=reshape(E,imax,jmax);
subplot(2,3,6);pcolor((abs(Ez(istart:iend,jstart:jend)').^2));colorbar;shading flat;daspect([1 1 1]);

end


W=([-3 2;-10 6]);



