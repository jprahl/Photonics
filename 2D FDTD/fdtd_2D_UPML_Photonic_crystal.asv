% FDTD simulation of photonic crystal 
%based on 1D code supplied by Dr.Gilbert Chang of NCKU
clear all %- deleted as runtime was being cleared

%***********************************************************************
%     Fundamental constants
%***********************************************************************
cc=3.e8;  %2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space

Emax=1; Emin=-Emax;  % plot Y value range for E field
Hmax=Emax/377.; Hmin=-Hmax; % plot Y value range for H field
               %frequency of source excitation
lambda=400.e-9; 
freq=cc/lambda;             %wavelength of source excitation
omega=2.0*pi*freq;

dtr = pi/180;               % What purpose do these serve?
rtd = 180/pi;
tpi = 2.0*pi;

lamdastart=300.e-9;
lamdaend=900.e-9;
iffmax=400;     %resolution of discretized frequency

                    % Band frequency? (highest) - see Spectral element method
                    % DFT replicates periodically on the freq axis with a
                    % period of 1/dt Hz, real part is symmetical about the
                    % middle freq w_N/2=pi/dt(rad/s)=>f_N/2=1/2dt(Hz)
                    % imaginary is anti-sym. because the DFT coeff in the
                    % freq range w_N/2<=w<=w_N are complex conjugates of
                    % those in the range 0<=w<=w_N/2 thus only the freq
                    % below w_N/2 are unique and those above are spurrious.
                    % Ergo - N real discrete time data are transformed into
                    % N/2 complex freq data and so the highest freq we can
                    % deal with using DFT is f_N/2=1/2dt (so is this N? a
                    % power of 2?

lamdaspan=(lamdaend-lamdastart)/iffmax;
omegaiff(1:iffmax)=0.;
lambdaiff(1:iffmax)=0.;
for iff=1:iffmax
   
	omegaiff(iff)=2*pi*cc/(lamdastart+iff*lamdaspan);
    lambdaiff(iff)=1e9*(lamdastart+iff*lamdaspan);  % wavelength in nm
end


%***********************************************************************
%     Grid parameters
%***********************************************************************

i_size=400;     %number of grid cells in x-direction
j_size=200;     %number of grid cells in y-direction

iebc=20;        %PML boundary
jebc=20;
ibbc=iebc+1;
jbbc=jebc+1;

ie=i_size+iebc*2;       %total grid size with PML 
je=j_size+jebc*2;
ib=ie+1;                % for Hy along x
jb=je+1;                % for Hx along y

dx=10.e-9  ;  %lambda/20.;             %space increment of 1-D lattice
dy=dx;
dt=dx/(2.0*cc);                   %time step - Courant factor S=0.5
omegadt=omega*dt;
xmax=ie*dx;         %x_max in meters
ymax=je*dy;

nmax=10000; %round(12.0e-9/dt);     %total number of time steps

iframe=(ie,je);  % plot range  USE [ie,je}?

%***********Define Source*************
isrc=ib/4; jsrc=jb/2;
period=1./freq/dt; % ndt=T_period
ndelay=3*(2.e-15/dt); %30;   %4*period;  % gaussian pulse center position
gwidth=(2.e-15/dt)^2/3 ;%ndelay^2/3; % gaussian pulse width

rmax=1.e-16;                    % ??????
orderbc=3;                      % ???????
delbc=iebc*dx;                  % length of PML
eta=sqrt(muz/(epsz));
sigmam=-log(rmax)*(orderbc+1)/(2*delbc*eta);        %????????
bcfactor=sigmam/(dx*(delbc^orderbc)*(orderbc+1));   %????????

kx=1.0;
ky=1.0;

%rbc=0.;             % for magic time step absorption boundary condition
%rbc2=0.;            % for S=0.5 absorption boundary condition
%rbc0=0;
%rbc1=0.;            % how were these supposed to relate to dt above?


%******Define Photonic crystal ***********
iL=50;                 %Position of photonic crystal
iR=ie-25;
ndelay2=1*(iR-iL+1);        

%***********************************************************************
%     Material parameters
%***********************************************************************

eps=1.0;
sig=0.0e-2;
eps2=2.0*2.0;  %n^2;

%***********************************************************************
%     Updating coefficients for space region with nonpermeable media
%***********************************************************************
%*********Initialize Grid Arrays***********
%	 skipped in original code?
      EZ=zeros(ie,je);
      DEZ=zeros(ie,je);         %?????
      HX=zeros(ib,jb);
      HY=zeros(ib,jb);		  
      BHX=zeros(ib,jb);         %?????????
      BHY=zeros(ib,jb);         %?????????
	  PEZ=zeros(ie,je);         %?????????

%=== sub_air
%scfact0=1./sqrt(muz/epsz);
scfact=1./sqrt(muz/(epsz));
%=== air_sub
%scfact=1./sqrt(muz/epsz);
%scfact0=1./sqrt(muz/(epsz*eps2));

epsr(1:ie,1:je)=1.0/epsz;       %Radial components of free space?
muer(1:ib,1:jb)=1.0/muz;        %Why is this necessary?

%**************Fill in PML*****************
% 2D code has axebc, ayebc, bxebc,...axhbc,...cyhbc =1.0 over 1:iebc/jebc
% AND axe, bxe,... ayh, ...cyh similarly defined for 1:ie/je
ca0=(1.0-(dt*sig)/(2.0*epsz*eps))/(1.0+(dt*sig)/(2.0*epsz*eps));
cb0=(dt/epsz/eps/dx)/(1.0+(dt*sig)/(2.0*epsz*eps));
da0=1.;
db0=dt/muz/dx;

ca2=(1.0-(dt*sig)/(2.0*epsz*eps2))/(1.0+(dt*sig)/(2.0*epsz*eps2));
cb2=(dt/(epsz*eps2)/dx)/(1.0+(dt*sig)/(2.0*epsz*eps2));


%global w_gap w_center SxR SxT SnormN ez hy sigmai sigmaiH da db ca cb ...
  %  w_left w_right nleftcount nrightcount;

  %  clear global
% set up criteria of photonic crystal - 
%WHY IS THIS BEFORE GRID IS INITIALIZED?

a_size(n_runtime)=a_initial_size + n_runtime; %  or 100.e-9/dx
b_size(n_runtime)=b_initial_size + n_runtime; %  or 50.e-9/dx
%n_a=2;
%n_b=1;
N_layer=15;
istart=iebc+50;
%iend=istart+(a_size+b_size)*N_layer;


% initialize the grid ca cb da db
ca(1:ib)=ca0;
cb(1:ib)=cb0;
da(1:ie)=da0;
db(1:ie)=db0;
% === dielectric slab  n_a=2

for nL=1:N_layer
ca(istart+1+(nL-1)*(a_size(n_runtime)+b_size(n_runtime))...
    :istart+(nL-1)*(a_size(n_runtime)+b_size(n_runtime))+a_size(n_runtime))=ca2;
cb(istart+1+(nL-1)*(a_size(n_runtime)+b_size(n_runtime)):istart+(nL-1)...
    *(a_size(n_runtime)+b_size(n_runtime))+a_size(n_runtime))=cb2;
end
% === Right substrate GaN n=2.5
%ca(iend+1:ib-iebc-1)=ca2;
%cb(iend+1:ib-iebc-1)=cb2;

%figure(3)
%plot(sqrt(muz/epsz)./cb)


%***********************************************************************
%     Field arrays
%***********************************************************************
sigmai(1:iebc)=0.;
sigmaiH(1:iebc)=0.;
ez(1:ib)=0.0;
hy(1:ie)=0.0;


Snorm(1:iffmax,1:2)=0.0;
SnormN(1:iffmax)=0.;
SxR(1:iffmax)=0.;
SxT(1:iffmax)=0.;

dSxR(1:iffmax-1)=0.0;
dSxT(1:iffmax-1)=0.0;

ezDFTReT(1:iffmax)=0.;  %Discrete Fourier Transform of field components
ezDFTImT(1:iffmax)=0.; 
hyDFTReT(1:iffmax)=0.; 
hyDFTImT(1:iffmax)=0.; 
ezDFTReR(1:iffmax)=0.;  
ezDFTImR(1:iffmax)=0.; 
hyDFTReR(1:iffmax)=0.; 
hyDFTImR(1:iffmax)=0.; 

kappa_x=1.0;      % I believe this is "kappa" - eq 7.91 - not wave vector
kappa_y=1.0;         

%*********Define UPML coefficients and location ***************
for i=2:iebc
  x1=(i-0.5)*dx;
  x2=(i-1.5)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  
  ca1=exp(-sigmax*dt/(epsz));
  cb1=(1-ca1)/(sigmax*dx);
  i1=ibbc+1-i;
  i2=ie-iebc+i;
  ca(i1)=ca1;
  cb(i1)=cb1;
  sigmai(i)=sigmax;
  ca(i2)=ca1;
  cb(i2)=cb1;
 % axebc(1,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m)) 
end

for j=2:iebc
  y1=(j-0.5)*dy;
  y2=(j-1.5)*dy;
  sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
  
  ca1=exp(-sigmay*dt/(epsz));
  cb1=(1-ca1)/(sigmay*dy);
  j1=jbbc+1-i;
  j2=je-jebc+i;
  ca(j1)=ca1;
  cb(j1)=cb1;
  sigmaj(j)=sigmay;
  ca(j2)=ca1;
  cb(j2)=cb1;
 % axebc(1,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m)) 
end
sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/(epsz));

cb1=(1-ca1)/(sigmax*dx);
sigmai(1)=sigmax;
ca(ibbc)=ca1;
cb(ibbc)=cb1;
%ca(ib-iebc)=ca1;
%cb(ib-iebc)=cb1;

for i=1:iebc
  x1=(i)*dx;
  x2=(i-1)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  sigmaxs=sigmax*(muz/(epsz));
  da1=exp(-sigmaxs*dt/muz);  
  % (1-sigmax*dt/2.0/epsz)/(1+sigmax*dt/2.0/epsz);
  db1=(1-da1)/(sigmaxs*dx);
  
  i1=ibbc-i;
  i2=ie-iebc+i;
  da(i1)=da1;
  db(i1)=db1;
  sigmaiH(i)=sigmaxs;
  da(i2)=da1;
  db(i2)=db1;
%  axhbc(i,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m))
end

for j=1:jebc
  y1=(j)*dy;
  y2=(y-1)*dy;
  sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
  sigmays=sigmay*(muz/(epsz));
  da1=exp(-sigmays*dt/muz);  
  % (1-sigmax*dt/2.0/epsz)/(1+sigmax*dt/2.0/epsz);
  db1=(1-da1)/(sigmays*dy);
  
  j1=jbbc-i;
  j2=je-jebc+i;
  da(j1)=da1;
  db(j1)=db1;
  sigmajH(j)=sigmays;
  da(j2)=da1;
  db(j2)=db1;
%  axhbc(i,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m))
end
%======== air end
% figure(1);
% %subplot(2,1,1);
% plot(sigmaiH);
% hold on;
% plot(sigmai*eta^2,'r');
% hold off;
% %stop

eta=sqrt(muz/epsz);
%eta=sqrt(muz/(epsz*eps2));
sigmam=-log(rmax)*(orderbc+1)/(2*delbc*eta);
bcfactor=sigmam/(dx*(delbc^orderbc)*(orderbc+1));


%for i=2:iebc
%  x1=(i-0.5)*dx;
%  x2=(i-1.5)*dx;
%  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
%  ca1=exp(-sigmax*dt/(epsz));
  %ca1=exp(-sigmax*dt/(epsz*eps2));
%  cb1=(1-ca1)/(sigmax*dx);
  %i1=ibbc+1-i;
%  i2=ie-iebc+i;
  %ca(i1)=ca1;
  %cb(i1)=cb1;
%  ca(i2)=ca1;
%  cb(i2)=cb1;
  %axebc(1,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m)) 
%end
sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/(epsz));
%ca1=exp(-sigmax*dt/(epsz*eps2));
cb1=(1-ca1)/(sigmax*dx);
%ca(ibbc)=ca1;
%cb(ibbc)=cb1;
ca(ib-iebc)=ca1;
cb(ib-iebc)=cb1;

%for i=1:iebc
%  x1=(i)*dx;
%  x2=(i-1)*dx;
%  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
%  sigmaxs=sigmax*(muz/(epsz));
  %sigmaxs=sigmax*(muz/(epsz*eps2));
%  da1=exp(-sigmaxs*dt/muz);  
  % (1-sigmax*dt/2.0/epsz)/(1+sigmax*dt/2.0/epsz);
%  db1=(1-da1)/(sigmaxs*dx);
%  i1=ibbc-i;
%  i2=ie-iebc+i;
%  da(i1)=da1;
%  db(i1)=db1;
%  da(i2)=da1;
%  db(i2)=db1;
  %axhbc(i,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m))
%end
%=============================
%***********************************************************************
%     Movie initialization
%***********************************************************************
 figure(1);
x=linspace(dx,ie*dx,ie);
y=linspace(dy,je*dy,je);
figure(1);

subplot(3,1,1);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(HX');   
shading flat;daspect([1 1 1]);
colorbar;
axis image;
axis off;
title(['Hx at time step = 0']);

subplot(3,1,2);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(HY');  
shading flat;daspect([1 1 1]);
colorbar;
axis image;
axis off;
title(['Ey at time step = 0']);

subplot(3,1,3);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(EZ');  
shading flat;daspect([1 1 1]);
colorbar;
axis image;
axis off;
title(['Hz at time step = 0']);

rect=get(gcf,'Position');
rect(1:2)=[0 0];

%subplot(2,1,1),plot(ez(1:ie)/scfact,'r'); %,axis([0 iframe Emin Emax]);
% ylabel('EZ');
% 
% subplot(2,1,2),plot(hy,'b') ;%,axis([0 iframe Hmin Hmax]);
% xlabel('x (dx)');ylabel('HY');
% 
% rect=get(gcf,'Position');
% rect(1:2)=[0 0];
% 
% M=moviein(nmax/100,gcf,rect);

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************

for n=1:nmax

%***********************************************************************
%     Update electric fields
%***********************************************************************
% incident wave source

%  gaussian
esource=exp(-(n-ndelay)^2/gwidth)*sin(omegadt*(n-ndelay));



ez(2:ie)=ca(2:ie).*ez(2:ie)-cb(2:ie).*(hy(2:ie)-hy(1:ie-1));


 

ez(iL)=ez(iL)-cb(iL)*scfact*(esource);



%***********************************************************************
%     Update magnetic fields
%***********************************************************************

hy(1:ie)=da(1:ie).*hy(1:ie)-db(1:ie).*(ez(2:ib)-ez(1:ie));
hy(iL-1)=hy(iL-1)-db(iL-1)*(esource);


%***********************************************************************
%     Fourier Transform parameters
%***********************************************************************
for iff=1:iffmax

	omegaf=2*pi*cc/(lamdastart+iff*lamdaspan);
    coswt=cos(omegaf*dt*n);
    sinwt=sin(omegaf*dt*n);
    Snorm(iff,1)=Snorm(iff,1)+coswt*esource;
    Snorm(iff,2)=Snorm(iff,2)+sinwt*esource;
 
    ezDFTReT(iff)=ezDFTReT(iff)+ coswt.*ez(ie-25) ;  
	ezDFTImT(iff)=ezDFTImT(iff)+ sinwt.*ez(ie-25);  
    hyDFTReT(iff)=hyDFTReT(iff)+ coswt.*(hy(ie-25)+hy(ie-26));  
	hyDFTImT(iff)=hyDFTImT(iff)+ sinwt.*(hy(ie-25)+hy(ie-26));  
  
    ezDFTReR(iff)=ezDFTReR(iff)+ coswt.* ez(25);  
	ezDFTImR(iff)=ezDFTImR(iff)+ sinwt.* ez(25); 
    hyDFTReR(iff)=hyDFTReR(iff)+ coswt.*(hy(25)+hy(24)); 
	hyDFTImR(iff)=hyDFTImR(iff)+ sinwt.*(hy(25)+hy(24)); 
end



%***********************************************************************
%     Visualize fields
%***********************************************************************

 if mod(n,100)==0;
 
 rtime=num2str(round(n*dt/1.0e-18));
 
%  subplot(2,1,1),plot(ez(1:ie),'r'); axis([0 iframe Emin Emax]);
%  title(['time = ',rtime,' fs']);
%  ylabel('EZ');
%  
%  subplot(2,1,2),plot(hy,'b');  axis([0 iframe Hmin Hmax]);
%  title(['time = ',rtime,' fs']);
%  xlabel('x (meters)');ylabel('HY');
 
% M(:,n/100)=getframe(gcf,rect);
 
end

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************
end











