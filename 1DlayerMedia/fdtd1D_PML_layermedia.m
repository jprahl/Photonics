%***********************************************************************
%     1-D FDTD code with simple radiation boundary conditions
%***********************************************************************
%
%     Program author: Susan C. Hagness
%                     Department of Electrical and Computer Engineering
%                     University of Wisconsin-Madison
%                     1415 Engineering Drive
%                     Madison, WI 53706-1691
%                     608-265-5739
%                     hagness@engr.wisc.edu
%
%     Date of this version:  February 2000
%
%     This MATLAB M-file implements the finite-difference time-domain
%     solution of Maxwell's curl equations over a one-dimensional space
%     lattice comprised of uniform grid cells.
%
%     To illustrate the algorithm, a sinusoidal wave (1GHz) propagating 
%     in a nonpermeable lossy medium (epsr=1.0, sigma=5.0e-3 S/m) is 
%     modeled.  The simplified finite difference system for nonpermeable
%     media (discussed in Section 3.6.6 of the text) is implemented.
%
%     The grid resolution (dx = 1.5 cm) is chosen to provide 20
%     samples per wavelength.  The Courant factor S=c*dt/dx is set to
%     the stability limit: S=1.  In 1-D, this is the "magic time step."
%
%     The computational domain is truncated using the simplest radiation
%     boundary condition for wave propagation in free space: 
%
%                      Ez(imax,n+1) = Ez(imax-1,n)
%
%     To execute this M-file, type "fdtd1D" at the MATLAB prompt.
%     This M-file displays the FDTD-computed Ez and Hy fields at every 
%     time step, and records those frames in a movie matrix, M, which is
%     played at the end of the simulation using the "movie" command.
%
%***********************************************************************

clear all;

%***********************************************************************
%     Fundamental constants
%***********************************************************************
cc=3.e8;  %2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space

               %frequency of source excitation
lambda=400.e-9 
freq=cc/lambda;             %wavelength of source excitation
omega=2.0*pi*freq;


lamdastart=300.e-9;
lamdaend=700.e-9;
iffmax=200;
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

i_size=400;  %number of grid cells in x-direction
iebc=20;
ibbc=iebc+1;
ie=i_size+iebc*2;
ib=ie+1;

dx=10.e-9  ;  %lambda/20.;             %space increment of 1-D lattice
dt=1.0*dx/cc;                   %time step
omegadt=omega*dt;

nmax=10000; %round(12.0e-9/dt);     %total number of time steps
isrc=ie/2;
period=1./freq/dt; % ndt=T_period
ndelay=3*(2.e-15/dt); %30;   %4*period;  % gaussian pulse center position
gwidth=(2.e-15/dt)^2/3 ;%ndelay^2/3; % gaussian pulse width
xmax=ie*dx; %x_max in meters
rbc=0.; % for magic time step absorption boundary condition
rbc2=0.; % for S=0.5 absorption boundary condition
rbc0=0;
rbc1=0.;
iframe=ie;  % plot range
Emax=1,Emin=-Emax;  % plot Y value range for E field
Hmax=Emax/377., Hmin=-Hmax; % plot Y value range for H field

iL=50;
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
%=== sub_air
scfact0=1./sqrt(muz/epsz);
scfact=1./sqrt(muz/(epsz));
%=== air_sub
%scfact=1./sqrt(muz/epsz);
%scfact0=1./sqrt(muz/(epsz*eps2));

ca0=(1.0-(dt*sig)/(2.0*epsz*eps))/(1.0+(dt*sig)/(2.0*epsz*eps));
cb0=(dt/epsz/eps/dx)/(1.0+(dt*sig)/(2.0*epsz*eps));
da0=1.;
db0=dt/muz/dx;

ca2=(1.0-(dt*sig)/(2.0*epsz*eps2))/(1.0+(dt*sig)/(2.0*epsz*eps2));
cb2=(dt/(epsz*eps2)/dx)/(1.0+(dt*sig)/(2.0*epsz*eps2));


maxn_runtime=10;
w_gap_target=250  % 250nm
w_gap(1:maxn_runtime)=0;
asize(1:maxn_runtime)=0;
asize(1:maxn_runtime)=0;
a_size=10; %  or 100.e-9/dx
b_size=10;
% run several times of iteration on a_size and b_size
for n_runtime=1:maxn_runtime
% set up criteria for a_size & b_size
asize(n_runtime)=a_size; %  or 100.e-9/dx
bsize(n_runtime)=b_size; %  or 50.e-9/dx
n_a=2;
n_b=1;
N_layer=15;
istart=iebc+50;
iend=istart+(a_size+b_size)*N_layer;


% initialize the grid ca cb da db
ca(1:ib)=ca0;
cb(1:ib)=cb0;
da(1:ie)=da0;
db(1:ie)=db0;
% === dielectric slab  n_a=2

for nL=1:N_layer
ca(istart+1+(nL-1)*(a_size+b_size):istart+(nL-1)*(a_size+b_size)+a_size)=ca2;
cb(istart+1+(nL-1)*(a_size+b_size):istart+(nL-1)*(a_size+b_size)+a_size)=cb2;
end
% === Right substrate GaN n=2.5
%ca(iend+1:ib-iebc-1)=ca2;
%cb(iend+1:ib-iebc-1)=cb2;

figure(3)
plot(cb)


%***********************************************************************
%     Field arrays
%***********************************************************************
sigmai(1:iebc)=0.
sigmaiH(1:iebc)=0.
ez(1:ib)=0.0;
hy(1:ie)=0.0;


Snorm(1:iffmax,1:2)=0.0;
SnormN(1:iffmax)=0.;
SxR(1:iffmax)=0.;
SxT(1:iffmax)=0.;

ezDFTReT(1:iffmax)=0.;  
ezDFTImT(1:iffmax)=0.; 
hyDFTReT(1:iffmax)=0.; 
hyDFTImT(1:iffmax)=0.; 
ezDFTReR(1:iffmax)=0.;  
ezDFTImR(1:iffmax)=0.; 
hyDFTReR(1:iffmax)=0.; 
hyDFTImR(1:iffmax)=0.; 

rmax=1.e-16;
orderbc=3;
delbc=iebc*dx;
%eta=sqrt(muz/epsz);
eta=sqrt(muz/(epsz));
sigmam=-log(rmax)*(orderbc+1)/(2*delbc*eta);
bcfactor=sigmam/(dx*(delbc^orderbc)*(orderbc+1));


for i=2:iebc
  x1=(i-0.5)*dx;
  x2=(i-1.5)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  %ca1=exp(-sigmax*dt/(epsz));
  ca1=exp(-sigmax*dt/(epsz));
  cb1=(1-ca1)/(sigmax*dx);
  i1=ibbc+1-i
  i2=ie-iebc+i
  ca(i1)=ca1;
  cb(i1)=cb1;
  sigmai(i)=sigmax;
%  ca(i2)=ca1;
%  cb(i2)=cb1;
 % axebc(1,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m)) 
end
sigmax=bcfactor*(0.5*dx)^(orderbc+1);
%ca1=exp(-sigmax*dt/(epsz));
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
  %sigmaxs=sigmax*(muz/(epsz));
  sigmaxs=sigmax*(muz/(epsz));
  da1=exp(-sigmaxs*dt/muz);  
  % (1-sigmax*dt/2.0/epsz)/(1+sigmax*dt/2.0/epsz);
  db1=(1-da1)/(sigmaxs*dx);
  i1=ibbc-i;
  i2=ie-iebc+i;
  da(i1)=da1;
  db(i1)=db1;
  sigmaiH(i)=sigmaxs;
%  da(i2)=da1;
%  db(i2)=db1;
%  axhbc(i,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m))
end
%======== air end
figure(4);
%subplot(2,1,1);
plot(sigmaiH);
hold on;
plot(sigmai*eta^2,'r');
hold off;
%stop

eta=sqrt(muz/epsz);
%eta=sqrt(muz/(epsz*eps2));
sigmam=-log(rmax)*(orderbc+1)/(2*delbc*eta);
bcfactor=sigmam/(dx*(delbc^orderbc)*(orderbc+1));


for i=2:iebc
  x1=(i-0.5)*dx;
  x2=(i-1.5)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  ca1=exp(-sigmax*dt/(epsz));
  %ca1=exp(-sigmax*dt/(epsz*eps2));
  cb1=(1-ca1)/(sigmax*dx);
  i1=ibbc+1-i
  i2=ie-iebc+i
%  ca(i1)=ca1;
%  cb(i1)=cb1;
  ca(i2)=ca1;
  cb(i2)=cb1;
 % axebc(1,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m)) 
end
sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=exp(-sigmax*dt/(epsz));
%ca1=exp(-sigmax*dt/(epsz*eps2));
cb1=(1-ca1)/(sigmax*dx);
%ca(ibbc)=ca1;
%cb(ibbc)=cb1;
ca(ib-iebc)=ca1;
cb(ib-iebc)=cb1;

for i=1:iebc
  x1=(i)*dx;
  x2=(i-1)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  sigmaxs=sigmax*(muz/(epsz));
  %sigmaxs=sigmax*(muz/(epsz*eps2));
  da1=exp(-sigmaxs*dt/muz);  
  % (1-sigmax*dt/2.0/epsz)/(1+sigmax*dt/2.0/epsz);
  db1=(1-da1)/(sigmaxs*dx);
  i1=ibbc-i;
  i2=ie-iebc+i;
%  da(i1)=da1;
%  db(i1)=db1;
  da(i2)=da1;
  db(i2)=db1;
%  axhbc(i,m)=(kapa-sigma*dt/2.D0/epsz/eps(m))/(kapa+sigma*dt/2.D0/epsz/eps(m))
end
%=============================
%***********************************************************************
%     Movie initialization
%***********************************************************************
figure(4);
x=linspace(dx,ie*dx,ie);

subplot(2,1,1),plot(ez(1:ie)/scfact,'r'); %,axis([0 iframe Emin Emax]);
ylabel('EZ');

subplot(2,1,2),plot(hy,'b') ;%,axis([0 iframe Hmin Hmax]);
xlabel('x (dx)');ylabel('HY');

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(nmax/100,gcf,rect);

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
%     Flourier Transtorm parameters
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
 
 subplot(2,1,1),plot(ez(1:ie),'r'); axis([0 iframe Emin Emax]);
 title(['time = ',rtime,' fs']);
 ylabel('EZ');
 
 subplot(2,1,2),plot(hy,'b');  axis([0 iframe Hmin Hmax]);
 title(['time = ',rtime,' fs']);
 xlabel('x (meters)');ylabel('HY');
 
 M(:,n/100)=getframe(gcf,rect);
 
end

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************
end
%stop

%************************************
% Poynting Vector
%************************************
for iff=1:iffmax
   
    SnormN(iff)=(Snorm(iff,1).^2+Snorm(iff,2).^2);

SxR(iff)=-(1./scfact)*0.5*(ezDFTReR(iff).*hyDFTReR(iff)+ezDFTImR(iff).*hyDFTImR(iff)); 
SxT(iff)=+(1./scfact)*0.5*(ezDFTReT(iff).*hyDFTReT(iff)+ezDFTImT(iff).*hyDFTImT(iff));
end
figure(5);
subplot(3,1,1);plot(lambdaiff(:),SnormN(:));
subplot(3,1,2);plot(lambdaiff(:),SxT(:)./SnormN(:));axis([lambdaiff(1) lambdaiff(iffmax) 0 2]) ;
subplot(3,1,3);plot(lambdaiff(:),SxR(:)./SnormN(:));axis([lambdaiff(1) lambdaiff(iffmax) 0 2])

T_criteria=0.01
w1(n_runtime)=0;w2(n_runtime)=0; % bandedge1
for iff=2:iffmax
    if ( SxT(iff)/SnormN(iff) >= T_criteria && SxT(iff-1)/SnormN(iff-1) <= T_criteria) 
        w1(n_runtime)=lambdaiff(iff)
    end    
    if ( SxT(iff)/SnormN(iff) <= T_criteria && SxT(iff-1)/SnormN(iff-1) >= T_criteria) 
        w2(n_runtime)=lambdaiff(iff-1)
    end  
end
w_center(n_runtime)=(w1(n_runtime)+w2(n_runtime))/2
w_gap(n_runtime)=w2(n_runtime)-w1(n_runtime)
 if ( w_gap(n_runtime) >= w_gap_target)
     stop
 end
 % rand generation to  have new a_size & b_size)
 a_size=asize(n_runtime)+fix(rand*5)
 b_size=bsize(n_runtime)+fix(rand*5)
 pause
 end  %  for n_runtime
     
     
fileoutput=0;
if (fileoutput==1)
file_id=fopen('transmission_lamda.csv','w');
 fprintf(file_id,'%6.5f \n',lambdaiff(:)/1.e-9);
 fclose(file_id);

 file_id=fopen('transmissionAgR_h5nmS05.csv','w');
 fprintf(file_id,'%6.5f \n',SxR(:)./SnormN(:));
 fclose(file_id);

file_id=fopen('transmissionAgT_h5nmS05.csv','w');
 fprintf(file_id,'%6.5f \n',SxT(:)./SnormN(:));
 fclose(file_id);
end

