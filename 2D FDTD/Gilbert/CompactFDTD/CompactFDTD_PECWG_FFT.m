%***********************************************************************
%     1-D scalar wave updating equation
%***********************************************************************
clear all;

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=1; %3.0e8; %2.99792458e8;            % speed of light in free space
muz=1; %4.0*pi*1.0e-7;          % permeability of free space
epsz=1.0/(cc*cc*muz);       % permittivity of free space

freq0=1;               % frequency of source excitation
lambda0=cc/freq0;              % wavelength of source excitation
omega0=2.0*pi*freq0;           % angular frequency of source excitation

%***********************************************************************
%     Grid parameters
%***********************************************************************

ib=50;                     % number of grid cells in x-direction
ie=ib+1;                    % ib+1
N_resolution=100;          % dx resolution =lambda/N_resolution
dx=lambda0/N_resolution;   % dx: space increment of 1-D lattice
             %each time step 
nmax=(2^10)*(2^6);             %total number of time steps
choice=1;              % incident wave source  1: sine wave 2: gaussian 3: square wave
iframe=ib;

dt=dx/cc*0.5
omegadt=omega0*dt;   % to reduce repeatly computation of omega * dt
period=1./freq0/dt;  % number of time step for one period
ndelay=2*period;    % gaussian pulse center position
gwidth=ndelay^2/20; % gaussian pulse width
Emax=2,Emin=-Emax;  % plot range for E field value



%***********************************************************************
%     Material parameters
%***********************************************************************

%***********************************************************************
%     Updating coefficients for space region with nonpermeable media
%***********************************************************************

scfact=1./sqrt(muz/epsz);  % scaling factor between E and H
eta0=sqrt(muz/epsz);
c0=(cc*dt/dx)^2;   % coefficient for updating equation C0
ca0=1;
cb0=dt/epsz/dx;
da0=1;
db0=dt/muz/dx;

beta=0.1:0.5:15.1;
nbeta=length(beta);

maxband=12;
EigenFreq=zeros(length(beta),maxband);
spectrumRecord=zeros(length(beta),nmax/2+1);

for i_beta=1:nbeta
 cbeta=beta(i_beta)/epsz*dt;
 dbeta=beta(i_beta)/muz*dt;
%***********************************************************************
%     Field arrays
%************
%***********************************************************
% initialization to zero field everywhere
%ez_old1(1:ie)=0.0;
%ez_old(1:ie)=0.0;
ey(1:ib)=0.0;
hx(1:ib)=0.0;
hz(1:ie)=0.0;


nobs=3;
iobs(1)=floor(ib/2+1);
iobs(2)=floor(ib/4+1);
iobs(3)=floor(ib/5+1);
isrc=floor(ib/2+3);
freq=0.01:0.01:5;
omega=2*pi*freq;

Eyfield1=zeros(nmax,1);
Eyfield2=zeros(nmax,1);
Eyfield3=zeros(nmax,1);
snorm(1:length(freq))=0;
    Ey1DFT(1:length(freq))=0;
    Ey2DFT(1:length(freq))=0;
    Ey3DFT(1:length(freq))=0;
    
    EDFT(1:length(freq),nobs)=0;
  

%***********************************************************************
%     Movie initialization
%***********************************************************************
figure(5);
plot(ey(1:iframe),'r'),axis([0 iframe Emin Emax]);
ylabel('EZ');

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************

for n=1:nmax

%***********************************************************************
%     Update electric fields
%***********************************************************************
% incident wave source


esource=sin(2*omegadt*n)*exp(-(n-ndelay)^2/gwidth);


for i=1:ib
ey(i)=ca0*ey(i)-cbeta*hx(i)+cb0*(hz(i)-hz(i+1));
end
ey(isrc)=ey(isrc)+esource;
%ez(iL+40)=0;
% point hard source excitation
%ez(isrc)=ez(isrc)+esource;
%ez(isrc+20)=ez(isrc+20)+0.2*esource;
% update eqn for Hy
for i=1:ib
    hx(i)=da0*hx(i)+dbeta*ey(i);
end
for i=2:ib
hz(i)=da0*hz(i)+db0*(ey(i-1)-ey(i));
end
hz(1)=da0*hz(1)-db0*ey(1);
hz(ie)=da0*hz(ie)+db0*ey(ib);



% fourier integral for spectrum analysis
Eyfield1(n)=ey(iobs(1));
Eyfield2(n)=ey(iobs(2));
Eyfield3(n)=ey(iobs(3));
%for nfreq=1:length(freq)
%    snorm(nfreq)=snorm(nfreq)+esource*exp(1i*omega(nfreq)*n*dt);
%    for k=1:nobs
%    EDFT(nfreq,k)=EDFT(nfreq,k)+ey(iobs(k))*exp(1i*omega(nfreq)*n*dt);
%    end
    %EtDFT(nfreq)=EtDFT(nfreq)+ey(iobst)*exp(1i*omega(nfreq)*n*dt);
%end
%***********************************************************************
%     Visualize fields
%***********************************************************************

% if (mod(n,1) ==0)
% rtime=num2str(n);
% 
% plot(ey,'b');
% %hold on;
% %plot(x_h,eta0*hy(1:iframe),'b'),axis([0 iframe Emin Emax]);
% %hold off;
% title(['time = ',rtime,' step']);
% ylabel('Ey');
% 
% getframe(gcf);
% 
% %pause
% end

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end
spectrum1=fft(Eyfield1,nmax)/nmax;
spectrum2=fft(Eyfield2,nmax)/nmax;
spectrum3=fft(Eyfield3,nmax)/nmax;
f=1/2/dt*linspace(0,1,nmax/2+1);
%Y = fft(y,NFFT)/L;
%f = Fs/2*linspace(0,1,NFFT/2+1);
if (i_beta ==1)
            nband=1;
            EigenFreq(i_beta,nband)=f(1);
else
nband=0;
end
nband
spectrum=abs(spectrum1)+abs(spectrum2)+abs(spectrum3);
    for nf=2:length(spectrum)/5
  
        
        if (  ( spectrum(nf) >= spectrum(nf-1)) && (spectrum(nf) >= spectrum(nf+1)) && (spectrum(nf) > 0.01));
            nband=nband+1;
            nband
            if (nband <= maxband)
            EigenFreq(i_beta,nband)=f(nf);
            end
        end;
    end
figure(2);hold on;
% Plot single-sided amplitude spectrum.
%plot(f,2*abs(spectrum(1:NFFT/2+1))) 
spectrumRecord(i_beta,1:nmax/2+1)=2*spectrum(1:nmax/2+1);
plot(f,2*spectrum(1:nmax/2+1)) ; axis([0 7 0 0.1])%hold off;

kxplot=beta(i_beta)*ones(1,maxband);
figure(3); hold on; plot(kxplot,EigenFreq(i_beta,:),'r.');axis([0 15 0 7]);hold off;
%pause
end  % for k=0:2*pi/a

