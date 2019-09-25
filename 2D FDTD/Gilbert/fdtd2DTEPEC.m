%***********************************************************************
%     2-D FDTD TE code with PML absorbing boundary conditions
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
%     solution of Maxwell's curl equations over a two-dimensional
%     Cartesian space lattice comprised of uniform square grid cells.
%
%     To illustrate the algorithm, a 6-cm-diameter metal cylindrical 
%     scatterer in free space is modeled. The source excitation is 
%     a Gaussian pulse with a carrier frequency of 5 GHz.
%
%     The grid resolution (dx = 3 mm) was chosen to provide 20 samples
%     per wavelength at the center frequency of the pulse (which in turn
%     provides approximately 10 samples per wavelength at the high end
%     of the excitation spectrum, around 10 GHz).
%
%     The computational domain is truncated using the perfectly matched
%     layer (PML) absorbing boundary conditions.  The formulation used 
%     in this code is based on the original split-field Berenger PML. The
%     PML regions are labeled as shown in the following diagram: 
%
%            ------------------------------------------
%           |                                   (ib,jb)|
%           |                                          |
%           |                                          |
%           |                                          |
%           |                MAIN GRID                 |
%           |                                          |
%           |                                          |
%           |(1,1)                                     |
%            ------------------------------------------
%
%     To execute this M-file, type "fdtd2D" at the MATLAB prompt.
%     This M-file displays the FDTD-computed Ex, Ey, and Hz fields at 
%     every 4th time step, and records those frames in a movie matrix, 
%     M, which is played at the end of the simulation using the "movie" 
%     command.
%
%***********************************************************************

clear

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space

freq=5.0e+9;                %center frequency of source excitation
lambda=cc/freq;             %center wavelength of source excitation
omega=2.0*pi*freq;          

%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=100;           %number of grid cells in x-direction
je=100;            %number of grid cells in y-direction

ib=ie+1;
jb=je+1;

isrc=ie/2;            %location of z-directed hard source
jsrc=je/2;          %location of z-directed hard source

dx=3.0e-3;        %space increment of square lattice
dt=dx/(2.0*cc);   %time step

nmax=100;         %total number of time steps

iebc=8;           %thickness of left and right PML region
jebc=8;           %thickness of front and back PML region


%***********************************************************************
%     Material parameters
%***********************************************************************

media=2;

eps=[1.0 1.0];
sig=[0.0 5.0e+7];
mur=[1.0 1.0];
sim=[0.0 0.0];

%***********************************************************************
%     Wave excitation
%***********************************************************************

rtau=160.0e-12;
tau=rtau/dt;
delay=3*tau;

source=zeros(1,nmax);
for n=1:nmax  %7.0*tau  
    sinenusoid=sin(omega*(n)*dt); % sin(omega*(n-delay)*dt);
    gaussianprofile=1. ; %exp(-((n-delay)^2/tau^2));
  source(n)=sinenusoid*gaussianprofile
end

%***********************************************************************
%     Field arrays
%***********************************************************************

ex=zeros(ie,jb);           %fields in main grid 
ey=zeros(ib,je);
hz=zeros(ie,je);



%***********************************************************************
%     Updating coefficients
%***********************************************************************

for i=1:media
  eaf  =dt*sig(i)/(2.0*epsz*eps(i));
  ca(i)=(1.0-eaf)/(1.0+eaf);
  cb(i)=dt/epsz/eps(i)/dx/(1.0+eaf);
  haf  =dt*sim(i)/(2.0*muz*mur(i));
  da(i)=(1.0-haf)/(1.0+haf);
  db(i)=dt/muz/mur(i)/dx/(1.0+haf);
end

%***********************************************************************
%     Geometry specification (main grid)
%***********************************************************************

%     Initialize entire main grid to free space

caex(1:ie,1:jb)=ca(1);     
cbex(1:ie,1:jb)=cb(1);

caey(1:ib,1:je)=ca(1);
cbey(1:ib,1:je)=cb(1);

dahz(1:ie,1:je)=da(1);
dbhz(1:ie,1:je)=db(1);

%     Add metal cylinder

diam=20;          % diameter of cylinder: 6 cm
rad=diam/2.0;     % radius of cylinder: 3 cm
icenter=4*ie/5;   % i-coordinate of cylinder's center
jcenter=je/2;     % j-coordinate of cylinder's center

geometryOn=0;
if (geometryOn == 1)
for i=1:ie
for j=1:je
  dist2=(i+0.5-icenter)^2 + (j-jcenter)^2;
  if dist2 <= rad^2 
     caex(i,j)=ca(2);
     cbex(i,j)=cb(2);
  end
  dist2=(i-icenter)^2 + (j+0.5-jcenter)^2;
  if dist2 <= rad^2 
     caey(i,j)=ca(2);
     cbey(i,j)=cb(2);
  end
end
end
end   % end if geometry


%***********************************************************************
%     Movie initialization
%***********************************************************************

subplot(3,1,1);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(ex');
shading flat;
caxis([-80.0 80.0]);
axis([1 ie 1 jb]);
colorbar;
axis image;
axis off;
title(['Ex at time step = 0']);

subplot(3,1,2);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(ey');
shading flat;
caxis([-80.0 80.0]);
axis([1 ib 1 je]);
colorbar;
axis image;
axis off;
title(['Ey at time step = 0']);

subplot(3,1,3);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(hz');
shading flat;
caxis([-0.2 0.2]);
axis([1 ie 1 je]);
colorbar;
axis image;
axis off;
title(['Hz at time step = 0']);

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(nmax/4,gcf,rect);

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************

for n=1:nmax

%***********************************************************************
%     Update electric fields (EX and EY) in main grid
%***********************************************************************

ex(:,2:je)=caex(:,2:je).*ex(:,2:je)+...
           cbex(:,2:je).*(hz(:,2:je)-hz(:,1:je-1));

ey(2:ie,:)=caey(2:ie,:).*ey(2:ie,:)+...
           cbey(2:ie,:).*(hz(1:ie-1,:)-hz(2:ie,:));


%***********************************************************************
%     Update magnetic fields (HZ) in main grid
%***********************************************************************

hz(1:ie,1:je)=dahz(1:ie,1:je).*hz(1:ie,1:je)+... 
              dbhz(1:ie,1:je).*(ex(1:ie,2:jb)-ex(1:ie,1:je)+...
                                ey(1:ie,1:je)-ey(2:ib,1:je));

hz(isrc,jsrc)=hz(isrc,jsrc)+source(n);


%***********************************************************************
%     Visualize fields
%***********************************************************************

if mod(n,4)==0;

timestep=int2str(n);

subplot(3,1,1);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(ex');%surf(ex'); %
shading flat;daspect([1 1 1]);
caxis([-80.0 80.0]);
axis([1 ie 1 jb]);
colorbar;;
%axis image;
%axis off;
title(['Ex at time step = ',timestep]);

subplot(3,1,2);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(ey');%surf(ey');
shading flat;daspect([1 1 1]);
caxis([-80.0 80.0]);
axis([1 ib 1 je]);
colorbar;
%axis image;
%axis off;
title(['Ey at time step = ',timestep]);

subplot(3,1,3);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(hz');%surf(hz');  %
shading flat;
caxis([-0.2 0.2]);daspect([1 1 1]);
axis([1 ie 1 je]);
colorbar;
%axis image;
%axis off;
title(['Hz at time step = ',timestep]);

nn=n/4;
M(:,nn)=getframe(gcf,rect);

end;

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end

movie(gcf,M,0,10,rect);

