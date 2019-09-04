%Total field scattered field 2-D  TE-z polarization
clear all;
%*********************************************************************
%   Constants

c=2.99792458e8;         %Speed of light in free space
mu_0=4.0*pi*1e-7;     %Permeability of free space
eps_0=1.0./(c.^2*mu_0);   %Permittivity of free space (in x-direction)
freq=1.0e14;            %frequency of source excitation
wavelength=c./freq;
angfreq=2.0*pi*freq;

%*********************************************************************
%       Grid parameters
ie=200;                 %Number of grid cells in x-direction
ib=ie+1;                %True start of grid, ie between each step

is1=20;               %Scattered field regime on left
is2=ie-is1-1;           %Scattered field regime on right

je=200;                 %Number of grid cells in y-direction
jb=je+1;

inc=2*ie;               % for inbetween grid? ? ? ? ? ? ?
incm=inc-1;

js1=20;                 %Scattered field regime on bottom
js2=je-js1-1;           %Scattered field regime on top

N=80;
dx=wavelength./N;   %N=20 space increment
dy=dx;


S=0.5;
dt=S.*dx./c;           %time step
angfreq_dt=angfreq.*dt; 

nmax=4000;               %number of time steps
period=1./freq/dt;      %why the division by dt?

ndelay=3*period;        %Center 
gwidth=ndelay^2/10;      %Gaussian pulse width  ???

Emax=1;Emin=-Emax;      %plot Y value range for E field
Hmax=Emax/377; Hmin=-Hmax; % plot Y value range for H field


%*************************************************************************
%       Material parameters
nmaterial=2;          %Number of different materials, 
                        %Here: freespace and metal

eps=[1.0,25.0];        %permittivity of material or relative permittivity
sig=[0.0,0.0e+7];     %electical conductivity
mu=[1.0,2.0];         %magnetic permiability
sim=[0.0,0.0];        %equivalent magnetic loss
%*************************************************************************
%   Coefficients for space region - nonpermeable media

scfact=1./sqrt(mu_0/eps_0);         
for i=1:nmaterial
    eaf =(dt*sig(i))/(2.0*eps_0*eps(i));
    ca0(i) =(1.0-eaf)/(1.0+eaf);
    cb0(i) =(dt/eps_0/eps(i)/dx)/(1.0+eaf);
    haf =dt*sig(i)/(2.0*mu_0*mu(i));
    da0(i) =(1-haf)/(1+haf);
    db0(i) =dt/(mu_0*mu(i))/dx/(1+haf);
end



ca(1:ib,1:jb)=ca0(1);

cb(1:ib,1:jb)=cb0(1);

da(1:ie,1:je)=da0(1);

db(1:ie,1:je)=db0(1);

%********************************************************************
%       Field arrays
ex(1:ie,1:jb)=0.0;
ey(1:ib,1:je)=0.0;
hz(1:ie,1:je)=0.0;

ex_inc(1:inc)=0.0;
ey_inc(1:inc)=0.0;
hz_inc(1:incm)=0.0;

ex_inc_low_m1=0;
ex_inc_low_m2=0;
ex_inc_high_m1=0;
ex_inc_high_m2=0;

ey_inc_low_m1=0;
ey_inc_low_m2=0;
ey_inc_high_m1=0;
ey_inc_high_m2=0;

hz_inc_low_m1=0;
hz_inc_low_m2=0;
hz_inc_high_m1=0;
hz_inc_high_m2=0;


%********************************************************************
%       Metal cylinder
diam=16;
rad=diam/2;
icenter=3*ie/5;
jcenter=3*je/5;

for i=1:ie
    for j=1:je
        dist=(i-icenter)^2+(j-jcenter)^2;
        if dist <= rad^2
            ca(i)=ca0(2);
            cb(i)=cb0(2);
            da(i)=da0(2);
            db(i)=db0(2);
        end
    end
end

%********************************************************************
% Movie initialization
x=linspace(dx,ie*dx,ie);
y=linspace(dy,je*dy,je);
rect=get(gcf,'Position');
rect(1:2)=[0 0];

%**********************************************************************
%       Time-stepping loop
choice=1;
for n=1:nmax

%***********************************************************************
%     Update field
%***********************************************************************
% incident wave source

if (choice ==1)
% 1. sinusoid
    esource=sin(angfreq_dt*n);
   elseif (choice ==2)
% 2. gaussian
    esource=1.*exp(-(n-ndelay)^2/gwidth);
   elseif (choice ==3)
% 3. square wave
   if (n >= 5 && n<=25)
    esource=1;
   else
    esource=0.;
    end
end

ex_inc(2:incm)=ex_inc(2:incm)-(dt/(eps_0*dx))*(hz_inc(2:incm)-hz_inc(1:incm-1));

% excite the incident grid
ex_inc(3)=ex_inc(3)+esource;
% simple radiative boundary condition (buffer twice , at both end point i=1, ib
ex_inc(1)=ex_inc_low_m2;
ex_inc_low_m2=ex_inc_low_m1;
ex_inc_low_m1=ex_inc(2);

ex_inc(inc)=ex_inc_high_m2;
ex_inc_high_m2=ex_inc_high_m1;
ex_inc_high_m1=ex_inc(incm);

% time up date of ex grid
ex(2:ie,2:je)=ca(2:ie,2:je).*ex(2:ie,2:je)+...
    cb(2:ie,2:je).*(hz(2:ie,2:je)-hz(2:ie,1:je-1));
ey(2:ie,2:je)=ca(2:ie,2:je).*ey(2:ie,2:je)-...
    cb(2:ie,2:je).*(hz(2:ie,2:je)-hz(1:ie-1,2:je));
%Incident field correction
    %left side
ey(is1,js1:js2)=ey(is1,js1:js2)+(dt/(eps_0*dx))*(hz_inc(is1-1));
    %right side
ey(is2,js1:js2)=ey(is2,js1:js2)-(dt/(eps_0*dx))*(hz_inc(is2));
    %bottom side
    for i=is1:is2-1
ex(i,js1)=ex(i,js1)+(dt/(eps_0*dx))*(hz_inc(i));
    end
    %top side
    for i=is1:is2-1
ex(i,js2+1)=ex(i,js2+1)-(dt/(eps_0*dx))*(hz_inc(i));
    end
% simple radiative boundary condition (buffer twice , at both end point i=1, ib
% ex(1,1:jb)=ex_low_m2;
% ex_low_m2=ex_low_m1;
% ex_low_m1=ex(2,1:jb);
% 
% ex(ib,1:jb)=ex_high_m2;
% ex_high_m2=ex_high_m1;
% ex_high_m1=ex(ib-1,1:jb);
% 
% ex(1:ib,1)=ex_low_m2;
% ex_low_m2=ex_low_m1;
% ex_low_m1=ex(1:ib,2);
% 
% ex(1:ib,jb)=ex_high_m2;
% ex_high_m2=ex_high_m1;
% ex_high_m1=ex(1:ib,jb-1);


%***********************************************************************
%     Magnetic field
%***********************************************************************
% 1-D incident wave reference grid updating
hz_inc(1:incm)=hz_inc(1:incm)-(dt/(mu_0*dx)).*(ex_inc(2:inc)-ex_inc(1:incm));

% hz grid time updating
hz(1:ie,1:je)=da(1:ie,1:je).*hz(1:ie,1:je)...
    +db(1:ie,1:je).*(ex(1:ie,2:jb)-ex(1:ie,1:je) ...
    +ey(2:ib,1:je)-ey(1:ib-1,1:je));

% hz total-scatter field correction
    %left side
hz(is1-1,js1:js2)=hz(is1-1,js1:js2)+(dt/(mu_0*dx))*(ex_inc(is1));
    %right side
hz(is2,js1:js2)=hz(is2,js1:js2)-(dt/(eps_0*dx)).*(ex_inc(is2));
    %bottom side
    for i=is1:is2-1
hz(i,js1-1)=hz(i,js1-1)-(dt/(eps_0*dx))*(ex_inc(i));
    end
    %top side
     for i=is1:is2-1
hz(i,js2+1)=hz(i,js2+1)+(dt/(eps_0*dx))*(ex_inc(i));
     end

%hz(ie/2,je/2)=esource/377;
% radiative boundary condition (buffer twice, at both end point i=2, ib+1
%-----------------------------
% hy(1,1:je)=hy_low_m2;
% hy_low_m2=hy_low_m1;
% hy_low_m1=hy(2,1:je);
% 
% hy(ie,1:jb)=hy_high_m2;
% hy_high_m2=hy_high_m1;
% hy_high_m1=hy(ie-1,1:jb);
% 
% hy(1:ie,1)=hy_low_m2;
% hy_low_m2=hy_low_m1;
% hy_low_m1=hy(1:ie,2);
% 
% hy(1:ie,je)=hy_low_m2;
% hy_low_m2=hy_low_m1;
% hy_low_m1=hy(1:ie,je-1);
% 
% hy_inc(1,1:je)=hy_inc_low_m2;
% hy_inc_low_m2=hy_inc_low_m1;
% hy_inc_low_m1=hy_inc(2,1:je);
% 
% hy_inc(ie,1:jb)=hy_inc_high_m2;
% hy_inc_high_m2=hy_inc_high_m1;
% hy_inc_high_m1=hy_inc(ie-1,1:jb);
% 
% hy_inc(1:ie,1)=hy_inc_low_m2;
% hy_inc_low_m2=hy_inc_low_m1;
% hy_inc_low_m1=hy_inc(1:ie,2);
% 
% hy_inc(1:ie,je)=hy_inc_low_m2;
% hy_inc_low_m2=hy_inc_low_m1;
% hy_inc_low_m1=hy_inc(1:ie,je-1);

%***********************************************************************
%     Visualize fields
%***********************************************************************

rtime=num2str((n));             %our time step

subplot(3,1,1),%pcolor(ex'),daspect([1 1 1]);axis([1 ie 1 je]);
plot(ex_inc);
set(gca,'units','points');
shading flat;
caxis([Emin, Emax]);
title(['Ex at time = ',rtime,' time step']);      %Labels time step process
%ylabel('EZ');

subplot(3,1,2),%pcolor(ey'); daspect([1 1 1]);axis([1 ie 1 je]);
plot(hz_inc);
set(gca,'units','points');
shading flat;
caxis([Emin, Emax]);
title(['Hy at time = ',rtime,' time step']);
%ylabel('HY');

subplot(3,1,3),pcolor(hz'); daspect([1 1 1]);axis([1 ie 1 je]);
set(gca,'units','points');
shading flat;
caxis([Hmin, Hmax]);
title(['Hy at time = ',rtime,' time step']);

M(:,n/1)=getframe(gcf,rect);


%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end

