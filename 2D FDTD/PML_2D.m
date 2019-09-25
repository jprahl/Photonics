% 2-D PML boundary
clear
%*********************************************************************
%   Constants

c=2.99792458e8;         %Speed of light in free space
mu0=4.0*pi*1e-7;     %Permeability of free space
eps0=1.0./(c.^2*mu0);   %Permittivity of free space (in x-direction)
freq=1.0e14;            %frequency of source excitation
wavelength=c./freq;
angfreq=2.0*pi*freq;

%*********************************************************************
%       Grid parameters
ie=200;                 %Number of grid cells in x-direction

ib=ie+1;                %outer grid - x

je=200;                 %Number of grid cells in y-direction

jb=je+1;                %outer grid - y
N=20;                   %Courant factor

dx=wavelength./N;      %N=20 space increment
dy=dx;

S=0.5;
dt=S.*dx./c;           %time step
angfreq_dt=angfreq.*dt; 

nmax=500;               %number of time steps
%isrc=ie./2;             %location of z-directed hard source?
period=1./freq/dt;      %why the division by dt?

iebc=10;                %Thickness of PML region
ibbc=iebc+1;
jebc=10;
jbbc=jebc+1;

ndelay=3*period;        %Center 
gwidth=ndelay^2/N;      %Gaussian pulse width  ???

Emax=2;Emin=-Emax;      %plot Y value range for E field
Hmax=Emax/377.; Hmin=-Hmax; % plot Y value range for H field - why the comma?


%*************************************************************************
%       Material parameters
nmaterial=2;          %Number of different materials, 
                        %Here: freespace and metal

eps=[1.0,3.0];        %permittivity of material or relative permittivity
sig=[0.0,2.0e+7];     %electical conductivity
mu=[1.0,2.0];         %magnetic permiability
sim=[0.0,0.0];        %equivalent magnetic loss
%********************************************************************
%       Field arrays
ex(1:ib,1:jb)=0.0;
hy(1:ie,1:je)=0.0;

%*************************************************************************
%   Coefficients for space region - nonpermeable media
      
for i=1:nmaterial
    eaf =(dt*sig(i))/(2.0*eps0*eps(i));
    ca0(i) =(1.0-eaf)/(1.0+eaf);
    cb0(i) =(dt/eps0/eps(i)/dx)/(1.0+eaf);
    haf =dt*sig(i)/(2.0*mu0*mu(i));
    da0(i) =(1-haf)/(1+haf);
    db0(i) =dt/(mu0*mu(i))/dx/(1);
end



ca(1:ib,1:jb)=ca0(1);

cb(1:ib,1:jb)=cb0(1);

da(1:ie,1:je)=da0(1);

db(1:ie,1:je)=db0(1);

% new approach ===========
rmax=1.e-16;
orderbc=3;
delbc=iebc*dx;

eta=sqrt(mu0/(eps0));
sigmam=-log(rmax)*(orderbc+1)/(2*delbc*eta);
bcfactor=sigmam/(dx*(delbc^orderbc)*(orderbc+1));

%Left and Right sides
for i=2:iebc
  x1=(i-0.5)*dx;
  x2=(i-1.5)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));

  ca1=(1.0-(dt*sigmax)/(2.0*eps0))/(1.0+(dt*sigmax)/(2.0*eps0));
  cb1=(dt/eps0/dx)/(1.0+(dt*sigmax)/(2.0*eps0));

  i1=ibbc+1-i;
  i2=ie-iebc+i;
  ca(i1)=ca1;
  cb(i1)=cb1;
  ca(i2)=ca1;
  cb(i2)=cb1;
end

 %  Top and Bottom sides
for j=2:jebc
  y1=(j-0.5)*dy;
  y2=(j-1.5)*dy;
  sigmay=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));

  ca1=(1.0-(dt*sigmay)/(2.0*eps0))/(1.0+(dt*sigmay)/(2.0*eps0));
  cb1=(dt/eps0/dy)/(1.0+(dt*sigmay)/(2.0*eps0));

  j1=jbbc+1-j;
  j2=je-jebc+j;
  ca(j1)=ca1;
  cb(j1)=cb1;
  ca(j2)=ca1;
  cb(j2)=cb1;
 end

sigmax=bcfactor*(0.5*dx)^(orderbc+1);
ca1=(1.0-(dt*sigmax)/(2.0*eps0))/(1.0+(dt*sigmax)/(2.0*eps0));
cb1=(dt/eps0/dx)/(1.0+(dt*sigmax)/(2.0*eps0));
ca(ibbc)=ca1;
cb(ibbc)=cb1;
ca(ib-iebc)=ca1;
cb(ib-iebc)=cb1;
ca(jbbc)=ca1;
cb(jbbc)=cb1;
ca(jb-jebc)=ca1;
cb(jb-jebc)=cb1;

for i=1:iebc
  x1=(i)*dx;
  x2=(i-1)*dx;
  sigmax=bcfactor*(x1^(orderbc+1)-x2^(orderbc+1));
  sigmaxs=sigmax*(mu0/(eps0*1));
  da1= (1-sigmaxs*dt/2.0/mu0)/(1+sigmaxs*dt/2.0/mu0);
  db1= (dt/mu0/dx)/(1+sigmaxs*dt/2.0/mu0);
  i1=ibbc-i;
  i2=ie-iebc+i;
  da(i1)=da1;
  db(i1)=db1;
  da(i2)=da1;
  db(i2)=db1;
end

for j=1:jebc
  y1=(j)*dy;
  y2=(j-1)*dy;
  sigmax=bcfactor*(y1^(orderbc+1)-y2^(orderbc+1));
  sigmaxs=sigmax*(mu0/(eps0*1));
  da1= (1-sigmaxs*dt/2.0/mu0)/(1+sigmaxs*dt/2.0/mu0);
  db1= (dt/mu0/dx)/(1+sigmaxs*dt/2.0/mu0);
  j1=jbbc-j;
  j2=je-jebc+j;
  da(j1)=da1;
  db(j1)=db1;
  da(j2)=da1;
  db(j2)=db1;
end


%********************************************************************
% Movie initialization
x=linspace(dx,ie*dx,ie);

rect=get(gcf,'Position');
rect(1:2)=[0 0];

%**********************************************************************
%       Time-stepping loop

for n=1:nmax

%***********************************************************************
%     Update field
%***********************************************************************
% incident wave source
choice=1;
if (choice ==1)
% 1. sinusoid
esource=sin(angfreq_dt*n);
elseif (choice ==2)
% 2. gaussian
esource=1.*exp(-(n-ndelay)^2/gwidth);
elseif (choice ==3)
% 3. square wave
if (n >= 5 & n<=25)
esource=1;
else
esource=0.;
end
end

ex(2:ie)=ca(2:ie).*ex(2:ie)+cb(2:ie).*(hy(1:ie-1)-hy(2:ie));


 
ex(100)=esource;

hy(1:ie)=da(1:ie).*hy(1:ie)-db(1:ie).*(ex(2:ib)-ex(1:ie));


%***********************************************************************
%     Visualize fields
%***********************************************************************

rtime=num2str((n));             %our time step

subplot(2,1,1),pcolor(ex'),axis([1 ie 1 je]);
set(gca,'units','points');
shading flat;
caxis([-2.0, 2.0]);
title(['Ex at time = ',rtime,' time step']);      %Labels time step process
%ylabel('EZ');

subplot(2,1,2),pcolor(hy'); axis([1 ie 1 je]);
set(gca,'units','points');
shading flat;
caxis([-2.0, 2.0]);
title(['Hy at time = ',rtime,' time step']);
%ylabel('HY');

M(:,n/1)=getframe(gcf,rect);


%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

end

