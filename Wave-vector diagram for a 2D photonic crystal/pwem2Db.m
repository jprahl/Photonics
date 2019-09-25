%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this program calculates and plots the wave-vector diagram (i.e. photonic bands at constant frequency)
%%% for a 2D photonic crystal consisting of cylinders with circular cross-section and
%%% infinite height, arranged in a triangular lattice; oblique propagation is implicit, so 
%%% the polarization states cannot be separated in E-pol and H-pol; 'omega'is taken as input; 
%%% Fourier coefficients for the expansion of dielectric constant are calculated analytically; 
%%% the materials considered here are dielectric and dispersionless;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% the package contains the following programs:
%%%     pwem2Db.m - main program
%%%     epsgg.m - routine for calculating the matrix of Fourier coefficients
%%%                 of dielectric function
%%%     bz_irr2.m - routine for 2D discretization of irreducible Brillouin zone polygon;
%%%     kvect2.m - routine for calculating diagonal matrices with elements
%%%                 (kx+Gx) and (ky+Gy), where G=(Gx,Gy) is a reciprocal
%%%                 lattice vector
%%%     oblic_eigs.m - routine for solving the eigenvalue problem for
%%%                     H-field

%%% Author: Cazimir-Gabriel Bostan, 2008 (Bucharest, Romania)
%%%	http://cgbostan.evonet.ro
%%%	cgbostan@yahoo.com

clear all
tic

omega=0.45; % normalized frequency "a/lambda"

r=0.43; % radius of cylindrical holes (normalized w.r.t. lattice constant "a")
na=1; nb=3.45; % refractive indices (cylinders-atoms, background) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
No1=7; No2=No1; 
N1=2*No1+1; N2=2*No2+1;
N=N1*N2; % total number of plane waves used in Fourier expansions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% primitive vectors of direct lattice (normalized w.r.t. lattice constant "a")
a1=[sqrt(3)/2, -1/2, 0]; a2=[sqrt(3)/2, 1/2, 0]; 
%%% area of primitive cell
ac=norm(cross(a1,a2)); 
%%% primitive vectors of direct lattice (normalized w.r.t. lattice constant "2*pi/a"): b1=[1/sqrt(3),-1]; b2=[1/sqrt(3),1]; 
b1=(1/ac)*[a2(2),-a2(1)]; b2=(1/ac)*[-a1(2), a1(1)]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% matrix of Fourier coefficients
eps1 = feval ('epsgg',r,na,nb,b1,b2,N1,N2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2D discretization of irreducible Brillouin zone polygon
Nr=20; % even number
[BZx,BZy]=feval('bz_irr2', Nr);
kx=[]; ky=[]; kz=[];

S=2.5; % point size for scatter plot

for j=1:length(BZx)
    %%% diagonal matrices with elements (kx+Gx) si (ky+Gy)
    [kGx, kGy] = feval('kvect2',BZx(j),BZy(j),b1,b2,N1,N2);
    [P, beta]=feval('oblic_eigs',omega,kGx,kGy,eps1,N);
    L=imag(beta)==0; 
    qp=sort(beta(L)); %%% keep only the propagative modes
    display(sprintf('Calculation for k[%d] is finished',j));
    for r=1:length(qp)
        kx(j,r)=BZx(j); ky(j,r)=BZy(j); kz(j,r)=qp(r);
    end  
end
M=length(BZx)*length(qp);
scatter3(reshape(kx,1,M), reshape(ky,1,M), reshape(kz,1,M), S,'r','filled'), view(65,10)
title(sprintf('Wavevector diagram for omega=%0.5g',omega)); 
xlabel('kx'); ylabel('ky'); zlabel('kz');

toc



