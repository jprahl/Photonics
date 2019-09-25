% this program calculates and plots the dispersion relations of elmag. waves
% in a 1D photonic crystal (PhC) for both "s" (TE) and "p" (TM) polarization states; 
% the 1D-PhC has a unit cell made of two different dielectric materials; 
% the propagation angle in one medium is taken as input.

% author: Cazimir-Gabriel Bostan, 2008 (web: http://cgbostan.evonet.ro; email: cgbostan@yahoo.com)

clear all
tic

% point-size used in the scatter plot
S=2.5;

% refractive indices (n1<n2)
n1 = 1.45;
n2 = 3.45;

% thicknesses of the two layers in the unit cell
d1 = n2/(n1+n2);
d2 = 1-d1;

M=300; N=M;

% L=d1+d2=1 is the normalized period and frec=L/lambda is the normalized
% frequency
freq = linspace(1e-3,1,N);

% ray angle w.r.t. normal inside medium "1" varies between 0 and pi/2
incid = linspace(1e-3,pi/2,M);

TE1d = []; TM1d = [];

for j=1:N
    for k=1:M
        f = freq(k);
        t1 = incid(j);
        % ray angle w.r.t. normal inside medium "2" comes from the Snell
        % law n1*sin(t1)=n2*sin(t2)
        t2 = asin(n1*sin(t1)/n2);
        k1 = 2*pi*n1*cos(t1);
        k2 = 2*pi*n2*cos(t2);
        % dispersion relations
        TE1d(j,k) = cos(k1*d1*f)*cos(k2*d2*f)-0.5*(k1/k2+k2/k1)*sin(k1*d1*f)*sin(k2*d2*f);
        TM1d(j,k) = cos(k1*d1*f)*cos(k2*d2*f)-0.5*((n2^2/n1^2)*(k1/k2)+(n1^2/n2^2)*(k2/k1))*sin(k1*d1*f)*sin(k2*d2*f);
    end
end

[x,y]=meshgrid(incid, freq);
% logical matrices, used to select points which belong to the forbidden bands 
L1=abs(TE1d')>=1;
L2=abs(TM1d')>=1;
yte=y(L1);ytm=y(L2);
% transform incidence angle in degrees
xte = x(L1)*180/pi; 
xtm = x(L2)*180/pi; 

scatter(xte,yte,S,'r','filled'); title('Forbidden bands - s pol.'); xlabel('Incidence angle (deg)'); ylabel('Normalized freq.');
figure
scatter(xtm,ytm,S,'b','filled'); title('Forbidden bands - p pol.'); xlabel('Incidence angle (deg)'); ylabel('Normalized freq.');

toc