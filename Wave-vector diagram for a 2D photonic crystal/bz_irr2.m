%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2D discretization of irreducible Brillouin zone (IBZ); here, example for triangular lattice;
%%% the coordinates should not be exactly "0", otherwise inverting the
%%% matrices does not work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [KX,KY]=bz_irr2(M)
%%% M is the number of discretization points
[KX,KY] = meshgrid(linspace(-0.5/sqrt(3),1e-4,M), linspace(1e-4,2/3,M));
%%% the vertices of triangular IBZ
KXv=[0, 0, -0.5/sqrt(3)]; KYv=[0, 2/3, 1/2]; 
%%% keep the points that are inside the IBZ polygon (triangle)
Z=inpolygon(KX,KY,KXv,KYv);
L=Z~=0; KX=KX(L); KY=KY(L);
%%%% plot the 2D IBZ
% plotmatrix(KX,KY)
        