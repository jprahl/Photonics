%      Program UPML

% TE code Hz, Ex, Ey, global error
% version 1.0 by Gilbert Chang 
% simplified variation from bitmapFDTD created in 1999 by Gilbert Chang
% add the metal model with All UPML
% add Total Scatter field

%	   logical gaussian,farfield
    

%      real*8 muz,Pi,cc,epsz,dtr,rtd,tpi,rmax,rnf,rns,kx,ky, &
%             sinegauss,sinusoid,gaussianProfile,dt_over_dx,sigmam, &
%			 sigma,bcfactor,DampEV 		
%	 real*8 Coeff000, Coeff01, Coeff02, Coeff03, Coeff04, Coeff05,dx,dt	
%	  real*8 lamdastart, lamdaend, lamda0

pi=3.141592654;
cc=3.0e+8;
muz=4.0*pi*1.0e-7
epsz=1.0/(muz*cc*cc);
dtr=pi/180.0D0;rtd=180.0D0/pi;tpi=2.0D0*pi;
iebc=16;jebc=16;
ibbc=iebc+1;jbbc=jebc+1;
rmax=1.e-14;orderbc=4.e0;mediabc=2;
rns=1.e0;rnf=1.e0;
lamda0=400.e-9;      
dx=lamda0/20;  
dt=dx/(2.D0*cc);nmax=500; %nmax=time step
dt_over_dx=dt/dx;
nfreq=1;
%gaussian=.true.,farfield=.true.,esrc_amp=1.0)

fcenter=cc/lamda0;
omegacenter=tpi*fcenter;
cj=0.0+1.0i;
ncycs=6;
nmin=1.0/fcenter/dt;
nenv=nmax-nmin;nmid=nmin;
i_size=100;j_size=100;MEDIA=2;
% ie,je  total gird size including PML
ie=i_size+2*iebc;je=j_size+2*jebc;
% ib for Hy x size, jb for Hx y size
ib=ie+1;jb=je+1;
%  problem parameter
isrc=ie/2;jsrc=je/2;
nobserver=3;
iobs1=isrc   ;jobs1=jsrc;
iobs2=isrc   ;jobs2=jsrc;
iobs3=isrc   ;jobs3=jsrc;
iobs4=isrc   ;jobs4=jsrc;  
%%.......Total Field/Scattered Field Parameters......
igaps=iebc+10; jgaps=jebc+10;
;igapsa=iebc+9; jgapsa=jebc+9;
	   maxdim=ib; inc=2*maxdim;
       incm = inc-1;
	  %        parameters for total/scattered boundary
       il = igaps + 1; jl = jgaps + 1;
       im = ib - igaps; jm = jb - jgaps;
       in = ie - igaps; jn = je - jgaps;
%        parameters for sa planes
       io = igapsa + 1; jo = jgapsa + 1;
       ip = ib - igapsa; jp = jb - jgapsa;
       iq = ip - 1; jq = jp - 1;


%.................Grid Arrays...............................
%	 
%............................................................
    EPS=  zeros (MEDIA,1);MUR=zeros(MEDIA,1);	 
     
       HZ=zeros(ie,je); BHZ=zeros(ie,je); EX=zeros(ie,jb); EY=zeros(ib,je);		  
	   DEX=zeros(ie,jb);DEY=zeros(ib,je);DEX_old=zeros(ie,jb);DEY_old=zeros(ib,je);
	   PEX=zeros(ie,jb);PEY=zeros(ib,je);PEX_old=zeros(ie,jb);PEY_old=zeros(ib,je);
       axe=zeros(ib,1);aye=zeros(jb,1);axh=zeros(ie,1);ayh=zeros(je,1);			
	   bxz=zeros(ib,1);bye=zeros(jb,1);bxh=zeros(ie,1);byh=zeros(je,1);
	   cxe=zeros(ib,1);cye=zeros(jb,1);cxh=zeros(ie,1);cyh=zeros(je,1);
	   epsrex=zeros(ie,jb);epsrey=zeros(ib,je);muer=zeros(ib,jb);
	   sigmazgeox=zeros(ie,jb);sigmazgeoy=zeros(ib,je);
%        arrays for the 1-d grid
 ezi=zeros(inc,1);hxi=zeros(incm,1);
  %        arrays required for rbc calculations
 freq=zeros(nfreq,1);omega=zeros(nfreq,1);
%        arrays for envelopes at sa plane
 ExDFT=zeros(ie,jb,nfreq,2);EyDFT=zeros(ib,je,nfreq,2);snorm=zeros(nfreq,2);rd=zeros(nfreq,2);

%.......General CA, DA, CB, DB Arrays.......
 axebc=zeros(ibbc,mediabc);ayebc=zeros(jebc,mediabc);
 bxebc=zeros(ibbc,mediabc);byebc=zeros(jebc,mediabc);
 cxebc=zeros(ibbc,mediabc);cyebc=zeros(jebc,mediabc);
 axhbc=zeros(ibbc,mediabc);ayhbc=zeros(jebc,mediabc);
 bxhbc=zeros(ibbc,mediabc);byhbc=zeros(jebc,mediabc);
 cxhbc=zeros(ibbc,mediabc);cyhbc=zeros(jebc,mediabc);

%           .... grid parameters....
      fil = il;
      fjl = jl;
      ril = fil - 0.5;
      rjl = fjl - 0.5;
      fim = im;
      fjm = jm;
      rim = fim + 0.5;
      rjm = fjm + 0.5;
%          .....wave-source case selection logic.....

	 phi=45.e0;
      cph = cos(phi * pi/180.e0);
      sph = sin(phi * pi/180.e0);
      if (phi <= 90.e0)
        wa = cph;
        wb = sph;
        xo = fil;
        yo = fjl;
      elseif (phi <= 180.e0) 
        wa = -sin(pi/180.e0 * (phi-90.e0));
        wb = cos(pi/180.e0 * (phi-90.e0));
        xo = fim;
        yo = fjl;
      elseif (phi <= 270.e0) 
        wa = -cos(pi/180.e0 * (phi-180.e0));
        wb = -sin(pi/180.e0 * (phi-180.e0));
        xo = fim;
        yo = fjm;
      else
        wa = sin(pi/180.e0 * (phi-270.e0));
        wb = -cos(pi/180.e0 * (phi-270.e0));
        xo = fil;
        yo = fjm;
    end

      txhz1 = wa*(fil-xo);
      txhz2 = wa*(fim-xo);
      tyhz1 = wb*(fjl-yo);
      tyhz2 = wb*(fjm-yo);
      txey1 = wa*(ril-xo) + 0.5e0;
      txey2 = wa*(rim-xo) + 0.5e0;
      tyex1 = wb*(rjl-yo) + 0.5e0;
      tyex2 = wb*(rjm-yo) + 0.5e0;
%
      es22 = 0.e0;
      es11 = 0.e0;
      ed22 = 0.e0;
      ed11 = 0.e0;

% Different Fourier Transform arrays

% Different Fourier Transform arrays



%.......Initialize Frequency and Angular Frequency Arrays.......
%      data freq/fcenter/
 lamdastart=500.e-9;
	  lamdaend=300.e-9;

    %  gwidth= gwidthtime / dt
	  wstart= cc/(lamdastart)*tpi ;% FT starting frequency
	  wspan= cc/(lamdaend)*tpi-cc/(lamdastart)*tpi;  %.254509e15
  for  i=1:nfreq
  omega(i)=wstart+(i-1)/(nfreq)*wspan;        %tpi*freq(i)
  freq(i)=omega(i)/tpi;
end
	 
	  
   
	  ndelay=floor(2.e-15/dt);   %floor(3.0*gwidth)
	  gwidth2=(ndelay)^2/9.0;	 
	          
%.......Initialize Update Coefficients.......
%.......Specify Electrical Properties......
	eps(1)=rns*rns;
	eps(2)=rnf*rnf;

	mur(1)=1.0;mur(2)=1.0;
     R = dt / 2.0 / epsz;
      RA = dt / epsz / dx;
      RB = dt / 2.0 / muz;
      RC = dt / muz / dx;

	  kx=1.e0;
	  ky=1.e0;
  
%.......Fill grid with ether......
      for  i=1:ib
      axe(i)=1.e0;
	  bxe(i)=1.e0;
	  cxe(i)=1.e0;
  end
      for j=1:jb
	  aye(j)=1.e0;
	  bye(j)=1.e0;
	  cye(j)=1.e0;
      end
	  for  i=1:ie
      axh(i)=1.e0;
	  bxh(i)=1.e0;
      cxh(i)=1.e0;
      end
	  for j=1:je
      ayh(j)=1.e0;
	  byh(j)=1.e0;
      cyh(j)=1.e0;
      end
      for  i=1:ie
	  for  j=1:jb

      epsrex(i,j)=1.e0/epsz ;
	  sigmazgeox(i,j)=0.e0;
      end
      end
	  for  i=1:ib
	  for  j=1:je
      epsrey(i,j)=1.e0/epsz;
	  sigmazgeoy(i,j)=0.e0 ;
      end
      end

    
	  for  i=1:ie
	  for  j=1:je

      muer(i,j)=1.e0/muz;
      end
      end
 
% Problem Geometry Specify
%
%...................GEOMETRY SPECIFICATION...........................
%
%        .....specify the random medium ....
%  
%      
		radius=20;
		geometal=1;
		if (geometal == 1) 
	 for  i=isrc-radius:isrc+radius
	 for  j=jsrc-radius:jsrc+radius
	  check=radius^2-((i-isrc)^2+(j-jsrc)^2);
	 if  (check >= 0)
	 sigmazgeox(i,j)=1.0;
	 sigmazgeox(i,j+1)=1.0;
	 sigmazgeoy(i,j)=1.0;
	 sigmazgeoy(i+1,j)=1.0;
     %epsrex(i,j)=1.e0/epsz/9.;
     %epsrex(i,j+1)=1.e0/epsz/9.;
     %epsrey(i,j)=1.e0/epsz/9.;
     %epsrey(i+1,j)=1.e0/epsz/9.;
 end
 end
end
end





% 
% -------- output geometry file -----------------
  
 % open(289,file='geometry.csv')
    
  for  i=iebc+1:ie-iebc
  for  j=jebc+1:je-jebc
		% sigmazgeox(i,j)
    end
end
%  close(289)
 % mur=
%
%
%.....................Fill PML Grid.........................
%


 %----------------------------
 	  kmax=1. ;
for m=1:mediabc
      delbc=(iebc)*dx;

	  eta=sqrt(muz*mur(m)/epsz/eps(m));
      sigmam=-log(rmax)*(orderbc+1.e0)/(2.e0*delbc*eta*eps(m));
      bcfactor=sigmam/(dx*(delbc^orderbc)*(orderbc+1.e0));
	  kfactor=(kmax-1.e0)/(dx*(delbc^orderbc)*(orderbc+1.e0)) ;
  
%
%       .........
%
        for i=1:iebc
          x1=dx*i;
          x2=dx*(i-1.0);
          sigma=bcfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0));
		  kapa=1.0+kfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0));

		
          axhbc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz);
          ayhbc(i,m)= axhbc(i,m);
		  bxhbc(i,m)=(kapa+sigma*dt/2.0/epsz);
		  byhbc(i,m)= bxhbc(i,m);
		  cxhbc(i,m)=(kapa-sigma*dt/2.0/epsz);
		  cyhbc(i,m)= cxhbc(i,m);
      end
          x1=dx*0.5e0;
          sigma=bcfactor*x1^(orderbc+1.0);
		  kapa=1.+kfactor*x1^(orderbc+1.0);
		 

          axebc(1,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz);
          ayebc(1,m)= axebc(1,m);
		  bxebc(1,m)=(kapa+sigma*dt/2.0/epsz);
		  byebc(1,m)= bxebc(1,m);
		  cxebc(1,m)=(kapa-sigma*dt/2.0/epsz);
		  cyebc(1,m)= cxebc(1,m);
        
        for  i=2:iebc
          x1=dx*(i-0.5);
          x2=dx*(i-1.5);
          sigma=bcfactor*(x1^(orderbc+1.)-x2^(orderbc+1.0));
		  kapa=1.0+kfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0));
	
          axebc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.e0/epsz);
          ayebc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.e0/epsz);
		  bxebc(i,m)=(kapa+sigma*dt/2.0/epsz);
		  byebc(i,m)=(kapa+sigma*dt/2.0/epsz);
		  cxebc(i,m)=(kapa-sigma*dt/2.0/epsz);
		  cyebc(i,m)=(kapa-sigma*dt/2.0/epsz);

      end
  end











 %----------------------------




%                  ....FRONT REGION.....
%
%.......Fill EZ and HX with PML Media in Front Region Everywhere.......
	  for i=1:iebc
		 i0=iebc+1-i;
	   	 i1=ibbc+1-i;
		 i2=ie-iebc+i;
		 i3=ie-iebc+i;
	   	 axh(i0) =axhbc(i,1);
		 axh(i2) =axhbc(i,1);
		 bxh(i0) =bxhbc(i,1);
		 bxh(i2) =bxhbc(i,1);
		 cxh(i0) =cxhbc(i,1);
		 cxh(i2) =cxhbc(i,1);

		 axe(i1) =axebc(i,1);
		 axe(i3) =axebc(i,1);
		 bxe(i1) =bxebc(i,1);
		 bxe(i3) =bxebc(i,1);
		 cxe(i1) =cxebc(i,1);
		 cxe(i3) =cxebc(i,1);

     end

	    for  j=1:jebc
	   	 j0=jebc+1-j;
		 j1=jbbc+1-j;
		 j2=je-jebc+j;
		 j3=je-jebc+j;
	   	 ayh(j0) =ayhbc(j,1);
		 ayh(j2) =ayhbc(j,1);
		 byh(j0) =byhbc(j,1);
		 byh(j2) =byhbc(j,1);
		 cyh(j0) =cyhbc(j,1);
		 cyh(j2) =cyhbc(j,1);
		 
		 aye(j1) =ayebc(j,1);
		 aye(j3) =ayebc(j,1);
		 bye(j1) =byebc(j,1);
		 bye(j3) =byebc(j,1);
		 cye(j1) =cyebc(j,1);
		 cye(j3) =cyebc(j,1);
     end
%open(145,file='axe.csv')
%for 155 i=1,ib
%155  write(145,199) 	  axe(i),',',axh(i),',',bxe(i),',',bxh(i),',',cxe(i),',',cxh(i)
%close(145)
%open(146,file='aye.csv')
%for 156 j=1,jb
%156  write(146,199) 	  aye(j),',',ayh(j),',',bye(j),',',byh(j),',',cye(j),',',cyh(j)
%close(146)    
%199 format(e14.6e3,a3,e14.6e3,a3,e14.6e3,a3,e14.6e3,a3,e14.6e3,a3,e14.6e3)
 
%.......Fill BC Grid Edges With PMC.......
       % axh(1) =-1.0
       % axh(ib)=-1.0
       % bxh(1) =-1.0
       % bxh(ib)=-1.0
	%	cxh(1) =-1.0
    %    cxh(ib)=-1.0
%		ayh(1) =-1.0
%        ayh(jb)=-1.0
%        byh(1) =-1.0
%        byh(jb)=-1.0
%		cyh(1) =-1.0
%        cyh(jb)=-1.0

%---------------------------------





%.......Initialize the Matrices.......
      for  i=1:ie
      for  j=1:je

	Hz(i,j)=0.0;
	BHz(i,j)=0.0;
end
end
      for  i=1:ie
      for  j=1:jb
	EX(i,j) = 0.0;
	DEX(i,j)=0.0;
	DEX_old(i,j)=0.0;
	PEX(i,j)=0.0;
	PEX_old(i,j)=0.0;
end
end
      for  i=1:ib
      for  j=1:je
	Ey(i,j) = 0.0;
	DEy(i,j)=0.e0;
	DEy_old(i,j)=0.0;
	PEy(i,j)=0.0;
	PEy_old(i,j)=0.0;

end
end

 %          .....initialize the matrices.....
      for  j=1:inc
      ezi(j) = 0.0;
      end
      for  j=1:incm
      hxi(j) = 0.0;
      end
 
  ExDFT(1:ie,1:jb,1:nfreq,1:2)=0.;
  EyDFT(1:ib,1:je,1:nfreq,1:2)=0.;
  snorm(1:nfreq,1:2)=0.;
  rd(1:nfreq,1:2)=0.;
  
% ============================
% metal plasma model parameter
% ============================
	  
 	  	 epsinfinite=8.926e0 ;   %6.60764  %*epsz
		 Wev= 11.585e0;
		 Wgap=Wev/27.2214e0/(2.41888e-17);
		 DampEV= 0.2030e0;
		 beta=DampEV/27.2214e0/(2.41888e-17);	%2.36136e14

Coeff000=epsinfinite/dt/dt+epsinfinite*beta/dt/2.e0+Wgap^2/2.e0;
Coeff01=(2.e0*epsinfinite/dt/dt)/Coeff000;
Coeff02=-(epsinfinite/dt/dt-epsinfinite*beta/dt/2.e0+Wgap^2/2.0)/Coeff000;
Coeff03=(1.e0/dt/dt+beta/dt/2.e0)/Coeff000;
Coeff04=(-2.e0/dt/dt)/Coeff000;
Coeff05=(1.e0/dt/dt-beta/dt/2.e0)/Coeff000;








% 	open(122,file='observer.csv')
%   	open(1222,file='observer2.csv')
%	open(124,file='globalerror.csv')
   
 

%***********************************************************************
%     Movie initialization
%***********************************************************************
% if ploting 3D image using surf(), change SurfOn=1
% if ploting 2D image using pcolor(), change SurfOn=0
 Emax=1.0;Emin=-Emax;
Hmax=Emax/377; Hmin=-Hmax;
SurfOn=0;
if (SurfOn ==1) 
subplot(3,1,1),surf(EX');
shading flat;
caxis([Emin Emax]);
axis([1 ie 1 jb Emin Emax]);
colorbar;
%axis image;
%axis off;
title(['EX at time step = 0']);

subplot(3,1,2),surf(EY');
shading flat;
caxis([Emin Emax]);
axis([1 ib 1 je Emin Emax]);
colorbar;
%axis image;
%axis off;
title(['Ey at time step = 0']);

subplot(3,1,3),  surf(HZ');
shading flat;
caxis([Hmin Hmax]);
axis([1 ie 1 je Hmin Hmax]);
colorbar;
%axis image;
%axis off;
title(['Hz at time step = 0']);

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(nmax/4,gcf,rect);

else
subplot(3,1,1),pcolor(EX');   
shading flat;
caxis([Emin Emax]);
axis([1 ie 1 jb]);
colorbar;
axis image;
axis off;
title(['Ex at time step = 0']);

subplot(3,1,2),pcolor(EY');  
shading flat;
caxis([Emin Emax]);
axis([1 ib 1 je]);
colorbar;
axis image;
axis off;
title(['Ey at time step = 0']);

subplot(3,1,3),pcolor(HZ');  
shading flat;
caxis([Hmin Hmax]);
axis([1 ie 1 je]);
colorbar;
axis image;
axis off;
title(['Hz at time step = 0']);

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(nmax/4,gcf,rect);
end  % 
%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************
	       for  n=1:nmax
		  
      t2=(n)*dt;
	  tt2=(n-ndelay)*dt;
%


%     Soft source
 
%      if (n.le. nlimit) then

%        sinusoid  =(10.e0-15.e0*cos(Pi*0.0125D0*(n))+6.e0*cos(tpi*0.0125D0*(n))	&
%		           -cos(3.e0*Pi*0.0125D0*(n)))/32.e0
%        gaussian=exp(-1.0*((n-ndelay)^2)/gwidth2)
%	  else
%	  sinusoid=0.e0
%	  endif

	   if ( n <= ndelay) 
	      gaussianProfile=exp(-1.0*((n-ndelay)^2)/gwidth2);
	   else
		   gaussianProfile= 1.e0;
       end
		  sinusoid=sin(omegacenter*tt2) ;
	   
        sinegauss=sinusoid*gaussianProfile;
		esource=sinegauss;

	 % 	write(122,*) esource,',', gaussianProfile,',',Ezi(10)

 
%          .....incident plane wave reference.....
      for  j=2:incm
      ezi(j) = ezi(j) + RA * (hxi(j-1)-hxi(j));
      end
      ezi(3) =ezi(3)+ esource;
      ezi(1) = es22;
      es22 = es11;
      es11 = ezi(2);
      ezi(inc) = ed22;
      ed22 = ed11;
      ed11 = ezi(incm);
	 

%	  Hz(isrc-30,jsrc)=Hz(isrc-30,jsrc)+ sinegauss

%

      
%.............EX and EY UPDATE...................................
%
%         .....MAIN for LOOPS.....
%
      
      for  i=1:ie
	  for  j=2:je


	  PEX_old2=PEX_old(i,j);
	  PEX_old(i,j)=PEX(i,j);
	  DEX_old2=DEX_old(i,j);
	  DEX_old(i,j)=DEX(i,j);
	
	   
        PEX(i,j)=aye(j)*PEX_old(i,j)+				...
                  (HZ(i,j)-HZ(i,j-1))*dt_over_dx/bye(j);
		
		%DEX(i,j)=  sigmazgeox(i,j)*(					...
				%	Coeff01*DEX_old(i,j)+Coeff02*DEX_old2+ ...
	            %    Coeff03*PEX(i,j)+Coeff04*PEX_old(i,j)+Coeff05*PEX_old2	  ...
		DEX(i,j)=		   	sigmazgeox(i,j)*(0.  ...
                   )+(1.e0-sigmazgeox(i,j))*PEX(i,j);
		DEX_old(i,j)=(1.e0-sigmazgeox(i,j))*PEX_old(i,j)+sigmazgeox(i,j)*0. ; %sigmazgeox(i,j)*DEX_old(i,j);
	
		EX(i,j)=EX(i,j)+				 ...
		            (bxh(i)*DEX(i,j)-cxh(i)*DEX_old(i,j))*epsrex(i,j);	
	
            end
        end


        for i=2:ie
	  for  j=1:je

		PEY_old2=PEY_old(i,j);
	    PEY_old(i,j)=PEY(i,j);
	    DEY_old2=DEY_old(i,j);
	    DEY_old(i,j)=DEY(i,j);
	  		
         PEY(i,j)=PEY_old(i,j)-					 ...
                  (HZ(i,j)-HZ(i-1,j))*dt_over_dx;

		 %DEY(i,j)= sigmazgeoy(i,j)*(					 ...
		           %Coeff01*DEY_old(i,j)+Coeff02*DEY_old2+ ...
	              %Coeff03*PEY(i,j)+Coeff04*PEY_old(i,j)+Coeff05*PEY_old2	...
			DEY(i,j)=	  sigmazgeoy(i,j)*(0.  ...
                  )+(1.e0-sigmazgeoy(i,j))*PEY(i,j);
		 DEY_old(i,j)=(1.e0-sigmazgeoy(i,j))*PEY_old(i,j)+sigmazgeoy(i,j)*0. ; %sigmazgeoy(i,j)*DEY_old(i,j);


		 EY(i,j)=axe(i)*EY(i,j)+        ...
		          (byh(j)*DEY(i,j)-cyh(j)*DEY_old(i,j))*epsrey(i,j)/bxe(i);
		

          end
      end



 %          .....ex wave source condition, front and back.....
      for i=il:in
        dis1 = wa*((i)+0.5-xo) + tyex1;
        dis2 = wa*((i)+0.5-xo) + tyex2;
        l1 = floor(dis1);
        l2 = floor(dis2);
        hzin1 = -((dis1-(l1))*hxi(l1+6) +		  ...
                ((l1)+1.0-dis1)*hxi(l1+5));
        hzin2 = -((dis2-(l2))*hxi(l2+6) +		  ...
                ((l2)+1.0-dis2)*hxi(l2+5));

			
        EX(i,jl) = EX(i,jl) - RA*hzin1;
        EX(i,jm) = EX(i,jm) + RA*hzin2;
      end
%

%
%          .....ey iteration.....
%
%          .....ey wave source condition, left and right.....
      for  j=jl:jn
        dis1 = wb*((j)+0.5-yo) + txey1;
        dis2 = wb*((j)+0.5-yo) + txey2;
        l1 = floor(dis1);
        l2 = floor(dis2);


        hzin1 = -((dis1-(l1))*hxi(l1+6) +	   ...
                ((l1)+1.0-dis1)*hxi(l1+5));
        hzin2 = -((dis2-(l2))*hxi(l2+6) +	   ...
                ((l2)+1.0-dis2)*hxi(l2+5));

				

        EY(il,j) = EY(il,j) + RA*hzin1;
        EY(im,j) = EY(im,j) - RA*hzin2;
    end


%
%................Hz  UPDATE...................................
%
%         .....MAIN for LOOPS.....
%
  

	  for  i=1:ie
	  for  j=1:je
	  BHzS=BHZ(i,j);
	  
	  BHZ(i,j)=axh(i)*BHZ(i,j)			     ...
              +(EX(i,j+1)-EX(i,j)-EY(i+1,j)+EY(i,j))*dt_over_dx/bxh(i);
	 
	  HZ(i,j)=ayh(j)*HZ(i,j)			     ...
             +(BHZ(i,j)-BHzS)*muer(i,j)/byh(j)	;	   % epsrfb(i,j)
	  
     end
     end




 %-------------
 %          .....incident plane wave reference.....
      for  j=1:incm
     hxi(j) = hxi(j) + RC*(ezi(j)-ezi(j+1));
 end
%


 %          .....hz wave source condition, left and right.....
      for  j=jl:jn
        dis1 = wb*((j)+0.5-yo) + txhz1;
        dis2 = wb*((j)+0.5-yo) + txhz2;
        l1 = floor(dis1);
        l2 = floor(dis2);
        eyin1 = ((dis1-(l1))*ezi(l1+7)+			   ...
               ((l1)+1.0-dis1)*ezi(l1+6)) * (-cph);
        eyin2 = ((dis2-(l2))*ezi(l2+7)+			   ...
               ((l2)+1.0-dis2)*ezi(l2+6)) * (-cph);

			   	

        HZ(igaps,j) = HZ(igaps,j) + RC*eyin1;

	

        HZ(im,j) = HZ(im,j) - RC*eyin2;
    end
%
%          .....hz wave source condition, front and back.....
      for  i=il:in
        dis1 = wa*((i)+0.5-xo) + tyhz1;
        dis2 = wa*((i)+0.5-xo) + tyhz2;
        l1 =floor(dis1);
        l2 =floor(dis2);
        exin1 = ((dis1-(l1))*ezi(l1+7)+			   ...
               ((l1)+1.0-dis1)*ezi(l1+6)) * sph;
        exin2 = ((dis2-(l2))*ezi(l2+7)+			   ...
               ((l2)+1.0-dis2)*ezi(l2+6)) * sph;
        HZ(i,jgaps) = HZ(i,jgaps) - RC*exin1;
        HZ(i,jm) = HZ(i,jm) + RC*exin2;
        end




 

%.......Next Time Step Iteration.......

   %...........  DFT...................

 % for i1=1:nfreq
 %       rd(i1,1) = sin(omega(i1)*t2);
 %       rd(i1,2) = cos(omega(i1)*t2);
 %       snorm(i1,1) = snorm(i1,1) + esource*rd(i1,1);
 %    snorm(i1,2) = snorm(i1,2) + esource*rd(i1,2);
 %end
%  for  i1=1:nfreq
 
%   for  i2=1:2
%		    for  i=iebc:ie-iebc
%            for  j=jebc:je-iebc
		
%		ExDFT(i,j,i1,i2)=ExDFT(i,j,i1,i2)+0.5*(EX(i,j)+EX(i,j+1))*rd(i1,i2);
%		EyDFT(i,j,i1,i2)=EyDFT(i,j,i1,i2)+0.5*(EY(i,j)+EY(i+1,j))*rd(i1,i2);
%    end
%end
%end
%end

 %

% ======================== file output ===================================






% -------- integrate over boundary and get pointing vector outward -----
 % next time step loop------------------------------
%***********************************************************************
%     Visualize fields
%***********************************************************************

if mod(n,4)==0;

timestep=int2str(n);

if (SurfOn ==1)
subplot(3,1,1),surf(EX');
shading flat;
caxis([Emin Emax]);
axis([1 ie 1 jb Hmin Hmax]);
colorbar;
title(['Ex at time step = ',timestep]);

subplot(3,1,2),surf(EY');
shading flat;
caxis([Emin Emax]);
axis([1 ib 1 je Hmin Hmax]);
colorbar;
title(['Ey at time step = ',timestep]);

subplot(3,1,3),surf(HZ');
shading flat;
caxis([Hmin Hmax]);
axis([1 ie 1 je -0.2 0.2]);
colorbar;
title(['Hz at time step = ',timestep]);

nn=n/4;
M(:,nn)=getframe(gcf,rect);
else
    
subplot(3,1,1),pcolor(EX');  
shading flat;
caxis([Emin Emax]);
axis([1 ie 1 jb ]);
colorbar;
axis image;
axis off;
title(['Ex at time step = ',timestep]);

subplot(3,1,2),pcolor(EY'); 
shading flat;
caxis([Emin Emax]);
axis([1 ib 1 je ]);
colorbar;
axis image;
axis off;
title(['Ey at time step = ',timestep]);

subplot(3,1,3),pcolor(HZ'); 
shading flat;
caxis([Hmin Hmax]);
axis([1 ie 1 je ]);
colorbar;
axis image;
axis off;
title(['Hz at time step = ',timestep]);

nn=n/4;
M(:,nn)=getframe(gcf,rect);    
    
end


end;

%***********************************************************************
%     END TIME-STEPPING LOOP
%***********************************************************************

  
end


%  	open(2, file='fieldHz.csv')
%	for 999 j=jebc+1,je-jebc
%	for 999 i=iebc+1,ie-iebc
%999	write(2,*) Hz(i,j)
%	close(2)
%   	open(21, file='fieldE2.csv')
%	for 9992 j=jebc+1,je-jebc
%	for 9992 i=iebc+1,ie-iebc
%9992	write(21,*) (0.5*(Ex(i,j)+Ex(i,j+1)))^2+(0.5*(Ey(i,j)+Ey(i+1,j)))^2
%	close(21)


%		      for  i1=1:nfreq
%        zz = sqrt(snorm(i1,1)^2+snorm(i1,2)^2);
%        snorm(i1,2) = -atan2(snorm(i1,1),snorm(i1,2));
%       snorm(i1,1) = zz;
%   end

  
%     for  i1=1:nfreq
%      for   j=jebc+1:je-jebc
%      for   i=iebc+1:ie-iebc

      
%        zz = sqrt(ExDFT(i,j,i1,1)^2+ExDFT(i,j,i1,2)^2)/snorm(i1,1);
%        ExDFT(i,j,i1,2) = -atan2(ExDFT(i,j,i1,1),ExDFT(i,j,i1,2))		 ...
%                         -snorm(i1,2);
%        ExDFT(i,j,i1,1) = zz;

%		zz = sqrt(EyDFT(i,j,i1,1)^2+EyDFT(i,j,i1,2)^2)/snorm(i1,1);
%        EyDFT(i,j,i1,2) = -atan2(EyDFT(i,j,i1,1),EyDFT(i,j,i1,2))		 ...
%                         -snorm(i1,2);
%        EyDFT(i,j,i1,1) = zz;
%    end
%end
%end

 %   	open(23, file='fieldE2DFT.csv')
%	for 9993 j=jebc+1,je-jebc
%	for 9993 i=iebc+1,ie-iebc
%9993	write(23,*) (ExDFT(i,j,1,1))^2+(EyDFT(i,j,1,1))^2
%	close(21)





%  stop
% end
