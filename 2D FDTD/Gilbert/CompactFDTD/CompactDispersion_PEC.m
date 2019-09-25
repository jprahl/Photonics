%---------------------------------------------------------------
%    Compact Dispersion by Su Yu-Lun
%---------------------------------------------------------------
clear all; clc;
c=3e8;
muz=4.0*pi*1.0e-7;
epsz=1.0/(c*c*muz);
epsr=1.0;      
sigm1=0.0;


lamdastart=400.e-9;  lamdaend=1000.e-9;
lamda0=500.D-9;
nfreq=50; 

nmax=10000;
i_size=50;
j_size=50;
dx=20e-9;
phi=135;

iebc=5; jebc=iebc;

NearToFar=0;  %!..................near to far field turn on/off ...............
corrdinatelamda=1;
%  problem parameter
for ProgramParameters=1:1
ie=i_size+2*iebc;
je=j_size+2*jebc;
isrc=round(ie/2);jsrc=round(je/2);
ib=ie-1;jb=je-1;
omegacenter=2*pi*c/lamda0;
dt=dx/2/c;
tpi=2*pi;      
end;
%.......initialize Frequency and Angular Frequency Arrays.......
%      data freq/fcenter/
      lamdaspan=(lamdaend-lamdastart);
if(nfreq  == 1.)
     wstart=omegacenter; %c/(lamdastart)*tpi !4.77204e15	   ! FT starting frequency
	  wspan=omegacenter; %c/(lamdaend)*tpi-c/(lamdastart)*tpi  !.254509e15
       for i=1:nfreq
	  omega(i)=omegacenter; % !   wstart+float(i-1)/float(nfreq)*wspan        !tpi*freq(i)
      freq(i)=omega(i)/tpi;
        end;
else
for i=1:nfreq
		if(corrdinatelamda == 1.)
			 lamda(i)=lamdastart+(i-1)*lamdaspan/(nfreq);
			 omega(i)=tpi*c/lamda(i);
		else  
		     wstart= c/(lamdastart)*tpi ;
             wspan= c/(lamdaend)*tpi-c/(lamdastart)*tpi;
			 omega(i)=wstart+(i-1)/(nfreq)*wspan;
        end  
			 freq(i)=omega(i)/tpi;
end;
end;       
      ndelay=(2.e-15/dt);   %!int(3.0*gwidth)
	  gwidth2=(ndelay)^2/9.0;     %!gwidth**2  

for setup=1:1

sig(1:ie,1:je)=0;
eps(1:ie,1:je)=0;

      Ex(1:ib,1:je)=0 ; 
      Ey(1:ie,1:jb)=0; 
      Hz(1:ib,1:jb)=0;  
      Ez(1:ie,1: je)=0; 
      Hx(1:ie,1: jb)=0;
      Hy(1:ib,1: je)=0;
      
%cpml
ProjExj1(1:ib,1:je)=0 ; 
ProjExj2(1:ib,1:je)=0 ; 
ProjEyi2(1:ie,1:jb)=0; 
ProjEyi1(1:ie,1:jb)=0; 
ProjHzj1(1:ib,1:jb)=0;  
ProjHzj2(1:ib,1:jb)=0;  
ProjHzi2(1:ib,1:jb)=0;  
ProjHzi1(1:ib,1:jb)=0;  

ProjEzj1(1:ie,1:je)=0 ; 
ProjEzj2(1:ie,1:je)=0 ; 
ProjEzi2(1:ie,1:je)=0; 
ProjEzi1(1:ie,1:je)=0; 
ProjHxj1(1:ie,1:jb)=0;  
ProjHxj2(1:ie,1:jb)=0;  
ProjHyi2(1:ib,1:je)=0;  
ProjHyi1(1:ib,1:je)=0;  


bwex(1:ie)=0;
cwex(1:ie)=0;
bwey(1:je)=0;
cwey(1:je)=0;
alphae(1:iebc)=0;
sigmae(1:iebc)=0;


bwhx(1:ie)=0;
cwhx(1:ie)=0;
bwhy(1:je)=0;
cwhy(1:je)=0;
alphah(1:iebc-1)=0;
sigmah(1:iebc-1)=0;

kapahx(1:ie)=1;
kapaex(1:ie)=1;
kapahy(1:je)=1;
kapaey(1:je)=1;

ehin_i(1:ie)=1;   %F_ey(1:jb,1)=0; 
hein_i(1:ib)=1;   %F_ex(1:ib,1)=0; 
ehin_j(1:je)=1;   % F_hx(1:ib,1)=0; ,F_hy(1:jb,1)=0; 
hein_j(1:jb)=1;

end;


%===============================
%                        TFSF 
%===============================
%......Total Field/Scattered Field Parameters......
igaps=iebc+14;jgaps=jebc+14;igapsa=iebc+13; jgapsa=jebc+13;
%igaps=iebc+54;jgaps=jebc+54;igapsa=iebc+13; jgapsa=jebc+13;
%.......Total Field/Scattered Field Parameters......
if (ie >= je)    maxdim=ie; 
else  maxdim=je; 
end;
 inc=4*maxdim; incm = inc-1;
 %parameters for total/scattered boundary
iL = igaps + 1; jL = jgaps + 1;im = ie - igaps; jm = je - jgaps;     in = ib - igaps;	  jn = jb - jgaps;
%!        parameters for sa planes
io = igapsa + 1; jo = jgapsa + 1;ip = ie - igapsa; jp = je - jgapsa;iq = ip - 1; jq = jp - 1;
%         .... grid parameters....
fil = (iL);      fjl = (jL);      ril = fil - 0.5;      rjl = fjl - 0.5;      fim = (im);
fjm = (jm);      rim = fim + 0.5;      rjm = fjm + 0.5;
for arraysTS=1:1
   %          .....initialize the matrices.....
  epsr1D(1:inc)=1.D0;
  ezi(1:inc) = 0.0;
  hxi(1:incm) = 0.0;
%       arrays for envelopes at sa plane
%TSi=iq-io+1;TSj=jp-jo+1;
TSi=iq;TSj=jp;
exef=zeros(TSi,nfreq,2);
eyer=zeros(TSj,nfreq,2);
exeb=zeros(TSi,nfreq,2);
eyel=zeros(TSj,nfreq,2);
hzef=zeros(TSi,nfreq,2);
hzer=zeros(TSj,nfreq,2);
hzeb=zeros(TSi,nfreq,2);
hzel=zeros(TSj,nfreq,2);
exief=zeros(TSi,nfreq,2);
eyier=zeros(TSj,nfreq,2);
exieb=zeros(TSi,nfreq,2);
eyiel=zeros(TSj,nfreq,2);
hzief=zeros(TSi,nfreq,2);
hzier=zeros(TSj,nfreq,2);
hzieb=zeros(TSi,nfreq,2);
hziel=zeros(TSj,nfreq,2);
 totalrcse=zeros(nfreq,1);
 totalrcseInc=zeros(nfreq,1);
 totalrcseExt=zeros(nfreq,1);
%        .....wave-source case selection logic.....
      cph = cos(phi * pi/180.D0);      sph = sin(phi * pi/180.D0);
      if (phi <= 90.D0)
        wa = cph;        wb = sph;        xo = fil;        yo = fjl;
      elseif (phi <= 180.D0)
        wa = -sin(pi/180.D0 * (phi-90.D0));        wb = cos(pi/180.D0 * (phi-90.D0));
        xo = fim;        yo = fjl;
      elseif (phi <= 270.D0)
        wa = -cos(pi/180.D0 * (phi-180.D0));        wb = -sin(pi/180.D0 * (phi-180.D0));
        xo = fim;        yo = fjm;
      else
        wa = sin(pi/180.D0 * (phi-270.D0));        wb = -cos(pi/180.D0 * (phi-270.D0));
        xo = fil;        yo = fjm;
      end;
      txhz1 = wa*(ril-xo);      txhz2 = wa*(rim-xo);      tyhz1 = wb*(rjl-yo);
      tyhz2 = wb*(rjm-yo);      txey1 = wa*(ril-xo) + 0.5D0;
      txey2 = wa*(rim-xo) + 0.5D0;      tyex1 = wb*(rjl-yo) + 0.5D0;
      tyex2 = wb*(rjm-yo) + 0.5D0;
	  ES33=0.D0;	      ES22 = 0.D0;
	  ES11 = 0.D0;	  ED33=0.D0;      ED22 = 0.D0;      ED11 = 0.D0;
	  Ex_old=0.D0;	  Ey_old=0.D0;	 PEy_old2=0.D0;% PEx_old2=0.D0;	  
end;
ra = dt / epsz / dx;     
rb = dt / 2.0 / muz;      
rc = dt / muz / dx; 


%===============================
%                          PML 
%===============================

ma=1;     m=3.d0;    mur=1;  rmax=1.e-14; epsr=1;
  eta=sqrt(muz*mur/epsz/epsr);
  SigmaMax=(0.8*(m+1)/(dx*(muz/epsz*epsr)^0.5)); 
AlphaMax=0.0003; 
 	  kmax=15.;
for PMLandValueinitial=1:1
  %!.......Fill grid with ether......
     FreespeceSigma=0;
      hzDa(1:ie,1:je)=1.0;
	  hzDb(1:ie,1:je)=(dt/(muz));	       
for i = 1:iebc   %=======  E  
      SigmaEbc(i) = SigmaMax * ( (iebc - i) / (iebc - 1.0) )^m ; 
      KapaEbc(i) = 1.0+(kmax-1.0)*((iebc - i) / (iebc - 1.0))^m ; 
      AlphaEbc(i) = AlphaMax*((i-1.0)/(iebc-1.0))^ma; 
      bwebc(i) = exp(-(SigmaEbc(i) / KapaEbc(i) + AlphaEbc(i))*dt/epsz) ; 
      if ( (SigmaEbc(i) == 0.0) && (AlphaEbc(i) == 0.0) && (i == iebc) )   
         cwebc(i) = 0.0; 
      else 
         cwebc(i) = SigmaEbc(i)*(bwebc(i)-1.0)/(SigmaEbc(i)+KapaEbc(i)*AlphaEbc(i)) / KapaEbc(i); 
      end 
  end 
for i = 1:iebc-1 %=======  H  
      SigmaHbc(i) = SigmaMax * ( (iebc - i - 0.5)/(iebc-1.0))^m; 
      KapaHbc(i) = 1.0+(kmax-1.0)*((iebc - i - 0.5) / (iebc - 1.0))^m ; 
      AlphaHbc(i) = AlphaMax*((i-0.5)/(iebc-1.0))^ma; 
      bwhbc(i) = exp(-(SigmaHbc(i) / KapaHbc(i) + AlphaHbc(i))*dt/epsz) ; 
      cwhbc(i) = SigmaHbc(i)*(bwhbc(i)-1.0)/ (SigmaHbc(i)+KapaHbc(i)*AlphaHbc(i)) / KapaHbc(i); 
end  


for i=1:iebc
		 el=iebc+1-i;   % iebc, iebc-1, ...1
		 er=ie-iebc+i;  % ie-iebc-1.........ie 
    
         kapaex(i)   =KapaEbc(i);
         kapaex(er) =KapaEbc(el);         
  		 bwex(i) =bwebc(i);
		 bwex(er) =bwebc(el);         
		 cwex(i) =cwebc(i);
		 cwex(er) =cwebc(el);
      end;
for i=1:iebc-1
		 el=iebc-i;   % iebc-1, iebc-1, ...1
		 er=ie-iebc+i;  % ie-iebc-1.........ib 

         kapahx(i)   =KapaHbc(i);
         kapahx(er) =KapaHbc(el);      
		 bwhx(i) =bwhbc(i);
		 bwhx(er) =bwhbc(el);         
		 cwhx(i) =cwhbc(i);
		 cwhx(er) =cwhbc(el);    
      end;

for j=1:jebc
		 el=jebc+1-j;   % iebc, iebc-1, ...1
		 er=je-jebc+j;  % ie-iebc-1.........ie 
    
         kapaey(j)   =KapaEbc(j);
         kapaey(er) =KapaEbc(el);         
  		 bwey(j) =bwebc(j);
		 bwey(er) =bwebc(el);         
		 cwey(j) =cwebc(j);
		 cwey(er) =cwebc(el);
      end;
for j=1:jebc-1
		 el=jebc-j;   % iebc-1, iebc-1, ...1
		 er=je-jebc+j;  % ie-iebc-1.........ib 

         kapahy(j)   =KapaHbc(j);
         kapahy(er) =KapaHbc(el);      
		 bwhy(j) =bwhbc(j);
		 bwhy(er) =bwhbc(el);         
		 cwhy(j) =cwhbc(j);
		 cwhy(er) =cwhbc(el);    
      end;      
   
end;  
%plot(kapae) 
%plot(ehCa)      plot(ehCb)         plot(heDa)       plot(heDb)   plot(a)

        
%---------------------------------------------------------------
%            1/k/dx
%---------------------------------------------------------------
   ehin_i(1:ie) = 1./(kapaex(1:ie)*dx);
   ehin_j(1:je) = 1./(kapaey(1:je)*dx);
   hein_i(1:ib) = 1./(kapahx(1:ib)*dx);
   hein_j(1:jb) = 1./(kapahy(1:jb)*dx);

      Ex_old(1:ib,1:je)=0 ; 
      Ey_old(1:ie,1:jb)=0; 
      DEx(1:ib,1:je)=0 ; 
      DEy(1:ie,1:jb)=0; 
      DEx_old(1:ib,1:je)=0; 
      DEy_old(1:ie,1:jb)=0; 
      DEx_old2(1:ib,1:je)=0; 
      DEy_old2(1:ie,1:jb)=0;
      
      sigmazgeox(1:ie,1:je)=0;
      sigmazgeoy(1:ie,1:je)=0;
      sigmazgeoz(1:ie,1:je)=0;

      Ez_old(1:ie,1:je)=0 ; 
      DEz(1:ie,1:je)=0 ; 
      DEz_old(1:ie,1:je)=0; 
      DEz_old2(1:ie,1:je)=0; 
 Ex0=Ex;
 Ey0=Ey;
 Ez0=Ez;     
 Hx0=Hx;
 Hy0=Hy;
 Hz0=Hz;
 
exCa(1:ib,1:je)=0;
exCb(1:ib,1:je)=0;
eyCa(1:ie,1:jb)=0;
eyCb(1:ie,1:jb)=0;
ezCa(1:ie,1:je)=0;
ezCb(1:ie,1:je)=0;
epsrxyz(1:ib,1:jb)=epsz;
epsrex(1:ib,1:je)=epsz;
epsrey(1:ie,1:jb)=epsz;
epsrez(1:ie,1:je)=epsz;

 
%=======Ag 250nm-1000nm======
epsInf= 2.3646;
Wgap= 8.7377           /27.2214D0/(2.41888D-17);
gama= 0.07489         /27.2214D0/(2.41888D-17);


cb0=dt/dx/epsz;
db0=dt/dx/muz;
dn1=(1-gama*dt/2)/(1+gama*dt/2)
dn2=epsz*Wgap^2*dt/2/(dt*gama/2+1)
betad=epsz*epsInf/dt+dn2/2
betad=1/betad
betab=(epsz*epsInf/dt-dn2/2)*betad


    geometry=1;
%Keff=0.01*omega/c;
% Slit
jinterface=jL+100*(1e-9)/dx;
radius= 50; %(nm)                ! width
thickness=0; %(nm)	        !  thickness    
distance=150; %nm
materials= 2;			% =1.D0 -->  Gold Au,  2.D0 --> PEC,   other--> Dielectric
radius = radius*(1e-9)/dx/2;  
thickness=thickness *(1e-9)/dx;distance=distance *(1e-9)/dx;

if geometry==1
   radius_a=800 ;  %nm
   radius_b=800 ;  %nm
   radius_a = radius_a*(1e-9)/dx/2    ;
   radius_b = radius_b*(1e-9)/dx/2    ;
for i=1:ib
for j=1:jb
if( i <= isrc-radius_a | i >= isrc+radius_a | j <= jsrc-radius_b | j >= jsrc+radius_b )
	  sigmazgeox(i,j)=materials;
	  sigmazgeox(i,j+1)=materials;
	  sigmazgeoy(i,j)=materials;
	  sigmazgeoy(i+1,j)=materials;
      sigmazgeoz(i,j)=materials;
      sigmazgeoz(i+1,j)=materials;
      sigmazgeoz(i,j+1)=materials;
      sigmazgeoz(i+1,j+1)=materials;
      epsrex(i,j)=epsz*1.0;
      epsrex(i,j+1)=epsz*1.0;
      epsrey(i,j)=epsz*1.0;
	  epsrey(i+1,j)=epsz*1.0;
      epsrez(i,j)=epsz*1.0;
	  epsrez(i+1,j)=epsz*1.0;
      epsrez(i,j+1)=epsz*1.0;
      epsrez(i+1,j+1)=epsz*1.0;
     end;
end;
end; %        !Rectangular waveguide
nobsmax=9;
iobs(1)=floor(isrc-radius_a+radius_a*0.5);
iobs(2)=floor(isrc-radius_a+radius_a*1);
iobs(3)=floor(isrc-radius_a+radius_a*1.5);

iobs(4)=floor(isrc-radius_a+radius_a*0.5);
iobs(5)=floor(isrc-radius_a+radius_a*1);
iobs(6)=floor(isrc-radius_a+radius_a*1.5);
iobs(7)=floor(isrc-radius_a+radius_a*0.5);
iobs(8)=floor(isrc-radius_a+radius_a*1);
iobs(9)=floor(isrc-radius_a+radius_a*1.5);

jobs(1)=floor(jsrc-radius_b+radius_b*0.5);
jobs(2)=floor(jsrc-radius_b+radius_b*0.5);
jobs(3)=floor(jsrc-radius_b+radius_b*0.5);
jobs(4)=floor(jsrc-radius_b+radius_b*1.0);
jobs(5)=floor(jsrc-radius_b+radius_b*1.0);
jobs(6)=floor(jsrc-radius_b+radius_b*1.0);
jobs(7)=floor(jsrc-radius_b+radius_b*1.5);
jobs(8)=floor(jsrc-radius_b+radius_b*1.5);
jobs(9)=floor(jsrc-radius_b+radius_b*1.5);
m=1;
n=1;
TEModes=1;
TMModes=0;
h2=((m)*pi/2./radius_a/dx)^2+((n)*pi/2./radius_b/dx)^2;
kc=sqrt(h2);
betamn=sqrt( (omegacenter/3e8)^2-kc^2)  ;  %   %%% lamda0  pulse=lamda
Keff=betamn;

if (TEModes== 1)   %%TEmode--------------------
for X=isrc-radius_a:isrc+radius_a
for Y=jsrc-radius_b:jsrc+radius_b
Ex0(X,Y)= omegacenter*muz/h2*(n*pi/2/radius_b)*cos(m*pi/2/radius_a*(X-isrc+radius_a))*sin(n*pi/2/(radius_b)*(Y-jsrc+radius_b));
Ey0(X,Y)=-omegacenter*muz/h2*(m*pi/2/radius_a)*sin(m*pi/2/radius_a*(X-isrc+radius_a))*cos(n*pi/2/(radius_b)*(Y-jsrc+radius_b));
Hx0(X,Y)=Keff/h2*(m*pi/2/radius_a)*sin(m*pi/2/radius_a*(X-isrc+radius_a))*cos(n*pi/2/(radius_b)*(Y-jsrc+radius_b));    
Hy0(X,Y)=Keff/h2*(n*pi/2/radius_b)*cos(m*pi/2/radius_a*(X-isrc+radius_a))*sin(n*pi/2/(radius_b)*(Y-jsrc+radius_b));    
Hz0(X,Y)=cos(m*pi/2/radius_a*(X-isrc+radius_a))*cos(n*pi/2/(radius_b)*(Y-jsrc+radius_b));

end  
end     
end    %%TEmode--------------------
if (TMModes == 1)  %%TMmode--------------------
for X=isrc-radius_a:isrc+radius_a
for Y=jsrc-radius_b:jsrc+radius_b
Ex0(X,Y)=sin(m*pi/2/radius_a*(X-isrc+radius_a))*sin(n*pi/2/(radius_b)*(Y-jsrc+radius_b));
Ey0(X,Y)=sin(m*pi/2/radius_a*(X-isrc+radius_a))*sin(n*pi/2/(radius_b)*(Y-jsrc+radius_b));
Ez0(X,Y)=sin(m*pi/2/radius_a*(X-isrc+radius_a))*sin(n*pi/2/(radius_b)*(Y-jsrc+radius_b));
Hx0(X,Y)=sin(m*pi/2/radius_a*(X-isrc+radius_a))*sin(n*pi/2/(radius_b)*(Y-jsrc+radius_b));
Hy0(X,Y)=sin(m*pi/2/radius_a*(X-isrc+radius_a))*sin(n*pi/2/(radius_b)*(Y-jsrc+radius_b));
end  
end     

end           %%TMmode--------------------

end
 Keff=betamn*dx;   

exCa(1:ib,1:je)=(1.0-FreespeceSigma*dt ./ (2.0*epsrex(1:ib,1:je))) ./(1.0 +FreespeceSigma* dt ./ (2.0*epsrex(1:ib,1:je)));
exCb(1:ib,1:je)=(dt./(epsrex(1:ib,1:je))) ./ (1.0 + FreespeceSigma*dt ./ (2.0*epsrex(1:ib,1:je)));
eyCa(1:ie,1:jb)=(1.0-FreespeceSigma*dt./ (2.0*epsrey(1:ie,1:jb))) ./(1.0 +FreespeceSigma* dt ./ (2.0*epsrey(1:ie,1:jb)));
eyCb(1:ie,1:jb)=(dt./(epsrey(1:ie,1:jb))) ./ (1.0 + FreespeceSigma*dt ./ (2.0*epsrey(1:ie,1:jb)));
ezCa(1:ie,1:je)=(1.0-FreespeceSigma*dt./ (2.0*epsrez(1:ie,1:je))) ./(1.0 +FreespeceSigma* dt ./ (2.0*epsrez(1:ie,1:je)));
ezCb(1:ie,1:je)=(dt./(epsrez(1:ie,1:je))) ./ (1.0 + FreespeceSigma*dt ./ (2.0*epsrez(1:ie,1:je)));
figure(99)
 pcolor((epsrex(1:ib,1:je)+sigmazgeox(1:ib,1:je))');shading flat;daspect([1 1 1]);colorbar; 
 %stop;%
 %pause;%stop;
 
 ExDFTobs(1:nobsmax,1:nfreq,1:2)=0.;
 EyDFTobs(1:nobsmax,1:nfreq,1:2)=0.;
 EzDFTobs(1:nobsmax,1:nfreq,1:2)=0.;
 HxDFTobs(1:nobsmax,1:nfreq,1:2)=0.;
 HyDFTobs(1:nobsmax,1:nfreq,1:2)=0.;
 HzDFTobs(1:nobsmax,1:nfreq,1:2)=0.;
 snorm(1:nfreq,1:2)=0.;
 %Hx=Hx0;Hy=Hy0;Hz=Hz0;
for  n=1:nmax
		  %n   
      t2=(n)*dt;tt2=(n-ndelay)*dt;



  gaussianProfile=exp(-1.0*((n-ndelay)^2)/gwidth2);

 sinusoid=sin(c/500e-9*2*pi*t2) ; 

        sinegauss=sinusoid*gaussianProfile;
		esource=sinegauss/377;          
        
%=======================
%        update  Ex Ey Ez
%=======================
      Ex_old2(1:ib,1:je)=Ex_old(1:ib,1:je); 
      Ex_old(1:ib,1:je)=Ex(1:ib,1:je); 
      DEx_old2(1:ib,1:je)=DEx_old(1:ib,1:je); 
      DEx_old(1:ib,1:je)=DEx(1:ib,1:je); 
  
for i=1:ib   
  for  j=2:jb
         if (sigmazgeox(i,j) == 1)
     Ex(i,j)=betab*Ex(i,j)+((Hz(i,j)-Hz(i,j-1))*ehin_j(j) -(dn1+1)*0.5*DEx(i,j))*betad +Keff*Hy(i,j)*exCb(i,j);
              DEx(i,j)=dn1*DEx_old(i,j)+dn2*(Ex(i,j)+Ex_old(i,j)) ;
              ProjExj1(i,j)=bwey(j)*ProjExj1(i,j)+cwey(j) *(Hz(i,j)-Hz(i,j-1))/dx;    
              Ex(i,j) = Ex(i,j) + exCb(i,j)*ProjExj1(i,j)*betad;
         elseif (sigmazgeox(i,j) == 2)       
     Ex(i,j)=0;
              ProjExj1(i,j)=bwey(j)*ProjExj1(i,j)+cwey(j) *(Hz(i,j)-Hz(i,j-1) )/dx;    
              Ex(i,j) = Ex(i,j) + exCb(i,j)*ProjExj1(i,j);
         else
     Ex(i,j)=Ex(i,j)+(Hz(i,j)-Hz(i,j-1)+Keff*Hy(i,j))*cb0  ;     
              %ProjExj1(i,j)=bwey(j)*ProjExj1(i,j)+cwey(j) *(Hz(i,j)-Hz(i,j-1))/dx;    
              %Ex(i,j) = Ex(i,j) + exCb(i,j)*ProjExj1(i,j);
         end             
   end;    
end;

      Ey_old2(1:ie,1:jb)=Ey_old(1:ie,1:jb); 
      Ey_old(1:ie,1:jb)=Ey(1:ie,1:jb); 
      DEy_old2(1:ie,1:jb)=DEy_old(1:ie,1:jb);  
      DEy_old(1:ie,1:jb)=DEy(1:ie,1:jb);        
for i=2:ib      
  for  j=1:jb 
         if (sigmazgeoy(i,j) == 1)
      Ey(i,j)=betab*Ey(i,j)+((Hz(i-1,j)-Hz(i,j))*ehin_i(i)  -(dn1+1)*0.5*DEy(i,j) )*betad -Keff*Hx(i,j)*eyCb(i,j);
              DEy(i,j)=dn1*DEy_old(i,j)+dn2*(Ey(i,j)+Ey_old(i,j));
              ProjEyi1(i,j) = bwex(i)*ProjEyi1(i,j) - cwex(i)*(Hz(i,j) - Hz(i-1,j))/dx;  
              Ey(i,j) = Ey(i,j) +  eyCb(i,j)*ProjEyi1(i,j)*betad;
         elseif (sigmazgeoy(i,j) == 2)    
      Ey(i,j)=0;
              ProjEyi1(i,j) = bwex(i)*ProjEyi1(i,j) - cwex(i)*(Hz(i,j) - Hz(i-1,j))/dx;  
              Ey(i,j) = Ey(i,j) +  eyCb(i,j)*ProjEyi1(i,j);
         else
      Ey(i,j)=Ey(i,j)+(-Hz(i,j)+Hz(i-1,j)-Keff*Hx(i,j) )*cb0;    
             % ProjEyi1(i,j) = bwex(i)*ProjEyi1(i,j) - cwex(i)*(Hz(i,j) - Hz(i-1,j))/dx;  
              %Ey(i,j) = Ey(i,j) +  eyCb(i,j)*ProjEyi1(i,j);
         end    
   end;    
end;

      Ez_old2(1:ib,1:jb)=Ez_old(1:ib,1:jb); 
      Ez_old(1:ib,1:jb)=Ez(1:ib,1:jb); 
      DEz_old2(1:ib,1:jb)=DEz_old(1:ib,1:jb);  
      DEz_old(1:ib,1:jb)=DEz(1:ib,1:jb); 
for i = 2:ib 
     for j = 2:jb 
          if (sigmazgeoz(i,j) == 1)      
      Ez(i,j)=betab*Ez(i,j)+( ((Hy(i,j) - Hy(i-1,j))*ehin_i(i) + (Hx(i,j-1) - Hx(i,j))*ehin_j(j)) -(dn1+1)*0.5*DEz(i,j) )*betad;
              DEz(i,j)=dn1*DEz_old(i,j)+dn2*(Ez(i,j)+Ez_old(i,j));
              ProjEzj1(i,j)=bwey(j)*ProjEzj1(i,j)+cwey(j) *(Hx(i,j-1) - Hx(i,j))/dx;    
              Ez(i,j) = Ez(i,j) +ezCb(i,j)*ProjEzj1(i,j)*betad;
              ProjEzi1(i,j) = bwex(i)*ProjEzi1(i,j) + cwex(i)*(Hy(i,j) - Hy(i-1,j))/dx;  
              Ez(i,j) = Ez(i,j) +ezCb(i,j)*ProjEzi1(i,j)*betad;
         elseif (sigmazgeoz(i,j) == 2)    
      Ez(i,j)=0;
              ProjEzj1(i,j)=bwey(j)*ProjEzj1(i,j)+cwey(j) *(Hx(i,j-1) - Hx(i,j))/dx;    
              Ez(i,j) = Ez(i,j) +ezCb(i,j)*ProjEzj1(i,j);
              ProjEzi1(i,j) = bwex(i)*ProjEzi1(i,j) + cwex(i)*(Hy(i,j) - Hy(i-1,j))/dx;  
              Ez(i,j) = Ez(i,j) +ezCb(i,j)*ProjEzi1(i,j);
          else
      Ez(i,j) = Ez(i,j) + cb0* ((Hy(i,j) - Hy(i-1,j) + (Hx(i,j-1) - Hx(i,j))));
             % ProjEzj1(i,j)=bwey(j)*ProjEzj1(i,j)+cwey(j) *(Hx(i,j-1) - Hx(i,j))/dx;    
              %Ez(i,j) = Ez(i,j) +ezCb(i,j)*ProjEzj1(i,j);
             % ProjEzi1(i,j) = bwex(i)*ProjEzi1(i,j) + cwex(i)*(Hy(i,j) - Hy(i-1,j))/dx;  
             % Ez(i,j) = Ez(i,j) +ezCb(i,j)*ProjEzi1(i,j);
          end
	end 
end 
 
 %  Ez=Ez+Ez0*esource; 
  for iff=1:nfreq
      rd(iff,1)=esource*cos(omega(iff)*t2);
      rd(iff,2)=esource*sin(omega(iff)*t2);
      snorm(iff,1)=snorm(iff,1)+rd(iff,1);
      snorm(iff,2)=snorm(iff,2)+rd(iff,2);
      for i1=1:2
           for nobs=1:nobsmax
              ExDFTobs(nobs,iff,i1)=ExDFTobs(nobs,iff,i1)+Ex(iobs(nobs),jobs(nobs))*rd(iff,i1);
              EyDFTobs(nobs,iff,i1)=EyDFTobs(nobs,iff,i1)+Ey(iobs(nobs),jobs(nobs))*rd(iff,i1);
              EzDFTobs(nobs,iff,i1)=EzDFTobs(nobs,iff,i1)+Ez(iobs(nobs),jobs(nobs))*rd(iff,i1);
              HxDFTobs(nobs,iff,i1)=HxDFTobs(nobs,iff,i1)+Hx(iobs(nobs),jobs(nobs))*rd(iff,i1);
              HyDFTobs(nobs,iff,i1)=HyDFTobs(nobs,iff,i1)+Hy(iobs(nobs),jobs(nobs))*rd(iff,i1);
              HzDFTobs(nobs,iff,i1)=HzDFTobs(nobs,iff,i1)+Hz(iobs(nobs),jobs(nobs))*rd(iff,i1);
           end
      end
  end
%=======================
%        update  Hx Hy Hz
%=======================
for i = 2:ib 
  for j = 1:jb 
     Hx(i,j)=Hx(i,j) + db0* (Ez(i,j) - Ez(i,j+1)  +Keff*Ey(i,j)) ;     
              %ProjHxj1(i,j) = bwhy(j)*ProjHxj1(i,j) + cwhy(j)*(Ez(i,j) - Ez(i,j+1))/dx;%  Hz, y direction
              %Hx(i,j) = Hx(i,j) + hzDb(i,j)*ProjHxj1(i,j);
	end 
end 
for i = 1:ib 
	 for j = 2:jb 
      Hy(i,j)=Hy(i,j) + db0* ( Ez(i+1,j) - Ez(i,j) - Keff*Ex(i,j));    
             % ProjHyi1(i,j) = bwhx(i)*ProjHyi1(i,j) + cwhx(i) *(Ez(i+1,j) - Ez(i,j))/dx;%   Hz, x direction
             % Hy(i,j) = Hy(i,j) +hzDb(i,j)*ProjHyi1(i,j);
	 end  
  end 
for i = 1:ib 
  for  j=1:jb	  
	  Hz(i,j)=Hz(i,j)+db0*(Ex(i,j+1)-Ex(i,j)-(Ey(i+1,j)-Ey(i,j)));
             % ProjHzi1(i,j) = bwhx(i)*ProjHzi1(i,j) - cwhx(i) *(Ey(i+1,j) - Ey(i,j))/dx;%   Hz, x direction
             % Hz(i,j) = Hz(i,j) +hzDb(i,j)*ProjHzi1(i,j);
             % ProjHzj1(i,j) = bwhy(j)*ProjHzj1(i,j) + cwhy(j)*(Ex(i,j+1) - Ex(i,j))/dx;%  Hz, y direction
             % Hz(i,j) = Hz(i,j) + hzDb(i,j)*ProjHzj1(i,j);
   end;
end

Hz=Hz+Hz0*esource; 

if mod(n,1)==0;
  n  
figure(9);
    timestep=int2str(n);

subplot(2,3,1);pcolor(Hx'); shading flat; daspect([1 1 1]);colorbar('horiz');%caxis([-0.001 0.001]);% title('Hx');
subplot(2,3,2);pcolor(Hy'); shading flat; daspect([1 1 1]);colorbar('horiz');%caxis([-0.001 0.001]);% title('Hy'); 
subplot(2,3,3);pcolor(Ez'); shading flat; daspect([1 1 1]);colorbar('horiz');%caxis([-0.5 0.5]);        %title('Ez');

subplot(2,3,4);pcolor(Ex'); shading flat; daspect([1 1 1]);colorbar('horiz');%caxis([-0.5 0.5]);%title('Ex');
subplot(2,3,5);pcolor(Ey'); shading flat; daspect([1 1 1]);colorbar('horiz');%caxis([-0.5 0.5]);%title('Ey');
subplot(2,3,6);pcolor(Hz'); shading flat; daspect([1 1 1]);colorbar('horiz');%caxis([-0.001 0.001]);%title('Hz');
 %pcolor( (Ex(1:ib,1:jb).^2+Ey(1:ib,1:jb).^2)'); shading flat;colorbar('horiz');
% title(['Hz at time step = ',timestep]);%colorbar;axis image;axis on;  % caxis([-0.01 0.01]);
%hold off;
% Hz(isrc,jsrc)

end


end;

for iff=1:nfreq
      
      snorm(iff,1)=(sqrt(snorm(iff,1)^2+snorm(iff,2)^2));

           for nobs=1:nobsmax
              ExDFTobs(nobs,iff,1)=sqrt(ExDFTobs(nobs,iff,1)^2+ExDFTobs(nobs,iff,2)^2)/snorm(iff,1);
              EyDFTobs(nobs,iff,1)=sqrt(EyDFTobs(nobs,iff,1)^2+EyDFTobs(nobs,iff,2)^2)/snorm(iff,1);
              EzDFTobs(nobs,iff,1)=sqrt(EzDFTobs(nobs,iff,1)^2+EzDFTobs(nobs,iff,2)^2)/snorm(iff,1);
              HxDFTobs(nobs,iff,1)=sqrt(HxDFTobs(nobs,iff,1)^2+HxDFTobs(nobs,iff,2)^2)/snorm(iff,1);
              HyDFTobs(nobs,iff,1)=sqrt(HyDFTobs(nobs,iff,1)^2+HyDFTobs(nobs,iff,2)^2)/snorm(iff,1);
              HzDFTobs(nobs,iff,1)=sqrt(HzDFTobs(nobs,iff,1)^2+HzDFTobs(nobs,iff,2)^2)/snorm(iff,1);
           end
end
ExDFT(1:nfreq)=0.;EyDFT(1:nfreq)=0.;EzDFT(1:nfreq)=0.;
HxDFT(1:nfreq)=0.;HyDFT(1:nfreq)=0.;HzDFT(1:nfreq)=0.;
for nobs=1:nobsmax
    for iff=1:nfreq
ExDFT(iff)= ExDFT(iff)+ExDFTobs(nobs,iff,1);
EyDFT(iff)= EyDFT(iff)+EyDFTobs(nobs,iff,1);
EzDFT(iff)= EzDFT(iff)+EzDFTobs(nobs,iff,1);
HxDFT(iff)= HxDFT(iff)+HxDFTobs(nobs,iff,1);
HyDFT(iff)= HyDFT(iff)+HyDFTobs(nobs,iff,1);
HzDFT(iff)= HzDFT(iff)+HzDFTobs(nobs,iff,1);
    end;end;

figure(32);subplot(3,2,1);plot(lamda,ExDFTobs(1,:,1));
subplot(3,2,2);plot(lamda,EyDFTobs(1,:,1));
subplot(3,2,3);plot(lamda,EzDFTobs(1,:,1));
subplot(3,2,4);plot(lamda,HxDFTobs(1,:,1));
subplot(3,2,5);plot(lamda,HyDFTobs(1,:,1));
subplot(3,2,6);plot(lamda,HzDFTobs(1,:,1));
figure(4);plot(lamda,ExDFT.*EyDFT.*HzDFT.*HxDFT.*HyDFT);


