% TE code Hz, Ex, Ey, global error
% All UPML TE code by Gilbert Chang 
 
		
clear all;
       pi=3.141592654;
       cc=3.000000000+8;
       muz=4.0*pi*1.0-7;
       epsz=1.0/(muz*cc*cc);
       dtr=pi/180.0;
       rtd=180.0/pi;
       tpi=2.0*pi;
  	
      iebc=16;jebc=16;
      rmax=1.e-14;orderbc=2.0;mediabc=2;
	       

	 dx=1.e-9;
	 dt=dx/(2.0*cc);
     nmax=500;
   	 dt_over_dx=dt/dx;
	
   	 i_size=100;
	 j_size=50;
	 i_size2=1000;
	 j_size2=500;      


	% total gird size including PML
 	 ie=i_size+2*iebc;je=j_size+2*jebc;
     ie2=i_size2+2*iebc;je2=j_size2+2*jebc;

     	 ib=ie+1;jb=je+1;
         ib2=ie2+1;jb2=je2+1;


 
	 isrc=ie/2;jsrc=je/2;
     isrc2=ie2/2;jsrc2=je2/2;

%.................Grid Arrays...............................
%	 


    







%............................................................
      real EPS(MEDIA),MUR(MEDIA),	 ...
     
       HZ(ie,je),BHZ(ie,je),EX(ie,jb),EY(ib,je),		   ...
	   DEX(ie,jb),DEY(ib,je), ...
       axe(ib),aye(jb),axh(ie),ayh(je),			...
	   bxe(ib),bye(jb),bxh(ie),byh(je), ...
	   cxe(ib),cye(jb),cxh(ie),cyh(je), ...
	   epsrex(ie,jb),epsrey(ib,je),muer(ib,jb)
 	  real  HZ2(ie2,je2),BHZ2(ie2,je2),EX2(ie2,jb2),EY2(ib2,je2),		   ...
	   DEX2(ie2,jb2),DEY2(ib2,je2), ...
       axe2(ib2),aye2(jb2),axh2(ie2),ayh2(je2),			...
	   bxe2(ib2),bye2(jb2),bxh2(ie2),byh2(je2), ...
	   cxe2(ib2),cye2(jb2),cxh2(ie2),cyh2(je2), ...
	   epsrex2(ie2,jb2),epsrey2(ib2,je2),muer2(ie2,je2)

  
%.......General CA, DA, CB, DB Arrays.......
      real axebc(ibbc,mediabc),ayebc(jebc,mediabc)
      real bxebc(ibbc,mediabc),byebc(jebc,mediabc)
	  real cxebc(ibbc,mediabc),cyebc(jebc,mediabc)
	  real axhbc(ibbc,mediabc),ayhbc(jebc,mediabc)
      real bxhbc(ibbc,mediabc),byhbc(jebc,mediabc)
	  real cxhbc(ibbc,mediabc),cyhbc(jebc,mediabc)

% Different Fourier Transform arrays
	   
       
%.......Initialize Update Coefficients.......
%.......Specify Electrical Properties......
	eps(1)=rns*rns
	eps(2)=rnf*rnf
DATA MUR/1.0,1.0/
     

	 write(*,*) isrc,jsrc,isrc2,jsrc2
	 

	  kx=1.0
	  ky=1.0

   
	  
%.......Fill grid with ether......
      do 820 i=1,ib
      axe(i)=1.0
	  bxe(i)=1.0
820	  cxe(i)=1.0
      do 821 j=1,jb
	  aye(j)=1.0
	  bye(j)=1.0
821	  cye(j)=1.0
	  do 822 i=1,ie
      axh(i)=1.0
	  bxh(i)=1.0
822   cxh(i)=1.0
	  do 823 j=1,je
      ayh(j)=1.0
	  byh(j)=1.0
823   cyh(j)=1.0

      do 830 i=1,ie
	  do 830 j=1,jb

      epsrex(i,j)=1.0/epsz 
 830 continue
	  do 832 i=1,ib
	  do 832 j=1,je

      epsrey(i,j)=1.0/epsz 
 832 continue

    write(*,*) 1.0/epsz
	  do 831 i=1,ie
	  do 831 j=1,je

 831     muer(i,j)=1.0/muz 

   write(*,*) 1.0/muz

   %.......Fill larger grid with ether......
      do 8202 i=1,ib2
      axe2(i)=1.0
	  bxe2(i)=1.0
8202  cxe2(i)=1.0
      do 8212 j=1,jb2
	  aye2(j)=1.0
	  bye2(j)=1.0
8212  cye2(j)=1.0
	  do 8222 i=1,ie2
      axh2(i)=1.0
	  bxh2(i)=1.0
8222  cxh2(i)=1.0
	  do 8232 j=1,je2
      ayh2(j)=1.0
	  byh2(j)=1.0
8232  cyh2(j)=1.0

      do 8302 i=1,ie2
	  do 8302 j=1,jb2

      epsrex2(i,j)=1.0/epsz 
 8302 continue
	  do 8322 i=1,ib2
	  do 8322 j=1,je2

      epsrey2(i,j)=1.0/epsz 
 8322 continue

    write(*,*) 1.0/epsz
	  do 8312 i=1,ie2
	  do 8312 j=1,je2

 8312     muer2(i,j)=1.0/muz 

   write(*,*) 1.0/muz


% Problem Geometry Specify
%
%...................GEOMETRY SPECIFICATION...........................
%
%        .....specify the random medium ....
%  
%      
% 
% -------- output geometry file -----------------
    
 % do 210 i=1,ie
 % do 210 j=1,je
 % epsr(i,j)=1./(epsz)
 % 210 continue
 % mur=
%
%
%.....................Fill PML Grid.........................
%


 %----------------------------
 	  kmax=1.
      do 8230 m=1,mediabc
      delbc=float(iebc)*dx

	  eta=sqrt(muz*mur(m)/epsz/eps(m))
      sigmam=-log(rmax)*(orderbc+1.0)/(2.0*delbc*eta*eps(m))

	   

      bcfactor=sigmam/(dx*(delbc**orderbc)*(orderbc+1.0))
	  kfactor=(Kmax-1.)/(dx*(delbc**orderbc)*(orderbc+1.0)) 

	  
%
%       .........
%
        do 824 i=1,iebc
          x1=dx*float(i)
          x2=dx*(float(i)-1.0)
          sigma=bcfactor*(x1**(orderbc+1.0)-x2**(orderbc+1.0))
		  kapa=1.+kfactor*(x1**(orderbc+1.0)-x2**(orderbc+1.0))

		
          axhbc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz)
          ayhbc(i,m)= axhbc(i,m)
		  bxhbc(i,m)=(kapa+sigma*dt/2.0/epsz)
		  byhbc(i,m)= bxhbc(i,m)
		  cxhbc(i,m)=(kapa-sigma*dt/2.0/epsz)
		  cyhbc(i,m)= cxhbc(i,m)
  824   continue
          x1=dx*0.50
          sigma=bcfactor*x1**(orderbc+1.0)
		  kapa=1.+kfactor*x1**(orderbc+1.0)
		 

          axebc(1,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz)
          ayebc(1,m)= axebc(1,m)
		  bxebc(1,m)=(kapa+sigma*dt/2.0/epsz)
		  byebc(1,m)= bxebc(1,m)
		  cxebc(1,m)=(kapa-sigma*dt/2.0/epsz)
		  cyebc(1,m)= cxebc(1,m)
        
        do 825 i=2,iebc
          x1=dx*(float(i)-0.50)
          x2=dx*(float(i)-1.50)
          sigma=bcfactor*(x1**(orderbc+1.0)-x2**(orderbc+1.0))
		  kapa=1.+kfactor*(x1**(orderbc+1.0)-x2**(orderbc+1.0))
	

          axebc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz)
          ayebc(i,m)=(kapa-sigma*dt/2.0/epsz)/(kapa+sigma*dt/2.0/epsz)
		  bxebc(i,m)=(kapa+sigma*dt/2.0/epsz)
		  byebc(i,m)=(kapa+sigma*dt/2.0/epsz)
		  cxebc(i,m)=(kapa-sigma*dt/2.0/epsz)
		  cyebc(i,m)=(kapa-sigma*dt/2.0/epsz)

  825   continue 
  8230 continue











 %----------------------------




%                  ....FRONT REGION.....
%
%.......Fill EZ and HX with PML Media in Front Region Everywhere.......
	  do 8260 i=1,iebc
		 i0=iebc+1-i
	   	 i1=ibbc+1-i
		 i2=ie-iebc+i
		 i3=ie-iebc+i
	   	 axh(i0) =axhbc(i,1)
		 axh(i2) =axhbc(i,1)
		 bxh(i0) =bxhbc(i,1)
		 bxh(i2) =bxhbc(i,1)
		 cxh(i0) =cxhbc(i,1)
		 cxh(i2) =cxhbc(i,1)

		 axe(i1) =axebc(i,1)
		 axe(i3) =axebc(i,1)
		 bxe(i1) =bxebc(i,1)
		 bxe(i3) =bxebc(i,1)
		 cxe(i1) =cxebc(i,1)
		 cxe(i3) =cxebc(i,1)

  8260 continue

	    do 8261 j=1,jebc
	   	 j0=jebc+1-j
		 j1=jbbc+1-j
		 j2=je-jebc+j
		 j3=je-jebc+j
	   	 ayh(j0) =ayhbc(j,1)
		 ayh(j2) =ayhbc(j,1)
		 byh(j0) =byhbc(j,1)
		 byh(j2) =byhbc(j,1)
		 cyh(j0) =cyhbc(j,1)
		 cyh(j2) =cyhbc(j,1)
		 
		 aye(j1) =ayebc(j,1)
		 aye(j3) =ayebc(j,1)
		 bye(j1) =byebc(j,1)
		 bye(j3) =byebc(j,1)
		 cye(j1) =cyebc(j,1)
		 cye(j3) =cyebc(j,1)

  8261 continue
%open(145,file='axe.csv')
%do 155 i=1,ib
%155  write(145,199) 	  axe(i),',',axh(i),',',bxe(i),',',bxh(i),',',cxe(i),',',cxh(i)
%close(145)
%open(146,file='aye.csv')
%do 156 j=1,jb
%156  write(146,199) 	  aye(j),',',ayh(j),',',bye(j),',',byh(j),',',cye(j),',',cyh(j)
%close(146)    
%199 format(e14.6e3,a3,e14.6e3,a3,e14.6e3,a3,e14.6e3,a3,e14.6e3,a3,e14.6e3)
  %%%%%%%---------- PML for larger grid
 	  do 82602 i=1,iebc
		 i0=iebc+1-i
	   	 i1=ibbc+1-i
		 i2=ie2-iebc+i
		 i3=ie2-iebc+i
	   	 axh2(i0) =axhbc(i,1)
		 axh2(i2) =axhbc(i,1)
		 bxh2(i0) =bxhbc(i,1)
		 bxh2(i2) =bxhbc(i,1)
		 cxh2(i0) =cxhbc(i,1)
		 cxh2(i2) =cxhbc(i,1)

		 axe2(i1) =axebc(i,1)
		 axe2(i3) =axebc(i,1)
		 bxe2(i1) =bxebc(i,1)
		 bxe2(i3) =bxebc(i,1)
		 cxe2(i1) =cxebc(i,1)
		 cxe2(i3) =cxebc(i,1)

  82602 continue

	    do 82612 j=1,jebc
	   	 j0=jebc+1-j
		 j1=jbbc+1-j
		 j2=je2-jebc+j
		 j3=je2-jebc+j
	   	 ayh2(j0) =ayhbc(j,1)
		 ayh2(j2) =ayhbc(j,1)
		 byh2(j0) =byhbc(j,1)
		 byh2(j2) =byhbc(j,1)
		 cyh2(j0) =cyhbc(j,1)
		 cyh2(j2) =cyhbc(j,1)
		 
		 aye2(j1) =ayebc(j,1)
		 aye2(j3) =ayebc(j,1)
		 bye2(j1) =byebc(j,1)
		 bye2(j3) =byebc(j,1)
		 cye2(j1) =cyebc(j,1)
		 cye2(j3) =cyebc(j,1)

  82612 continue

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
      do 805 i=1,ie
      do 805 j=1,je

	Hz(1:ie,1:je)=0.0
	BHz(i,j)=0.0
805   continue

      
	Ex(1:ie,1:jb) = 0.0
	DEx(1:ie,1:jb)=0.0

    
	Ey(1:ib,1:je) = 0.0
	DEy(1:ib,1:je)=0.0
 

%.......Initialize the larger grid Matrices.......
    

	Hz2(1:ie2,1:je2)=0.0
	BHz2(1:ie2,1:je2)=0.0

    
	Ex2(1:ie2,1:jb2) = 0.0
	DEx2(1:ie2,1:jb2)=0.0

     
	 Ey2(1:ib2,1:je2) = 0.0
	DEy2(1:ib2,1:je2)=0.0
  

 
   
 
 
%...............TIME-STEPPING LOOP...........................
%
	nlimit=40  %1.e-9/2.5e-11


	       DO 500 N=1,nmax
      t2=float(n)*dt
%
%
%
%................Hz  UPDATE...................................
%
%         .....MAIN DO LOOPS.....
%
  

	  for i=1:ie
	  for j=1:je
	  BHzS=BHz(i,j)
	  
	  BHz(i,j)=axh(i)*BHz(i,j)			  ...
              +(EX(i,j+1)-EX(i,j)-EY(i+1,j)+EY(i,j))*dt_over_dx/bxh(i)
	 
	  HZ(i,j)=ayh(j)*HZ(i,j)			     ...
             +(BHZ(i,j)-BHzS)*muer(i,j)/byh(j)		   
	  
      end
      end

 %----------- Hz update for larger grid ---------------
 
      for i=1:ie2
	  for j=1:je2
	  BHS2=BHz2(i,j)
      BHz2(i,j)=axh2(i)*BHz2(i,j)			  ...
              -(EY2(i+1,j)-EY2(i,j)-EX2(i,j+1)+EX2(i,j))*dt_over_dx/bxh2(i)
	  HZ2(i,j)=ayh2(j)*HZ2(i,j)			     ...
             +(BHZ2(i,j)-BHS2)*muer2(i,j)/byh2(j)		   
      end
      end


%     Soft source
 
      if (n <= nlimit) t

        esource  =(10.0-15.0*cos(Pi*0.05*(n))+6.0*cos(tpi*0.05*(n))	...
		           -cos(3.0*Pi*0.05*(n)))/32.0;
      else
	    esource=0.0;
      end
	    

	 

	  Hz(isrc,jsrc)= esource;
	  Hz2(isrc2,jsrc2)= esource;

      
%.............EX and EY UPDATE...................................
%
%         .....MAIN DO LOOPS.....
%
      
      for i=1:ie
	  for j=2:je

	  	 DEXS=DEX(i,j)
	  
        DEX(i,j)=aye(j)*DEX(i,j)+				...
                  (HZ(i,j)-HZ(i,j-1))*dt_over_dx/bye(j)
	
		EX(i,j)=EX(i,j)+				  ...
		            (bxh(i)*DEX(i,j)-cxh(i)*DEXS)*epsrex(i,j)			   
      end
      end



      
      for i=2:ie
	  for j=1:je
	  	DEYS=DEY(i,j)
	
         DEY(i,j)=DEY(i,j)-					 ...
                  (HZ(i,j)-HZ(i-1,j))*dt_over_dx
		
		 EY(i,j)=axe(i)*EY(i,j)+        ...
		          (byh(j)*DEY(i,j)-cyh(j)*DEYS)*epsrey(i,j)/bxe(i)
      end
      end
 
       
% ------------ Ex Ey update for larger grid size
   
      for i=1:ie2
	  for j=2:je2

		 DEXS2=DEX2(i,j)
	
	   
        DEX2(i,j)=aye2(j)*DEX2(i,j)+				...
                  (HZ2(i,j)-HZ2(i,j-1))*dt_over_dx/bye2(j)
	
		EX2(i,j)=EX2(i,j)+				  ...
		            (bxh2(i)*DEX2(i,j)-cxh2(i)*DEXS2)*epsrex2(i,j)	
      end
      end
 
      for j=2:je2
      for i=1:ie2
		
		DEYS2=DEY2(i,j)
	
         DEY2(i,j)=DEY2(i,j)-					 ...
                  (HZ2(i,j)-HZ2(i-1,j))*dt_over_dx
		
		 EY2(i,j)=axe2(i)*EY2(i,j)+        ...
		          (byh2(j)*DEY2(i,j)-cyh2(j)*DEYS2)*epsrey2(i,j)/bxe2(i)
      end
      end

 
%.......Next Time Step Iteration.......


% ======================== file output ===================================





 % calculate global error
   for i=1:100
       for j=1:50
   Globalerror(n)=Globalerror(n)+(HZ(isrc-50+i,jsrc-25+j)-HZ2(isrc2-50+i,jsrc2-25+j)).^2;
       end
   end
	


 

% -------- integrate over boundary and get pointing vector outward -----
 % next time step loop------------------------------


  
 end

  
 