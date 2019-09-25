 % TM code Ez, Hx, Hy
 % All UPML 2D TM code by Gilbert Chang
 	
clear all;	
       pi=3.141592654;
       cc=3.000000000+8;
       muz=4.0*pi*1.0-7;
       epsz=1.0/(muz*cc*cc);
       dtr=pi/180.0;                    %
       rtd=180.0/pi;                    %
       tpi=2.0*pi;
  	
      iebc=16;jebc=16;                  %PML boundary
      rmax=1.e-14;                      %Reflection factor? eq7.59 R(0)
      orderbc=2.0;     % m - large m yields a sigma distribution that is relatively flat near the PML surface
                       %deeper within the PML sigma increases more rapidly for small m 
                       % between 3 and for is optimal in many dimulations pg293
	       

	 dx=1.e-9;
	 dt=dx/(2.0*cc);
     nmax=500;              % number of time steps - spectral resolution
   	 dt_over_dx=dt/dx;
	
   	 i_size=100;
	 j_size=70;
	


	% ie,je  total gird size including PML
 	 ib=i_size+2*iebc;jb=j_size+2*jebc;
     
	% ib for Hy x size, jb for Hx y size
     	 ie=ib+1;je=jb+1;
    % source position 
        isrc=ib/4;jsrc=jb/2;

	

%...... Initialize Grid Arrays............
%	 
      EZ=zeros(ie,je);
      DEZ=zeros(ie,je);
      HX=zeros(ib,jb);
      HY=zeros(ib,jb);		  
      BHX=zeros(ib,jb);
      BHY=zeros(ib,jb); 
	  PEZ=zeros(ie,je);
      PEZ_old=zeros(ie,je);
      DEZ_old=zeros(ie,je);
   
 	  
%.......General Update Coefficients.......
      axebc(1:iebc)=1.0;
      ayebc(1:jebc)=1.0;
      bxebc(1:iebc)=1.0;
      byebc(1:jebc)=1.0;
	  cxebc(1:iebc)=1.0;
      cyebc(1:jebc)=1.0;
	  axhbc(1:iebc)=1.0;
      ayhbc(1:jebc)=1.0;
      bxhbc(1:iebc)=1.0;
      byhbc(1:jebc)=1.0;
	  cxhbc(1:iebc)=1.0;
      cyhbc(1:jebc)=1.0;

% Different Fourier Transform arrays
	   
       
%.......Initialize Update Coefficients.......

      
      axe(1:ie)=1.0;
	  bxe(1:ie)=1.0;
	  cxe(1:ie)=1.0;
      
	  aye(1:je)=1.0;
	  bye(1:je)=1.0;
	  cye(1:je)=1.0;
	  
      axh(1:ib)=1.0;
	  bxh(1:ib)=1.0;
      cxh(1:ib)=1.0;
	  
      ayh(1:jb)=1.0;
	  byh(1:jb)=1.0;
      cyh(1:jb)=1.0;

      
      epsr(1:ie,1:je)=1.0/epsz;
	  
      muer(1:ib,1:jb)=1.0/muz ;
      

%
%.....................Fill PML Grid.........................


      
      delbc=(iebc)*dx;
	  eta=sqrt(muz/epsz);       % eq 7.7a
      % polynomial grading - eq 7.62 - sigma max
      sigmam=-log(rmax)*(orderbc+1.0)/(2.0*delbc*eta);  
      %delbc --> d
      bcfactor=sigmam/(dx*(delbc^orderbc)*(orderbc+1.0)) ;
      kx=1.0;
	  ky=1.0;
%
%       .....Fill General UPML coefficient.....
%
        for i=1:iebc;
          x1=dx*(i);
          x2=dx*((i)-1.0);
          sigma=bcfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0));
          %what is the purpose of "D" below?
          axebc(i)=(kx-sigma*dt/2.D0/epsz)/(kx+sigma*dt/2.D0/epsz); %Why the difference from 7.91? or 7.55
          ayebc(i)=(ky-sigma*dt/2.D0/epsz)/(ky+sigma*dt/2.D0/epsz);
		  bxebc(i)=(kx+sigma*dt/2.D0/epsz);
		  byebc(i)=(ky+sigma*dt/2.D0/epsz);
		  cxebc(i)=(kx-sigma*dt/2.D0/epsz);
		  cyebc(i)=(ky-sigma*dt/2.D0/epsz);
        end
          x1=dx*0.5;
          sigma=bcfactor*x1^(orderbc+1.0);

		  axhbc(1)=(kx-sigma*dt/2.0/epsz)/(kx+sigma*dt/2.0/epsz);
          ayhbc(1)=(ky-sigma*dt/2.0/epsz)/(ky+sigma*dt/2.0/epsz);
		  bxhbc(1)=(kx+sigma*dt/2.0/epsz);
		  byhbc(1)=(ky+sigma*dt/2.0/epsz);
		  cxhbc(1)=(kx-sigma*dt/2.0/epsz);
		  cyhbc(1)=(ky-sigma*dt/2.0/epsz);
        
          for i=2:iebc;
          x1=dx*(i-0.5);
          x2=dx*(i-1.5);
          sigma=bcfactor*(x1^(orderbc+1.0)-x2^(orderbc+1.0));

		  axhbc(i)=(kx-sigma*dt/2.0/epsz)/(kx+sigma*dt/2.0/epsz);
          ayhbc(i)=(ky-sigma*dt/2.0/epsz)/(ky+sigma*dt/2.0/epsz);
		  bxhbc(i)=(kx+sigma*dt/2.0/epsz);
		  byhbc(i)=(ky+sigma*dt/2.0/epsz);
		  cxhbc(i)=(kx-sigma*dt/2.0/epsz);
		  cyhbc(i)=(ky-sigma*dt/2.0/epsz);

          end

 
%
%.......Fill EZ, HY and HX with PML Media .......
	  for i=1:iebc;
		 i0=iebc+1-i;
	   	 i1=iebc+1-i;
		 i2=ie-iebc+i;
		 i3=ib-iebc+i;
	   	 axe(i0) =axebc(i);
		 axe(i2) =axebc(i);
		 bxe(i0) =bxebc(i);
		 bxe(i2) =bxebc(i);
		 cxe(i0) =cxebc(i);
		 cxe(i2) =cxebc(i);

		 axh(i1) =axhbc(i);
		 axh(i3) =axhbc(i);
		 bxh(i1) =bxhbc(i);
		 bxh(i3) =bxhbc(i);
		 cxh(i1) =cxhbc(i);
		 cxh(i3) =cxhbc(i);
       end

	    for j=1:jebc;
	   	 j0=jebc+1-j;
		 j1=jebc+1-j;
		 j2=je-jebc+j;
		 j3=jb-jebc+j;
	   	 aye(j0) =ayebc(j);
		 aye(j2) =ayebc(j);
		 bye(j0) =byebc(j);
		 bye(j2) =byebc(j);
		 cye(j0) =cyebc(j);
		 cye(j2) =cyebc(j);
	 
		 ayh(j1) =ayhbc(j);
		 ayh(j3) =ayhbc(j);
		 byh(j1) =byhbc(j);
		 byh(j3) =byhbc(j);
		 cyh(j1) =cyhbc(j);
		 cyh(j3) =cyhbc(j);
        end




%.......Initialize the Matrices.......
      
%***********************************************************************
%     Movie initialization
%***********************************************************************

figure(1);

subplot(3,1,1);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(HX');   
shading flat;daspect([1 1 1]);
colorbar;
axis image;
axis off;
title(['Hx at time step = 0']);

subplot(3,1,2);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(HY');  
shading flat;daspect([1 1 1]);
colorbar;
axis image;
axis off;
title(['Ey at time step = 0']);

subplot(3,1,3);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(EZ');  
shading flat;daspect([1 1 1]);
colorbar;
axis image;
axis off;
title(['Hz at time step = 0']);

rect=get(gcf,'Position');
rect(1:2)=[0 0];



 
 
%...............TIME-STEPPING LOOP...........................
%
	nlimit=40;  


 for n=1:nmax;      %Begin time-stepping loop
		 
%................Ez UPDATE...................................

      for j=2:jb;
      for i=2:ib;
	
	  PEZ_old(i,j)=PEZ(i,j);
	 
	  DEZ_old(i,j)=DEZ(i,j);
    	DEZ(i,j)=axe(i)*DEZ(i,j)			  ...
              +(HY(i,j)-HY(i-1,j)-HX(i,j)+HX(i,j-1))*dt_over_dx/bxe(i);
		        
		EZ(i,j)=aye(j)*EZ(i,j)			     ...
             +(DEZ(i,j)-DEZ_old(i,j))*epsr(i,j)/bye(j);
      end
      end

%     Soft source
      if (n <= nlimit) 

	   esource  =(10.0-15.0*cos(pi*0.05*n)+6.0*cos(tpi*0.05*n)	...
		           -cos(3.0*pi*0.05*n))/32.0;
	  else
	  esource=0.0;
      end
    
	  EZ(isrc,jsrc)= esource;
	 
     
%.............HX and HY UPDATE...................................

      for j=1:jb;
      for i=1:ib;
		BHXS=BHX(i,j);
        BHX(i,j)=ayh(j)*BHX(i,j)-				...
                  (EZ(i,j+1)-EZ(i,j))*dt_over_dx/byh(j);
		HX(i,j)=HX(i,j)+				  ...
		            (bxe(i)*BHX(i,j)-cxe(i)*BHXS)*muer(i,j)		;	   
      end
      end
     for j=1:jb
      for i=1:ib
		BHYS=BHY(i,j);
        BHY(i,j)=BHY(i,j)+					 ...
                  (EZ(i+1,j)-EZ(i,j))*dt_over_dx;
		HY(i,j)=axh(i)*HY(i,j)+        ...
		          (bye(j)*BHY(i,j)-cye(j)*BHYS)*muer(i,j)/bxh(i);
      end
     end

 

%.......Next Time Step Iteration.......

%***********************************************************************
%     Visualize fields
%***********************************************************************

if mod(n,10)==0;

timestep=int2str(n);

   
subplot(3,1,1);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(HX');  
shading flat;daspect([1 1 1]);
colorbar;

title(['Hx at time step = ',timestep]);

subplot(3,1,2);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(HY'); 
shading flat;daspect([1 1 1]);
colorbar;

title(['Hy at time step = ',timestep]);

subplot(3,1,3);
set(gca,'units','points'); %without this setting, the axes will shrink in 2009a version
pcolor(EZ'); 
shading flat;daspect([1 1 1]);
colorbar;

title(['Ez at time step = ',timestep]);
getframe(gcf,rect);
    

end;


  
         

 

 
 end