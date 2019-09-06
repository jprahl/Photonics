
% True FDTD band structure of 2-D periodic material
% TM mode
clear all;

%parameters
cc = 3e8;
a = 1e-6;           % lattice constant in meters

jb=50; % # of mesh elements at space discretization stage
ib=50;
dy = a/jb;   % size of mesh cell
dx = a/ib;
dt = 0.5*dx/cc;    %1/cc/sqrt(1/dy^2+1/dx^2);      %Courant condition

f0=cc/a;            %central freq of Gaussian pulse
omega0=2*pi*f0;
lambda0 = a;
period_T=1/f0;
ndelay = 3*period_T/dt;
Nlambda = lambda0/dx;
gausswidth = ndelay/Nlambda;%100*period_T;
gausswidth2=(0.5*period_T/dt)^2;

nfreq=2000;
freqspan=2*omega0/nfreq;
omega(1:nfreq)=freqspan:freqspan:2*omega0;
rd(1:nfreq,1:2)=0;
snorm(1:nfreq,1:2)=0;
EzDFT(1:nfreq,1:2)=0;

mib = ib*4;           %grid to define material
mjb = jb*4;
nx = 11;  ny = 11;        %Resolution of grid (NEEDS TO BE ODD???)

dxm = a/mib;
dym = a/mjb;

isrc1=round(ib/5);
jsrc1=round(jb/5);

iobs1=round(ib/10);
jobs1=round(jb/10);
iobs2=round(ib/5*2);
jobs2=round(jb/5*2);
iobs3=round(ib/5*3);
jobs3=round(jb/5*3);
iobs4=round(ib/5*4);
jobs4=round(jb/5*4);
%Where is N_slant defined? - see eq 13.14/15 from Tavlove
eps1 = 1;           %Material parameters, switching between will create a 
eps2 = 9;           %hole/rod structure

%l1 = 0.5;           % thickness of layers defined as fractions of period
%l2 = 1-l1;

%Circle
r = 0.3*a;          %Radius of holes/rods of lattice
icenter=ib/2;
jcenter=jb/2;

%Define triangle by vertices
vertex1=[icenter,jcenter+r/dy];
vertex2=[icenter+((r/dy)*sind(60)),jcenter-((r/dy)*cosd(60))+4];
vertex3=[icenter-((r/dy)*sind(60)),jcenter-((r/dy)*cosd(60))-4];
%NOTE - ASSUMES EQUILATERAL TRIANGLE- ALTER ANGLES FOR COMPLEX SHAPES

slope1=(vertex1(2)-vertex3(2))/(vertex1(1)-vertex3(1));
slope2=(vertex2(2)-vertex1(2))/(vertex2(1)-vertex1(1));
slope3=(vertex2(2)-vertex3(2))/(vertex2(1)-vertex3(1));

j1=vertex1(2)-slope1*vertex1(1);
j2=vertex2(2)-slope2*vertex2(1);
j3=vertex3(2)-slope3*vertex3(1);

cond1=(j<=(i*slope1)+j1);
cond2=(j<=(i*slope2)+j2);
cond3=(j>=(i*slope3)+j3);

material(1:ib,1:jb)=0;

%    for i=1:ib
%        for j=1:jb
%            for mx=1:nx
%                for my=1:ny
%               ii=i+(mx-(nx+1)/2)/nx;
%               jj=j+(my-(ny+1)/2)/ny;
%                if (jj<=(ii*slope1)+j1) && (jj<=(ii*slope2)+j2) && (jj>=(ii*slope3)+j3)
% 
%               material(i,j)=material(i,j)+eps2;
%                else
%               material(i,j)=material(i,j)+eps1;
%                end  % end if
%                end   % end my
%            end  % end mx
%            material(i,j)=1/(material(i,j)/nx/ny);
%        end  % end i
%    end     % end j

   for i=1:ib
       for j=1:jb
           for mx=1:nx
               for my=1:ny
              ii=i+(mx-(nx+1)/2)/nx;
              jj=j+(my-(ny+1)/2)/ny;
               if ((ii-icenter)^2+(jj-jcenter)^2 <= (r/dx)^2)
              material(i,j)=material(i,j)+eps2;
               else
              material(i,j)=material(i,j)+eps1;
               end  % end if
               end   % end my
           end  % end mx
           material(i,j)=1/(material(i,j)/nx/ny);
       end  % end i
   end     % end j

%    for i=1:ib
%        for j=1:jb
%           
%           material(i,j)= real(eps1+(eps2-eps1)*...
%               (1-sqrt((((i-icenter)^2+(j-jcenter)^2)-(r/dx))/nx)));
%               
%            material(i,j)=1/(material(i,j)/nx/ny);
%        end  % end i
%    end     % end j

figure(2);pcolor(1./material);shading flat;colorbar;daspect([1 1 1]); 
%pause

t_max = 2^14;       %number of Field steps
nmax=t_max;                    %for Fourier image using FFT, # should be multiple of 2

numK = 0;  %counts wave vecctors for PBG formation
%numKy = 0;

%wEigen=[zeros,0];

dt1 =  cc*dt/dy; %coefficents at field components
%dt2 = dt/dy/muz;  %NOTE - IF dy does not = dx, MODIFICATION NEEDED!


%ax = axes;
hold on;



%Begin loop
%Define high symmetry points for square lattice Gamma, X, M
%For now let a1=a2(lattice) ==> g1=g2 (reciprocal lattice) 
precis=10;
Kx(1:precis+1)=0:pi/a/precis:pi/a;
Ky(1:precis+1)=zeros(1,precis+1);
Kx(precis+2:precis+precis+1)=pi/a;
Ky(precis+2:precis+precis+1)=pi/a/precis:pi/a/precis:pi/a;
Kx(precis+2+precis:precis+precis+1+precis)=...
    pi/a-pi/a/precis:-pi/a/precis:0;      %step modified for side
Ky(precis+2+precis:precis+precis+1+precis)=...
    pi/a-pi/a/precis:-pi/a/precis:0;

%K = [Kx;Ky]';
simulation_data=fopen('band_diagram_2D.txt','w');
for k=1:length(Kx)
   
    numK=numK+1;
                
    Ez=zeros(ib,jb);      %Initialize field arrays
    Hx=zeros(ib,jb);      %different in other code - ie and ib
    Hy=zeros(ib,jb);      %due to uncertain number of layers in current
                              %Should this be changed back to seperate
                              %discretization points for E and H?
    
    %calculate coefficient which defines phase rotation
    %when the field passes through the boundary with Bloch B.C.
    %rotate_x = exp(1i*Kx*a);
    %rotate_y = exp(1i*Ky*a);
    
    clear Field;     %receives values of one of the field components
                    %at some point in space at various Fields
    Field(1:nmax)=0;                
    for nt=1:nmax
       % dt:((t_max)*dt)-dt  %calculates field variation in Field
       esource=exp(-(nt-ndelay)^2/gausswidth2)*sin(omega0*(nt-ndelay)*dt);

%View field intensity during simulation
%     if(mod(nt,100)==0)
%         nt
%         esource
%         Ez(iobs1,jobs1)
%     end
     %initial condition (excitation) cmplx hrmnic fnc in a unit ccell
         
       for i=1:ib 
       Hx(i,1)=Hx(i,1)-dt1*(Ez(i,1)-Ez(i,jb)*exp(-1i*Ky(k)*a)); %def boundary cconditions
      %Hx(i,1)=Hx(i,1)-dt1*(Ez(i,1)-Ez(i,jb));
       end                                  %NOTE - BC ASSUME MU = 0
       
       for j=1:jb                            % WHAT ABOUT THE HALF STEP?
       Hy(1,j)=Hy(1,j)+dt1*(Ez(1,j)-Ez(ib,j)*exp(-1i*Kx(k)*a));
       %Hy(1,j)=Hy(1,j)+dt1*(Ez(1,j)-Ez(ib,j));
       end;
        
     %calculate magnetic component
        for j=2:jb
            for i=1:ib  
            
                Hx(i,j)=Hx(i,j)-dt1*(Ez(i,j)-Ez(i,j-1));
           
            end
        end
        
        for j=1:jb
            for i=2:ib  
            
                Hy(i,j)=Hy(i,j)+dt1*(Ez(i,j)-Ez(i-1,j));
           
            end
        end
     %calculate electric component 
     %note mu not used = 1
     
      for i=1:ib-1 
        Ez(i,jb)=Ez(i,jb)+...
            dt1*material(i,jb)*(Hy(i+1,jb)-Hy(i,jb))+ ...
            dt1*material(i,jb)*(Hx(i,jb)-Hx(i,1)*exp(1i*Ky(k)*a)); 
            %dt1*material(i,j)*(Hx(i,jb)-Hx(i,1));
         %def boundary conditions
      end
       for j=1:jb-1
        Ez(ib,j)=Ez(ib,j)+...
            dt1*material(ib,j)*(-Hx(ib,j+1)+Hx(ib,j)...
            +dt1*material(ib,j)*(Hy(1,j)*exp(1i*Kx(k)*a)...
            -Hy(ib,j)));
           %+dt1*material(i,j)*(Hy(1,j) ... 
       end;
        Ez(ib,jb)=Ez(ib,jb)+...
            dt1*material(ib,jb)*(Hx(ib,jb)-...
            Hx(ib,1)*exp(1i*Ky(k)*a))...
            +dt1*material(ib,jb)*(Hy(1,jb)*exp(1i*Kx(k)*a)...
           -Hy(ib,jb));
           
           %Hx(ib,1)) ...  
           % 
           %  +dt1*material(ib,jb)*(Hy(1,jb)
%        Ez(ib,jb)=Ez(ib,jb)-...
%             dt1/material(ib,jb)*(Hx(ib,jb)-...
%             Hx(ib,1)*exp(-1i*Kx(k)*a)*exp(-1i*Ky(k)*a))...
%             +dt1/material(ib,jb)*(Hy(1,jb)*exp(-1i*Kx(k)*a)*exp(-1i*Ky(k)*a)-Hy(ib,jb));
     for j=1:jb-1
         for i=1:ib-1
         
             Ez(i,j)=Ez(i,j)+dt1*material(i,j)*(Hy(i+1,j)-Hy(i,j)...
                 +Hx(i,j)-Hx(i,j+1));
          
         end
     end
       
       Ez(isrc1,jsrc1)=Ez(isrc1,jsrc1)+esource; 
       
     %Intensity of the Electric field at a point of a unit cell a 
     %suitable distance from excitation point
     %Read fields from several points and average them
     
        Field(nt)=Ez(iobs1,jobs1);
        Field02(nt)=(Ez(iobs2,jobs2));
        Field03(nt)=(Ez(iobs3,jobs3));
        Field04(nt)=(Ez(iobs4,jobs4));
        Field2(nt)=esource;
        
        
         for iff=1:nfreq
             rd(iff,1)=cos(omega(iff)*nt*dt);
             rd(iff,2)=sin(omega(iff)*nt*dt);
         snorm(iff,1)=snorm(iff,1)+esource*rd(iff,1);
         snorm(iff,2)=snorm(iff,2)+esource*rd(iff,2);
         EzDFT(iff,1)=EzDFT(iff,1)+Ez(iobs1,jobs1)*rd(iff,1);
         EzDFT(iff,2)=EzDFT(iff,2)+Ez(iobs1,jobs1)*rd(iff,2);
         end
        if(mod(nt/dt,100) ==0)
            figure(3);pcolor(Ez);shading flat;colorbar;daspect([1 1 1]);
          %  pause
         end
    end
    
    %Fourier image of Field-dependent pulse response of structure
    fourier=fft(Field,nmax)/nmax;
    fourier2=fft(Field2,nmax)/nmax;
    f=1/(2*dt)*linspace(0,1,nmax/2+1);
    figure(3);plot(2*pi*f,2*abs(fourier2(1:nmax/2+1)));
    figure(4);plot(2*pi*f,2*abs(fourier(1:nmax/2+1)));
    figure(5);plot(2*pi*f,abs(fourier(1:nmax/2+1))./abs(fourier2(1:nmax/2+1)));
    snorm(1:nfreq,1)=sqrt(snorm(1:nfreq,1).^2+snorm(1:nfreq,2).^2);
    EzDFT(1:nfreq,1)=sqrt(EzDFT(1:nfreq,1).^2+EzDFT(1:nfreq,2).^2)./snorm(1:nfreq,1);
    figure(6);plot(omega,snorm(1:nfreq,1));
    figure(7);plot(omega,log(EzDFT(1:nfreq,1)));
    figure(8);plot(omega,(EzDFT(1:nfreq,1)));
    %pause;
    %fourier=fourier./abs(fourier2)
    %Absolute frequency values
    
    %Discrete Fourier transform

    
    wcount=1;       %counter of eigen-freq for current K 
    
    %Analyze spectrum of pulse response, 1st element seperately
   % if (fourier(1)/fourier(2)>1.0001)
     %  wEigen(numK, wcount)=f(1);
     %  wcount=wcount+1;
   % end
   
    %other elements, if the current element is larger it is considered as
    %an eigen-state and corresponding frequenccy is recorded into the array
    %of eigen-states of the Phc
    wcountmax=length(fourier)/2/20;
    for u=2:wcountmax
   % while(wcount <=8 )
        if(abs(fourier(u))/abs(fourier(u-1))>1.0001)&&...
                (abs(fourier(u))/abs(fourier(u+1))>1.0001)&&...
                ((abs(fourier(u))/abs(fourier2(u)))>0.05)
            wEigen(numK, wcount)=f(u);
            numK
            wcount
            f(u)
            wcount=wcount+1;
           
     %   end
    end
    end
    
    %Display the eigen-states
    %If 16 or more eigen-states are found only 16 are shown
    figure(1);
    hold on
    %if(wcount-1>=16)
        plot(ones(1,8)*(k-1)*a/2/pi,...
            abs(wEigen(numK,1:8)*(a/cc)),'o','LineWidth',3);
    %else
    %    plot(ones(1,wcount-1)*sqrt(Kx(k)^2+Ky(k)^2)*a/2/pi,...
    %        abs(wEigen(numK,1:wcount-1)*(a/2/pi/cc)),'o','LineWidth',3);
    %end
    
    %Axes
    %ylim([0 0.2]);
    %xlim([0 3.15E6]);
    %ylim([0 2]);
    %xlim([0 2]);
    xlabel('K,m^{-1}');
    ylabel('Frequenccy, Hz');
    %set(axes,'XGrid','on');
   % set(axes,'YGrid','on');
    drawnow;
 
 for mm=1:8
 fprintf(simulation_data, '%13.9f %13.9f\n',sqrt(Kx(k)^2+Ky(k)^2)*a/2/pi, abs(wEigen(numK,mm)*(a/2/pi/cc)));
 end
 
 toc
%figure(4);plot(f(1:2^14)*(a/2/pi/cc),abs(fourier2))
% tt=0:dt:((t_max)*dt)-dt;
%figure(5);plot(tt,exp(-(tt-0.5e-12).^2/gausswidth2).*sin(5e15*tt));

   
end             %Kx

fclose(simulation_data);  
    