function PBGBand(ea,eb,R,PCType,Keach,TEorTM)

%function PBGBand(ea,eb,R,PCType,Keach,TEorTM)
%--------------------------------------------------------------
%| This is a program to calculate the Photonic Bands of two   |
%| dimension Photonic Crystal with circular inclusions.       |
%| It calculates both TE and TM modes (E and H polarization)  |
%--------------------------------------------------------------
%Parameters:
%ea: The dielectric constant of the circular inclusions.
%eb: The dielectric constant of background.
%R: The radius of dielectric columns
%PCType  =1: Square lattice
%        =2: Triangular lattice
%        =3: Honeycomb
%Keach: The number of k vectors in each wave vector branch.
%TEorTM: =0: TE modes
%        =1: TM modes
%--------------------------------------------------------------

tic;

disp('--------------------------------------------------')
if (TEorTM==0)
    disp('Plane wave expansion method for PC bands: TE modes');
else;
    disp('Plane wave expansion method for PC bands: TM modes');
end    
disp('--------------------------------------------------')
if (PCType==1)
    disp('Square lattice');
end
if (PCType==2)
    disp('Triangular lattice');
end
if (PCType==3)
    disp('Honeycomb lattice');
end

%Control parameters
Ktype=3; % The number of band parts, such as X->M, T->X, ...
NumberK=Ktype*Keach; %The total number of K vector;

NEIG=20; %NEIG: The number cut of the Eigvalue.

%Initial parameters
a=1;   %Lattice constance.

a1=a*[1,0];
if PCType==1
   a2=a*[0,1];
end
if PCType==2
   a2=a*[0.5,sqrt(3)/2];
end
%a1,a2 are the basic vectors of lacctice cell.

ac=abs(a1(1)*a2(2)-a1(2)*a2(1)); 
%ac: Area of lattice cell.
b1=2*pi/ac*[a2(2),-a2(1)];
b2=2*pi/ac*[-a1(2),a1(1)];
%b1, b2 are vectors in reciprocal space.

f=pi*R*R/ac;
%f: The filling fraction, i.e. the fraction of
%   the total volume occupied by the rods.

MaxDimForG=10;  % The max Potive Number of the reciprocal lattice, G
DimForG=2*MaxDimForG+1; 
NPW=DimForG*DimForG; %NPW: The number of Plane Waves

%
%            ^ 
%            | Y
%      O  O  O  O  O  -->(MaxDimForG,MaxDimForG)for this point!
%      O  O  O  O  O 
%    --O--O--O--O--O--> X
%      O  O  O  O  O
%      O  O  O  O  O 
%            |   
%            |
%initial the G matrix
disp('--------------------------------------------------')
disp('Dielectric constant FT--- BEGIN')

gtemp=-MaxDimForG:MaxDimForG;
gtemp1=repmat(gtemp,DimForG,1);
Gx=b1(1)*gtemp1+b2(1)*gtemp1';
Gy=b1(2)*gtemp1+b2(2)*gtemp1';
Gx=Gx(:)';
Gy=Gy(:)';

disp(strcat('The number of plane waves is--',num2str(NPW)));

Gx_m=repmat(Gx,NPW,1);
Gx_n=Gx_m';
Gy_m=repmat(Gy,NPW,1);
Gy_n=Gy_m'; 

%calculate the Matrix coefficience.
ek0=f/ea+(1-f)/eb; 
ekc=(1/ea-1/eb)*f*2;  

%Calculate the ek matrix, the coefficence of Fourier transform of ek.
GR_mat=sqrt((Gx_m-Gx_n).*(Gx_m-Gx_n)+(Gy_m-Gy_n).*(Gy_m-Gy_n))*R;
if PCType==1|PCType==2
   %eliminate the division on zero in the calculatation of ek
   na=find(GR_mat==0);
   GR_mat(na)=1;   
   ek_mat=ekc*besselj(1,GR_mat)./GR_mat; 
   ek_mat(na)=ek0;
end
if PCType==3
   %eliminate the division on zero in the calculatation of ek
   na=find(GR_mat==0);
   GR_mat(na)=1;   
   ek_mat=cos((Gx_m-Gx_n).*a/2+(Gy_m-Gy_n).*a*sqrt(3)/6).*ekc.*besselj(1,GR_mat)./GR_mat; 
   ek_mat(na)=ek0;
end
%toc
%tic
%Calculated points:
Point=zeros(Ktype+1,2);
if PCType==1 %Square lattice, or rectangular lattice
   Point(1,:)=[0,0]; %Gama Point
   Point(2,:)=[b1(1)/2,0]; %X Point
   Point(3,:)=[(b1(1)+b2(1))/2,(b1(2)+b2(2))/2]; %M point
   Point(4,:)=[0,0]; %Gama Point
end
if PCType==2|PCType==3 %Triangular lattice and Honeycomb lattice.
   Point(1,:)=[0,0]; %Gama Point
   Point(2,:)=[b2(2)*sqrt(3)/6,b2(2)/2]; %K Point
   Point(3,:)=[0,b2(2)/2]; % M point
   Point(4,:)=[0,0]; %Gama Point
end

%These three are for the K vectors, for the different case.
K1=[];
K2=[];
for ktnum=1:Ktype
   K1temp=linspace(Point(ktnum,1),Point(ktnum+1,1),Keach+1);
   K2temp=linspace(Point(ktnum,2),Point(ktnum+1,2),Keach+1);
   K1=[K1,K1temp(1:Keach)];
   K2=[K2,K2temp(1:Keach)];
end
disp('Dielectric constant FT--- END')
disp('--------------------------------------------------')

disp('Eigen value calculations--- BEGIN')

eigval=[]; %Initial the eigvalue matrix.

for knum=1:NumberK
   disp(strcat('---K vector No.',num2str(knum),'---',num2str(NumberK))) 
   
   kx=K1(knum);
   ky=K2(knum);
   %tic
   
   %Now begin to calculate the matrix H:      
   if (TEorTM==0)
       %TE part
       KGmn_mat=(kx-Gx_m).*(kx-Gx_n)+(ky-Gy_m).*(ky-Gy_n);
       H=KGmn_mat.*ek_mat;
   else
       %TM part
       %Place your codes below for Task 1.
   
   end       
       
   %Find the eigenvalues
   eigvalue=sort(eig(H));
   eigval=[eigval,eigvalue(1:NEIG)];   
end

eigval=[eigval,eigval(:,1)];  
eigval=sqrt(eigval)*a*0.5/pi;
eigval=real(eigval); %Complete All the things

disp('Eigen value calculation0s--- END')
disp('--------------------------------------------------')

%Plot the figures
%x=linspace(0,10,NumberK+1);
for m=1:Ktype
   D(m)=sqrt((Point(m+1,1)-Point(m,1))^2+(Point(m+1,2)-Point(m,2))^2);
   xtemp(m,:)=linspace(0,D(m),Keach+1);
end
x=xtemp(1,1:Keach);
Dtotal=0;
for m=2:Ktype
   Dtotal=Dtotal+D(m-1);
   x=[x,xtemp(m,1:Keach)+Dtotal];
end
x=[x,xtemp(Ktype,Keach+1)+Dtotal];
x=x/max(x);

MaxB=0.8;
x1=x(Keach+1);
x2=x(Keach*2+1);

figure;
clf;

h=plot(x,eigval,'b-',[x1 x1],[0 MaxB],'k:',[x2 x2],[0 MaxB],'k:');
set(h,'LineWidth',2.0);
if (TEorTM==0)
    legend('TE modes',4);
else
    legend('TM modes',4);
end    

axis([0 1 0 MaxB]);
h=ylabel('Normalized frequency (a/\lambda)');
set(h,'FontSize',14);
if (PCType==1)
   titletext=strcat('Square Lattice (ea=',num2str(ea),', eb=',num2str(eb),', R=',num2str(R),')');
   text(x(1)-0.02,-0.03, '\Gamma','FontSize',14)
   text(x1-0.02,-0.03, 'X','FontSize',14)
   text(x2-0.02,-0.03, 'M','FontSize',14)
   text(x(Keach*Ktype+1)-0.02,-0.03, '\Gamma','FontSize',14)
   
end
if (PCType==2)
   titletext=strcat('Triangular Lattice (ea=',num2str(ea),', eb=',num2str(eb),', R=',num2str(R),')');
   text(x(1)-0.02,-0.03, '\Gamma','FontSize',14)
   text(x1-0.02,-0.03, 'K','FontSize',14)
   text(x2-0.02,-0.03, 'M','FontSize',14)
   text(x(Keach*Ktype+1)-0.02,-0.03, '\Gamma','FontSize',14)
end
if (PCType==3)
   titletext=strcat('Honeycomb Lattice (ea=',num2str(ea),', eb=',num2str(eb),', R=',num2str(R),')');
   text(x(1)-0.02,-0.03, '\Gamma','FontSize',14)
   text(x1-0.02,-0.03, 'K','FontSize',14)
   text(x2-0.02,-0.03, 'M','FontSize',14)
   text(x(Keach*Ktype+1)-0.02,-0.03, '\Gamma','FontSize',14)
end
h=title(titletext);
set(h,'FontSize',14);
set(gca,'xtick',[]);

%Save the data
if (TEorTM==0)
    save datate.mat x ea eb R PCType MaxB Keach eigval;
else
    save datatm.mat x ea eb R PCType MaxB Keach eigval;
end    

toc;



