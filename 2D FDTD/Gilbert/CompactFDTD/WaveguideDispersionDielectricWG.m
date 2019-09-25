clear all;
%k=0:0.05:3;
epsz=1;muz=1;cc=1;
epsr=9;
w0=0.05:0.1:20;
k=w0;
figure(4);
plot(w0,k);
k=w0.*sqrt(epsr);
figure(4); hold on;
plot(w0,k,'k');
%pause;
wgwidth=0.5
 	wgwidth2=wgwidth/2.
		
	epsnwg=epsz*epsr
	%F(x)=Tan(x)-Sqrt(a-x**2)/x
	%G(x)=1/(Cos(x)**2)-(-x**2/Sqrt(a-x**2)-sqrt(a-x**2))/x**2
    ncount=0;
    for ncount=1:length(w0)
        w=w0(ncount)
	a0=(w*wgwidth2).^2*muz*(epsr-epsz)  
    if (w < 0.4)
     xx=(w*wgwidth2)^2*muz*epsz-(w/cc*wgwidth2)^0.5;   
    %    xx=w/cc*wgwidth2
    %elseif (w>=5 || w<=15)
    else   
	xx=(w*wgwidth2)^2*muz*epsr-(0.5*(w/cc*wgwidth2+sqrt(epsr)*w/cc*wgwidth2))^0.5;
    %xx=(w*wgwidth2)^2*muz*epsr-(sqrt(epsr)*w/cc*wgwidth2)^0.5;
    end
    %else
    % xx=sqrt(epsr)*w/cc*wgwidth2*0.7
    %end
        %*wgwidth/wavelength(m)/rnf/2.
	%Fx=epsz/epsr*xx*tan(xx)-sqrt(a0-xx^2)

%xx2=fsolve(@(x)epsz/epsr.*x.*tan(x)-sqrt(a0-x.^2)./x,xx);
xx2=lsqnonlin(@(x)epsz/epsr.*x.*tan(x)-sqrt(a0-x.^2)./x,xx,0*pi/2,1*pi/2);
%	for n=1:400
%		Fx=epsz/epsr*xx*tan(xx)-sqrt(a0-xx.^2)./xx
%        Gx=epsz/epsr*(tan(xx)+xx./(cos(xx).^2))+(xx./sqrt(a0-xx.^2))
		
%		xx=xx-Fx./Gx
        
  %  end
	
%	kk(ncount)=xx
%	Kair(ncount)=(a0/(wgwidth2^2)-kk(ncount)^2)^0.5
%    Kz(ncount)=1/wgwidth2*((w*wgwidth2)^2*muz*epsr-kk(ncount)^2)^0.5
    
    kk2(ncount)=xx2
	Kair2(ncount)=(a0/(wgwidth2^2)-kk2(ncount)^2)^0.5
    Kz2(ncount)=1/wgwidth2*((w*wgwidth2)^2*muz*epsr-kk2(ncount)^2)^0.5
	%pause
    end
	
figure(4);hold on;
%plot(w0,real(Kz),'r.');
plot(w0,real(Kz2),'b.');

