var1=30;
var2=30;
x=linspace(-15,15,var1);
y=linspace(-15,15,var2);
sizepoints=var1*var2;
[X,Y]=ndgrid(x,y);
tspan=0:pi/2:200*pi;

[t,x]=ode45(@integrator,tspan,[X(:);Y(:)]);
sizexvector=size(x);
x=reshape(x,sizexvector(1),sizexvector(2));



axisx = x(1:end,1:sizexvector(2)/2)
axisy = x(1:end,(sizexvector(2)/2)+1:sizexvector(2))
%
plot(axisx(1:end,:),axisy(1:end,:),'k.','Markersize',1.5)
xlim([-15,15])
ylim([-15 15])
title('KAM Curves')
xlabel('x')
ylabel('y')

hold on 

%%

function dxdt=integrator(t,x)
    N=numel(x)/2;

    alpha=-0.015;
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1 = sin(4*t).*xf+(cos(4*t)+0.5).*yf+alpha*(xf.^2-yf.^2);
    dx2 = (cos(4*t)-0.5).*xf-sin(4*t).*yf+alpha*(-2*xf.*yf);
    
    dx1(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    dx2(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;
    

 end