
%% code for the instanteous streamlines at time t=0 
alpha=0.005;
title('The streamlines of the velocity field, Okubo Weiss region and Kam curves for t=0 and alpha=0.005')
xspan = -10:1:10; yspan = -10:1:10;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;

u1=1.5*y+alpha*(x.^4-6*(x.^2).*(y.^2)+y.^4);
v=0.5*x+alpha*(-4*y.*x.^3+4*x.*y.^3);
streamline(x,y,u1,v,startx,starty,[0.04,500])
hold on

u2=-1.5*y-alpha*(x.^4-6*(x.^2).*(y.^2)+y.^4);
v2=-0.5*x-alpha*(-4*y.*x.^3+4*x.*y.^3);
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.04,500])
hold on 

xlabel('x')
ylabel('y')
% code for okubo Weiss region 

r = -10:0.03:10;  %Okubo Weiss elliptic criterion for t=0
[x,y] = meshgrid(r,r);

condition=Q(x,y)<0;
condition = double(condition); % convert to double for plotting
condition(condition == 0) = NaN; %set the 0s to NaN so they are not plotted
surf(x,y,condition,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','y')
view(0,90)
hold on 

%%
var1=60;
var2=60;
x=linspace(-15,15,var1);
y=linspace(-15,15,var2);
sizepoints=var1*var2;
[X,Y]=ndgrid(x,y);
tspan=0:pi/2:300*pi;

[t,x]=ode45(@integrator,tspan,[X(:);Y(:)]);
sizexvector=size(x)
x=reshape(x,sizexvector(1),sizexvector(2));



axisx = x(1:end,1:sizexvector(2)/2)
axisy = x(1:end,(sizexvector(2)/2)+1:sizexvector(2))
%%
plot(axisx(400:end,:),axisy(400:end,:),'k.','MarkerSize',0.1')
xlim([-10,10])
ylim([-10,10])
hold on 

%%

function dxdt=integrator(t,x)
    N=numel(x)/2;

    alpha=0.005;
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1 = sin(4*t).*xf+(cos(4*t)+1/2).*yf+alpha*(xf.^4-6*(xf.^2).*(yf.^2)+yf.^4);
    dx2 = (cos(4*t)-1/2).*xf-sin(4*t).*yf+alpha*(-4*yf.*xf.^3+4*xf.*yf.^3);
    
    dx1(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    dx2(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;
    

 end
 function Q=Q(x,y)
alpha=0.005;
a=4*alpha*(x.^3)-12*alpha*x.*(y.^2);
b=0.5-12*alpha*(x.^2).*y+4*alpha*y.^3;
c=1.5-12*alpha*y.*x.^2+4*alpha*y.^3;
Q=a.^2+b.*c
 end
 
 
