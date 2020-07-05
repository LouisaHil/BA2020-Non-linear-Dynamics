
%% code for the instanteous streamlines at time t=0 
alpha=-0.1;
title('The streamlines of the velocity field, Okubo Weiss region,Poincare Map for t=0 ,alpha=-0.1, C=-4,w=-4')
xspan = -5:0.5:5; yspan = -5:0.5:5;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;

u1=3*y+alpha*(x.^4-6*(x.^2).*(y.^2)+y.^4);
v=-1*x+alpha*(-4*y.*x.^3+4*x.*y.^3);
streamline(x,y,u1,v,startx,starty,[0.09,100])
hold on

u2=-3*y-alpha*(x.^4-6*(x.^2).*(y.^2)+y.^4);
v2=1*x-alpha*(-4*y.*x.^3+4*x.*y.^3);
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.09,100])
hold on 

xlabel('x')
ylabel('y')
% code for okubo Weiss region 

r = -15:0.05:15;  %Okubo Weiss elliptic criterion for t=0
[x,y] = meshgrid(r,r);

condition=Q(x,y)<0;
condition = double(condition); % convert to double for plotting
condition(condition == 0) = NaN; %set the 0s to NaN so they are not plotted
surf(x,y,condition,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','y')
view(0,90)
hold on 


function dxdt=integrator(t,x)
    N=numel(x)/2;

    alpha=-0.1;
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1 = sin(4*t).*xf+(cos(4*t)+2).*yf+alpha*(xf.^4-6*(xf.^2).*(yf.^2)+yf.^4);
    dx2 = (cos(4*t)-2).*xf-sin(4*t).*yf+alpha*(-4*yf.*xf.^3+4*xf.*yf.^3);
    
    dx1(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    dx2(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;
    

 end
 function Q=Q(x,y)
alpha=-0.1;
a=4*alpha*(x.^3)-12*alpha*x.*(y.^2);
b=-1-12*alpha*(x.^2).*y+4*alpha*y.^3;
c=3-12*alpha*y.*x.^2+4*alpha*y.^3;
Q=a.^2+b.*c
 end
 

