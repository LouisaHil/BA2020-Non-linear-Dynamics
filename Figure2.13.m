
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
%%
C=-4;
w=-4;
A=[0 1+0.5*(C-w);1-0.5*(C-w) 0];
[eiv,lamda]=eig(A)
%technically we would need to transform it back to the original frame x so
%do the computation x=M*eiv to find the correct eigenavlues. 
%but here we are only looking at time t=pi/3 wich is the period of the flow, 
%Therefore M is the Unity matrix and we observe no change. 


eivstable= eiv(:,2)
eivunstable=eiv(:,1)


xaxisS=abs(eivstable(1,1));
yaxisS=abs(eivstable(2,1));

xaxisU=abs(eivunstable(1,1));
yaxisU=abs(eivunstable(2,1));

x0=-2:0.002:2
k=-2:1:2
%k=0
Nk=length(k)
Nx=length(x0)
coeffS=yaxisS./xaxisS;
coeffU=yaxisU./xaxisU;
% for i=1:Nk
%     
% y0f=coeffS.*x0+k(i)
% y0b=-coeffU.*x0+k(i)
% plot(x0,y0f,'r.')
% hold on 
% plot(x0,y0b,'k.')
% hold on 
% end 
% xlim([-20 20])
% ylim([-20 20])
%%

Nk=length(k)
Nx=length(x0)
tspan=[0 pi/2] %forward in time
tspanb=[pi/2 0]; %backward in time 
for i=1:Nk 
    
y0f=coeffS.*x0+k(i);
y0b=-coeffU.*x0+k(i);

[t,x]=ode45(@integrator,tspan,[x0(:);y0f(:)]);
[tb,xb]=ode45(@integrator,tspanb,[x0(:);y0b(:)]);


sizexvector=size(x);
sizexvectorb=size(xb);
x=reshape(x,sizexvector(1),sizexvector(2));
xb=reshape(xb,sizexvectorb(1),sizexvectorb(2));

axisx = x(sizexvector(1),1:sizexvector(2)/2);
axisy = x(sizexvector(1),(sizexvector(2)/2)+1:sizexvector(2));
axisxb = xb(sizexvectorb(1),1:sizexvectorb(2)/2);
axisyb = xb(sizexvectorb(1),(sizexvectorb(2)/2)+1:sizexvectorb(2));

plot(axisx,axisy,'k.')
hold on 
plot(axisxb,axisyb,'r.')
hold on 
xlim([-5 5])
ylim([-5 5])
hold on 
i
end
%%


%%


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
 

