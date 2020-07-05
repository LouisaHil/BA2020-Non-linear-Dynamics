%% streamlines
alpha=0.005;
title('The streamlines of the velocity field for t=0, Okubo Weiss region and Poincare Map')
xspan = -20:3:20; yspan = -20:3:20;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;

u1=(7/2*y)+alpha*(x.*(x.^2-3*y.^2));
v=-3/2*x+alpha*(-y.*(3*x.^2-y.^2));
streamline(x,y,u1,v,startx,starty,[0.006,800])
hold on

u2=(-7/2)*y-alpha*(x.*(x.^2-3*y.^2));
v2=3/2*x-alpha*(-y.*(3*x.^2-y.^2));
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.006,800])
hold on 

xlabel('x')
ylabel('y')
hold on 


r = -15:0.05:15;  %Okubo Weiss elliptic criterion for t=0
[x,y] = meshgrid(r,r);

condition=Q(x,y)<0;
condition = double(condition); % convert to double for plotting
condition(condition == 0) = NaN %set the 0s to NaN so they are not plotted
surf(x,y,condition,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','y')
xlim([-20 20 ])
ylim([-20 20])
view(0,90)
hold on 
%% plot eiv of linearized system 
C=-4;
w=-4;
A=[0 1+0.5*(C-w);1-0.5*(C-w) 0];
[eiv,lamda]=eig(A)
%technically we would need to transform it back to the original frame x so
%do the computation x=M*eiv to find the correct eigenavlues. 
%but here we are only looking at time t=pi/3 wich is the period of the flow, 
%Therefore M is the Unity matrix and we observe no change. 

%%
eivstable= eiv(:,2)
eivunstable=eiv(:,1)

%%
xaxisS=abs(eivstable(1,1));
yaxisS=abs(eivstable(2,1));

xaxisU=abs(eivunstable(1,1));
yaxisU=abs(eivunstable(2,1));

x0=-15:0.5:15;
k=-20:10:20 ;

Nk=length(k);
Nx=length(x0);
coeffS=yaxisS./xaxisS;
coeffU=yaxisU./xaxisU;
%%
% for i=1:Nk
    
y0f=coeffS.*x0
y0b=-coeffU.*x0
plot(x0,y0f,'r.')
hold on 
plot(x0,y0b,'k.')
hold on 
% end 
xlim([-20 20])
ylim([-20 20])
title('Points in direction of stable(red) and unstable(black) Eigenvectors')
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

sizexvector=size(x)
sizexvectorb=size(xb)
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
xlim([-20 20 ])
ylim([-20 20])
hold on 
end

%%
%%


function Q=Q(x,y)
alpha=0.005;
a=3*alpha*x.^2-3*y.^2*alpha;
b=-3/2-6*alpha*y.*x;
c=7/2-6*alpha*y.*x;
Q=a.^2+b.*c
end 

function dxdt=integrator(t,x)
    N=numel(x)/2;

    alpha=0.005;
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1 = sin(6*t).*xf+(cos(6*t)+5/2).*yf+alpha*xf.*(xf.^2-3*yf.^2);
    dx2 = (cos(6*t)-5/2).*xf-sin(6*t).*yf+alpha*(-yf.*(3*xf.^2-yf.^2));
    
    dx1(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    dx2(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;
    

 end
