%% streamlines
alpha=-0.01;
title('The streamlines of the velocity field for t=0, Okubo Weiss region and KAM Map for alpha=-0.01')
xspan = -20:3:20; yspan = -20:3:20;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;

u1=(3*y)+alpha*(x.*(x.^2-3*y.^2));
v=-x+alpha*(-y.*(3*x.^2-y.^2));
streamline(x,y,u1,v,startx,starty,[0.006,800])
hold on

u2=(-3)*y-alpha*(x.*(x.^2-3*y.^2));
v2=x-alpha*(-y.*(3*x.^2-y.^2));
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.006,800])
hold on 

xlabel('x')
ylabel('y')
hold on 


%% plot OW
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

%%
var1=20;
var2=20;
x=linspace(-15,15,var1);
y=linspace(-15,15,var2);
sizepoints=var1*var2;
[X,Y]=ndgrid(x,y);
tspan=0:pi/5:90*pi;

[t,x]=ode45(@integrator,tspan,[X(:);Y(:)]);
sizexvector=size(x)
x=reshape(x,sizexvector(1),sizexvector(2));



axisx = x(1:end,1:sizexvector(2)/2)
axisy = x(1:end,(sizexvector(2)/2)+1:sizexvector(2))

plot(axisx(300:end,:),axisy(300:end,:),'k.','MarkerSize',0.5')

hold on
xlim([-20 20 ])
ylim([-20 20])

% title('KAM Curves')
% xlabel('x')
% ylabel('y')


%%



function Q=Q(x,y)
alpha=-0.01;
a=3*alpha*x.^2-3*y.^2*alpha;
b=-1-6*alpha*y.*x;
c=3-6*alpha*y.*x;
Q=a.^2+b.*c
end 

function dxdt=integrator(t,x)
    N=numel(x)/2;

    alpha=-0.01;
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1 = sin(10*t).*xf+(cos(10*t)+2).*yf+alpha*xf.*(xf.^2-3*yf.^2);
    dx2 = (cos(10*t)-2).*xf-sin(10*t).*yf+alpha*(-yf.*(3*xf.^2-yf.^2));
    
    dx1(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    dx2(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;
    

 end