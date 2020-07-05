

%%
%%TRAJECTORIES FOR DIFFERENT INITIAL CONDITIONS
%%change initial conditions for different trajectories 
x0=-1.8
y0=1.8



tspan=linspace(0,300,1000001);

[t,x]=ode45(@(t,x) u(t,x), tspan,[x0;y0]);
title('Trajectories of particles near the origin for t=[0,300]')
plot(x(:,1),x(:,2),'k.','MarkerSize',1') %$This is plotting the  fluid particle 
 hold on 
 xlabel('x')
ylabel('y')
 xlim([-5 5])
ylim([-5 5])
%trajectory for the time interval [t_0,t_1]=[0,300]
%title('Fluid particle trajectory [to,t1]=[0,300] and x0=(-50,50)')


%%
function dxdt=u(t,x)

dxdt=zeros(2,1);
 alpha=sin(log(t+1));
 dxdt(1,:) = sin(10*t).*x(1)+(cos(10*t)+2).*x(2)+alpha.*x(1).*(x(1).^2-3*x(2).^2);
 dxdt(2,:) = (cos(10*t)-2).*x(1)-sin(10*t).*x(2)+alpha.*(-x(2).*(3*x(1).^2-x(2).^2));

 
end 
function dxdt=integrator(t,x)
    N=numel(x)/2;

    alpha=sin(log(t+1));
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1 = sin(10*t).*xf+(cos(10*t)+2).*yf+alpha.*xf.*(xf.^2-3*yf.^2);
    dx2 = (cos(10*t)-2).*xf-sin(10*t).*yf+alpha.*(-yf.*(3*xf.^2-yf.^2));
    
    dx1(xf>10^6 | xf<-10^6 | yf>10^6 | yf<-10^6)=0;
    dx2(xf>10^6 | xf<-10^6| yf>10^6 | yf<-10^6)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;
    

 end

function Q=Q(x,y)
t=0.01
alpha=sin(log(t+1));
a=sin(10*t)+3*alpha*x.^2-3*y.^2*alpha;
b=(cos(10*t)-2)-6*alpha*y.*x;
c=(cos(10*t)+2)-6*alpha*y.*x;
Q=a.^2+b.*c
end 
