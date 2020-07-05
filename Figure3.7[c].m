
 %%
%%TRAJECTORIES FOR DIFFERENT INITIAL CONDITIONS
%change intial conditions for different trajectories 
x0=0.2
y0=0.1



tspan=linspace(0,400,10000001);

[t,x]=ode45(@(t,x) u(t,x), tspan,[x0;y0]);
title('Trajectories of particles near the origin for t=[0,300]')
plot(x(:,1),x(:,2),'y.','MarkerSize',1') %$This is plotting the  fluid particle 
 hold on 
 xlabel('x')
ylabel('y')
 xlim([-5 5])
ylim([-5 5])
%trajectory for the time interval [t_0,t_1]=[0,300]
%title('Fluid particle trajectory [to,t1]=[0,300] and x0=(-50,50)')




%%
plot(axisx(end,:),axisy(end,:),'k.','MarkerSize',2') %$This is plotting the  fluid particle 
 hold on
xlabel('x')
ylabel('y')

%%
function dxdt=u(t,x)
dxdt=zeros(2,1);
alpha=sin(log(t+1));
dxdt(1)= sin(4*t)*x(1)+(cos(4*t)+2)*x(2)+alpha.*(x(1).^2-x(2).^2);
dxdt(2) = (cos(4*t)-2)*x(1)-sin(4*t)*x(2)+alpha.*(-2*x(1).*x(2));
end 

function dxdt=integrator(t,x)
    N=numel(x)/2;
    a=log(t+1);

    alpha=sin(a);
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1 = sin(4*t).*xf+(cos(4*t)+2).*yf+alpha*(xf.^2-yf.^2);
    dx2 = (cos(4*t)-2).*xf-sin(4*t).*yf+alpha*(-2*xf.*yf);
    
    dx1(xf>10^10 | xf<-10^10 | yf>10^10 | yf<-10^10)=0;
    dx2(xf>10^10 | xf<-10^10 | yf>10^10 | yf<-10^10)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;
    

 end
function Q=Q(x,y)
t=0.1;
alpha=sin(log(t+1));
a=sin(4*t)+2*alpha*x;
b=(cos(4*t)-2)+alpha*(-2*y);
c=(cos(4*t)+2)-alpha*2*y;
Q=a.^2+b.*c;
 end