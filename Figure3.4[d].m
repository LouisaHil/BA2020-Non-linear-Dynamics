%%TRAJECTORIES FOR DIFFERENT INITIAL CONDITIONS
%% change the initial positions for more trajectories 
x0=0.1 
y0=0.2
tspan=linspace(0,300,40001)


[t,x]=ode45(@u, tspan,[x0;y0]);
plot(x(:,1),x(:,2),'k.','MarkerSize',1') %$This is plotting the  fluid particle 
 hold on
xlabel('x')
ylabel('y')
plot(x(:,1),x(:,2),'k.','MarkerSize',1') %$This is plotting the  fluid particle 
 hold on

xlabel('x')
ylabel('y')

%%
function dxdt=u(t,x)
dxdt=zeros(2,1);
 alpha=0.005;
 dxdt(1)= sin(4*t)*x(1)+(cos(4*t)+2)*x(2)+alpha*(x(1).^2-x(2).^2)+sin(t^2);


 dxdt(2) = (cos(4*t)-2)*x(1)-sin(4*t)*x(2)+alpha*(-2*x(1).*x(2))+sin(t^2);
end 
function dxdt=integrator(t,x)
    N=numel(x)/2;

    alpha=0.005;
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1 = sin(4*t).*xf+(cos(4*t)+2).*yf+alpha*(xf.^2-yf.^2)+sin(t.^2);
    dx2 = (cos(4*t)-2).*xf-sin(4*t).*yf+alpha*(-2*xf.*yf)+sin(t.^2);
    
    dx1(xf>10^12 | xf<-10^12 | yf>10^12 | yf<-10^12)=0;
    dx2(xf>10^12 | xf<-10^12 | yf>10^12 | yf<-10^12)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;
    

end

 function Q=Q(x,y)
alpha=0.005;
t=0.1;
a=sin(4*t)+2*alpha*x;
b=(cos(4*t)-2)+alpha*(-2*y);
c=(cos(4*t)+2)-alpha*2*y;
Q=a.^2+b.*c;
 end