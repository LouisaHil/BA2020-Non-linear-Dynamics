x0=-4:0.4:4
y0=-4:0.4:4
Nx=length(x0)
tspan=[0 pi/2]
[t,x]=ode45(@integrator,tspan,[x0(:);y0(:)]);
[m,n]=size(x)

title('The streamlines of the velocity field for the Pendulum ')
x1=x(1,1:n/2);
y1=x(1,n/2+1:end)
[xnew,ynew] = meshgrid(x1,y1);
startx=xnew;
starty=ynew;

u1=ynew;
v=-sin(xnew);
streamline(xnew,ynew,u1,v,startx,starty,[0.1,100])
hold on
u2=-ynew;
v2=sin(xnew);
streamline(xnew,ynew,u2,v2,startx,starty,[0.1,100])
hold on
xlabel('x')
ylabel('y')


xlim([-4 4 ])
ylim([-4 4])
%%
function dvdt=integrator(t,x) 
    N=numel(x)/2;
    dvdt=zeros(2*N,1);
    dx1=x(N+1:2*N,1); %v=xdot=solution of ode45(2)
    dx2=-sin(x(1:N,1)); %vdot=sin(x)
    
    dvdt=[dx1;dx2];
    
end
