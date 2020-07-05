figure(1)
[t,x]=ode45(@(t,x) u(t,x), [0 200],[2;0]);
plot(x(:,1),x(:,2)) %$This is plotting the  fluid particle 
%trajectory for the time interval [t_0,t_1]=[0,200]
title('Fluid particle trajectory [to,t1]=[0,200] and x0=(2,0)')
xlabel('x')
ylabel('y')

[x,y] = meshgrid(-5:1:5,-5:1:5);
startx=x;
starty=y;
u1=-1.5*y;
v=-0.5*x;
%Instantaneous streamlines 
figure(2)
title('The instanteneous streamlines of the velocity field for t=0')
streamline(x,y,u1,v,startx,starty,[0.1,2000])
hold on 
u2=1.5*y;
v2=0.5*x;
streamline(x,y,u2,v2,startx,starty,[0.1,2000])
xlabel('x')
ylabel('y')

%%
function dxdt=u(t,x)
    dxdt=zeros(2,1);
    dxdt(1)=sin(4*t)*x(1)+(cos(4*t)+0.5)*x(2);
    dxdt(2)=(cos(4*t)-0.5)*x(1)-sin(4*t)*x(2);
end