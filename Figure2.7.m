title('The streamlines of the velocity field for t=0') %This is an autonomous streamline
%function from matlab
xspan = -15:2:15; yspan = -15:2:15;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;

u1=1.5*y+alpha*(x.^2-y.^2); % when plugging t=0 into u(x,t)
v=0.5*x+alpha*(-2*x.*y);
streamline(x,y,u1,v,startx,starty,[0.1,1000])
hold on

u2=-1.5*y-alpha*(x.^2-y.^2);
v2=-0.5*x-alpha*(-2*x.*y);
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.1,2000])
hold on 

xlabel('x')
ylabel('y')
