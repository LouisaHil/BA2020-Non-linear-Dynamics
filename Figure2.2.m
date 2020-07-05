alpha=-0.1;
figure()
r = -20:0.05:20;  %Okubo Weiss elliptic criterion for t=0
[x,y] = meshgrid(r,r);
condition=(y+15).*(y-5)+x.^2<0; % we have implemented directly 
%the condition that we have computed analytically. 
condition = double(condition); % convert to double for plotting
condition(condition == 0) = NaN;%set the 0s to NaN so they are not plotted
surf(x,y,condition,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','y')
view(0,90)

hold on 

title('The streamlines of the velocity field for t=0') %This is an autonomous streamline
%function from matlab
xspan = -20:2:20; yspan = -20:2:20;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;

u1=3*y+alpha*(x.^2-y.^2); % when plugging t=0 into u(x,t)
v=-x+alpha*(-2*x.*y);
streamline(x,y,u1,v,startx,starty,[0.1,1000])
hold on

u2=-3*y-alpha*(x.^2-y.^2);
v2=x-alpha*(-2*x.*y);
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.1,1000])
hold on 

xlabel('x')
ylabel('y')
