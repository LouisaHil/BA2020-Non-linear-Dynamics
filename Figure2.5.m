
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

%%
C=-4;
w=-4;
A=[0 1+0.5*(C-w);1-0.5*(C-w) 0];
[eiv,lamda]=eig(A)
eiv

eivstable= eiv(2,:)
eivunstable=eiv(1,:)


xaxisS=abs(eivstable(1,1));
yaxisS=abs(eivstable(1,2));

xaxisU=abs(eivunstable(1,1));
yaxisU=abs(eivunstable(1,2));

x0=-15:0.005:15

k=-20:10:20 

Nk=length(k)
Nx=length(x0)
coeffS=yaxisS./xaxisS;
coeffU=yaxisU./xaxisU;
% for i=1:Nk
%     
% y0f=coeffS.*x0+k(i)
% y0b=-coeffU.*x0+k(i)
% % plot(x0,y0f,'r.')
% % hold on 
% % plot(x0,y0b,'k.')
% % hold on 
% end 

%% 

Nk=length(k)
Nx=length(x0)
tspan=[0 pi/2] %forward in time
tspanb=[pi/2 0]; %backward in time 
%run the ode45 solver and find x and y for each initial condition, Nk parallel lines in the direction of the Eigenvectors
for i=1:Nk 
    
 y0f=coeffS.*x0+k(i);
y0b=-coeffU.*x0+k(i);
[t,x]=ode45(@integrator,tspan,[x0(:);y0f(:)]);
[tb,xb]=ode45(@integrator,tspanb,[x0(:);y0b(:)]);

sizexvector=size(x)
sizexbvector=size(xb)
x=reshape(x,sizexvector(1),sizexvector(2));
xb=reshape(xb,sizexbvector(1),sizexbvector(2));

%extract the x %coordinate and the y coordinate from the vector x 
axisx = x(sizexvector(1),2:sizexvector(2)/2); 
axisy = x(sizexvector(1),(sizexvector(2)/2)+2:sizexvector(2));
axisxb = xb(sizexbvector(1),1:sizexbvector(2)/2);
axisyb = xb(sizexbvector(1),(sizexbvector(2)/2)+1:sizexbvector(2));

plot(axisx,axisy,'k.')
hold on 
plot(axisxb,axisyb,'r.')
hold on 
xlim([-20 20 ])
ylim([-20 20])
hold on 
i
end 

%%
alpha=-0.1;




r = -20:0.05:20;  %Okubo Weiss elliptic criterion for t=0
[x,y] = meshgrid(r,r);
condition=(y+15).*(y-5)+x.^2<0;
condition = double(condition); % convert to double for plotting
condition(condition == 0) = NaN;%set the 0s to NaN so they are not plotted
surf(x,y,condition,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','y')
view(0,90)

hold on 

title('The streamlines of the velocity field for t=0, OW region, Poincare map of a fixed saddle point')
xspan = -20:3:20; yspan = -20:3:20;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;

u1=3*y+alpha*(x.^2-y.^2);
v=-x+alpha*(-2*x.*y);
streamline(x,y,u1,v,startx,starty,[0.1,1000])
hold on
u2=-3*y-alpha*(x.^2-y.^2);
v2=x-alpha*(-2*x.*y);
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.1,2000])
hold on 

xlabel('x')
ylabel('y')

%%
function dxdt=integrator(t,x)
    N=numel(x)/2;
    

    alpha=-0.1;
    dxdt=zeros(2*N,1);
    
    xf=x(1:N,1); %create a vector where the frist half entries are the x coordinates
    yf=x(N+1:2*N,1);%create a vector where the second half entries are the y coordinates
    
    dx1=sin(4*t)*xf+(cos(4*t)+2)*yf+alpha*(xf.^2-yf.^2);
    dx2=(cos(4*t)-2)*xf-sin(4*t)*yf+alpha*(-2*xf.*yf);
    
    dx1(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    dx2(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;

    dxdt=[dx1;dx2];
 end


