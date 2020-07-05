% code for the instanteous streamlines at time t 
t=0.01;
alpha=sin(log(t+1));

title('The streamlines  and Okubo Weiss region of the velocity field for t=0.01')
xspan = -50:5:50; yspan = -50:5:50;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;

u1=sin(10*t).*x+(cos(10*t)+2).*y+alpha.*x.*(x.^2-3*y.^2);
v=(cos(10*t)-2).*x-sin(10*t).*y+alpha.*(-y.*(3*x.^2-y.^2));
streamline(x,y,u1,v,startx,starty,[0.02,400])
hold on

u2=-sin(10*t).*x-(cos(10*t)+2).*y-alpha.*x.*(x.^2-3*y.^2);
v2=-(cos(10*t)-2).*x+sin(10*t).*y-alpha.*(-y.*(3*x.^2-y.^2));
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.02,400])
hold on 

xlabel('x')
ylabel('y')
% try out okubo weiss 
% code for okubo Weiss region 

r = -200:0.5:200;  %Okubo Weiss elliptic criterion for t=0
[x,y] = meshgrid(r,r);

condition=Q(x,y)<0;
condition = double(condition); % convert to double for plotting
condition(condition == 0) = NaN; %set the 0s to NaN so they are not plotted
surf(x,y,condition,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','y')
view(0,90)
hold on 

%%
var=60;

%var2=50;
xi=linspace(-20,20,var);
yi=linspace(-20,20,var);
sizepoints=var*var;
[X,Y]=ndgrid(xi,yi);
tspan=linspace(0,0.06,4001);

rho=0.1;

[m, n]=size(X); 
Nrad = 4; %nb of points next to the actual point so new initial conditions
xt=zeros(m,n,Nrad);
yt=zeros(m,n,Nrad);

for k=1:Nrad
    xt(:,:,k) = X + rho*cos( (k-1)*pi/2 ); % we put in a trigonometric function because we don't have simultaneously y and x changing
    yt(:,:,k) = Y + rho*sin( (k-1)*pi/2 );%rho.x is delat1 and rho.y is delta2. 
end
%%
[t,x]=ode45(@integrator,tspan,[xt(:);yt(:)]); 
%%
% check if for t_n and tn-1,2 we don't have stagnating values that have the same location 
[Nt,Nx]=size(x)
xvector=x(:,1:end/2);
yvector=x(:,end/2+1:end);
i=(Nt-10);
condition= i<=(Nt-1);
for condition=1
     if xvector(i,:)==xvector(i+1,:) & yvector(i,:)==yvector(i+1,:)
        endval=i;
     end
     endval=Nt;
     i=i+1;
end 
endval

 %
xvector=x(endval,1:end/2);
yvector=x(endval,end/2+1:end);
xvector = reshape(xvector,m, n, Nrad);
yvector = reshape(yvector, m,n, Nrad);


    


%
 % %% Compute gradient of flow  and eigenvalues of gradient FTLE computation 

t0=tspan(1);
t1=tspan(end);


 F11=(xvector(:,:,1)-xvector(:,:,3))/(2*rho);
 F12=(xvector(:,:,2)-xvector(:,:,4))/(2*rho);
 F21=(yvector(:,:,1)-yvector(:,:,3))/(2*rho);
 F22=(yvector(:,:,2)-yvector(:,:,4))/(2*rho);
 
 % we have var different initial conditions and for each initial condition a 2*2 Matrix. 

 C11=F11.*F11+F21.*F21;
 C12=F11.*F12+F21.*F22;
 C22=F12.*F12+F22.*F22;
 % prove symetric and positve definite 
 detC=C11.*C22-C12.^2;
 traceC=C11+C22; 


lambda1 = 0.5*traceC-sqrt((0.5*traceC).^2-detC);
lambda2 = 0.5*traceC+sqrt((0.5*traceC).^2-detC);
%
if lambda1>lambda2 % in this case lambda2 is always bigger than lambda1 
    lambda=lambda1;
else 
    lambda=lambda2;
end



FTLEf=1/(2*(t1-t0))*log(lambda); %forward in time

FTLEb=1/(2*(t0-t1))*log(lambda); %backward in time 


%
figure(1)
title("FTlE forward")
pcolor(X,Y,FTLEf)

shading interp
hold on 

title("FTlE forward")
colorbar

colormap('jet');
figure (2)

pcolor(X,Y,FTLEb)

shading interp
hold on 
title("FTlE backward")

colorbar
colormap('jet');


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
