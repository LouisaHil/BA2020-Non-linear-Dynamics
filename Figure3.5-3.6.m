%%
% code for the instanteous streamlines at time t 

t=0.1
title('The streamlines and Okubo Weiss elliptic region of the velocity field for t=0.1')
xspan = -60:6:60; yspan = -60:6:60;
[x,y] = meshgrid(xspan,yspan);
startx=x;
starty=y;
alpha=sin(log(t+1));
u1=sin(4*t)*x+(cos(4*t)+2)*y+alpha*(x.^2-y.^2);
v=(cos(4*t)-2)*x-sin(4*t)*y+alpha*(-2*x.*y);
streamline(x,y,u1,v,startx,starty,[0.1,400])
hold on

u2=-sin(4*t)*x-(cos(4*t)+2)*y-alpha*(x.^2-y.^2);
v2=-(cos(4*t)-2)*x+sin(4*t)*y-alpha*(-2*x.*y);
startx=-x;
starty=-y;
streamline(x,y,u2,v2,startx,starty,[0.1,400])
hold on 

xlabel('x')
ylabel('y')
% try out okubo weiss 
% code for okubo Weiss region 

r = -60:0.09:60;  %Okubo Weiss elliptic criterion for t=0
[x,y] = meshgrid(r,r);

condition=Q(x,y)<0;
condition = double(condition); % convert to double for plotting
condition(condition == 0) = NaN; %set the 0s to NaN so they are not plotted
surf(x,y,condition,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','y')
view(0,90)
hold on 
%%
var=40;

%var2=50;
xi=linspace(-60,60,var);
yi=linspace(-60,60,var);
sizepoints=var*var;
[X,Y]=ndgrid(xi,yi);
tspan=linspace(0,0.1,4001)

rho=0.1;

[m, n]=size(X); 
Nrad = 4; %nb of points next to the actual point so new initial conditions
xt=zeros(m,n,Nrad);
yt=zeros(m,n,Nrad);

for k=1:Nrad
    xt(:,:,k) = X + rho*cos( (k-1)*pi/2 ); % we put in a trigonometric function because we don't have simultaneously y and x changing
    yt(:,:,k) = Y + rho*sin( (k-1)*pi/2 );%rho.x is delat1 and rho.y is delta2. 
end

[t,x]=ode45(@integrator,tspan,[xt(:);yt(:)]); 
[Nt,Nx]=size(x)

%%check if for t_End and t_End-1,2,3 we don't have stagnating values 
[Nt,Nx]=size(x)
xvector=x(:,1:end/2);
yvector=x(:,end/2+1:end);



i=Nt-100;

condition= i<=(Nt-1);
    for condition=1
        if xvector(i,:)==xvector(i+1,:) & yvector(i,:)==yvector(i+1,:)
        enval=Nt-i;
        else
        endval=Nt;
        i=i+1;
        end
    end 

endval
if endval<Nt
    xvector=x(1:endval,1:Nx/2);
    yvector=x(1:endval,Nx/2+1:end);
end 
% check that we don't have coordinates that are outsde of the domain, hence been put to 0 
for j=1:Nx/2
if xvector(endval,j)~=0 & yvector(endval,j)~=0
    disp(' All particles are within our region of interest. No error')
else 
    disp('Error, particles diverge to infinity')
    break
end 
end


 
 
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
 
 % we have #var different initial conditions and for each initial condition a 2*2 Matrix. 

 C11=F11.*F11+F21.*F21;
 C12=F11.*F12+F21.*F22;
 C22=F12.*F12+F22.*F22;
 % prove symetric and positve definite 
 detC=C11.*C22-C12.^2;
 traceC=C11+C22; 


lambda1 = 0.5*traceC-sqrt((0.5*traceC).^2-detC);
lambda2 = 0.5*traceC+sqrt((0.5*traceC).^2-detC);

%
if lambda1>lambda2 
    lambda=lambda1;
    disp('lambda1 is the bigger eigenvalue')
else 
    lambda=lambda2;
    disp('lambda2 is the bigger eigenvalue')
end
%
 %%The eigenvalue of the Cauchy strain tensor needs to be bigger than 1
 %%for it it be able to be the maximum eigenvalue of relative stretching of
 %%infinitesimal perturbations.
 

if lambda<1
    disp('Error')
end 


FTLEf=1/(2*(t1-t0))*log(lambda); %forward in time
condf=FTLEf>0;
if condf==0
    disp('Error')
end 
    
FTLEb=1/(2*(t0-t1))*log(lambda); %backward in time 
condb=FTLEb<0;
if condb==0
    disp('Error')
end 

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