

var=70;

%var2=50;
xi=linspace(-200,200,var);
yi=linspace(-200,200,var);
sizepoints=var*var;
[X,Y]=ndgrid(xi,yi);
tspan=linspace(0,0.1,4001);

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
%%TRAJECTORIES FOR DIFFERENT INITIAL CONDITIONS
x0=0
y0=100
tspan=linspace(0,0.1,40001)


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