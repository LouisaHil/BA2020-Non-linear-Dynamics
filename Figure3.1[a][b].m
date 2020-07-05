var=100;

%var2=50;
xi=linspace(-pi,pi,var);
yi=linspace(-pi,pi,var);
sizepoints=var*var;
[X,Y]=ndgrid(xi,yi);
tspan=0:20;

rho=0.1;

[m, n]=size(X); %xi being initial conditions 200 points with each 2 coordinates
Nrad = 4; %nb of points next to the actual point so new initial conditions
xt=zeros(m,n,Nrad);
yt=zeros(m,n,Nrad);

for k=1:Nrad
    xt(:,:,k) = X + rho*cos( (k-1)*pi/2 ); % we put in a
    % trigonometric function because we don't have simultaneously y and x changing
    yt(:,:,k) = Y + rho*sin( (k-1)*pi/2 );%rho.x is delat1 and rho.y is delta2. 
end
%%
[t,x]=ode45(@integrator,tspan,[xt(:);yt(:)]); %just for my first initial point[xi(1),yi(1]
%[a,b]=size(x);
%%
xvector=x(end,1:end/2);
yvector=x(end,end/2+1:end);
xvector = reshape(xvector,m, n, Nrad);
yvector = reshape(yvector, m,n, Nrad);


%% Compute gradient of flow  and eigenvalues of gradient FTLE computation 

t0=tspan(1);
t1=tspan(end);

 gradF_ic=zeros(2,2);
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
%%
if lambda1>lambda2 % in this case lambda2 is always bigger than lambda1 
    lambda=lambda1;
else 
    lambda=lambda2;
end



FTLEf=1/(2*(t1-t0))*log(lambda); %forward in time

FTLEb=1/(2*(t0-t1))*log(lambda); %backward in time 

plot (FTLEf,X,'k.')
pcolor(X,Y,FTLEf)
shading interp
hold on 
plot(FTLEb,X,'r.')

colorbar

colormap('jet');
title('FTLE in the flow of a simple pendulum')
 xlabel('v=xdot')
 ylabel('vdot=sin(x)')
%%

function dvdt=integrator(t,x) 
    N=numel(x)/2;
    dvdt=zeros(2*N,1);
    dx1=x(N+1:2*N,1); %xdot=v=solution of ode45(2)
    dx2=-sin(x(1:N,1)); %vdotx=sin(x)
    
    dvdt=[dx1;dx2];
    
end


