%%MATLAB CODE. DENSITY EQUATION 
deltai=0.1
x0=-15:deltai:15
y0=x0;
Nx=length(x0);
tspan=[0 pi/2];
[t,x]=ode45(@integrator,tspan,[x0(:);y0(:)]);
sizexvector=size(x);
x=reshape(x,sizexvector(1),sizexvector(2));
xi = x(sizexvector(1),1:sizexvector(2)/2);
yi = x(sizexvector(1),(sizexvector(2)/2)+1:sizexvector(2));


%% looking for t=pi/2
x=transpose([x0 y0]);
t=pi/2
u=integrator(t,x);

[m,n]=size(u)

ui=u(1:m/2,1)
vi=u(m/2+1:m,1)

%density integrator 
%% initial conditions a0 b0 and aN need to be known
[N,M]=size(ui)
%%
aN=2
a0=3
b0=1
%%the initial conditions need to put in 
k(N,1)=(1-(aN/a0)^(-1))*ui(N,1)/deltai;

for j=(N-1):-1:1
    k(j,1)=k(j+1,1)*ui(j+1,1)/ui(j,1);
        for i=1:N
            rhoi(i)=a0*b0*(k(j,1)*deltai*(1\ui(i,1)-1/vi(i,1))-k(j,1)^2*deltai^2/(ui(i,1)*vi(i,1)))^(-i);
        end 
end
rhoi
%%
plot(rhoi(:),'.')

%%
function dxdt=integrator(t,x)
    N=numel(x)/2;
   

    alpha=-0.1;
    dxdt=zeros(2*N,1);
    xf=x(1:N,1);
    yf=x(N+1:2*N,1);
    
    dx1=sin(4*t)*xf+(cos(4*t)+2)*yf+alpha*(xf.^2-yf.^2);
    dx2=(cos(4*t)-2)*xf-sin(4*t)*yf+alpha*(-2*xf.*yf);
    
    dx1(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    dx2(xf>150 | xf<-150 | yf>150 | yf<-150)=0;
    
    dxdt(1:N,1) = dx1;
    dxdt(N+1:2*N,1) = dx2;

    dxdt=[dx1;dx2];
 end


