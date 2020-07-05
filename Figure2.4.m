%% plot eiv of linearized system 
C=-4; %look at the velocity field and determine the constants
w=-4;
A=[0 1+0.5*(C-w);1-0.5*(C-w) 0]; %plugged into the transformation %matrix
[eiv,lamda]=eig(A)
%technically we would need to transform it back to the original frame x so we would need to 
%do the computation x=M*eiv to find the correct eigenvalues. 
%but here we are only looking at time t=pi/2 which is the period of the flow, 
%Therefore M is the Unity matrix and we observe no change so we don't need to do the transformation. 

%%
eivstable= eiv(:,2)
eivunstable=eiv(:,1)

%%
xaxisS=abs(eivstable(1,1));
yaxisS=abs(eivstable(2,1));

xaxisU=abs(eivunstable(1,1));
yaxisU=abs(eivunstable(2,1));

x0=-15:0.005:15
k=0

Nk=length(k)
Nx=length(x0)
coeffS=yaxisS./xaxisS; %in this case we have coeff=1 but this is not %always the case
coeffU=yaxisU./xaxisU;
for i=1:Nk
    
    y0f=coeffS.*x0+k(1)
    y0b=-coeffU.*x0+k(1)
    plot(x0,y0f,'r.')
    hold on 
    plot(x0,y0b,'k.')
    hold on 
end 
xlim([-20 20])
ylim([-20 20])


