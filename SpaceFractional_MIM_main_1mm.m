clc
clear
alpha=1.6;%1.6
lambda=0.5;
mu=0.35;
beta=1;
v=0.05763;%0.2854;
D=0.06;%0.8563;
%%
a=1;
b=1-a;
xL=0;
xR=20;
xN=200;
dx=(xR-xL)/xN;
tL=0;
tM=400;
tN=400;
dt=(tM-tL)/tN;
x=xL:dx:xR;
t=tL:dt:tM;
%%
k1=-v*dt/dx;
k2=D*dt*a/(dx)^alpha;
k3=D*dt*b/(dx)^alpha;
C=zeros(xN+1,tN+1);
Cim=zeros(xN+1,tN+1);
M=zeros(xN+1,xN+1);
rhs=zeros(xN+1,1);
C(1,1)=1;%Instantaneous source
g=zeros(xN+1,1);
g(1)=1;
for i=1:length(g)
    g(i+1)=-g(i)*(alpha-i+1)/i;
end
M1=zeros(xN+1,xN+1);
M2=zeros(xN+1,xN+1);
M3=zeros(xN+1,xN+1);
for i=2:xN+1
    for j=1:xN+1
        if i==j
            M1(i,j)=1+k1;
        elseif i==j-1
            M1(i,j)=-k1;
        else
            M1(i,j)=0;
        end
    end
end
for i=2:xN
    M2_trans=zeros(1,xN+1);
    M3_trans=zeros(1,xN+1);
    for j=1:i+1
        M2_trans(i-j+2)=-k2*g(j);
    end
    f=1;
    for j=i-1:xN+1
        M3_trans(j)=-k3*g(f);
        f=f+1;
    end
    M2(i,:)=M2_trans;
    M3(i,:)=M3_trans;
end
M=M1+M2+M3;
M(1,:)=0;M(xN+1,:)=0;
M(1,1)=1;M(xN+1,end)=1; M(xN+1,end-1)=-1;
rhs=zeros(xN+1,1);
for j=1:tN
    for i=1:xN+1
        Cim(i,j+1)=Cim(i,j)+dt*lambda*C(i,j)-dt*mu*Cim(i,j);
    end
    for i=1:xN+1
        rhs(i)=C(i,j)-beta*(Cim(i,j+1)-Cim(i,j));
    end
     if j>1
     rhs(1)=1;
     rhs(xN+1)=0;
     end
    C(:,j+1)=M\rhs;  
end
%C=C+Cim;
hold on
T=155
t_plot=ceil(T/tM*tN);
x_plot=x(25:end)';
y_plot=C(25:end,t_plot)/max(C(25:end,t_plot));
data=[x_plot y_plot];
plot(x(25:end),C(25:end,t_plot)/max(C(25:end,t_plot)));