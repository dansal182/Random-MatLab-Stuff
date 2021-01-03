%met new n dimension para sitemaF(x)=0
function[Indicador,r]=Newton_ND(x0,MaxNumIter,Tol)
Indicador=0; r=Inf; %0 cuando no hay convergencia

for i=1:MaxNumIter
    x=x0-inv(Fp(x0))*F(x0);
    if norm(x-x0)/norm(x)<=Tol
        Indicador=1;
        r=x;
        break
    end
    x0=x;
end
end

function[y]=F(u)
m=length(u);
n=m+1;
h=1/n;
z1=-2/(h*h);
z2=1/(h*h);
y=zeros(m,1);
y(1)=z1*u(1)+u(1)^2+z2*u(2)-f(h);
y(n-1)=z2*u(n-2)+z1*u(n-1)+u(n-1)^2-f((n-1)*h)-z2;
for i=2:n-2
    y(i)=z2*u(i-1)+z1*u(i)+u(i)^2+z2*u(i+1)-f(i*h);
end
end

%la derivada
function[J]=Fp(x)
    h=eps^(1/3);
    n=length(x);
    J=zeros(n,n);
    for j=1:n
        ej=zeros(n,1);
        ej(j)=1;
        J(:,j)=(F(x+h*ej)-F(x-h*ej))/(2*h);
    end
end

function [y]=f(x)
    y=x^4-2;
end
