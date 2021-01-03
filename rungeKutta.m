function [t,y] = rungeKutta(fname,a,b,y0,n)
%erunge kutta para aproximar solucion de una edo con valor inicial

t=linspace(a,b,n)';
m=length(y0);
y=zeros(n,m);
y(1,:)=y0';

h=t(2)-t(1);
for k=1:n-1
    s1=feval(fname,t(k), y(k,:))';
    s2=feval(fname,t(k)+(h/2),y(k,:)'+(h/2)*s1')';
    s3=feval(fname,t(k)+(h/2),y(k,:)'+(h/2)*s2')';
    s4=feval(fname,t(k)+h, y(k,:)'+h*s3')';
    y(k+1,:)=y(k,:)+(h/6)*(s1+2*s2+2*s3+s4);
end
end
