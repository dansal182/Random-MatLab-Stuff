%script disparo
close all

fname='ftar1';
a=0; b=1.5;
ya=-1; yb=3;
n=16;

sa=2; sa2=4;

%met numerico

[t,z]=rungeKutta(fname,a,b,[sa, ya]',n);

[t1,z2]=rungeKutta(fname,a,b,[sa2, ya]',n);

v=z(n,2)-yb;

v2=z2(n,2)-yb;

%biseccion

alfa=2;  beta=4;

falfa=v; fbeta=v2;

tol=1.e-08; maxiter=50; iter=0;

while(abs(beta-alfa)>tol && iter <maxiter && falfa*fbeta<0)
    c=(alfa+beta)/2;
    [t,zz]=rungeKutta(fname,a,b,[c ya]',n);
    fc=zz(n,2)-yb;
    if (falfa*fc<0)
        beta=c;
        fbeta=fc;
    else
        alfa=c;
        falfa=fc;
    end
    iter=iter+1;
   
end

if(fc==0)
    raiz=c;
else
    raiz=(alfa+beta)/2;
    [t,zz]=rungeKutta(fname,a,b,[raiz ya],n);
end


%sol num
yy=zz(:,2);

%sol. exacta
s=linspace(a,b,n);
r=3*sin(pi/3*s)-cos(pi/3*s);

%graficamos
subplot(2,1,1)
plot(t,yy,'b',a,ya,'dr',b,yb,'dr','Linewidth',3)
title('Método del disparo');
subplot(2,1,2)
plot(s,r,'m',a,ya,'dr',b,yb,'dr','Linewidth',3)
title('Solución exacta')
hold on
figure(2)
plot(t,yy,'-*b','Linewidth',3)
hold on
plot(s,r,'-m','Linewidth',3)
legend('sol num','sol tray','solanal','sol')



