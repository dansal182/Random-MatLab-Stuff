% script: galerkin.m
% Resolver por medio de Galerkin
%9y''= -pi^*2y
% y(0)=-1
% y(1.5)=3
close all
n=20; % la particion tiene n+1 puntos
np=n+1;
h=(1.5)/np;

ya=-1;  yb=3;

%Diagonal de la matriz
alfa=((-2*pi*pi/(3*9))*h) + (2/h);
beta=((-pi*pi/(6*9))*h) - (1/h);

e=ones(n,1);
M=spdiags([beta*e alfa*e beta*e],-1:1,n,n);

d=zeros(n,1);
d(1)=-ya*beta;
d(n)=-yb*beta;

c=M\d;
cs=[ya , c',yb];

t=linspace(0,1.5,np+1);
s=linspace(0,1.5,np+1);
r=3*sin(pi/3*s)-cos(pi/3*s);

plot(t,cs,'-*r','Linewidth',3)
title('metodo de galerkin vs sol analitica')
hold on 
%plot(s,r,'c','Linewidth',3)
%legend('sol num','sol anal')

