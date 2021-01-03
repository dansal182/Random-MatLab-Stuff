function[wp]=ftar1(t,w)
   %EDo para problema de valores iniciales
%valores en la frontera
%w(2)==y' ; wp(1)==s'
wp=zeros(2,1);
wp(1)=-(pi*pi*w(2))/9;
wp(2)=w(1);
end