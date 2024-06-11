function [U] =  analyticIC2D(x,y)
t=0;
% h=cos(-t + x).*cos(-t + y) + 2;
% q=sin(x-t);
% p=sin(y-t);

h=ones(size(x))*5;
h(x < 0)=10;
q=zeros(size(x));
p=zeros(size(y));

U(:,:,1)=h; U(:,:,2)=q; U(:,:,3)=p;

return