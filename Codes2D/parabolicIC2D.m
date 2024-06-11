function [U] =  parabolicIC2D(x,y)
h=0.33-(8<= x & x<=12).*(0.2-0.05*(x-10).^2);
% q=0.18*ones(size(x));
q=zeros(size(x));
p=zeros(size(x));

U(:,:,1)=h; U(:,:,2)=q; U(:,:,3)=p;
return