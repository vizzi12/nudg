function [U] =  guusIC2D(x,y)
g=9.81;
h=ones(size(x));
q=3*sqrt(g)*ones(size(x));
% q=zeros(size(x));
p=zeros(size(x));

U(:,:,1)=h; U(:,:,2)=q; U(:,:,3)=p;
return