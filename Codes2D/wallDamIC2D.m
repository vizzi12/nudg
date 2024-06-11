function [U] =  wallDamIC2D(x,y)
fd = @(p) drectangle(p,0,95,0,200);
fd2 = @(p) drectangle(p,95,105,87.5,162.5);

fd = @(p) dunion(fd(p),fd2(p));

h=zeros(size(x));
% % h(x<105) = 10;
% % h(~(x<105)) =5;
h(:) = 5*(fd([x(:),y(:)]) <= 0 )+5; 

h(:,mean(h)<9)=5;
h(:,mean(h)>9)=10;

U(:,:,1)=h; U(:,:,2)=zeros(size(x)); U(:,:,3)=zeros(size(x));
return