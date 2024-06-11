function [U] =  circularDamIC2D(x,y)
  
R=2.5;
h_init = @(x,y) 2*(x.^2+y.^2<=R^2)+0.5; 

U(:,:,1) = h_init(x,y); 
U(:,:,2)= zeros(size(x)); U(:,:,3) = zeros(size(x));
return