function S=parabolicSource2D(S)
g=9.81;
s1 = @(x,y,t,U) 0.*x;
s2 = @(x,y,t,U) -g*U(:,:,1).*(8<= x & x<=12).*(-0.10*x + 1.00);
s3 = @(x,y,t,U) 0.*x;

S.source = @(x,y,t,U) cat(3,s1(x,y,t,U),s2(x,y,t,U),s3(x,y,t,U));

end