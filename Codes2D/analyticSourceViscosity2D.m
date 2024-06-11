function S=analyticSourceViscosity2D(S)
g=9.81;
c=1;
s1 = @(x,y,t,U) -sin(-x+t).*cos(-y+t)-sin(-y+t).*cos(-x+t)+cos(-x+t)+cos(-y+t)+2*cos(-x+t).*cos(-y+t);
s2 = @(x,y,t,U) -cos(-x+t)-1./(cos(-x+t).*cos(-y+t)+2).^2.*sin(-x+t).^3.*cos(-y+t)-2./(cos(-x+t).*cos(-y+t)+2).*sin(-x+t).*cos(-x+t)+g*(cos(-x+t).*cos(-y+t)+2).*sin(-x+t).*cos(-y+t)-1./(cos(-x+t).*cos(-y+t)+2).^2.*sin(-y+t).^2.*sin(-x+t).*cos(-x+t)-1./(cos(-x+t).*cos(-y+t)+2).*cos(-y+t).*sin(-x+t)-sin(-x+t);
s3 = @(x,y,t,U) -cos(-y+t)-1./(cos(-x+t).*cos(-y+t)+2).^2.*sin(-y+t).*sin(-x+t).^2.*cos(-y+t)-1./(cos(-x+t).*cos(-y+t)+2).*sin(-y+t).*cos(-x+t)-1./(cos(-x+t).*cos(-y+t)+2).^2.*sin(-y+t).^3.*cos(-x+t)-2./(cos(-x+t).*cos(-y+t)+2).*sin(-y+t).*cos(-y+t)+g*(cos(-x+t).*cos(-y+t)+2).*sin(-y+t).*cos(-x+t)-sin(-y+t);

S.source = @(x,y,t,U) cat(3,s1(x,y,t,U),s2(x,y,t,U),s3(x,y,t,U));

end