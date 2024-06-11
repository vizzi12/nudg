clear all
close all
clc

MCELLS=1;
CHALEN=0;
hL=1.0;
uL=-5.0;
hR=1.0;
uR=5.0;
gate=0;
t=20;

vR=0;
vL=0;

[h,u,v]=SWRPExact(hL,uL,vL,hR,uR,vR,MCELLS,CHALEN,gate,t);

figure
plot(h)
figure
plot(u)