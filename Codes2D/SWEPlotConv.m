clear all
close all
clc
NN=4:9;
hh=[1 0.5 0.25 0.125 0.0625 0.03125];

err=zeros(length(hh),length(NN));
t=0.1;
i=1;
for N=NN
    for h=hh
        load(append('test',num2str(i),'.mat'));
        x=points(:,1);
        y=points(:,2);


        hf=cos(x-t).*cos(y-t)+2;
     
       err(i)= max(abs(hf(:)-Density(:)));
       i=i+1;
    end
end

% err(2:3,:)=[];

for i=1:6
    semilogy(NN,err(i,:),'-*','LineWidth',2);
    hold on
end
xlabel('$N$', 'Interpreter', 'Latex', 'FontSize', 15)
ylabel('$L_\infty$', 'Interpreter', 'Latex','FontSize', 15)
legend('$h=1$','$h=0.5$','$h=0.25$','$h=0.125$','$h=0.0625$','$h=0.03125$','Interpreter', 'Latex','FontSize', 15,'Location','eastoutside')
legend('boxoff')  
box off
axis square
grid on
% hh(2:3)=[];
figure
for i=1:6
    loglog(hh,err(:,i));
    hold on
end
% loglog(hh,1e-4*hh.^(5.5),'--')
% err2=[NN;err']
% fprintf(fopen('errav.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', err2);

% 
% err3=[hh;err]
% fprintf(fopen('erravh.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', err3);


% log10(err(3,2:end)./err(3,1:end-1))./log10(hh(2:end)./hh(1:end-1))

log10(err(:,end)./err(:,1))./log10(hh(end)./hh(1))

%%
clc

load('errtest.dat')
load('errcurv.dat')
load('errav.dat')
load('errv.dat')

load('errtesth.dat')
load('errcurvh.dat')
load('erravh.dat')
load('errvh.dat')

errtest(:,2:end)=errtest(:,2:end)';
errcurv(:,2:end)=errcurv(:,2:end)';
errav(:,2:end)=errav(:,2:end)';
errv(:,2:end)=errv(:,2:end)';

errtesth(:,2:end)=errtesth(:,2:end)';
errcurvh(:,2:end)=errcurvh(:,2:end)';
erravh(:,2:end)=erravh(:,2:end)';
errvh(:,2:end)=errvh(:,2:end)';


% log(errtest(1:3,end)./errtest(1:3,3))./log(hh(end)./hh(2))
% log(errtest(4,end-1)./errtest(4,3))./log(hh(end-1)./hh(2))
% log(errtest(5,end-2)./errtest(5,3))./log(hh(end-2)./hh(2))
% log(errtest(6,end-3)./errtest(6,3))./log(hh(end-3)./hh(2))

% log(errcurv(1:3,end)./errcurv(1:3,3))./log(hh(end)./hh(2))
% log(errcurv(4,end-1)./errcurv(4,3))./log(hh(end-1)./hh(2))
% log(errcurv(5,end-2)./errcurv(5,3))./log(hh(end-2)./hh(2))
% log(errcurv(6,end-3)./errcurv(6,3))./log(hh(end-3)./hh(2))

log(errav(1,end)./errav(1,3))./log(hh(end)./hh(2))
log(errav(2,end-1)./errav(2,3))./log(hh(end-1)./hh(2))
log(errav(3,end-2)./errav(3,3))./log(hh(end-2)./hh(2))
log(errav(4,end-3)./errav(4,3))./log(hh(end-3)./hh(2))

% log(errav(:,end-3)./errav(:,3))./log(hh(end-3)./hh(2))
% 
% log(errv(:,end-3)./errv(:,3))./log(hh(end-3)./hh(2))

log(errv(1:3,end)./errv(1:3,3))./log(hh(end)./hh(2))
log(errv(4,end-1)./errv(4,3))./log(hh(end-1)./hh(2))
log(errv(5,end-2)./errv(5,3))./log(hh(end-2)./hh(2))
log(errv(6,end-3)./errv(6,3))./log(hh(end-3)./hh(2))

% fprintf(fopen('errtest.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', errtest');
% fprintf(fopen('errtesth.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', errtesth');
% 
% fprintf(fopen('errcurv.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', errcurv');
% fprintf(fopen('errcurvh.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', errcurvh');

% fprintf(fopen('errav.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', errav');
% fprintf(fopen('erravh.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', erravh');

% fprintf(fopen('errv.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', errv');
% fprintf(fopen('errvh.dat','w'), '%.15f %.15f %.15f %.15f %.15f %.15f %.15f\n', errvh');

