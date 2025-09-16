clear all;
close all;
clc

tic;

alpha=0.3;
beta=0.8;
rho=1/0.8-1;
delta=0.2;
B=(rho+delta)/alpha;
maxit=1000;
kss=(1/(alpha*beta*B)+(delta-1)/(alpha*B))^(1/(alpha-1));
nk=501;
tol=10^-6;


% Grid Generation
gridk=zeros(1,nk);
for i=1:1:nk
    gridk(i)=(1+0.25/250*(i-251))*kss;
end


% Numerical Solution

% Initial Value Function

vfun(1,1:nk)=0.0;

% Value function iteration

it=1;
diff=10^10;

while it<=maxit && diff>tol
   it=it+1;
    for kc=1:nk
   	vfun(it,kc)=-10^10;
      for kcc=1:nk
         vhelp=-10^10;
         cons=B*gridk(kc)^alpha+(1-delta)*gridk(kc)-gridk(kcc);
         if cons<=0.0
            vhelp=-10^12;
         else
            vhelp=log(cons)+beta*vfun(it-1,kcc);
         end
         if vhelp>=vfun(it,kc)
            vfun(it,kc)=vhelp;
            gfun(it,kc)=gridk(kcc);
         end
      end
    end
   diff=max(abs(vfun(it,:)-vfun(it-1,:)));
end

% Plot some stuff

pl=plot(gridk,[ vfun(1,:)' vfun(2,:)' vfun(3,:)' vfun(11,:)' vfun(it,:)']);
xla=xlabel('Capital Stock k Today');
yla=ylabel('Value Function');
tit=title('Value Function: True and Approximated');
le=legend('v0','v1','v2','v10','Converged v');
axis([gridk(1) gridk(nk) 0 2 ]);
%gtext('V_{0}');
%gtext('V_{1}');
%gtext('V_{2}');
%gtext('V_{10}');
%gtext('True Value Function');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'Fontsize',12,'Fontweight','bold');

%print -dpsc2 valuepic.eps;
%print -dpsc2 valuepic.ps;

figure;

pl=plot(gridk,[ gfun(2,:)' gfun(3,:)' gfun(11,:)' gfun(it,:)']);
xla=xlabel('Capital Stock k Today');
yla=ylabel('Policy Function');
tit=title('Policy Function: True and Approximated');
le=legend('g1','g2','g10','Converged g');
axis([gridk(1) gridk(nk) gridk(1) gridk(nk)]);
%gtext('g_{1}');
%gtext('g_{2}');
%gtext('g_{10}');
%gtext('True Policy Function');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'Fontsize',12,'Fontweight','bold');


k0=0.9*kss;
mink=(1+0.25/250*(1-251))*kss;
maxk=(1+0.25/250*(501-251))*kss;
ctcolumn=zeros(101,1);
ktcolumn=zeros(102,1);
ktcolumn(1)=k0;
for i=1:101
    ktcolumn(i+1)=gfun(it,int32((ktcolumn(i)-mink)/(0.25/250*kss)+1));
    ctcolumn(i)=B*ktcolumn(i)^alpha+(1-delta)*ktcolumn(i)-ktcolumn(i+1);
end
optimalvalue=zeros(101,2);
optimalvalue(:,1)=ctcolumn;
optimalvalue(:,2)=ktcolumn(2:102);

t=linspace(0,20,21);
figure;
plot(t,optimalvalue(1:21,1),'-o');
xlabel('Time');
ylabel('c(t)')
title('c(t) over time')
grid on;

figure;
plot(t,optimalvalue(1:21,2),'-o');
xlabel('Time');
ylabel('k(t)')
title('k(t) over time')
grid on;



toc
