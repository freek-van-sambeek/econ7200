clear all;
close all;
clc

tic;

alpha=0.3;
beta=0.6;
B=(alpha*beta)^(-1);
maxit=1000;
mink=0.8;
maxk=1.2;
nk=5;
tol=10^-6;


% Grid Generation

gridk=mink:((maxk-mink)/(nk-1)):maxk;

% Analytical solution

a0=(log(B*(1-alpha*beta))+alpha*beta/(1-alpha*beta)*log(alpha*beta*B))/(1-beta);
a1=alpha/(1.0-alpha*beta);

for kc=1:nk
	valanal(kc)=a0+a1*log(gridk(kc));
    polanal(kc)=alpha*beta*B*gridk(kc)^alpha;
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
         cons=B*gridk(kc)^alpha-gridk(kcc);
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

pl=plot(gridk,[valanal(:) vfun(1,:)' vfun(2,:)' vfun(3,:)' vfun(11,:)' vfun(it,:)']);
xla=xlabel('Capital Stock k Today');
yla=ylabel('Value Function');
tit=title('Value Function: True and Approximated');
le=legend('Analytical','v0','v1','v2','v10','Converged v');
axis([gridk(1) gridk(nk) 1 5 ]);
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

pl=plot(gridk,[polanal(:) gfun(2,:)' gfun(3,:)' gfun(11,:)' gfun(it,:)']);
xla=xlabel('Capital Stock k Today');
yla=ylabel('Policy Function');
tit=title('Policy Function: True and Approximated');
le=legend('Analytical','g1','g2','g10','Converged g');
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




valuediff=valanal-vfun(it,:);

toc;

