clear all;
close all;
clc

tic;

alpha=0.3;
beta=0.6;
B=(alpha*beta)^(-1);
delta=1;
maxit=1000;
kss=(alpha*beta*B)^(1/(1-alpha));
nk=501;
tol=10^-6;


% Grid Generation
gridk=zeros(1,nk);
for i=1:1:nk
    gridk(i)=(1+0.25/250*(i-251))*kss;
end


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

k0=0.9*kss;
mink=(1+0.25/250*(1-251))*kss;
maxk=(1+0.25/250*(501-251))*kss;
ctcolumn=zeros(101,1);
ktcolumn=zeros(102,1);
ktcolumn(1,1)=k0;
for i=1:101
    ktcolumn(i+1)=gfun(it,int32((ktcolumn(i)-mink)/(0.25/250*kss)+1));
    ctcolumn(i)=B*ktcolumn(i)^alpha+(1-delta)*ktcolumn(i)-ktcolumn(i+1);
end
optimalvalue=zeros(101,2);
optimalvalue(:,1)=ctcolumn;
optimalvalue(:,2)=ktcolumn(2:102);

t=linspace(0,10,11);
figure;
plot(t,optimalvalue(1:11,1),'-o');
xlabel('Time');
ylabel('c(t)')
title('c(t) over time')
grid on;

figure;
plot(t,optimalvalue(1:11,2),'-o');
xlabel('Time');
ylabel('k(t)')
title('k(t) over time')
grid on;



toc
