function [BEAR1W,BEAR1V,BEAR2W,BEAR2V,BEAR3W,BEAR3V,BEAR3Z]=forcetransmatrix(onetwo)

% warning('off')

ns=1.9;
water_v=4.63;

index=20;  L=1.6;
delta=L/index;
start=delta/2;

w=1:1:300;

BEAR1W=zeros(index,length(w));
BEAR1V=zeros(index,length(w));
BEAR1Z=zeros(index,length(w));

BEAR2W=zeros(index,length(w));
BEAR2V=zeros(index,length(w));
BEAR2Z=zeros(index,length(w));

BEAR3W=zeros(index,length(w));
BEAR3V=zeros(index,length(w));
BEAR3Z=zeros(index,length(w));

for jj=0:index-1
    
    a=L-(start+delta*jj);      % 0.08----1.6
    attack_ang=atan(ns*2*pi*(start+delta*jj)/water_v);

    [bear1w,bear1v,bear1z,bear2w,bear2v,bear2z,bear3w,bear3v,bear3z]=response(a,attack_ang,onetwo);
    
    BEAR1W(jj+1,:)=bear1w;
    BEAR1V(jj+1,:)=bear1v;
    BEAR1Z(jj+1,:)=bear1z;
    
    BEAR2W(jj+1,:)=bear2w;
    BEAR2V(jj+1,:)=bear2v;
    BEAR2Z(jj+1,:)=bear2z;
    
    BEAR3W(jj+1,:)=bear3w;
    BEAR3V(jj+1,:)=bear3v;
    BEAR3Z(jj+1,:)=bear3z;
    
end
end


% figure(1)

% semilogy(w,abs(bear1w),'--b','LineWidth',2)
% hold on
% semilogy(w,abs(bear1v),'--g','LineWidth',2)
% semilogy(w,bear1z,'r','LineWidth',2)

% figure(2)
% 
% semilogy(w,abs(bear2w),'--b','LineWidth',2)
% hold on
% semilogy(w,abs(bear2v),'--g','LineWidth',2)
% semilogy(w,bear2z,'r','LineWidth',2)

% figure(3)
% 
% semilogy(w,abs(bear3w),'--b','LineWidth',2)
% hold on
% semilogy(w,abs(bear3v),'--g','LineWidth',2)
% semilogy(w,abs(bear3z),'r','LineWidth',2)
% set(gca,'XTick',0:10:300); 

