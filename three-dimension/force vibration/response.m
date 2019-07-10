function [bear1w,bear1v,bear1z,bear2w,bear2v,bear2z,bear3w,bear3v,bear3z]=response(a,attack_ang,onetwo)
w=1:300;   % Ô²ÆµÂÊ
% onetwo=1;
% FF=[0;0;1;0;zeros(48,1)];
% a=0.4;   attack_ang=0.5;
COFFC=zeros(52,length(w));

for ii=1:length(w)
    
    if onetwo==1
        
        
        [matrix,fu,fv]=first_condition_matrix_space(2*pi*w(ii),a,attack_ang);
        FF=[0;0;fu;0;0;-fv;0;0;zeros(42,1)];
        
    else
        [matrix,f3u,f3w]=condition_matrix_space(2*pi*w(ii),a,attack_ang);
        FF=[zeros(20,1);0;0;f3u;0;0;f3w;0;0;zeros(22,1)];
        
        
    end

    COFFC(:,ii)=matrix\FF;
        
end


E=2.1e11*(1+0.01i);  coe=0.9; G=E/2/(1+0.3); density=7500;

L1=0.3;   L2=5;   L3=10;
d1=0.4;  d2=0.36;   d3=0.32;   %¸÷¿ç¶Î³¤¶È¡¢Ö±¾¶   

spring1=1.8e8;  %ÖáÏµÖá³Ð1,2Ö§³Å¸Õ¶È
spring2=1.8e8;

k0u=2.3e8;  k0v=8.8e7;  k0w=8.8e7;

Area2=0.25*pi*d2^2;  
Area3=0.25*pi*d3^2;  
hspring1=spring1/(coe*G*Area2);

hspring2=spring2/(coe*G*Area3);


L01=0:0.1:L1;
L01=L01';

L02=0:0.2:L2;
L02=L02';

L03=0:0.2:L3;
L03=L03';

Area1=0.25*pi*d1^2;  
MI1=pi/64*d1^4;

Area3=0.25*pi*d3^2;  
MI3=pi/64*d3^4;
    
Area2=0.25*pi*d2^2;  
MI2=pi/64*d2^4;

hpros1v=zeros(length(L01),length(w));
hpros1w=zeros(length(L01),length(w));
zpros1=zeros(length(L01),length(w));

hpros2v=zeros(length(L02),length(w));
hpros2w=zeros(length(L02),length(w));
zpros2=zeros(length(L02),length(w));

hpros3v=zeros(length(L03),length(w));
hpros3w=zeros(length(L03),length(w));
zpros3=zeros(length(L03),length(w));

% ½°Ò¶1
% lpro1=1.6*0.5; b=0.03;h=0.3;
% l1=0:0.1:lpro1;
% l1=l1';
% hpro1u=zeros(length(l1),length(w));
% hpro1v=zeros(length(l1),length(w));
% zpro1=zeros(length(l1),length(w));

for nn=1:length(w)

coffC0v=COFFC(1:4,nn);
coffC0w=COFFC(5:8,nn);
coffD0=COFFC(9:10,nn);
% coffH0=COFFC(11:12,nn);


[hengmatrix,zongmatrix,~,heng2,zong2,~]=transfer_matrix(hspring1,hspring2,L1,L2,d1,d2,d3,2*pi*w(nn));

coffC02v=heng2*coffC0v;
coffC02w=heng2*coffC0w;
coffD02=zong2*coffD0;


coffC03v=hengmatrix*coffC0v;
coffC03w=hengmatrix*coffC0w;
coffD03=zongmatrix*coffD0;

%µÚÒ»¿ç¶Î
[~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area1,coe,G,MI1,2*pi*w(nn),density);
hpros1v(:,nn)=[cosh(lambda*L01),sinh(lambda*L01),cos(hlambda*L01),sin(hlambda*L01)]*coffC0v;
hpros1w(:,nn)=[cosh(lambda*L01),sinh(lambda*L01),cos(hlambda*L01),sin(hlambda*L01)]*coffC0w;

zpros1(:,nn)=-[cos(mue*L01),sin(mue*L01)]*coffD02;

%µÚ¶þ¿ç¶Î
[~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area2,coe,G,MI2,2*pi*w(nn),density);
hpros2v(:,nn)=[cosh(lambda*L02),sinh(lambda*L02),cos(hlambda*L02),sin(hlambda*L02)]*coffC02v;
hpros2w(:,nn)=[cosh(lambda*L02),sinh(lambda*L02),cos(hlambda*L02),sin(hlambda*L02)]*coffC02w;

zpros2(:,nn)=-[cos(mue*L02),sin(mue*L02)]*coffD02;


%µÚÈý¿ç¶Î
[~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area3,coe,G,MI3,2*pi*w(nn),density);
hpros3v(:,nn)=[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)]*coffC03v;
hpros3w(:,nn)=[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)]*coffC03w;

zpros3(:,nn)=-[cos(mue*L03),sin(mue*L03)]*coffD03;

end

%½°Ò¶
% density=9900;
% E=9.2e10*(1+0.001i); 
% coe=5/6; G=E/2/(1+0.32); 
% for nn=1:length(w)
% 
% 
% %½°Ò¶
% CoffC1u=COFFC(13:16,nn);
% CoffC1v=COFFC(17:20,nn);
% CoffD1=COFFC(21:22,nn);
% 
% Area=b*h;  
% MI=b^3*h/12;
% [~,lambda,~,hlambda,~,~,~,~]=solution_par(E,Area,coe,G,MI,2*pi*w(nn),density);
% 
% hpro1u(:,nn)=[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*CoffC1u;
% 
% vhpro1u(:,nn)=[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*CoffC1u*2*pi*w(nn)*i;
% 
% MI=b*h^3/12;
% [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*w(nn),density);
% 
% hpro1v(:,nn)=[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*CoffC1v;
% 
% zpro1(:,nn)=[cos(mue*l1),sin(mue*l1)]*CoffD1; 
% 
% 
% end


% figure(1)
% semilogy(w,abs(hpros1w(end,:)),'--r','LineWidth',2)
% hold on
% semilogy(w,abs(hpros1v(end,:)),'r','LineWidth',2)
% semilogy(w,abs(zpros1(1,:)),'g','LineWidth',2)
% hold on
bear1w=hpros1w(end,:)*spring1;
bear1v=hpros1v(end,:)*spring1;
bear1z=zpros1(end,:);


% figure(2)
% semilogy(w,abs(hpros2w(end,:)),'--r','LineWidth',2)
% hold on
% semilogy(w,abs(hpros2v(end,:)),'r','LineWidth',2)

% semilogy(w,abs(zpros2(end,:)),'r','LineWidth',2)

bear2w=hpros2w(end,:)*spring2;
bear2v=hpros2v(end,:)*spring2;
bear2z=zpros2(end,:);


% figure(3)
% semilogy(w,abs(hpros3w(end,:)),'--k','LineWidth',2)
% hold on
% semilogy(w,abs(hpros3v(end,:)),'r','LineWidth',2)
% semilogy(w,abs(zpros3(end,:)),'--r','LineWidth',2)
% hold on
% set(gca,'XTick',0:10:300); 

bear3w=hpros3w(end,:)*k0w;
bear3v=hpros3v(end,:)*k0v;
bear3z=zpros3(end,:)*k0u;
end


% figure(4)
% semilogy(w,abs(vhpro1u(1,:)),'--g','LineWidth',2);
% hold on



    