clear all
clc
format long

w=1:1:300;   % 频率hz

COFFC=zeros(18,length(w));


for ii=1:length(w)
    
     [matrix,f]=xiu_matrix_planer(2*pi*w(ii));
%      FF=[0;0;0;0;0;0;0;0;0;0;-f;0;0;0;0;0;0];     %加载在桨叶2
    FF=[0;0;0;0;-f;0;0;0;0;0;0;0;0;0;0;0;0];     %加载在桨叶1
    COFFC(:,ii)=matrix\FF;
        
end


E=2.1e11*(1+0.001i);  coe=0.9; G=E/2/(1+0.3);  density=7500;
L1=0.48; L2=1.8;   L3=1.8;
d1=0.16;  d2=0.14;   d3=0.11;   %各跨段长度、直径   

k0x=7.66e6;  k0y=1e8;   %  主轴右端轴向、径向刚度（推力轴承）

spring1=1e8;  %轴系轴承1,2支撑刚度
spring2=1e8;

L01=0:0.12:L1;
L01=L01';

L02=0:0.3:L2;
L02=L02';

L03=0:0.3:L3;
L03=L03';

Area1=0.25*pi*d1^2;  
MI1=pi/64*d1^4;

Area3=0.25*pi*d3^2;  
MI3=pi/64*d3^4;
    
Area2=0.25*pi*d2^2;  
MI2=pi/64*d2^4;


hspring1=spring1/(coe*G*Area2);

hspring2=spring2/(coe*G*Area3);

hpros1=zeros(length(L01),length(w));
zpros1=zeros(length(L01),length(w));

hpros2=zeros(length(L02),length(w));
zpros2=zeros(length(L02),length(w));

hpros3=zeros(length(L03),length(w));
zpros3=zeros(length(L03),length(w));

for nn=1:length(w)

    CoffC1=COFFC(1:4,nn);
    CoffD1=COFFC(5:6,nn);
       
    CoffC01=COFFC(13:16,nn);
    CoffD01=COFFC(17:18,nn);
    
    
    [hengmatrix,zongmatrix,~,heng2,zong2,~]=transfer_matrix(hspring1,hspring2,L1,L2,d1,d2,d3,2*pi*w(nn),density);
    
    
    CoffC02=heng2*CoffC01;
    CoffD02=zong2*CoffD01;
    
    CoffC03=hengmatrix*CoffC01;
    CoffD03=zongmatrix*CoffD01;

    %第一跨段
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area1,coe,G,MI1,2*pi*w(nn),density);

    hpros1(:,nn)=-[cosh(lambda*L01),sinh(lambda*L01),cos(hlambda*L01),sin(hlambda*L01)]* CoffC01;

    zpros1(:,nn)=-[cos(mue*L01),sin(mue*L01)]* CoffD01;
    
    
     %第二跨段
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area2,coe,G,MI2,2*pi*w(nn),density);
    hpros2(:,nn)=-[cosh(lambda*L02),sinh(lambda*L02),cos(hlambda*L02),sin(hlambda*L02)]* CoffC02;

    zpros2(:,nn)=-[cos(mue*L02),sin(mue*L02)]* CoffD02;
    
    %第三跨段

    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area3,coe,G,MI3,2*pi*w(nn),density);
    hpros3(:,nn)=-[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)]* CoffC03;

    zpros3(:,nn)=-[cos(mue*L03),sin(mue*L03)]* CoffD03;

end


%第一跨段
figure(1)
semilogy(w,abs(hpros1(5,:)*spring1),'--b','LineWidth',2)
% hold on
% semilogy(w,abs(zpros1(5,:)),'r','LineWidth',2)
set(gca,'XTick',0:20:300); 
set(gca,'fontsize',16)
% legend('横向','纵向');
xlabel('频率/Hz')
ylabel('幅值/N')
legend('径向轴承反力');


%第二跨段
figure(2)
semilogy(w,abs(hpros2(7,:)*spring2),'--b','LineWidth',2)
% hold on
% semilogy(w,abs(zpros2(7,:)),'r','LineWidth',2)
set(gca,'XTick',0:20:300); 
set(gca,'fontsize',16)
xlabel('频率/Hz')
ylabel('幅值/N')

%第三跨段
figure(3)
semilogy(w,abs(hpros3(7,:)*k0y),'--b','LineWidth',2)
hold on
semilogy(w,abs(zpros3(7,:)*k0x),'r','LineWidth',2)
set(gca,'XTick',0:20:300); 
set(gca,'fontsize',16)
legend('基座反力','推力轴承反力');
xlabel('频率/Hz')
ylabel('幅值/N')

% bear1=abs(hpros1(5,:)*spring1)';
% 
% bear2=abs(hpros2(5,:)*spring2)';
% 
% tuibearh=abs(hpros3(7,:)*k0y)';
% 
% tuibearz=abs(zpros3(7,:)*k0x)';
% save('bear1.mat','bear1')
% save('bear2.mat','bear2')
% save('tuibearh.mat','tuibearh')
% save('tuibearz.mat','tuibearz')
% save('w.mat','w')


    