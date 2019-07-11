w=1:300;   % 频率Hz
onetwo=1;
a=0.4;  %  加载点

parameters2;

COFFC=zeros(60,length(w));

for ii=1:length(w)
    
    [matrix,ffu,ff2]=condition_matrix_space(2*pi*w(ii),a,onetwo);
    
    if onetwo==1
        FF=[0;0;ffu;0;0;ff2;0;0;zeros(50,1)];
    else
        FF=[zeros(24,1);0;0;ffu;0;0;ff2;0;0;zeros(26,1)];
    end
    COFFC(:,ii)=matrix\FF;
        
end

l01=0:0.01:L01;
l01=l01';

l02=0:0.01:L02;
l02=l02';

l03=0:0.01:L03;
l03=l03';

l1=0:0.01:L1;
l1=l1';


Area1=0.25*pi*d1^2;  
MI1=pi/64*d1^4;
    
Area2=0.25*pi*d2^2;  
MI2=pi/64*d2^4;

Area3=0.25*pi*d3^2;  
MI3=pi/64*d3^4;

hpro1u=zeros(length(l1),length(w));

hpros1v=zeros(length(l01),length(w));
hpros1w=zeros(length(l01),length(w));
zpros1=zeros(length(l01),length(w));

hpros2v=zeros(length(l02),length(w));
hpros2w=zeros(length(l02),length(w));
zpros2=zeros(length(l02),length(w));

hpros3v=zeros(length(l03),length(w));
hpros3w=zeros(length(l03),length(w));
zpros3=zeros(length(l03),length(w));


for nn=1:length(w)

coffC0v=COFFC(1:4,nn);
coffC0w=COFFC(5:8,nn);
coffD0=COFFC(9:10,nn);

coffC1u=COFFC(13:16,nn);

[hengmatrix,zongmatrix,~,heng2,zong2,~]=transfer_matrix(spring1,spring2,L01,L02,d1,d2,d3,2*pi*w(nn),density);
% [hengmatrixv,zongmatrix,~,heng2v,zong2,~]=transfer_matrix(spring1v,spring2v,L01,L02,d1,d2,d3,2*w(nn),density);

% [hengmatrixw,~,~,heng2w,~,~]=transfer_matrix(spring1w,spring2w,L01,L02,d1,d2,d3,2*w(nn),density);


coffC02v=heng2*coffC0v;
coffC02w=heng2*coffC0w;
coffD02=zong2*coffD0;

coffC03v=hengmatrix*coffC0v;
coffC03w=hengmatrix*coffC0w;
coffD03=zongmatrix*coffD0;

%第一跨段
[~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area1,coe,G,MI1,2*pi*w(nn),density);
hpros1v(:,nn)=[cosh(lambda*l01),sinh(lambda*l01),cos(hlambda*l01),sin(hlambda*l01)]*coffC0v;
hpros1w(:,nn)=[cosh(lambda*l01),sinh(lambda*l01),cos(hlambda*l01),sin(hlambda*l01)]*coffC0w;
zpros1(:,nn)=-[cos(mue*l01),sin(mue*l01)]*coffD02;

%第二跨段
[~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area2,coe,G,MI2,2*pi*w(nn),density);
hpros2v(:,nn)=[cosh(lambda*l02),sinh(lambda*l02),cos(hlambda*l02),sin(hlambda*l02)]*coffC02v;
hpros2w(:,nn)=[cosh(lambda*l02),sinh(lambda*l02),cos(hlambda*l02),sin(hlambda*l02)]*coffC02w;
zpros2(:,nn)=-[cos(mue*l02),sin(mue*l02)]*coffD02;


%第三跨段
[~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area3,coe,G,MI3,2*pi*w(nn),density);
hpros3v(:,nn)=[cosh(lambda*l03),sinh(lambda*l03),cos(hlambda*l03),sin(hlambda*l03)]*coffC03v;
hpros3w(:,nn)=[cosh(lambda*l03),sinh(lambda*l03),cos(hlambda*l03),sin(hlambda*l03)]*coffC03w;
zpros3(:,nn)=-[cos(mue*l03),sin(mue*l03)]*coffD03;


% [~,lambda,~,hlambda,~,~,~,~]=solution_par(pro_E,b*h,5/6,pro_E/2/(1+ue),h*b^3/12,2*pi*w(nn),pro_density);
% 
% hpro1u(:,nn)=2*pi*w(nn)*i*[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*coffC1u;

end

% pro1=hpro1u(1,:);
% 
% figure(1)
% semilogy(w,abs(pro1),'--r','LineWidth',2)

figure(1)
semilogy(w,abs(hpros1w(end,:)),'--r','LineWidth',2)
% hold on
% semilogy(w,abs(hpros1v(end,:)),'r','LineWidth',2)
% semilogy(w,abs(zpros1(1,:)),'g','LineWidth',2)
% hold on


figure(2)
semilogy(w,abs(hpros2w(end,:)),'--r','LineWidth',2)
% hold on
% semilogy(w,abs(hpros2v(end,:)),'r','LineWidth',2)
% semilogy(w,abs(zpros2(end,:)),'r','LineWidth',2)


figure(3)
semilogy(w,abs(hpros3w(end,:)),'--k','LineWidth',2)
hold on
% semilogy(w,abs(hpros3v(end,:)),'r','LineWidth',2)
semilogy(w,abs(zpros3(end,:)),'--r','LineWidth',2)
% hold on
% set(gca,'XTick',0:10:300); 

