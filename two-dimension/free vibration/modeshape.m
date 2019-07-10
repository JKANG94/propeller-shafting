% clear,clc
% load('result.mat')

for nn=1:6

    subplot(2,3,nn)
    
    %主轴参数
    E=2.1e11; 
    
    coe=9/10; G=E/2/(1+0.3); 
    density=7500;
    
    L1=0.48;   L2=1.8;   L3=1.8;
    d1=0.16;  d2=0.14;   d3=0.11;   %各跨段长度、直径   
    
    spring1=1e8;  %轴系轴承1,2支撑刚度
    spring2=1e8;
    
    Area2=0.25*pi*d2^2;  
    Area3=0.25*pi*d3^2;  
    hspring1=spring1/(coe*G*Area2);

    hspring2=spring2/(coe*G*Area3);
   
    %桨叶参数
    lpro1=0.62;b=0.04; h=0.2;
    lpro2=0.62;
    
    [hengmatrix,zongmatrix,~,heng2,zong2,~]=transfer_matrix(hspring1,hspring2,L1,L2,d1,d2,d3,2*pi*result_fre(nn),density);
    
    matrix=xiu_matrix_planer(2*pi*result_fre(nn));
    
    [~,ss,vv]=svd(matrix);
    v=vv(:,18);
    
    CoffC1=v(1:4);
    CoffD1=v(5:6);
    
    CoffC2=v(7:10);
    CoffD2=v(11:12);
    
    CoffC01=v(13:16);
    CoffD01=v(17:18);
    
    CoffC02=heng2*CoffC01;
    CoffD02=zong2*CoffD01;
    
    CoffC03=hengmatrix*CoffC01;
    CoffD03=zongmatrix*CoffD01;
    
    hold on
    axis equal
    
    %主轴
    % 第一跨段
    L01=0:0.02:L1;
    L01=L01';
    y01=zeros(length(L01),1);

    Area=0.25*pi*d1^2;  
    MI=pi/64*d1^4;
    
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),density);
    
    hpros1=[cosh(lambda*L01),sinh(lambda*L01),cos(hlambda*L01),sin(hlambda*L01)]*CoffC01;
    
    zpros1=-[cos(mue*L01),sin(mue*L01)]*CoffD01;
 
    e = [(1:length(L01)-1)' (2:length(L01))'];
    c01 = [L01+zpros1 y01+hpros1];
    fig01 = patch('Faces',e,'Vertices',c01, 'FaceColor','w','LineWidth',3);
    set(fig01,'EdgeColor','interp','FaceVertexCData',abs(zpros1),'CDataMapping','scaled');
    
    %第二跨段
    Area=0.25*pi*d2^2;  
    MI=pi/64*d2^4;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),density);
    
    L02=linspace(0,L2,40);
    L02=L02';
    y02=zeros(length(L02),1);
    hpros2=[cosh(lambda*L02),sinh(lambda*L02),cos(hlambda*L02),sin(hlambda*L02)]*CoffC02;
    
    zpros2=-[cos(mue*L02),sin(mue*L02)]*CoffD02;
    
    e = [(1:length(L02)-1)' (2:length(L02))'];
    c02 = [L1+L02+zpros2 y02+hpros2];
    fig02 = patch('Faces',e,'Vertices',c02, 'FaceColor','w','LineWidth',3);
    set(fig02,'EdgeColor','interp','FaceVertexCData',abs(zpros2),'CDataMapping','scaled');
    
    %第三跨段
    Area=0.25*pi*d3^2;  
    MI=pi/64*d3^4;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),density);
    
    L03=linspace(0,L3,40);
    L03=L03';
    y03=zeros(length(L03),1);
    hpros3=[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)]*CoffC03;
    
    zpros3=-[cos(mue*L03),sin(mue*L03)]*CoffD03;
    
    e = [(1:length(L03)-1)' (2:length(L03))'];
    c03 = [L1+L2+L03+zpros3 y03+hpros3];
    fig03 = patch('Faces',e,'Vertices',c03, 'FaceColor','w','LineWidth',3);
    set(fig03,'EdgeColor','interp','FaceVertexCData',abs(zpros3),'CDataMapping','scaled');
    
    
    %桨叶物理参数
    E=9.2e10; 
    
    coe=5/6; G=E/2/(1+0.32); 
    
    %桨叶1
    l1=0:0.02:lpro1;
    l1=l1';
    y1=flipud(l1);
    x1=zeros(length(l1),1);
    
    Area=b*h;  
    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),density);
    
    hpro1=-[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*CoffC1;
    
    zpro1=[cos(mue*l1),sin(mue*l1)]*CoffD1;
    
    e = [(1:length(l1)-1)' (2:length(l1))'];
    c1 = [x1+hpro1 y1+zpro1];
    fig1 = patch('Faces',e,'Vertices',c1, 'FaceColor','w','LineWidth',3);
    set(fig1,'EdgeColor','interp','FaceVertexCData',abs(hpro1),'CDataMapping','scaled');
    
%    桨叶2
    l2=0:0.02:lpro2;
    l2=l2';
%     y2=-flipud(l2);
    y2=-l2;
    x2=zeros(length(l2),1);
    hpro2=-[cosh(lambda*l2),sinh(lambda*l2),cos(hlambda*l2),sin(hlambda*l2)]*CoffC2;
    
    zpro2=[cos(mue*l2),sin(mue*l2)]*CoffD2;
 
    e = [(1:length(l2)-1)' (2:length(l2))'];
    c2 = [x2+hpro2 y2+zpro2];
    fig2 = patch('Faces',e,'Vertices',c2, 'FaceColor','w','LineWidth',3);
    set(fig2,'EdgeColor','interp','FaceVertexCData',abs(hpro2),'CDataMapping','scaled');
    
    
    
%     plot(x1,y1,'--b',x1+hpro1,y1+zpro1,'r','LineWidth',2)
%     plot(x2,y2,'--b',x2+hpro2,y2+zpro2,'r','LineWidth',2)
%     plot(L01,y01,'--b',L01+zpros1,y01+hpros1,'r','LineWidth',2)
%     plot(L02+L1,y02,'--b',L02+L1+zpros2,y02+hpros2,'r','LineWidth',2)
%     plot(L03+L1+L2,y03,'--b',L03+L1+L2+zpros3,y03+hpros3,'r','LineWidth',2)

%     legend('原结构','模态振型');
    text(1.5,0.2,[num2str(result_fre(nn)),'Hz'],'Fontsize',16);
    set(gca,'fontsize',16)
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    view(2)

    
end