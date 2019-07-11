%参数初始化
parameters2;

Area2=0.25*pi*d2^2;
Area3=0.25*pi*d3^2;

MI2=pi/64*d2^4;
MI3=pi/64*d3^4;

pro_Aera=b*h;

for nn=4:4
    
    [hengmatrix,zongmatrix,niumatrix,heng2,zong2,niu2]=transfer_matrix(spring1,spring2,L01,L02,d1,d2,d3,2*pi*result_fre(nn),density);
    matrix=condition_matrix_space(2*pi*result_fre(nn));
    [~,ss,vv]=svd(matrix);
    v=vv(:,60);
    
    CoffC01v=v(1:4);
    CoffC01w=v(5:8);
    CoffD01=v(9:10);
    CoffH01=v(11:12);
    
    CoffC1u=v(13:16);
    CoffC1v=v(17:20);
    CoffD1=v(21:22);
    
    CoffC2u=v(25:28);
    CoffC2v=v(29:32);
    CoffD2=v(33:34);
    
    CoffC3u=v(37:40);
    CoffC3w=v(41:44);
    CoffD3=v(45:46);
    
    CoffC4u=v(49:52);
    CoffC4w=v(53:56);
    CoffD4=v(57:58);
    
    CoffC02v=heng2*CoffC01v;
    CoffC02w=heng2*CoffC01w;
    CoffD02=zong2*CoffD01;

    
    CoffC03v=hengmatrix*CoffC01v;
    CoffC03w=hengmatrix*CoffC01w;
    CoffD03=zongmatrix*CoffD01;

    
    hold on

%     xmin=-0.18;xmax=0.18; ymin=-0.18;ymax=2;zmin=-0.18;zmax=0.18;
% 	axis([xmin xmax ymin ymax zmin zmax]);
    
    %主轴第一跨段
    Area=0.25*pi*d1^2;
    MI=pi/64*d1^4;
    
    l01=0:0.01:L01;
    l01=l01';
    
    x01=zeros(length(l01),1);
    y01=l01;
    z01=zeros(length(l01),1);
    
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),density);
    
    hpros1v=[cosh(lambda*l01),sinh(lambda*l01),cos(hlambda*l01),sin(hlambda*l01)]*CoffC01v;
    hpros1w=-[cosh(lambda*l01),sinh(lambda*l01),cos(hlambda*l01),sin(hlambda*l01)]*CoffC01w;
    zpros1=[cos(mue*l01),sin(mue*l01)]*CoffD01;
%     plot3(x01,y01,z01,'--b',x01+hpros1w,y01+zpros1,z01+hpros1v,'r','LineWidth',2);
    
    e = [(1:length(x01)-1)' (2:length(x01))'];
    c01 = [x01+hpros1v y01+zpros1 z01+hpros1w];
    fig01 = patch('Faces',e,'Vertices',c01, 'FaceColor','w','LineWidth',4);
    set(fig01,'EdgeColor','interp','FaceVertexCData',abs(hpros1v),'CDataMapping','scaled');
    plot3(x01,y01,z01,'--b','LineWidth',2);
    
    %第二跨段

    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area2,coe,G,MI2,2*pi*result_fre(nn),density);
    
    l02=0:0.01:L02;
    l02=l02';
    
    x02=zeros(length(l02),1);
    z02=zeros(length(l02),1);
    
    hpros2v=[cosh(lambda*l02),sinh(lambda*l02),cos(hlambda*l02),sin(hlambda*l02)]*CoffC02v;
    hpros2w=-[cosh(lambda*l02),sinh(lambda*l02),cos(hlambda*l02),sin(hlambda*l02)]*CoffC02w;
    
    zpros2=[cos(mue*l02),sin(mue*l02)]*CoffD02;
    
%     plot3(x02,L01+l02,z02,'--b',x02+hpros2w,L01+l02+zpros2,z02+hpros2v,'r','LineWidth',2)
    
    e = [(1:length(x02)-1)' (2:length(x02))'];
    c02 = [x02+hpros2v L01+l02+zpros2 z02+hpros2w];
    fig02 = patch('Faces',e,'Vertices',c02, 'FaceColor','w','LineWidth',4);
    set(fig02,'EdgeColor','interp','FaceVertexCData',abs(hpros2v),'CDataMapping','scaled');
    plot3(x02,L01+l02,z02,'--b','LineWidth',2)
    
    %第三跨段
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area3,coe,G,MI3,2*pi*result_fre(nn),density);
    

    l03=0:0.01:L03;
    l03=l03';
    x03=zeros(length(l03),1);
    z03=zeros(length(l03),1);
    
    hpros3v=[cosh(lambda*l03),sinh(lambda*l03),cos(hlambda*l03),sin(hlambda*l03)]*CoffC03v;
    hpros3w=-[cosh(lambda*l03),sinh(lambda*l03),cos(hlambda*l03),sin(hlambda*l03)]*CoffC03w;
    
    zpros3=[cos(mue*l03),sin(mue*l03)]*CoffD03;
    
%     plot3(x03,L01+L02+l03,z03,'--b',x03+hpros3w,L01+L02+l03+zpros3,z03+hpros3v,'r','LineWidth',2)
    e = [(1:length(x03)-1)' (2:length(x03))'];
    c03 = [x03+hpros3v L01+L02+l03+zpros3 z03+hpros3w];
    fig03 = patch('Faces',e,'Vertices',c03, 'FaceColor','w','LineWidth',4);
    set(fig03,'EdgeColor','interp','FaceVertexCData',abs(hpros3v),'CDataMapping','scaled');
    plot3(x03,L01+L02+l03,z03,'--b','LineWidth',2)
     
    %桨叶
    coe=5/6; G=pro_E/2/(1+ue);
    
    %桨叶1
    l1=0:0.01:L1;
    l1=l1';
    x1=-flipud(l1);
    y1=zeros(length(l1),1);
    z1=zeros(length(l1),1);
    

    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,~,~]=solution_par(pro_E,pro_Aera,coe,G,MI,2*pi*result_fre(nn),pro_density);
    hpro1u=-[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*CoffC1u;
    
    MI=b*h^3/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(pro_E,pro_Aera,coe,G,MI,2*pi*result_fre(nn),pro_density);
    hpro1v=[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*CoffC1v;
    
    zpro1=[cos(mue*l1),sin(mue*l1)]*CoffD1;
    
%     plot3(x1,y1,z1,'--b',x1+zpro1,y1+hpro1u,z1+hpro1v,'r','LineWidth',2)

   
    plot3(x1,y1,z1,'--b','LineWidth',2)    

    
    % 桨叶2
    l2=0:0.01:L2;
    l2=l2';
    x2=l2;
    y2=zeros(length(l2),1);
    z2=zeros(length(l2),1);
    
    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,~,~]=solution_par(pro_E,pro_Aera,coe,G,MI,2*pi*result_fre(nn),pro_density);
    hpro2u=[cosh(lambda*l2),sinh(lambda*l2),cos(hlambda*l2),sin(hlambda*l2)]*CoffC2u;
    
    MI=b*h^3/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(pro_E,pro_Aera,coe,G,MI,2*pi*result_fre(nn),pro_density);
    hpro2v=[cosh(lambda*l2),sinh(lambda*l2),cos(hlambda*l2),sin(hlambda*l2)]*CoffC2v;
    
    zpro2=[cos(mue*l2),sin(mue*l2)]*CoffD2;
    
%     plot3(x2,y2,z2,'--b',x2+zpro2,y2+hpro2u,z2+hpro2v,'r','LineWidth',2)
   
    plot3(x2,y2,z2,'--b','LineWidth',2)
    
    % 桨叶3
    l3=0:0.01:L3;
    l3=l3';
    z3=flipud(l3);
    x3=zeros(length(l3),1);
    y3=zeros(length(l3),1);
    
    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,~,~]=solution_par(pro_E,pro_Aera,coe,G,MI,2*pi*result_fre(nn),pro_density);
    hpro3u=[cosh(lambda*l3),sinh(lambda*l3),cos(hlambda*l3),sin(hlambda*l3)]*CoffC3u;
    
    MI=b*h^3/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(pro_E,pro_Aera,coe,G,MI,2*pi*result_fre(nn),pro_density);
    hpro3w=[cosh(lambda*l3),sinh(lambda*l3),cos(hlambda*l3),sin(hlambda*l3)]*CoffC3w;
    
    zpro3=[cos(mue*l3),sin(mue*l3)]*CoffD3;
%     plot3(x3,y3,z3,'--b',x3+hpro3w,y3+hpro3u,z3+zpro3,'r','LineWidth',2)
  
    plot3(x3,y3,z3,'--b','LineWidth',2)
    
    
    % 桨叶4
    l4=0:0.01:L4;
    l4=l4';
    z4=-l4;
    x4=zeros(length(l4),1);
    y4=zeros(length(l4),1);
    
    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,~,~]=solution_par(pro_E,pro_Aera,coe,G,MI,2*pi*result_fre(nn),pro_density);
    hpro4u=[cosh(lambda*l4),sinh(lambda*l4),cos(hlambda*l4),sin(hlambda*l4)]*CoffC4u;
    
    MI=b*h^3/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(pro_E,pro_Aera,coe,G,MI,2*pi*result_fre(nn),pro_density);
    hpro4w=[cosh(lambda*l4),sinh(lambda*l4),cos(hlambda*l4),sin(hlambda*l4)]*CoffC4w;
    
    zpro4=[cos(mue*l4),sin(mue*l4)]*CoffD4;

    
%     hpro1u=flipud(hpro2u);
     e = [(1:length(x1)-1)' (2:length(x1))'];
    c12 = [x1+zpro1 y1+hpro3u z1+0*hpro1v];
    fig1 = patch('Faces',e,'Vertices',c12, 'FaceColor','w','LineWidth',4);
    set(fig1,'EdgeColor','interp','FaceVertexCData',abs(hpro3u),'CDataMapping','scaled');
    
    
    
     c22 = [x2+zpro2 y2+hpro4u z2+0*hpro2v];
    fig2 = patch('Faces',e,'Vertices',c22, 'FaceColor','w','LineWidth',4);
    set(fig2,'EdgeColor','interp','FaceVertexCData',abs(hpro4u),'CDataMapping','scaled');
    
    
      c32 = [x3+hpro1v y3+hpro1u z3+0*zpro3];
    fig3 = patch('Faces',e,'Vertices',c32, 'FaceColor','w','LineWidth',4);
    set(fig3,'EdgeColor','interp','FaceVertexCData',abs(hpro1u),'CDataMapping','scaled');
    
    
    %     plot3(x4,y4,z4,'--b',x4+hpro4w,y4+hpro4u,z4+zpro4,'r','LineWidth',2)
    c42 = [x4+hpro2v y4+hpro2u z4+0*zpro4];
    fig4 = patch('Faces',e,'Vertices',c42, 'FaceColor','w','LineWidth',4);
    set(fig4,'EdgeColor','interp','FaceVertexCData',abs(hpro2u),'CDataMapping','scaled');
    plot3(x4,y4,z4,'--b','LineWidth',2)

%     xlabel('横向模态位移/m')
%     ylabel('纵向模态位移/m')
%     zlabel('垂向模态位移/m')
%     legend('固有结构','模态振型');
    text(1,1,1,[num2str(result_fre(nn)),'Hz'],'FontSize',15);
    view(3)
    
end