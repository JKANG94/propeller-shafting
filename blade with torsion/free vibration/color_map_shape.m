% 轴系参数
coe=0.9;    E=2.1e11; G=E/2/(1+0.3); shaft_density=7500;

L1=0.1;   L2=0.6;   L3=1.3;
d1=0.017;  d2=0.017;   d3=0.017;     %各跨段长度、直径

%轴系轴承1,2支撑刚度
spring1=1.8e8;
spring2=1.8e8;

%桨叶参数
lpro1=0.18; lpro2=0.18;lpro3=0.18;lpro4=0.18;
b=0.01; h=0.06;
pro_density=7500;

%

Area2=0.25*pi*d2^2;
Area3=0.25*pi*d3^2;
hspring1=spring1/(coe*G*Area2);
hspring2=spring2/(coe*G*Area3);

for nn=2:2
    
    [hengmatrix,zongmatrix,niumatrix,heng2,zong2,niu2]=transfer_matrix(hspring1,hspring2,L1,L2,d1,d2,d3,2*pi*result_fre(nn));
    
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
    CoffH02=niu2*CoffH01;
    
    CoffC03v=hengmatrix*CoffC01v;
    CoffC03w=hengmatrix*CoffC01w;
    CoffD03=zongmatrix*CoffD01;
    CoffH03=niumatrix*CoffH01;
    
    
    hold on
    
        xmin=-0.18;xmax=0.18; ymin=-0.18;ymax=2;zmin=-0.18;zmax=0.18;
        axis([xmin xmax ymin ymax zmin zmax]);
    
    % 主轴
    
    %  第一跨段
    Area=0.25*pi*d1^2;
    
    MI=pi/64*d1^4;
    
    l01=0:0.01:L1;
    l01=l01';
    
    x01=zeros(length(l01),1);
    y01=l01;
    z01=zeros(length(l01),1);
    
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),shaft_density);
    
    hpros1v=[cosh(lambda*l01),sinh(lambda*l01),cos(hlambda*l01),sin(hlambda*l01)]*CoffC01v;
    
    hpros1w=[cosh(lambda*l01),sinh(lambda*l01),cos(hlambda*l01),sin(hlambda*l01)]*CoffC01w;
    
    zpros1=-[cos(mue*l01),sin(mue*l01)]*CoffD01;
    
    e = [(1:length(x01)-1)' (2:length(x01))'];
    c01 = [x01+hpros1w y01+zpros1 z01+hpros1v];
    fig01 = patch('Faces',e,'Vertices',c01, 'FaceColor','w','LineWidth',3);
    set(fig01,'EdgeColor','interp','FaceVertexCData',abs(hpros1w),'CDataMapping','scaled');
    plot3(x01,y01,z01,'--b','LineWidth',2);
    
    %     第二跨段
    Area=0.25*pi*d2^2;
    MI=pi/64*d2^4;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),shaft_density);
   
    
    L02=0:0.01:L2;
    L02=L02';
    
    x02=zeros(length(L02),1);
    z02=zeros(length(L02),1);
    
    hpros2v=[cosh(lambda*L02),sinh(lambda*L02),cos(hlambda*L02),sin(hlambda*L02)]*CoffC02v;
    hpros2w=[cosh(lambda*L02),sinh(lambda*L02),cos(hlambda*L02),sin(hlambda*L02)]*CoffC02w;
    
    zpros2=-[cos(mue*L02),sin(mue*L02)]*CoffD02;
    
    
    e = [(1:length(x02)-1)' (2:length(x02))'];
    c02 = [x02+hpros2w L1+L02+zpros2 z02+hpros2v];
    fig02 = patch('Faces',e,'Vertices',c02, 'FaceColor','w','LineWidth',3);
    set(fig02,'EdgeColor','interp','FaceVertexCData',abs(hpros2w),'CDataMapping','scaled');
    plot3(x02,L1+L02,z02,'--b','LineWidth',2)
    
    %      第三跨段
    Area=0.25*pi*d3^2;
    MI=pi/64*d3^4;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),shaft_density);
    
    L03=linspace(0,L3,40);
    L03=L03';
    x03=zeros(length(L03),1);
    z03=zeros(length(L03),1);
    
    hpros3v=[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)]*CoffC03v;
    hpros3w=[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)]*CoffC03w;
    zpros3=-[cos(mue*L03),sin(mue*L03)]*CoffD03;
    
    e = [(1:length(x03)-1)' (2:length(x03))'];
    c03 = [x03+hpros3w L1+L2+L03+zpros3 z03+hpros3v];
    fig03 = patch('Faces',e,'Vertices',c03, 'FaceColor','w','LineWidth',3);
    set(fig03,'EdgeColor','interp','FaceVertexCData',abs(hpros3w),'CDataMapping','scaled');
    plot3(x03,L1+L2+L03,z03,'--b','LineWidth',2)
    
    
    % 桨叶
    
    E=2.1e11;
    coe=5/6;   G=E/2/(1+0.32);
    
    %桨叶1
    l1=0:0.01:lpro1;
    l1=l1';
    x1=-flipud(l1);
    y1=zeros(length(l1),1);
    z1=zeros(length(l1),1);
    
    Area=b*h;
    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,~,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),pro_density);
    
    hpro1u=[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*CoffC1u;
    
    MI=b*h^3/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),pro_density);
    
    hpro1v=-[cosh(lambda*l1),sinh(lambda*l1),cos(hlambda*l1),sin(hlambda*l1)]*CoffC1v;
    
    zpro1=[cos(mue*l1),sin(mue*l1)]*CoffD1;
    
    
    e = [(1:length(x1)-1)' (2:length(x1))'];
    c12 = [x1+zpro1 y1+hpro1u z1+hpro1v];
    fig1 = patch('Faces',e,'Vertices',c12, 'FaceColor','w','LineWidth',3);
    set(fig1,'EdgeColor','interp','FaceVertexCData',abs(hpro1u),'CDataMapping','scaled');
    %     plot3(x1,y1,z1,'--b','LineWidth',2)
    
    
    % 桨叶2
    l2=0:0.02:lpro2;
    l2=l2';
    x2=l2;
    y2=zeros(length(l2),1);
    z2=zeros(length(l2),1);
    
    Area=b*h;
    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,~,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),pro_density);
    
    hpro2u=[cosh(lambda*l2),sinh(lambda*l2),cos(hlambda*l2),sin(hlambda*l2)]*CoffC2u;
    
    MI=b*h^3/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),pro_density);
    
    hpro2v=[cosh(lambda*l2),sinh(lambda*l2),cos(hlambda*l2),sin(hlambda*l2)]*CoffC2v;
    
    zpro2=[cos(mue*l2),sin(mue*l2)]*CoffD2;
    
    c22 = [x2+zpro2 y2+hpro2u z2+hpro2v];
    fig2 = patch('Faces',e,'Vertices',c22, 'FaceColor','w','LineWidth',3);
    set(fig2,'EdgeColor','interp','FaceVertexCData',abs(hpro2u),'CDataMapping','scaled');
    %     plot3(x2,y2,z2,'--b','LineWidth',2)
    
    
    % 桨叶3
    l3=0:0.02:lpro3;
    l3=l3';
    z3=flipud(l3);
    x3=zeros(length(l3),1);
    y3=zeros(length(l3),1);
    
    Area=b*h;
    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,~,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),pro_density);
    
    
    hpro3u=-[cosh(lambda*l3),sinh(lambda*l3),cos(hlambda*l3),sin(hlambda*l3)]*CoffC3u;
    
    MI=b*h^3/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),pro_density);
    
    hpro3w=[cosh(lambda*l3),sinh(lambda*l3),cos(hlambda*l3),sin(hlambda*l3)]*CoffC3w;
    
    zpro3=[cos(mue*l3),sin(mue*l3)]*CoffD3;
    
    c32 = [x3+hpro3w y3+hpro3u z3+zpro3];
    fig3 = patch('Faces',e,'Vertices',c32, 'FaceColor','w','LineWidth',3);
    set(fig3,'EdgeColor','interp','FaceVertexCData',abs(hpro3w),'CDataMapping','scaled');
    
    %     plot3(x3,y3,z3,'--b','LineWidth',2)
    % 桨叶4
    l4=0:0.02:lpro4;
    l4=l4';
    z4=-l4;
    x4=zeros(length(l4),1);
    y4=zeros(length(l4),1);
    
    Area=b*h;
    MI=b^3*h/12;
    [~,lambda,~,hlambda,~,~,~,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),pro_density);
    
    hpro4u=-[cosh(lambda*l4),sinh(lambda*l4),cos(hlambda*l4),sin(hlambda*l4)]*CoffC4u;
    
    MI=b*h^3/12;
    [~,lambda,~,hlambda,~,~,mue,~]=solution_par(E,Area,coe,G,MI,2*pi*result_fre(nn),pro_density);
    
    hpro4w=[cosh(lambda*l4),sinh(lambda*l4),cos(hlambda*l4),sin(hlambda*l4)]*CoffC4w;
    
    zpro4=[cos(mue*l4),sin(mue*l4)]*CoffD4;
    
    c42 = [x4+hpro4w y4+hpro4u z4+zpro4];
    fig4 = patch('Faces',e,'Vertices',c42, 'FaceColor','w','LineWidth',3);
    set(fig4,'EdgeColor','interp','FaceVertexCData',abs(hpro4w),'CDataMapping','scaled');
    %     plot3(x4,y4,z4,'--b','LineWidth',2)
    
    %     xlabel('横向模态位移/m')
    %     ylabel('纵向模态位移/m')
    %     zlabel('垂向模态位移/m')
    %     legend('固有结构','模态振型');
    %     text(0.5,1,0.5,[num2str(result_fre(nn)),'Hz']);
    view(3)
end
