function [result_matrix,ffu,ff2]=condition_matrix_space(w,a,one_two)

% 返回桨-轴系统给定频率值处平衡条件矩阵
parameters2;   %参数初始化

%%此处参数固定

% k1u=2e7; k1v=2.5e8; k1w=3e7;     %桨叶与轴各向弹簧刚度   N/m
% k2u=2e7; k2v=2.5e8; k2w=3e7;
% k3u=2e7; k3v=2.5e8; k3w=3e7;
% k4u=2e7; k4v=2.5e8; k4w=3e7;
% 
% kr1u=4e7; kr1v=1e7; kr1w=2e7;  %扭转刚度 Nm/rad
% kr2u=4e7; kr2v=1e7; kr2w=2e7;
% kr3u=4e7; kr3v=1e7; kr3w=2e7;
% kr4u=4e7; kr4v=1e7; kr4w=2e7;

k1u=2e8; k1v=2e9; k1w=3e9;     %桨叶与轴各向弹簧刚度   N/m
k2u=2e8; k2v=2e9; k2w=3e9;
k3u=2e8; k3v=2e9; k3w=3e9;
k4u=2e8; k4v=2e9; k4w=3e9;

kr1u=1e9; kr1v=1e8; kr1w=2e9;  %扭转刚度
kr2u=1e9; kr2v=1e8; kr2w=2e9;
kr3u=1e9; kr3v=1e8; kr3w=2e9;
kr4u=1e9; kr4v=1e8; kr4w=2e9;

% k1u=2e8; k1v=2e9; k1w=3e9;     %桨叶与轴各向弹簧刚度   N/m
% k2u=2e8; k2v=2e9; k2w=3e9;
% k3u=2e8; k3v=3e9; k3w=2e9;
% k4u=2e8; k4v=3e9; k4w=2e9;
% 
% kr1u=1e9; kr1v=1e8; kr1w=2e9;  %扭转刚度
% kr2u=1e9; kr2v=1e8; kr2w=2e9;
% kr3u=1e9; kr3v=2e9; kr3w=1e8;
% kr4u=1e9; kr4v=2e9; kr4w=1e8;

%%固定参数结束

sumku=k1u+k2u+k3u+k4u;  sumkv=k1v+k2v+k3v+k4v;    sumkw=k1w+k2w+k3w+k4w;

sumkru=kr1u+kr2u+kr3u+kr4u;   sumkrv=kr1v+kr2v+kr3v+kr4v;  sumkrw=kr1w+kr2w+kr3w+kr4w;


% 第一跨段两向弯曲+纵向+扭转

%第一跨段横向振动---xy平面
Area1=0.25*pi*d1^2;   J1=pi*d1^4/32;
MI1=pi/64*d1^4;

[q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area1,coe,G,MI1,w,density);

Sc0v=[-m*w^2+sumkv,coe*G*Area1*belta,-m*w^2+sumkv,-coe*G*Area1*yita];
Swc0v=[E*MI1*q*lambda,w^2*Imw*q-sumkrw*q,E*MI1*hq*hlambda,-w^2*Imw*hq+sumkrw*hq];

q0v=q;  hq0v=hq;


%第一跨段横向振动--xz平面， 圆形截面参数、惯量相同

Sc0w=[-m*w^2+sumkw,coe*G*Area1*belta,-m*w^2+sumkw,-coe*G*Area1*yita];

Swc0w=[E*MI1*q*lambda,w^2*Imv*q-sumkrv*q,E*MI1*hq*hlambda,-w^2*Imv*hq+sumkrv*hq];

q0w=q;  hq0w=hq;

%第一跨段纵向振动扭转振动

Sd0=[w^2*m-sumku,Area1*E*mue];

Snh0=[w^2*Imu-sumkru,G*J1*gamma];

% 第一跨段结束

% 第三跨段

MI3=pi/64*d3^4;
Area3=0.25*pi*d3^2;
[q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area3,coe,G,MI3,w,density);

% xy平面-v向边界条件

B0v_3=[coe*G*Area3*[belta*sinh(lambda*L03),belta*cosh(lambda*L03),yita*sin(hlambda*L03),-yita*cos(hlambda*L03)]-...
    k0v*[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)];...
    E*MI3*[q*lambda*cosh(lambda*L03),q*lambda*sinh(lambda*L03),hq*hlambda*cos(hlambda*L03),hq*hlambda*sin(hlambda*L03)]];

% xz平面-w向边界条件

B0w_3=[coe*G*Area3*[belta*sinh(lambda*L03),belta*cosh(lambda*L03),yita*sin(hlambda*L03),-yita*cos(hlambda*L03)]-...
    k0w*[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)];...
    E*MI3*[q*lambda*cosh(lambda*L03),q*lambda*sinh(lambda*L03),hq*hlambda*cos(hlambda*L03),hq*hlambda*sin(hlambda*L03)]];

% 纵向边界条件
Guzong_3=[-Area3*E*mue*sin(mue*L03)+k0u*cos(mue*L03),Area3*E*mue*cos(mue*L03)+k0u*sin(mue*L03)];

% 扭转边界条件
Sngu_3=[cos(gamma*L03),sin(gamma*L03)];

% 最右侧到第一跨段传递

% 传递矩阵
[hengmatrix,zongmatrix,niumatrix,~,~,~]=transfer_matrix(spring1,spring2,L01,L02,d1,d2,d3,w,density);

% [hengmatrixv,zongmatrix,niumatrix,~,~,~]=transfer_matrix(spring1v,spring2v,L01,L02,d1,d2,d3,w,density);
% 
% [hengmatrixw,~,~,~,~,~]=transfer_matrix(spring1w,spring2w,L01,L02,d1,d2,d3,w,density);


% 等效到第一跨段
B0v=B0v_3*hengmatrix;
B0w=B0w_3*hengmatrix;
Guzong=Guzong_3*zongmatrix;
Sngu0=Sngu_3*niumatrix;


% 轴第一跨段边界条件
mid_shaft_0=[zeros(1,4),zeros(1,4),Sd0,zeros(1,2);...
    zeros(1,4),zeros(1,4),Guzong,zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2),Snh0;....
    zeros(1,4),zeros(1,4),zeros(1,2),Sngu0;...
    Sc0v,zeros(1,4),zeros(1,2),zeros(1,2);...
    Swc0v,zeros(1,4),zeros(1,2),zeros(1,2);
    B0v,zeros(2,4),zeros(2,2),zeros(2,2);...
    zeros(1,4),Sc0w,zeros(1,2),zeros(1,2);...
    zeros(1,4),Swc0w,zeros(1,2),zeros(1,2);...
    zeros(2,4),B0w,zeros(2,2),zeros(2,2)];


%第一根桨叶

ii=1; L=L1;   %调用函数得到第i根桨叶条件矩阵

if one_two==1

    [propeller_1,mid_shaft_1,ffu,ff2]=first_propeller1_matrix_space(pro_E,pro_density,ue,k1u,k1v,k1w,kr1u,kr1v,kr1w,b,h,L,q0w,hq0w,w,ii,n,a,ratio);

else
    [propeller_1,mid_shaft_1]=propeller1_matrix_space(pro_E,pro_density,ue,k1u,k1v,k1w,kr1u,kr1v,kr1w,b,h,L,q0w,hq0w,w,ii,n,ratio);
end

%第二根桨叶

ii=2;  L=L2;   %调用函数得到第i根桨叶条件矩阵
[propeller_2,mid_shaft_2]=propeller2_matrix_space(pro_E,pro_density,ue,k2u,k2v,k2w,kr2u,kr2v,kr2w,b,h,L,q0w,hq0w,w,ii,n,ratio);


%第三根桨叶
ii=3;   L=L3;   %调用函数得到第i根桨叶条件矩阵
if one_two==1

    [propeller_3,mid_shaft_3]=first_propeller3_matrix_space(pro_E,pro_density,ue,k3u,k3v,k3w,kr3u,kr3v,kr3w,b,h,L,q0v,hq0v,w,ii,n,ratio);
else
    [propeller_3,mid_shaft_3,ffu,ff2]=propeller3_matrix_space(pro_E,pro_density,ue,k3u,k3v,k3w,kr3u,kr3v,kr3w,b,h,L,q0v,hq0v,w,ii,n,a,ratio);
end
%第四根桨叶

ii=4;   L=L4;   %调用函数得到第i根桨叶条件矩阵
[propeller_4,mid_shaft_4]=propeller4_matrix_space(pro_E,pro_density,ue,k4u,k4v,k4w,kr4u,kr4v,kr4w,b,h,L,q0v,hq0v,w,ii,n,ratio);

shaft=[mid_shaft_0,mid_shaft_1,mid_shaft_2,mid_shaft_3,mid_shaft_4];

result_matrix=[propeller_1;propeller_2;propeller_3;propeller_4;shaft];
end