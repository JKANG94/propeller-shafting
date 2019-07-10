function [result_matrix,fu,fv]=first_condition_matrix_space(w,a,attack_ang)

n=4;   %桨叶数
k1u=2e8; k1v=2e9; k1w=3e9;     %各向弹簧刚度   
k2u=2e8; k2v=2e9; k2w=3e9;    
k3u=2e8; k3v=2e9; k3w=3e9;    
k4u=2e8; k4v=2e9; k4w=3e9;    

kr1u=1e9; kr1v=1e8; kr1w=2e9;
kr2u=1e9; kr2v=1e8; kr2w=2e9;
kr3u=1e9; kr3v=1e8; kr3w=2e9;
kr4u=1e9; kr4v=1e8; kr4w=2e9;

L1=1.6;  L2=1.6;  L3=1.6; L4=1.6;  %各桨叶长度
b=0.07;  h=0.3;  %桨叶截面

sumku=k1u+k2u+k3u+k4u;  sumkv=k1v+k2v+k3v+k4v;    sumkw=k1w+k2w+k3w+k4w;

sumkru=kr1u+kr2u+kr3u+kr4u;   sumkrv=kr1v+kr2v+kr3v+kr4v;  sumkrw=kr1w+kr2w+kr3w+kr4w;

%推进轴  物理参数

E=2.1e11*(1+0.01i);  coe=0.9; G=E/2/(1+0.3);  density=7500;
m=1860;                                  
 
L01=0.3;   L02=5;   L03=10;
d1=0.4;  d2=0.36;   d3=0.32;      %各跨段长度、直径   


Imu=73;  Imv=85;  Imw=85;         %桨毂三向转动惯量
k0u=2.3e8;  k0v=8.8e7;  k0w=8.8e7;       %轴后端三向移动弹簧刚度

spring1=1.8e8;   %轴承支撑刚度,两向相同
spring2=1.8e8;

Area2=0.25*pi*d2^2;  
Area3=0.25*pi*d3^2;  
hspring1=spring1/(coe*G*Area2);

hspring2=spring2/(coe*G*Area3);

%第一跨段横向振动---v
Area=0.25*pi*d1^2;   J=pi*d1^4/32; 
MI=pi/64*d1^4;  
[q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area,coe,G,MI,w,density);    

   
Sc0v=[-m*w^2+sumkv,coe*G*Area*belta,-m*w^2+sumkv,-coe*G*Area*yita];   
Swc0v=[E*MI*q*lambda,w^2*Imw*q-sumkrw*q,E*MI*hq*hlambda,-w^2*Imw*hq+sumkrw*hq];

q0v=q;  hq0v=hq;


%第一跨段 横向振动2--w 圆形截面参数惯量相同

% MI=pi/64*d1^4;
% [q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area,coe,G,MI,w);    

Sc0w=[-m*w^2+sumkw,coe*G*Area*belta,-m*w^2+sumkw,-coe*G*Area*yita]; 

Swc0w=[E*MI*q*lambda,w^2*Imv*q-sumkrv*q,E*MI*hq*hlambda,-w^2*Imv*hq+sumkrv*hq];  

q0w=q;  hq0w=hq;


%第一跨段纵向振动扭转振动

Sd0=[w^2*m-sumku,Area*E*mue];

Snh0=[w^2*Imu-sumkru,G*J*gamma];


%第三跨段
Area=0.25*pi*d3^2;  
MI=pi/64*d3^4;

[q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area,coe,G,MI,w,density);

%v向边界
B0v_3=[coe*G*Area*[belta*sinh(lambda*L03),belta*cosh(lambda*L03),yita*sin(hlambda*L03),-yita*cos(hlambda*L03)]-...
    k0v*[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)];...
    E*MI*[q*lambda*cosh(lambda*L03),q*lambda*sinh(lambda*L03),hq*hlambda*cos(hlambda*L03),hq*hlambda*sin(hlambda*L03)]];  %主轴右端边界条件

%w向边界
B0w_3=[coe*G*Area*[belta*sinh(lambda*L03),belta*cosh(lambda*L03),yita*sin(hlambda*L03),-yita*cos(hlambda*L03)]-...
    k0w*[cosh(lambda*L03),sinh(lambda*L03),cos(hlambda*L03),sin(hlambda*L03)];...
    E*MI*[q*lambda*cosh(lambda*L03),q*lambda*sinh(lambda*L03),hq*hlambda*cos(hlambda*L03),hq*hlambda*sin(hlambda*L03)]];

%纵向边界
Guzong_3=[-Area*E*mue*sin(mue*L03)+k0u*cos(mue*L03),Area*E*mue*cos(mue*L03)+k0u*sin(mue*L03)];
%扭转边界
Sngu_3=[cos(gamma*L03),sin(gamma*L03)];

%传递矩阵
[hengmatrix,zongmatrix,niumatrix,~,~,~]=transfer_matrix(hspring1,hspring2,L01,L02,d1,d2,d3,w);

B0v=B0v_3*hengmatrix;
B0w=B0w_3*hengmatrix;
Guzong=Guzong_3*zongmatrix;
Sngu0=Sngu_3*niumatrix;

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
[propeller_1,mid_shaft_1,fu,fv]=first_propeller1_matrix_space(k1u,k1v,k1w,kr1u,kr1v,b,h,L,q0w,hq0w,w,ii,n,a,attack_ang);


%第二根桨叶

ii=2;  L=L2;   %调用函数得到第i根桨叶条件矩阵
[propeller_2,mid_shaft_2]=propeller2_matrix_space(k2u,k2v,k2w,kr2u,kr2v,b,h,L,q0w,hq0w,w,ii,n);


%第三根桨叶
ii=3;   L=L3;   %调用函数得到第i根桨叶条件矩阵
[propeller_3,mid_shaft_3]=first_propeller3_matrix_space(k3u,k3v,k3w,kr3u,kr3w,b,h,L,q0v,hq0v,w,ii,n);


%第四根桨叶

ii=4;   L=L4;   %调用函数得到第i根桨叶条件矩阵
[propeller_4,mid_shaft_4]=propeller4_matrix_space(k4u,k4v,k4w,kr4u,kr4w,b,h,L,q0v,hq0v,w,ii,n);

shaft=[mid_shaft_0,mid_shaft_1,mid_shaft_2,mid_shaft_3,mid_shaft_4];

result_matrix=[propeller_1;propeller_2;propeller_3;propeller_4;shaft];



