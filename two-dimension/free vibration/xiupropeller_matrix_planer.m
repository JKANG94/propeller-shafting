function [propeller,mid_shaft_1,mid_shaft_2]=xiupropeller_matrix_planer(spring1_x,spring1_y,spring2_x, spring2_y,kr1z,kr2z,q0,hq0,w)

%平面内交叉梁，2根桨叶

%桨叶1
E=9.2e10; 

L=0.62;    %弹性模量、铁木辛柯剪切系数、切变弹性模量
coe=5/6; G=E/2/(1+0.32);  %矩形截面  泊松比0.32 
density=8900;  % 铜

b=0.04; h=0.2;
Area=b*h;  
MI=b^3*h/12;

[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);  % 基本参数

Boundary1=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];


pc1=coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)]-...
    spring1_x*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];


pd0=[spring1_x,0];


pwc1=E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]+...
    kr1z*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];



pwc0=[0,-kr1z*q0,0,kr1z*hq0];

%纵向振动
pc0=[-spring1_y,0,-spring1_y,0];


pd1=Area*E*[-mue*sin(mue*L), mue*cos(mue*L)]+spring1_y*[cos(mue*L),sin(mue*L)];


free1=Area*E*[0,mue];


%主轴与桨叶1关联的参数
swc1=kr1z*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];


sd1=-spring1_y*[cos(mue*L),sin(mue*L)]; 


sc1=spring1_x*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];



pro1=[pc1,zeros(1,2),zeros(1,4),zeros(1,2),zeros(1,4),pd0;...
    pwc1,zeros(1,2),zeros(1,4),zeros(1,2),pwc0,zeros(1,2);...
    zeros(1,4),pd1,zeros(1,4),zeros(1,2),pc0,zeros(1,2);...
    zeros(1,4),free1,zeros(1,4),zeros(1,2),zeros(1,4),zeros(1,2);...
    Boundary1,zeros(2,14)];


mid_shaft_1=[zeros(1,4),sd1;swc1,zeros(1,2);sc1,zeros(1,2);...
    zeros(1,4),zeros(1,2);zeros(2,4),zeros(2,2)];


%桨叶2
E=9.2e10; 

L=0.62;    %弹性模量、铁木辛柯剪切系数、切变弹性模量
coe=5/6; G=E/2/(1+0.32); 
density=8900;

b=0.04; h=0.2;
Area=b*h;  
MI=b^3*h/12;

[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);


Boundary2=[coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)];...
    E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]];


pc2=[spring2_x,coe*G*Area*belta,spring2_x,-coe*G*Area*yita];

pd0=[-spring2_x,0];


pwc2=E*MI*[q*lambda,0,hq*hlambda,0]-kr2z*[0,q,0,-hq];


pwc0=[0,kr2z*q0,0,-kr2z*hq0];


pc0=[spring2_y,0,spring2_y,0];


pd2=[0, Area*E*mue]-spring2_y*[1,0];


free2=Area*E*[-mue*sin(mue*L),mue*cos(mue*L)];


%主轴与桨叶2关联的参数
swc2=kr2z*[0,q,0,-hq];  

sd2=-spring2_y*[1,0];   

sc2=spring2_x*[1,0,1,0];

pro2=[zeros(1,4),zeros(1,2),pc2,zeros(1,2),zeros(1,4),pd0;...
    zeros(1,4),zeros(1,2),pwc2,zeros(1,2),pwc0,zeros(1,2);...
    zeros(1,4),zeros(1,2),zeros(1,4),pd2,pc0,zeros(1,2);...
    zeros(1,4),zeros(1,2),zeros(1,4),free2,zeros(1,4),zeros(1,2);...
    zeros(2,4),zeros(2,2),Boundary2,zeros(2,8)];

mid_shaft_2=[zeros(1,4),sd2;swc2,zeros(1,2);sc2,zeros(1,2);...
    zeros(1,4),zeros(1,2);zeros(2,4),zeros(2,2)];

propeller=[pro1;pro2];





