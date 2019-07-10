function [hengmatrix,zongmatrix,niumatrix,heng2,zong2,niu2]=transfer_matrix(hspring1,hspring2,L1,L2,d1,d2,d3,w,density)


E=2.1e11;  coe=0.9; G=E/2/(1+0.3);   % 轴系参数，各跨段相同，圆形截面
% density=7500;
%横向振动

%第一个轴承支撑
MI1=pi/64*d1^4;
MI2=pi/64*d2^4;
Area1=0.25*pi*d1^2;  
Area2=0.25*pi*d2^2;  
J1=pi*d1^4/32;
J2=pi*d2^4/32;
J3=pi*d3^4/32;

ei=MI1/MI2;  kga=Area1/Area2;       %ei=E1*MI1/E1*MI1   kga=coe1*G1*Area1/coe2*G2*Area2    %轴系圆柱截面

[q,lambda,hq,hlambda,~,~,mue,gamma]=solution_par(E,Area1,coe,G,MI1,w,density);

epsilion1=ei*q*lambda;   epsilion2=ei*hq*hlambda;   epsilion3=kga*(q-lambda);   epsilion4=kga*(hq+hlambda);

laL=lambda*L1; hlaL=hlambda*L1;    %间断点位置

T_left_1=[cosh(laL),sinh(laL),cos(hlaL),sin(hlaL); ...
    q*sinh(laL),q*cosh(laL),hq*sin(hlaL),-hq*cos(hlaL); ...
    epsilion1*cosh(laL),epsilion1*sinh(laL),epsilion2*cos(hlaL),epsilion2*sin(hlaL);...
    epsilion3*sinh(laL),epsilion3*cosh(laL),epsilion4*sin(hlaL),-epsilion4*cos(hlaL)];

N_left_1=[cos(mue*L1),sin(mue*L1);-Area1*E*mue*sin(mue*L1),Area1*E*mue*cos(mue*L1)];

Z_left_1=[cos(gamma*L1),sin(gamma*L1);-J1*G*gamma*sin(gamma*L1),J1*G*gamma*cos(gamma*L1)];

[q,lambda,hq,hlambda,~,~,mue,gamma]=solution_par(E,Area2,coe,G,MI2,w,density);

T_right_2=[1,0,1,0; 0,q,0,-hq; q*lambda,0,hq*hlambda,0; hspring1,q-lambda,hspring1,-hq-hlambda];   %左右传递矩阵
N_right_2=[1,0;0,Area2*E*mue];
Z_right_2=[1,0;0,J2*G*gamma];

%第二个轴承支撑
MI3=pi/64*d3^4;
Area3=0.25*pi*d3^2;  
ei=MI2/MI3;  kga=Area2/Area3; 

[q,lambda,hq,hlambda,~,~,mue,gamma]=solution_par(E,Area2,coe,G,MI2,w,density);
epsilion1=ei*q*lambda;   epsilion2=ei*hq*hlambda;   epsilion3=kga*(q-lambda);   epsilion4=kga*(hq+hlambda);

laL=lambda*L2; hlaL=hlambda*L2;    %间断点位置

T_left_2=[cosh(laL),sinh(laL),cos(hlaL),sin(hlaL); ...
    q*sinh(laL),q*cosh(laL),hq*sin(hlaL),-hq*cos(hlaL); ...
    epsilion1*cosh(laL),epsilion1*sinh(laL),epsilion2*cos(hlaL),epsilion2*sin(hlaL);...
    epsilion3*sinh(laL),epsilion3*cosh(laL),epsilion4*sin(hlaL),-epsilion4*cos(hlaL)];
N_left_2=[cos(mue*L2),sin(mue*L2);-Area2*E*mue*sin(mue*L2),Area2*E*mue*cos(mue*L2)];
Z_left_2=[cos(gamma*L2),sin(gamma*L2);-J2*G*gamma*sin(gamma*L2),J2*G*gamma*cos(gamma*L2)];

[q,lambda,hq,hlambda,~,~,mue,gamma]=solution_par(E,Area3,coe,G,MI3,w,density);

T_right_3=[1,0,1,0; 0,q,0,-hq; q*lambda,0,hq*hlambda,0; hspring2,q-lambda,hspring2,-hq-hlambda];   %左右传递矩阵

N_right_3=[1,0;0,Area3*E*mue];
Z_right_3=[1,0;0,J3*G*gamma];


heng2=inv(T_right_2)*T_left_1;
zong2=inv(N_right_2)*N_left_1;
niu2=inv(Z_right_2)*Z_left_1;

hengmatrix=inv(T_right_3)*T_left_2*(inv(T_right_2)*T_left_1);  %横向振动传递矩阵

zongmatrix=inv(N_right_3)*N_left_2*(inv(N_right_2)*N_left_1);   %纵向振动传递矩阵

niumatrix=inv(Z_right_3)*Z_left_2*(inv(Z_right_2)*Z_left_1);    %扭转振动传递矩阵
end