function [result_matrix,f]=xiu_matrix_planer(w)

%ƽ���ڽ�������2����Ҷ

E=2.1e11*(1+0.001i);  

coe=0.9; G=E/2/(1+0.3);   density=7500;  % ��ϵ�������������ͬ��Բ�ν���

%�ɱ����

spring1_x=4e10;  spring1_y=5e10;  spring2_x=4e10;  spring2_y=5e10;   kr1z=5e10; kr2z=5e10;  %��Ҷ-��ϵ�����ӵ��ɸն�

m=223.68; Imz=2.98;  %���������ת������

k0x=7.66e6;  k0y=1e8;   %  �����Ҷ����򡢾���նȣ�������У�

spring1=1e8;    %��ϵ�������1,2֧�Ÿն�
spring2=1e8;

L1=0.48;   L2=1.8;   L3=1.8;
d1=0.16;  d2=0.14;   d3=0.11;   %����γ��ȡ�ֱ��   

%�ɱ����

sumkx=spring1_x+spring2_x;
sumky=spring1_y+spring2_y;
sumkrz=kr1z+kr2z;

%����

%��һ���
Area=0.25*pi*d1^2;  
MI=pi/64*d1^4;

Area2=0.25*pi*d2^2;  
Area3=0.25*pi*d3^2;  
hspring1=spring1/(coe*G*Area2);

hspring2=spring2/(coe*G*Area3);

[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);   %��һ���

% ������
sc0=[-m*w^2+sumky,coe*G*Area*belta,-m*w^2+sumky,-coe*G*Area*yita];  

swc0=[E*MI*q*lambda,(w^2*Imz-sumkrz)*q,E*MI*hq*hlambda,-(w^2*Imz-sumkrz)*hq];

%������
sd0=[w^2*m-sumkx,Area*E*mue];

q0=q;  hq0=hq;


%������Σ��Ҷ�֧�ű߽�
Area=0.25*pi*d3^2;  
MI=pi/64*d3^4;

[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);

Boundary_3=[coe*G*Area*[belta*sinh(lambda*L3),belta*cosh(lambda*L3),yita*sin(hlambda*L3),-yita*cos(hlambda*L3)]-...
    k0y*[cosh(lambda*L3),sinh(lambda*L3),cos(hlambda*L3),sin(hlambda*L3)];
E*MI*[q*lambda*cosh(lambda*L3),q*lambda*sinh(lambda*L3),hq*hlambda*cos(hlambda*L3),hq*hlambda*sin(hlambda*L3)]];


gu3=[-Area*E*mue*sin(mue*L3)+k0x*cos(mue*L3),Area*E*mue*cos(mue*L3)+k0x*sin(mue*L3)];


[hengmatrix,zongmatrix,~,~,~,~]=transfer_matrix(hspring1,hspring2,L1,L2,d1,d2,d3,w,density);   % ���ݾ���


Boundary0=Boundary_3*hengmatrix;   %�Ҷ�������Ч�����
gu=gu3*zongmatrix;

shaft0=[sc0,zeros(1,2);swc0,zeros(1,2);zeros(1,4),sd0;...
    zeros(1,4),gu;Boundary0,zeros(2,2)];

[propeller,mid_shaft_1,mid_shaft_2,f]=xiupropeller_matrix_planer(spring1_x,spring1_y,spring2_x, spring2_y,kr1z,kr2z,q0,hq0,w);
shaft=[mid_shaft_1,mid_shaft_2,shaft0];

result_matrix=[propeller;shaft];


end