function [propeller,mid_shaft_1,mid_shaft_2,f]=xiupropeller_matrix_planer(spring1_x,spring1_y,spring2_x, spring2_y,kr1z,kr2z,q0,hq0,w)

%��άģ�ͣ�������Ҷ

%��Ҷ1
E=9.2e10*(1+0.001i); %����ģ��

L=0.62;    
coe=5/6; G=E/2/(1+0.32);  %��ľ���¼���ϵ�������ɱ�0.32  �б䵯��ģ��
density=8900;    % ͭ�ܶ�

b=0.04; h=0.2;    %���ν���  
a=0; Force=1;   %���ص����Ҷ�Ҿ���,������ֵ

Area=b*h;  
MI=b^3*h/12;

[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);  % ��������

%������

if a~=0

    [T_2,T_1]=proper_trans(E,G,coe,MI,Area,q,lambda,hq,hlambda,belta,yita,a);

    trans1=T_2\T_1;
    Boundary_pro1=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

    Boundary1=Boundary_pro1*trans1;

    f=Boundary_pro1*(T_2\[zeros(3,1);Force]);
else
    Boundary1=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];
    f=[Force;0];

end


pc1=coe*G*Area*[belta*sinh(lambda*(L-a)),belta*cosh(lambda*(L-a)),yita*sin(hlambda*(L-a)),-yita*cos(hlambda*(L-a))]-...
    spring1_x*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];

pd0=[spring1_x,0];

pwc1=E*MI*[q*lambda*cosh(lambda*(L-a)),q*lambda*sinh(lambda*(L-a)),hq*hlambda*cos(hlambda*(L-a)),hq*hlambda*sin(hlambda*(L-a))]+...
    kr1z*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];

pwc0=[0,-kr1z*q0,0,kr1z*hq0];

%������
pc0=[-spring1_y,0,-spring1_y,0];

pd1=Area*E*[-mue*sin(mue*L), mue*cos(mue*L)]+spring1_y*[cos(mue*L),sin(mue*L)];

free1=Area*E*[0,mue];


%�����뽰Ҷ1�����Ĳ���
swc1=kr1z*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];

sd1=-spring1_y*[cos(mue*L),sin(mue*L)]; 

sc1=spring1_x*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];

pro1=[pc1,zeros(1,2),zeros(1,4),zeros(1,2),zeros(1,4),pd0;...
    pwc1,zeros(1,2),zeros(1,4),zeros(1,2),pwc0,zeros(1,2);...
    zeros(1,4),pd1,zeros(1,4),zeros(1,2),pc0,zeros(1,2);...
    zeros(1,4),free1,zeros(1,4),zeros(1,2),zeros(1,4),zeros(1,2);...
    Boundary1,zeros(2,14)];


mid_shaft_1=[zeros(1,4),sd1;swc1,zeros(1,2);sc1,zeros(1,2);...
    zeros(1,4),zeros(1,2);zeros(2,4),zeros(2,2)];


%��Ҷ2

%������,���ڽ�Ҷ2�ϼ���
% trans=T_1\T_2;    %�������ഫ�ݾ���

% Boundary_pro2=[coe*G*Area*[belta*sinh(lambda*a),belta*cosh(lambda*a),yita*sin(hlambda*a),-yita*cos(hlambda*a)];...
%     E*MI*[q*lambda*cosh(lambda*a),q*lambda*sinh(lambda*a),hq*hlambda*cos(hlambda*a),hq*hlambda*sin(hlambda*a)]];

Boundary2=[coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)];...
    E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]];


% Boundary2=Boundary_pro2*trans;   % Ҷ�ұ߽�--->Ҷ��

% f=Boundary_pro2*(T_1\[zeros(3,1);Force]);     % ��λ����

pc2=[spring2_x,coe*G*Area*belta,spring2_x,-coe*G*Area*yita];

pd0=[-spring2_x,0];


pwc2=E*MI*[q*lambda,0,hq*hlambda,0]-kr2z*[0,q,0,-hq];


pwc0=[0,kr2z*q0,0,-kr2z*hq0];


pc0=[spring2_y,0,spring2_y,0];


pd2=[0, Area*E*mue]-spring2_y*[1,0];


free2=Area*E*[-mue*sin(mue*L),mue*cos(mue*L)];


%��Ҷ2�������໥Ӱ��Ĳ���
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





