function [propeller,mid_shaft_i]=propeller1_matrix_space(E,density,ue,ku,kv,kw,kru,krv,b,h,L,q0w,hq0w,w,ii,n)

% ���ؽ�Ҷ1�߽��������������������������
coe=5/6; G=E/2/(1+ue);       %��Ҷ����ģ������ľ���¼���ϵ�����б䵯��ģ��

Area=b*h;

% ����1---��(u����xzƽ��
MI=h*b^3/12;

[q,lambda,hq,hlambda,belta,yita,~,~]=solution_par(E,Area,coe,G,MI,w,density);

Boundary_u=[0,belta,0,-yita;q*lambda,0,hq*hlambda,0];

Pcu=coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)]-...
    ku*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

Pd0=[ku,0];

Pwcu=E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]+...
    krv*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];

Pwc0w=[0,-krv*q0w,0,krv*hq0w];    %��ϵ����q0z,hq0z

% ����ϵ��ϲ���
Swcu=krv*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];   %����w����xzƽ��

Scu=ku*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

% ����2---��v����yzƽ��

MI=b*h^3/12;   %������Ծؼ����뷽���й�
[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);

Boundary_v=[0,belta,0,-yita;q*lambda,0,hq*hlambda,0];

Pcv=coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)]-...
    kv*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

Pc0v=[kv,0,kv,0];

Pwcv=E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]+...
    kru*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];

Pwh0u=[-kru,0];

% ������
Pc0w=[-kw,0,-kw,0];
Pdw=Area*E*[-mue*sin(mue*L), mue*cos(mue*L)]+kw*[cos(mue*L),sin(mue*L)];
Fzong=[0,mue];

% ��Ҷ-->���� �м����
Scv=-kv*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];    %����V����

Sdw=-kw*[cos(mue*L),sin(mue*L)];   %����W����

Sncv=kru*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];


% ��Ҷ�߽�����
propeller=[zeros(1,4),zeros(1,4),Pd0,zeros(1,2),zeros(1,(ii-1)*10),Pcu,zeros(1,4),zeros(1,2),zeros(1,(n-ii)*10);...
    zeros(1,4),Pwc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),Pwcu,zeros(1,4),zeros(1,2),zeros(1,(n-ii)*10);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*10),Boundary_u,zeros(2,4),zeros(2,2),zeros(2,(n-ii)*10);...
    Pc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),Pcv,zeros(1,2),zeros(1,(n-ii)*10);...
    zeros(1,4),zeros(1,4),zeros(1,2),Pwh0u,zeros(1,(ii-1)*10),zeros(1,4),Pwcv,zeros(1,2),zeros(1,(n-ii)*10);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*10),zeros(2,4),Boundary_v,zeros(2,2),zeros(2,(n-ii)*10);...
    zeros(1,4),Pc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),zeros(1,4),Pdw,zeros(1,(n-ii)*10);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),zeros(1,4),Fzong,zeros(1,(n-ii)*10)];

% ��ϵ�߽�������ؽ�Ҷ
mid_shaft_i=[Scu,zeros(1,4),zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2);...
    zeros(1,4),Sncv,zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2);...
    zeros(1,4),Scv,zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2);...
    zeros(2,4),zeros(2,4),zeros(2,2);...
    zeros(1,4),zeros(1,4),Sdw;...
    Swcu,zeros(1,4),zeros(1,2);...
    zeros(2,4),zeros(2,4),zeros(2,2)];

end