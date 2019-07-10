function [propeller,mid_shaft_i]=propeller3_matrix_space(E,density,ue,ku,kv,kw,kru,krw,b,h,L,q0v,hq0v,w,ii,n)

coe=5/6; G=E/2/(1+ue);       %弹性模量、铁木辛柯剪切系数、切变弹性模量
Area=b*h;

%%横向1---振动(u方向）xy平面

MI=h*b^3/12;

[q,lambda,hq,hlambda,belta,yita,~,~]=solution_par(E,Area,coe,G,MI,w,density);

Boundary_u=[0,belta,0,-yita;q*lambda,0,hq*hlambda,0];

Pcu=coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)]-...
    ku*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

Pd0=[ku,0];

Pwcu=E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]+...
    krw*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];

Pwc0v=[0,-krw*q0v,0,krw*hq0v];    %传参q0z,hq0z


% 轴系耦合参数
Swcu=krw*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];   %主轴V方向

Scu=ku*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

%%横向2---振动（w方向）yz平面

MI=b*h^3/12;
[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);   %截面惯性矩计算与方向有关


Boundary_w=[0,belta,0,-yita;q*lambda,0,hq*hlambda,0];

Pcw=coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)]-...
    kw*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

Pc0w=[kw,0,kw,0];

Pwcw=E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]+...
    kru*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];

Pwh0u=[-kru,0];

%%纵向振动
Pc0v=[-kv,0,-kv,0];
Pdv=Area*E*[-mue*sin(mue*L), mue*cos(mue*L)]+kv*[cos(mue*L),sin(mue*L)];
Fzong=[0,mue];

% 桨叶-->主轴 中间矩阵
Sdv=-kv*[cos(mue*L),sin(mue*L)];   %主轴V方向

Scw=-kw*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];   %主轴w方向

Sncw=kru*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];


%%桨叶边界条件
propeller=[zeros(1,4),zeros(1,4),Pd0,zeros(1,2),zeros(1,(ii-1)*10),Pcu,zeros(1,4),zeros(1,2),zeros(1,(n-ii)*10);...
    Pwc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),Pwcu,zeros(1,4),zeros(1,2),zeros(1,(n-ii)*10);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*10),Boundary_u,zeros(2,4),zeros(2,2),zeros(2,(n-ii)*10);...
    zeros(1,4),Pc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),Pcw,zeros(1,2),zeros(1,(n-ii)*10);...
    zeros(1,4),zeros(1,4),zeros(1,2),Pwh0u,zeros(1,(ii-1)*10),zeros(1,4),Pwcw,zeros(1,2),zeros(1,(n-ii)*10);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*10),zeros(2,4),Boundary_w,zeros(2,2),zeros(2,(n-ii)*10);...
    Pc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),zeros(1,4),Pdv,zeros(1,(n-ii)*10);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),zeros(1,4),Fzong,zeros(1,(n-ii)*10)];

%不考虑桨叶扭转的条件矩阵,主轴有扭转
mid_shaft_i=[Scu,zeros(1,4),zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2);...
    zeros(1,4),Sncw,zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2);...
    zeros(1,4),zeros(1,4),Sdv;...
    Swcu,zeros(1,4),zeros(1,2);...
    zeros(2,4),zeros(2,4),zeros(2,2);...
    zeros(1,4),Scw,zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2);...
    zeros(2,4),zeros(2,4),zeros(2,2)];

end