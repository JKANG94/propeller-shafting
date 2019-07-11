function [propeller,mid_shaft_i]=propeller1_matrix_space(E,density,ue,ku,kv,kw,kru,krv,krw,b,h,L,q0w,hq0w,w,ii,n,ratio)


coe=5/6; G=E/2/(1+ue);       %弹性模量、铁木辛柯剪切系数、切变弹性模量
Area=b*h;

% ratio=0.299;
J=ratio*h*b^3;

%横向1---振动(u方向）
MI=h*b^3/12;

[q,lambda,hq,hlambda,belta,yita,~,~]=solution_par(E,Area,coe,G,MI,w,density);


Boundary_u=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

Pcu=coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)]-...
    ku*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

Pd0=[ku,0];

Pwcu=E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]+...
    krv*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];

Pwc0w=[0,-krv*q0w,0,krv*hq0w];    %传参q0z,hq0z

% 轴系相关参数
Swcu=krv*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];   %主轴w方向

Scu=ku*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

% 横向2---振动（v方向）

MI=b*h^3/12;
[q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area,coe,G,MI,w,density);   %截面惯性矩计算与方向有关

Boundary_v=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

Pcv=coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)]-...
    kv*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];

Pc0v=[kv,0,kv,0];

Pwcv=E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]+...
    kru*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];

Pwh0u=[-kru,0];

%纵向振动
Pc0w=[-kw,0,-kw,0];
Pdw=Area*E*[-mue*sin(mue*L), mue*cos(mue*L)]+kw*[cos(mue*L),sin(mue*L)];
Fzong=Area*E*[0,mue];

%扭转振动
PNc0v=-krw*[0,q0w,0,-hq0w];   %q0v=q0w;;;  hq0v=hq0w
Ph1=G*J*[-gamma*sin(gamma*L),gamma*cos(gamma*L)]+krw*[cos(gamma*L),sin(gamma*L)];
FN1=G*J*[0,gamma];

% 桨叶-->主轴 中间矩阵
Scv=-kv*[cosh(lambda*L),sinh(lambda*L),cos(hlambda*L),sin(hlambda*L)];    %主轴V方向
Swh1=krw*[cos(gamma*L),sin(gamma*L)];

Sdw=-kw*[cos(mue*L),sin(mue*L)];   %主轴W方向

Sncv=kru*[q*sinh(lambda*L),q*cosh(lambda*L),hq*sin(hlambda*L),-hq*cos(hlambda*L)];


propeller=[zeros(1,4),zeros(1,4),Pd0,zeros(1,2),zeros(1,(ii-1)*12),Pcu,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(1,4),Pwc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),Pwcu,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*12),Boundary_u,zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(n-ii)*12);...
    Pc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),Pcv,zeros(1,2),zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(1,4),zeros(1,4),zeros(1,2),Pwh0u,zeros(1,(ii-1)*12),zeros(1,4),Pwcv,zeros(1,2),zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*12),zeros(2,4),Boundary_v,zeros(2,2),zeros(2,2),zeros(2,(n-ii)*12);...
    zeros(1,4),Pc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),zeros(1,4),Pdw,zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),zeros(1,4),Fzong,zeros(1,2),zeros(1,(n-ii)*12);...
    PNc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),zeros(1,4),zeros(1,2),Ph1,zeros(1,(n-ii)*12);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),zeros(1,4),zeros(1,2),FN1,zeros(1,(n-ii)*12)];


mid_shaft_i=[Scu,zeros(1,4),zeros(1,2),zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2);...
    zeros(1,4),Sncv,zeros(1,2),zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2);...
    zeros(1,4),Scv,zeros(1,2),zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2),Swh1;...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2);...
    zeros(1,4),zeros(1,4),Sdw,zeros(1,2);...
    Swcu,zeros(1,4),zeros(1,2),zeros(1,2);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2)];

end