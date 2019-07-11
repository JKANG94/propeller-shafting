function [propeller,mid_shaft_i]=propeller2_matrix_space(E,density,ue,ku,kv,kw,kru,krv,krw,b,h,L,q0w,hq0w,w,ii,n,ratio)
% 返回桨叶2边界条件及其与主轴耦合条件矩阵
coe=5/6; G=E/2/(1+ue);       %弹性模量、铁木辛柯剪切系数、切变弹性模量
Area=b*h;

J=ratio*h*b^3;

% 横向1---振动(u方向）  xz平面

MI=h*b^3/12;

[q,lambda,hq,hlambda,belta,yita,~,~]=solution_par(E,Area,coe,G,MI,w,density);

Boundary_u=[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L);...
    q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)];

% Boundary_u=[coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)];...
%     E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]];

Pcu=coe*G*Area*[0,belta,0,-yita]+ku*[1,0,1,0];

Pd0=[-ku,0];


Pwcu=E*MI*[q*lambda,0,hq*hlambda,0]-krv*[0,q,0,-hq];

Pwc0w=[0,krv*q0w,0,-krv*hq0w];    %传参q0z,hq0z

% 桨叶-->主轴 耦合中间矩阵

Swcu=krv*[0,q,0,-hq];   %主轴w方向

Scu=ku*[1,0,1,0];

% 横向2---振动（v方向）yz平面

MI=b*h^3/12;
[q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area,coe,G,MI,w,density);   %截面惯性矩计算与方向有关

Boundary_v=[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L);...
    q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)];

% Boundary_v=[coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)];...
%     E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]];


Pcv=coe*G*Area*[0,belta,0,-yita]+kv*[1,0,1,0];

Pc0v=[-kv,0,-kv,0];

Pwcv=E*MI*[q*lambda,0,hq*hlambda,0]-kru*[0,q,0,-hq];


Pwh0u=[kru,0];

%纵向振动

Pc0w=[kw,0,kw,0];

Pdw=Area*E*[0, mue]-kw*[1,0];

Fzong=[-mue*sin(mue*L),mue*cos(mue*L)];
% Fzong=Area*E*[-mue*sin(mue*L),mue*cos(mue*L)];

%扭转振动
PNc0v=krw*[0,q0w,0,-hq0w];   %q0v=q0w;;;  hq0v=hq0w
Ph2=[0,G*J*gamma]-krw*[1,0];
FN2=[-gamma*sin(gamma*L),gamma*cos(gamma*L)];
% FN2=G*J*[-gamma*sin(gamma*L),gamma*cos(gamma*L)];

% 桨叶-->主轴 耦合中间矩阵
%主轴V方向
Scv=-kv*[1,0,1,0];
Swh2=krw*[1,0];

%主轴W方向
Sdw=-kw*[1,0];

Sncv=kru*[0,q,0,-hq];


% 桨叶边界条件
propeller=[zeros(1,4),zeros(1,4),Pd0,zeros(1,2),zeros(1,(ii-1)*12),Pcu,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(1,4),Pwc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),Pwcu,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*12),Boundary_u,zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(n-ii)*12);...
    Pc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),Pcv,zeros(1,2),zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(1,4),zeros(1,4),zeros(1,2),Pwh0u,zeros(1,(ii-1)*12),zeros(1,4),Pwcv,zeros(1,2),zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*12),zeros(2,4),Boundary_v,zeros(2,2),zeros(2,2),zeros(2,(n-ii)*12);...
    zeros(1,4),Pc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),zeros(1,4),Pdw,zeros(1,2),zeros(1,(n-ii)*12);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),zeros(1,4),Fzong,zeros(1,2),zeros(1,(n-ii)*12);...
    PNc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),zeros(1,4),zeros(1,2),Ph2,zeros(1,(n-ii)*12);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*12),zeros(1,4),zeros(1,4),zeros(1,2),FN2,zeros(1,(n-ii)*12)];


mid_shaft_i=[Scu,zeros(1,4),zeros(1,2),zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2);...
    zeros(1,4),Sncv,zeros(1,2),zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2);...
    zeros(1,4),Scv,zeros(1,2),zeros(1,2);...
    zeros(1,4),zeros(1,4),zeros(1,2),Swh2;...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2);...
    zeros(1,4),zeros(1,4),Sdw,zeros(1,2);...
    Swcu,zeros(1,4),zeros(1,2),zeros(1,2);...
    zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2)];

end