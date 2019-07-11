function [propeller,mid_shaft_i,fu,fv]=first_propeller1_matrix_space(E,density,ue,ku,kv,kw,kru,krv,krw,b,h,L,q0w,hq0w,w,ii,n,a,ratio)


coe=5/6; G=E/2/(1+ue);       %桨叶铁木辛柯剪切系数、切变弹性模量,矩形截面
Area=b*h;

J=ratio*h*b^3;

% Forceu=1*cos(attack_ang);  %距离叶梢距离，单位幅值脉动力载荷 水动力攻角
Forceu=1;
%横向1---振动(u方向）

MI=h*b^3/12;

[q,lambda,hq,hlambda,belta,yita,~,~]=solution_par(E,Area,coe,G,MI,w,density);

% 加载导致的间断-传递矩阵法
if a~=0    % 加载点距离叶梢
    
    [T_2,T_1]=proper_trans(E,G,coe,MI,Area,q,lambda,hq,hlambda,belta,yita,a);
    
    trans1=T_2\T_1;
    Boundaryu_pro1=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];
    
    Boundary_u=Boundaryu_pro1*trans1;
    
    fu=-Boundaryu_pro1*(T_2\[zeros(3,1);Forceu]);
    
else
    Boundary_u=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];
    fu=-[Forceu;0];
    
end

Pcu=coe*G*Area*[belta*sinh(lambda*(L-a)),belta*cosh(lambda*(L-a)),yita*sin(hlambda*(L-a)),-yita*cos(hlambda*(L-a))]-...
    ku*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];

Pd0=[ku,0];

Pwcu=E*MI*[q*lambda*cosh(lambda*(L-a)),q*lambda*sinh(lambda*(L-a)),hq*hlambda*cos(hlambda*(L-a)),hq*hlambda*sin(hlambda*(L-a))]+...
    krv*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];

Pwc0w=[0,-krv*q0w,0,krv*hq0w];    %传参q0z,hq0z

% 轴系相关参数
Swcu=krv*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];   %主轴w方向

Scu=ku*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];

% 横向2---振动（v方向）

MI=b*h^3/12;
[q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area,coe,G,MI,w,density);   %截面惯性矩计算与方向有关

% Forcev=1*sin(attack_ang);
Forcev=1;

if a~=0
    
    [T_2,T_1]=proper_trans(E,G,coe,MI,Area,q,lambda,hq,hlambda,belta,yita,a);
    
    trans1=T_2\T_1;
    Boundaryv_pro1=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];
    
    Boundary_v=Boundaryv_pro1*trans1;
    
    fv=-Boundaryv_pro1*(T_2\[zeros(3,1);Forcev]);
else
    Boundary_v=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];
    fv=-[Forcev;0];
    
end

Pcv=coe*G*Area*[belta*sinh(lambda*(L-a)),belta*cosh(lambda*(L-a)),yita*sin(hlambda*(L-a)),-yita*cos(hlambda*(L-a))]-...
    kv*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];

Pc0v=[kv,0,kv,0];

Pwcv=E*MI*[q*lambda*cosh(lambda*(L-a)),q*lambda*sinh(lambda*(L-a)),hq*hlambda*cos(hlambda*(L-a)),hq*hlambda*sin(hlambda*(L-a))]+...
    kru*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];

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
Scv=-kv*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];    %主轴V方向
Swh1=krw*[cos(gamma*L),sin(gamma*L)];

Sdw=-kw*[cos(mue*L),sin(mue*L)];   %主轴W方向

Sncv=kru*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];


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