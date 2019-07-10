function [propeller,mid_shaft_i,fu,fw]=propeller3_matrix_space(ku,kv,kw,kru,krw,b,h,L,q0v,hq0v,w,ii,n,a,attack_ang)

E=4.2e11*(1+0.01i);  coe=5/6; G=E/2/(1+0.32);       %弹性模量、铁木辛柯剪切系数、切变弹性模量
Area=b*h;  density=9900;

%横向1---振动(u方向）
Forceu=1*cos(attack_ang);
% Forceu=1;
MI=h*b^3/12;

[q,lambda,hq,hlambda,belta,yita,~,~]=solution_par(E,Area,coe,G,MI,w,density);

if a~=0

    [T_2,T_1]=proper_trans(E,G,coe,MI,Area,q,lambda,hq,hlambda,belta,yita,a);

    trans1=T_2\T_1;
    Boundaryu_pro1=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

    Boundary_u=Boundaryu_pro1*trans1;

    fu=-Boundaryu_pro1*(T_2\[zeros(3,1);Forceu]);
      
else
    Boundary_u=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];
    fu=-[Forceu;0];

end

% Boundary_u=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

Pcu=coe*G*Area*[belta*sinh(lambda*(L-a)),belta*cosh(lambda*(L-a)),yita*sin(hlambda*(L-a)),-yita*cos(hlambda*(L-a))]-...
    ku*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];

Pd0=[ku,0];  

Pwcu=E*MI*[q*lambda*cosh(lambda*(L-a)),q*lambda*sinh(lambda*(L-a)),hq*hlambda*cos(hlambda*(L-a)),hq*hlambda*sin(hlambda*(L-a))]+...
    krw*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];

Pwc0v=[0,-krw*q0v,0,krw*hq0v];    %传参q0z,hq0z
 

% 轴系相关参数
Swcu=krw*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];   %主轴V方向

Scu=ku*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];

% 横向2---振动（w方向）
Forcew=1*sin(attack_ang);
% Forcew=1;
MI=b*h^3/12;
[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);   %截面惯性矩计算与方向有关

if a~=0

    [T_2,T_1]=proper_trans(E,G,coe,MI,Area,q,lambda,hq,hlambda,belta,yita,a);

    trans1=T_2\T_1;
    Boundaryw_pro1=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

    Boundary_w=Boundaryw_pro1*trans1;

    fw=-Boundaryw_pro1*(T_2\[zeros(3,1);Forcew]);
      
else
    Boundary_w=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];
    fw=-[Forcew;0];

end


% Boundary_w=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

Pcw=coe*G*Area*[belta*sinh(lambda*(L-a)),belta*cosh(lambda*(L-a)),yita*sin(hlambda*(L-a)),-yita*cos(hlambda*(L-a))]-...
    kw*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];

Pc0w=[kw,0,kw,0]; 

Pwcw=E*MI*[q*lambda*cosh(lambda*(L-a)),q*lambda*sinh(lambda*(L-a)),hq*hlambda*cos(hlambda*(L-a)),hq*hlambda*sin(hlambda*(L-a))]+...
    kru*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))];

Pwh0u=[-kru,0];

%纵向振动
Pc0v=[-kv,0,-kv,0];
Pdv=Area*E*[-mue*sin(mue*L), mue*cos(mue*L)]+kv*[cos(mue*L),sin(mue*L)];
Fzong=Area*E*[0,mue];

% 桨叶-->主轴 中间矩阵
Sdv=-kv*[cos(mue*L),sin(mue*L)];   %主轴V方向

Scw=-kw*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];   %主轴w方向

Sncw=kru*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))]; 


% Swhi3=[kr1w,0];
   
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