function [propeller,mid_shaft_i]=propeller1_matrix_space(ku,kv,kw,kru,krv,b,h,L,q0w,hq0w,w,ii,n)


E=4.2e11*(1+0.01i);  coe=5/6; G=E/2/(1+0.32);       %弹性模量、铁木辛柯剪切系数、切变弹性模量
Area=b*h;  density=9900;
a=0;
%横向1---振动(u方向）
MI=h*b^3/12;

[q,lambda,hq,hlambda,belta,yita,~,~]=solution_par(E,Area,coe,G,MI,w,density);


Boundary_u=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

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
[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);   %截面惯性矩计算与方向有关

Boundary_v=[coe*G*Area*[0,belta,0,-yita];E*MI*[q*lambda,0,hq*hlambda,0]];

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

% 桨叶-->主轴 中间矩阵
Scv=-kv*[cosh(lambda*(L-a)),sinh(lambda*(L-a)),cos(hlambda*(L-a)),sin(hlambda*(L-a))];    %主轴V方向

Sdw=-kw*[cos(mue*L),sin(mue*L)];   %主轴W方向

Sncv=kru*[q*sinh(lambda*(L-a)),q*cosh(lambda*(L-a)),hq*sin(hlambda*(L-a)),-hq*cos(hlambda*(L-a))]; 


% Swhi3=[kr1w,0];
   
propeller=[zeros(1,4),zeros(1,4),Pd0,zeros(1,2),zeros(1,(ii-1)*10),Pcu,zeros(1,4),zeros(1,2),zeros(1,(n-ii)*10);...
zeros(1,4),Pwc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),Pwcu,zeros(1,4),zeros(1,2),zeros(1,(n-ii)*10);...
zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*10),Boundary_u,zeros(2,4),zeros(2,2),zeros(2,(n-ii)*10);...
Pc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),Pcv,zeros(1,2),zeros(1,(n-ii)*10);...
zeros(1,4),zeros(1,4),zeros(1,2),Pwh0u,zeros(1,(ii-1)*10),zeros(1,4),Pwcv,zeros(1,2),zeros(1,(n-ii)*10);...
zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*10),zeros(2,4),Boundary_v,zeros(2,2),zeros(2,(n-ii)*10);...
zeros(1,4),Pc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),zeros(1,4),Pdw,zeros(1,(n-ii)*10);...
zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),zeros(1,4),Fzong,zeros(1,(n-ii)*10)];
    
    %不考虑桨叶扭转的条件矩阵,主轴有扭转
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