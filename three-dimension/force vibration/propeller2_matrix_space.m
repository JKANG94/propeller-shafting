function [propeller,mid_shaft_i]=propeller2_matrix_space(ku,kv,kw,kru,krv,b,h,L,q0w,hq0w,w,ii,n)


E=4.2e11*(1+0.01i);  coe=5/6; G=E/2/(1+0.32);       %����ģ������ľ���¼���ϵ�����б䵯��ģ��
Area=b*h;   density=9900;

%����1---��(u����

MI=h*b^3/12;

[q,lambda,hq,hlambda,belta,yita,~,~]=solution_par(E,Area,coe,G,MI,w,density);

Boundary_u=[coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)];...
    E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]];


Pcu=coe*G*Area*[0,belta,0,-yita]+ku*[1,0,1,0];

Pd0=[-ku,0];  


Pwcu=E*MI*[q*lambda,0,hq*hlambda,0]-krv*[0,q,0,-hq];

Pwc0w=[0,krv*q0w,0,-krv*hq0w];    %����q0z,hq0z
 
% ��Ҷ-->���� �м����

Swcu=krv*[0,q,0,-hq];   %����w����

Scu=ku*[1,0,1,0];


% ����2---�񶯣�v����

MI=b*h^3/12;
[q,lambda,hq,hlambda,belta,yita,mue,~]=solution_par(E,Area,coe,G,MI,w,density);   %������Ծؼ����뷽���й�

Boundary_v=[coe*G*Area*[belta*sinh(lambda*L),belta*cosh(lambda*L),yita*sin(hlambda*L),-yita*cos(hlambda*L)];...
    E*MI*[q*lambda*cosh(lambda*L),q*lambda*sinh(lambda*L),hq*hlambda*cos(hlambda*L),hq*hlambda*sin(hlambda*L)]];


Pcv=coe*G*Area*[0,belta,0,-yita]+kv*[1,0,1,0];

Pc0v=[-kv,0,-kv,0]; 


Pwcv=E*MI*[q*lambda,0,hq*hlambda,0]-kru*[0,q,0,-hq];


Pwh0u=[kru,0];

%������

Pc0w=[kw,0,kw,0];

Pdw=Area*E*[0, mue]-kw*[1,0];

Fzong=Area*E*[-mue*sin(mue*L),mue*cos(mue*L)];

% ��Ҷ-->���� �м����
%����V����
Scv=-kv*[1,0,1,0];
%����W����
Sdw=-kw*[1,0]; 

Sncv=kru*[0,q,0,-hq]; 

   
propeller=[zeros(1,4),zeros(1,4),Pd0,zeros(1,2),zeros(1,(ii-1)*10),Pcu,zeros(1,4),zeros(1,2),zeros(1,(n-ii)*10);...
zeros(1,4),Pwc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),Pwcu,zeros(1,4),zeros(1,2),zeros(1,(n-ii)*10);...
zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*10),Boundary_u,zeros(2,4),zeros(2,2),zeros(2,(n-ii)*10);...
Pc0v,zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),Pcv,zeros(1,2),zeros(1,(n-ii)*10);...
zeros(1,4),zeros(1,4),zeros(1,2),Pwh0u,zeros(1,(ii-1)*10),zeros(1,4),Pwcv,zeros(1,2),zeros(1,(n-ii)*10);...
zeros(2,4),zeros(2,4),zeros(2,2),zeros(2,2),zeros(2,(ii-1)*10),zeros(2,4),Boundary_v,zeros(2,2),zeros(2,(n-ii)*10);...
zeros(1,4),Pc0w,zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),zeros(1,4),Pdw,zeros(1,(n-ii)*10);...
zeros(1,4),zeros(1,4),zeros(1,2),zeros(1,2),zeros(1,(ii-1)*10),zeros(1,4),zeros(1,4),Fzong,zeros(1,(n-ii)*10)];
    
%�����ǽ�ҶŤת����������,������Ťת
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