function [T_2,T_1]=proper_trans(E,G,coe,MI,Area,q,lambda,hq,hlambda,belta,yita,a)

% ╪сть

laL=lambda*a; hlaL=hlambda*a;

T_2=[cosh(laL),sinh(laL),cos(hlaL),sin(hlaL); ...
    q*sinh(laL),q*cosh(laL),hq*sin(hlaL),-hq*cos(hlaL); ...
    E*MI*[q*lambda*cosh(laL),q*lambda*sinh(laL),hq*hlambda*cos(hlaL),hq*hlambda*sin(hlaL)];...
    coe*G*Area*[belta*sinh(laL),belta*cosh(laL),yita*sin(hlaL),-yita*cos(hlaL)]];

T_1=[1,0,1,0;0,q,0,-hq;...
    E*MI*q*lambda,0,E*MI*hq*hlambda,0;...
    0,coe*G*Area*belta,0,-coe*G*Area*yita];



end