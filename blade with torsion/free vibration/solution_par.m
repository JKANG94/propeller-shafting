
function [q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area,coe,G,MI,w,density)                                  


%梁横向振动中间参数
sigma=density*w^2/E;   tau=density*w^2/(coe*G) ;   alpha=density*Area*w^2/(E*MI);   

lambda=sqrt(sqrt((sigma-tau)^2/4+alpha)-(sigma+tau)/2);      hlambda=sqrt(sqrt((sigma-tau)^2/4+alpha)+(sigma+tau)/2);

q=w^2*density/(coe*G*lambda)+lambda;           hq=w^2*density/(coe*G*hlambda)-hlambda;

belta=q-lambda;       yita=hq+hlambda;

% 纵向振动中间参数
mue=w*sqrt(density/E);

% 扭转振动中间参数
gamma=w*sqrt(density/G);
end