
function [q,lambda,hq,hlambda,belta,yita,mue,gamma]=solution_par(E,Area,coe,G,MI,w,density)                                  


%���������м����
sigma=density*w^2/E;   tau=density*w^2/(coe*G) ;   alpha=density*Area*w^2/(E*MI);   

lambda=sqrt(sqrt((sigma-tau)^2/4+alpha)-(sigma+tau)/2);      hlambda=sqrt(sqrt((sigma-tau)^2/4+alpha)+(sigma+tau)/2);

q=w^2*density/(coe*G*lambda)+lambda;           hq=w^2*density/(coe*G*hlambda)-hlambda;

belta=q-lambda;       yita=hq+hlambda;

% �������м����
mue=w*sqrt(density/E);

% Ťת���м����
gamma=w*sqrt(density/G);
end