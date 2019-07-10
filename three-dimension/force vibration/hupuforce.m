Nb=4;
NR=20;

% F1spe=force_skew(Nb);
load('F1spe.mat')

onetwo=1;

[BEAR1W,BEAR1V,BEAR2W,BEAR2V,BEAR3W,BEAR3V,BEAR3Z]=forcetransmatrix(onetwo);

onetwo=0;


[BEAR1W_3,BEAR1V_3,BEAR2W_3,BEAR2V_3,BEAR3W_3,BEAR3V_3,BEAR3Z_3]=forcetransmatrix(onetwo);

totalv1=[BEAR1V;-BEAR1V_3;-BEAR1V;BEAR1V_3];
totalw1=[BEAR1W;-BEAR1W_3;-BEAR1W;BEAR1W_3];

totalv2=[BEAR2V;-BEAR2V_3;-BEAR2V;BEAR2V_3];
totalw2=[BEAR2W;-BEAR2W_3;-BEAR2W;BEAR2W_3];

totalv3=[BEAR3V;-BEAR3V_3;-BEAR3V;BEAR3V_3];
totalw3=[BEAR3W;-BEAR3W_3;-BEAR3W;BEAR3W_3];
totalz3=[BEAR3Z;BEAR3Z_3;BEAR3Z;BEAR3Z_3];

RFR1_2W=zeros(300,1);
RFR1_3W=zeros(300,1);
RFR1_3WZ=zeros(300,1);

RFR2_3W=zeros(300,1);
RFR2_3WZ=zeros(300,1);

RFR3_3WZ=zeros(300,1);

for i=1:300
    
    forcespectrum=zeros(Nb*NR,Nb*NR);
    for j=0:Nb*NR-1
        forcespectrum(j+1,:)=F1spe((Nb*NR)*j+1:(Nb*NR)*(j+1),i);
    end
    
    %%%

    conju1_2w=totalw1(:,i)*totalw2(:,i)';
    conju1_3w=totalw1(:,i)*totalw3(:,i)';
    conju1_3wz=totalw1(:,i)*totalz3(:,i)';
    
   
    conju2_3w=totalw2(:,i)*totalw3(:,i)';
    conju2_3wz=totalw2(:,i)*totalz3(:,i)';
    

    conju3_3wz=totalw3(:,i)*totalz3(:,i)';
    
    %%%
    
    
    RFR1_2W(i)=sum(sum(forcespectrum.*conju1_2w));
    RFR1_3W(i)=sum(sum(forcespectrum.*conju1_3w));
    RFR1_3WZ(i)=sum(sum(forcespectrum.*conju1_3wz));
    
    
    RFR2_3W(i)=sum(sum(forcespectrum.*conju2_3w));
    RFR2_3WZ(i)=sum(sum(forcespectrum.*conju2_3wz));
    
    
    RFR3_3WZ(i)=sum(sum(forcespectrum.*conju3_3wz));

end

w=1:300;

figure(1) 
% semilogy(w,abs(RFR1W),'--b','LineWidth',2);       
% hold on
semilogy(w,abs(RFR1_2W),'r','LineWidth',2);


figure(2) 
% semilogy(w,abs(RFR2W),'--b','LineWidth',2);
% hold on
semilogy(w,abs(RFR1_3W),'r','LineWidth',2);

figure(3)
% semilogy(w,abs(RFR3W),'--b','LineWidth',2);
% hold on
semilogy(w,abs(RFR1_3WZ),'r','LineWidth',2);
% semilogy(w,abs(RFR3Z),'--b','LineWidth',2);
% hold on
