function F1spe=force_skew(Nb)

% clear
% clc
%--------------------------------------------------------------------------
% �������������ϵ����������
%--------------------------------------------------------------------------

scale=1;
% Nb=4;   %��Ҷ�� ????????????????????????????????????????????????
Geob=load ('propgeo1.txt');   %�����������ļ�����������ÿһ��Ҷ���������ߵ��λ�ã�xyz����ÿ��������ҳ�C
NR=size(Geob,1);             %ÿһƬ��Ҷ�ķֶ���
dR=scale*(sqrt(Geob(2,2)^2+Geob(2,3)^2)-sqrt(Geob(1,2)^2+Geob(1,3)^2));
RT=scale*sqrt(Geob(NR,2)^2+Geob(NR,3)^2)+dR/2;                   %Ҷ�Ұ뾶
RH=0;     %��챰뾶
R=RT-RH; % ���İ뾶
thita0=90*pi/180;            %ҶƬ����ʼ�Ƕ�
%--------------------------------------------------------------------------
%��������h
%--------------------------------------------------------------------------
ns=1.9;               %��ת��-rpm
M=4*0.0254; %դ��ĳߴ磬��Ӣ��
lbda=2.5;   %
U=32.41;     %�����ٶ�
u=0.03*U;   %�����ȣ���
A=0.04;     %�������Ļ��ֳ߶�
bt=0.15;  %�ҳ���һ��
Vt=sqrt(U^2+(2*pi*ns*(RH+0.7*R))^2);  %
%--------------------------------------------------------------------------

numw=300;

Tf1=zeros(1,numw); Tf2=zeros(1,numw);  Co=zeros(1,numw);
% for k=1:numw
%     Tf1(1,k)=0;
%     Tf2(1,k)=0;
%     Co(1,k)=0;
% end

%--------------------------------------------------------------------------
Polb(:,1)=0*Geob(:,1);             %�������������������rake
Skew(:,1)=0*pi*Geob(:,5)/180;      %�����������Ĳ�б��skew
Polb(:,2)=scale*abs(Geob(:,2)+Geob(:,3)*1i); %��������İ뾶
Polb(:,3)=angle(Geob(:,2)+Geob(:,3)*1i)+thita0;   %���������thita��
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%����ҶƬ��������н�ҶҶ�������ĵ�����꣬�ܵĵ���Nb*NR
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% ���ü������������������λ��
%--------------------------------------------------------------------------
Pcor=zeros(Nb*NR,3);
Cbla=zeros(Nb*NR,1);
Ccor=zeros(Nb*NR,3);
thita=zeros(Nb*NR,1);
F1spe=zeros(Nb*NR*Nb*NR,numw);
Rn=zeros(1,64*1024);

for  m=1:Nb
    
    Thita_a=2*pi*(m-1)/Nb;             %���ÿƬ���ӵļн�
    
    for n=1:NR
        
        Pcor((m-1)*NR+n,1)=Polb(n,1);      %�������������ĵļ�����z,r,thita
        Pcor((m-1)*NR+n,2)=Polb(n,2);
        Pcor((m-1)*NR+n,3)=Polb(n,3)+Skew(n)+Thita_a;
        
        Cbla((m-1)*NR+n,1)=Geob(n,4);      %�������������ҳ�
        Ccor((m-1)*NR+n,1)=Polb(n,1);       %�������������ĵĵѿ�������
        Ccor((m-1)*NR+n,2)= Pcor((m-1)*NR+n,2)*cos(Pcor((m-1)*NR+n,3));
        Ccor((m-1)*NR+n,3)= Pcor((m-1)*NR+n,2)*sin(Pcor((m-1)*NR+n,3));
        thita((m-1)*NR+n,1)=Pcor((m-1)*NR+n,3);
        
    end
end

%--------------------------------------------------------------------------
%���һҶ�����������н�Ҷ������֮��Ļ���
% Nb*DivR
%--------------------------------------------------------------------------

for mm=1:Nb*NR
    
    ya = Ccor(mm,2);
    za = Ccor(mm,3);
    ra = Pcor(mm,2);
    oi = Pcor(mm,3);
    %       mm
    for nn=1:Nb*NR
        
        yb = Ccor(nn,2);
        zb = Ccor(nn,3);
        rb = Pcor(nn,2);
        oj = Pcor(nn,3);
        %--------------------------------------------------------------------------
        % ��������Ҷ�任�ĸ��ֲ���������Ƶ��һ����ò�Ҫ����64kHz
        %--------------------------------------------------------------------------
        NP=64;
        Fs = NP*1024;                    % Sampling frequency
        T = 1/Fs;                     % Sample time
        L = NP*1024;                   % Length of signal
        t = (0:L-1)*T;                % Time vector
        NFFT = 2^nextpow2(L); % Next power of 2 from length of y
%         f = Fs/2*linspace(0,1,NFFT/2+1);
        f=2:601;
        
        %--------------------------------------------------------------------------
        % �������֮����ٶ���غ���Rnm
        %--------------------------------------------------------------------------
        for ii=1:length(t)
            ojtao=oj+2*pi*ns*t(ii);
            dx=U*t(ii);                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dy=abs(rb*sin(ojtao)-ra*sin(oi));
            dz=abs(rb*cos(ojtao)-ra*cos(oi));
            r=sqrt(dx^2+dy^2+dz^2);
            
            if r == 0;
                Rxx = u^2;
                Ryy = u^2;
                Rzz = u^2;
                Rxy=0;
                Ryz=0;
                Rxz=0;
            else
                Rxx = u^2*(dx^2/(2*r*A)+(1-r/(2*A)))*exp(-r/A);
                Ryy = u^2*(dy^2/(2*r*A)+(1-r/(2*A)))*exp(-r/A);
                Rzz = u^2*(dz^2/(2*r*A)+(1-r/(2*A)))*exp(-r/A);
                
                
                Rxy = u^2*(dx*dy/(2*A*r))*exp(-r/A);
                Rxz = u^2*(dx*dz/(2*A*r))*exp(-r/A);
                Ryz = u^2*(dy*dz/(2*A*r))*exp(-r/A);
            end
            Ryx = Rxy;
            Rzx = Rxz;
            Rzy = Ryz;
            
            Rxo = Rxy*cos(ojtao) - Rxz*sin(ojtao);
            Rox = Ryx*cos(oi) - Rzx*sin(oi);
            Roo = Ryy*cos(oi)*cos(ojtao) + Rzz*sin(oi)*sin(ojtao) - Ryz*cos(oi)*sin(ojtao)- Rzy*sin(oi)*cos(ojtao);
            %     qa=atan(U/(2*pi*ns*0.7*R));
            %     qb=atan(U/(2*pi*ns*0.7*R));
            qa=atan(U/(2*pi*ns*ra));
            qb=atan(U/(2*pi*ns*rb));
            Rn(ii) = Rxx*cos(qa)*cos(qb) + Rxo*cos(qa)*sin(qb) + Rox*sin(qa)*cos(qb) +Roo*sin(qa)*sin(qb);
            %
            %   if mm==1
            %      if nn==1
            %       rr(ii,1)=dy;
            %       rr(ii,2)=dz;
            %       rr(ii,3)=r;
            %       rr(ii,4)=Rxx;
            %       rr(ii,5)=Ryy;
            %      end
            %  end
            
            
        end
        %--------------------------------------------------------------------------
        % Ȼ��������Ҷ�任�õ��ٶ����ܶ�
        %--------------------------------------------------------------------------
        
        Gn= 2*fft(Rn,NFFT)/L;
        
        
        Ri=sqrt(ya^2+za^2);
        Rj=sqrt(yb^2+zb^2);
        betai=atan(U/(2*pi*ns*0.7*R));
        betaj=atan(U/(2*pi*ns*0.7*R));
        betaT=atan(U/(2*pi*ns*0.7*R));
        Vi=sqrt(U^2+(2*pi*ns*Ri)^2);
        Vj=sqrt(U^2+(2*pi*ns*Rj)^2);
        q=sqrt((yb-ya)^2+(zb-za)^2);
        %--------------------------------------------------------------------------
        % ����sevik�ķ�����õ�������
        %--------------------------------------------------------------------------
%         for k=1:numw
%             ki=2*pi*k*bt/Vi;
%             H0=besselh(0,2,ki);
%             H1=besselh(1,2,ki);
%             Lei=(2*1i/pi)/(ki*(H1+1i*H0));      %���Ķ�άsears����
%             
%             kj=2*pi*k*bt/Vj;
%             H0=besselh(0,2,kj);
%             H1=besselh(1,2,kj);
%             Lej=(2*1i/pi)/(kj*(H1+1i*H0));    %���Ķ�άsears����
%             
%             omega1=2*pi*(k-1);
%             HG11(k)=u^2*2*((lbda*U/M-1i*omega1)/((lbda*U/M)^2+omega1^2))*exp(-lbda*q/M);
%             HG22(k)=u^2*2*((lbda*U/M-1i*omega1)/((lbda*U/M)^2+omega1^2))*exp(-lbda*q/M)*(1-0.5*lbda*U/M*((lbda*U/M-i*omega1)/((lbda*U/M)^2+omega1^2)));
%             QQ=omega1*bt/Vt;
%             Hsquare=(2*pi*1000*Vt*bt)^2*dR^2/(1+2*pi*QQ);
%             spec2=Hsquare*cos(betaT)^2*cos(betaT)^2*HG11(k)+0.25*Hsquare*sin(2*betaT)^2*sin(2*betaT)^2*HG22(k)*0.1;
%             F1spe_R((mm-1)*Nb*NR+nn,k)=real(spec2/(cos(betai)*cos(betaj)));
%             F1spe_I((mm-1)*Nb*NR+nn,k)=imag(spec2/(cos(betai)*cos(betaj)));
%             F2spe_R((mm-1)*Nb*NR+nn,k)=real(spec2/(sin(betai)*sin(betaj)));
%             F2spe_I((mm-1)*Nb*NR+nn,k)=imag(spec2/(sin(betai)*sin(betaj)));
%             Tf2(1,k)=Tf2(1,k)+ spec2;
        
        
        
        
        %    Vi=Vt;
        %    Vj=Vt;
        %    betai=betaT;
        %    betaj=betaT;
        %--------------------------------------------------------------------------
        % ����jiang�ķ�����õ�������
        %--------------------------------------------------------------------------
        for k=1:numw
            ki=2*pi*f(k)*bt/Vi;
            H0=besselh(0,2,ki);
            H1=besselh(1,2,ki);
            Lei=(2*1i/pi)/(ki*(H1+1i*H0));      %���Ķ�άsears����
%             �Ե�i����������������ά����
%             if  TDflag==1
%                 Ls=L3d(ki/bt,bigc,deltaR/(2*bt));
%             else
%                 Ls=1;
%             end
%             Lei=Lei*Ls;
            
            kj=2*pi*f(k)*bt/Vj;
            H0=besselh(0,2,kj);
            H1=besselh(1,2,kj);
            Lej=(2*1i/pi)/(kj*(H1+1i*H0));    %���Ķ�άsears����
%             �Ե�j����������������ά����
%             if  TDflag==1
%                 Ls=L3d(kj/bt,bigc,deltaR/(2*bt));
%             else
%                 Ls=1;
%             end
%             Lej=Lej*Ls;
            
            Hsquare=(2*pi*1000*Vi*bt)*(2*pi*1000*Vj*bt)*dR^2*conj(Lei)*Lej;
            spec1=Hsquare*Gn(k)*0.1;
            Tf1(1,k)=Tf1(1,k)+ spec1;
            %spec1=Hsquare*Gn(k);       %��Ϊ����Ԫ����ʱ������ֱ�Ӹ����˽�Ҷ��λ�ã����Բ��ó��ԣ���
            
            %�����к��е�˳��д�����׵ľ���
            
            %Fspe_R((mm-1)*Nb*DivR+nn,k)=abs(spec1/(cos(betai)*cos(betaj)));
            %Vspe((mm-1)*Nb*DivR+nn,k)=abs(core);
            F1spe((mm-1)*Nb*NR+nn,k)=spec1;

        end
    end
end

% ����õ������׾���ʵ�����鲿д���ļ�����abaqus��ȡ
%-------------------------------------------------------------------------


% fid3=fopen('PSd.txt','w');  % ���ļ�
% for k1=1:numw
%     fprintf(fid3,'%e\n',Tf1(1,k1));
% end
% fclose(fid3);
% for k=1:numw
%     filename1=strcat('veol_R',num2str(k),'.txt');
%     fid1=fopen(filename1,'a+');  % ���ļ�
%     indicator=0;
%     for (mm=1:Nb*DivR)
%         for (nn=1:Nb*DivR)
%             indicator=indicator+1;
%             fprintf(fid1,'%e\n',Vspe(indicator,k));
%         end
%     end
%     fclose(fid1);
% end

%-------------------------------------------------------------------------
%

% plot(10*log10(real(Tf1)));
% hold on
% plot(10*log10(real(Tf2)),'red');
% xlim([0 numw])
% Co=abs(Co');




