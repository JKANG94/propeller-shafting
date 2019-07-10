clear all
clc
w=2:0.01:500;     %  ԲƵ��

egn=zeros(length(w),1);   % ɨƵ�������������ֵ

parfor ii=1:length(w)
    
    matrix=condition_matrix_space(w(ii));
    
    egn(ii)=det(matrix);
end

result_fre=zeros(50,1);   % ģ̬Ƶ��
kk=1;

for jj=2:length(w)
    
    if egn(jj)>0 && egn(jj-1)<0
        (w(jj-1)+w(jj))/4/pi
        result_fre(kk)=(w(jj-1)+w(jj))/4/pi;
        kk=kk+1;
    elseif egn(jj)<0 && egn(jj-1)>0
        (w(jj-1)+w(jj))/4/pi
        result_fre(kk)=(w(jj-1)+w(jj))/4/pi;
        kk=kk+1;
    end
end