function Answer= Fc(c,d,SL,r)
%��ȡ��SL��r�����Ŷ������붩��������ֵc�Ĳ�ֵ
%c:����������ֵ
%d������������
%SL :����ˮƽ
%r :�����ۼ���  
%Answer         ���Ŷ������붩��������ֵc�Ĳ�ֵ
n=length(d);
B=0;
D(1)=d(1);
for ts=1:n-1               %ts--�۲�ĵڼ�����   n-��������
    if (d(ts)+sum(B))<c
        B=0;                                                    %�������ڶ���Ļ�����Ϊ0
    else
        B=syl([d(ts) B],c);%%�������Ķ����������ã�����ڵ���ĵ�������ǰһ���������--ʣ����
    end
    B_zs=B*((1+r).^(0:(length(B)-1)))';                      %���ݲ����ۼ��������ʣ������ֵ
    D(ts+1)=d(ts+1)+B_zs;
end
% [f,xi]=ksdensity(D);
% plot(xi,f,'k-','LineWidth',1,'MarkerSize',3);
% cdfplot(D);
% hold on
C=prctile(D,SL*100);  %�ɡ��ڶ�������������ʷֲ��õ������Ų���
Answer=C-c;
