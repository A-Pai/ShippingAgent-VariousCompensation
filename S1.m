%�˳����о��Ӵ󲹳��µ��������������ı�ͯģ��--�������������
%������ֻ�յ�һ�۸������췢�ͣ��Ƴٵ��ڶ������ṩ�⳥�������Ų�������---����
%������ȷ�����Ų��ݣ��������Ľ���֮�Ƚϣ���������������
%��������ֵ���ȶ���--��ͬ���������
%���ñ�ͯ-������(������ŵ���ַ�)


clc
clear
close all

%% ��������
sj=10;     %����������Ĵ���--�ɡ������ݸ�sj��
n=5000;    %�۲������
rate=1;       %�˼�
fare=0.6;   %���ݵ�λ�۸�
p=0.2;      %ÿ�Ӻ�һ��Ļ�׼������
r=0.5;       %ÿ�Ӻ�һ��Ĳ����ۼ���      ������=p*(1+r)^(k-1)     kΪ�Ӻ������
D1=0.2E4;     %������ݵ����ֵ
D2=10.2E4;     %������ݵ����ֵ,D1��D2��ȡֵ�Ǹ��ݵ����ķֲ�
Gap=1000;         %����ÿ�ε����ֵ
cs=(D2-D1)/Gap;    %����Ĵ���cs
c=D1+Gap:Gap:D2; %����Ĳ������У�������Ĵ�����Ӧ
m=10;v=0.5;           %ÿ�쵽������������lognrnd��������ֵ�ͷ����ʵ�Ƕ��ڵ���̬�ֱ�ģ�
Co=fare;                      %��������ĳɱ�--�����ÿ��λ�Ļ���
Cu=rate-fare+p;          %����Ӧ��ĳɱ�--�����ÿ��λ�Ļ���
SL=Cu/(Cu+Co);            %��ͯ����ķ���ˮƽ

rng(1);      % �����������
d= lognrnd(m,v,sj*2,n);      %���ɵ����Ļ�������������,��sj*2��(�У���ǰsj��������⣬��sj��������֤

%% �������͵����������Ų���

for rseed=1:sj    %��ͬ���������������--���ݲ�ͬ
    %% �����������еĲ��ݣ�һ���Ěi�ࣩȥ���飬�õ����Ų���
    for cr=1:cs                %c(cr)-��������cs-ѭ���Ĵ���
        B=0;                     %�����������
        for ts=1:n               %ts--�۲�ĵڼ�����   n-��������
            if (d(rseed,ts)+sum(B))<c(cr)
                revenue1(cr,ts)=(d(rseed,ts)+sum(B))*rate-c(cr)*fare;    %��ʵ�����ͳ�ȥ�Ļ�����Ϊ׼
                 revenue1_bb(cr,ts)=(d(rseed,ts)+sum(B))*rate-c(cr)*fare;  
                B=0;    %�õ���������һ�ڵĻ�������һ��Ҫ���ڼ���revenue1���ĺ��
            else
                %%�������ڶ���Ļ�����Ϊ��ǰһ��ĵ�����+������ǰһ��Ļ�����-��ǰ�Ĳ���
                B=syl([d(rseed,ts) B],c(cr));%%�������Ķ����������ã�����ڵ���ĵ�������ǰһ���������,syl()�Ǻ���-ʣ����
                bc=B*(p*(1+r).^(0:(length(B)-1)))';       %bc--����Ĳ����ܶ��������һ�ڵĻ������ڱ��ڲ���
                revenue1(cr,ts)=c(cr)*(rate-fare)-bc;
                revenue1_bb(cr,ts)=c(cr)*(rate-fare)-sum(B)*p;
            end
            if mod(ts,100)==0&mod(cr,100)==0
                fprintf('������⣬��ǰ�ǵ�%d-%d-%d�Σ���%d-%d-%d��\n', ts,cr,rseed,n,cs,sj)
            end
        end
        revenue1_m(rseed,cr)=mean(revenue1(cr,:));                   %��ǰ�������µĶ��������ƽ��ֵ
        revenue1_m_bb(rseed,cr)=mean(revenue1_bb(cr,:));                   %��ǰ�������µĶ��������ƽ��ֵ        
    end
    [revenue1_m_max,id2]=max(revenue1_m,[],2);   %id2--�����⣬�ҵ�sj������������ݵ�������ֵ��λ��
    [revenue1_m_bb_max,id2_bb]=max(revenue1_m_bb,[],2);   %id2--�����⣬�ҵ�sj������������ݵ�������ֵ��λ��    
    
    
    
    
    
    %%  ������(������ŵ���ַ���----���ǰ���ȥ��
    
    ddcs=1;  %��������
    a=D1;      %D1�����Եġ�������Ϊ��Ѳ��ݡ��Ľ�Сֵ
    b=D2;      %D2�����Եġ�������Ϊ��Ѳ��ݡ��Ľϴ�ֵ
    c_dd(rseed,ddcs)=(a+b)/2;      %�������
    Distance(rseed,ddcs)=Fc(c_dd(rseed,ddcs),d(rseed,:),SL,r);   %��ȡ��ǰ�������c_dd(rseed,ddcs)����ǰ��������d(rseed,:)������ˮƽSL�µ�C-c
    threshold=0.1;    %ԭʼ��ֵ
    times=1;          %������ֵ�Ĵ���
    while abs(Distance(rseed,ddcs))>threshold         %����ǰ�ġ�C-c�����ϴ�
        if  Fc(c_dd(rseed,ddcs),d(rseed,:),SL,r)*Fc(b,d(rseed,:),SL,r)<0
            a=c_dd(rseed,ddcs);
        else
            b=c_dd(rseed,ddcs);
        end
        c_dd(rseed,ddcs+1)=(a+b)/2;
        Distance(rseed,ddcs+1)=Fc(c_dd(rseed,ddcs+1),d(rseed,:),SL,r);
        ddcs=ddcs+1;
        if ddcs>20*times       %�µġ���ֵ�����ٴ�����20��
            threshold=threshold+0.1;      %��������ֵ��5����λ
            times=times+1;
        end
    end
%     % '���Խ�','��ͯ��'�ӽ����ͼ
%         figure
%         yx=(find(c_dd(rseed,:)~=0));       %c_dd�е�ǰ������Ч��yx�����У���Щ��Ϊ0��
%         h1=plot(c_dd(rseed,yx),'k-o','LineWidth',1,'MarkerSize',3);
%         hold on
%         h2=plot(Distance(rseed,yx),'k:o','LineWidth',1,'MarkerSize',3);
%         grid on
%         legend([h1 h2],'������','f(c)','Location','best')
%         title('����������')
    
    C_dd(rseed)=c_dd(rseed,max(find(c_dd(rseed,:)>0))) ; %ÿ�ε���ʵ������Ų���ֵ
    DistanceGg(rseed)=Distance(rseed,max(find(Distance(rseed,:)~=0))) ;     %���ԡ����ء��ġ�C-c��
    
    
    
    %%%��������ڵ�ǰ���Ų����µ�����
    B=0;
    for ts=1:n              %ts--�۲�ĵڼ�����   n-��������
        if (d(rseed,ts)+sum(B))<C_dd(rseed)
            revenue2(ts)=(d(rseed,ts)+sum(B))*rate-C_dd(rseed)*fare;
            B=0;
        else
            B=syl([d(rseed,ts) B],C_dd(rseed));%%�������Ķ����������ã�����ڵ���ĵ�������ǰһ���������
            bc=sum(B.*(p*(1+r).^(0:(length(B)-1))));       %bc--����Ĳ����ܶ��������һ�ڵĻ������ڱ��ڲ���
            revenue2(ts)=C_dd(rseed)*(rate-fare)-bc;
        end
        
    end
    revenue2_m_max(rseed)=mean(revenue2);     %��ǰ������ݣ����������Ų����µ�ƽ������
    %         figure
    %         h1=plot(c_dd(rseed,find(c_dd(rseed,:)>0)),'k-o','LineWidth',1,'MarkerSize',3);
    %         hold on
    %         h2=plot(2:length(find(C(rseed,:)>0))+1,C(rseed,find(C(rseed,:)>0)),'r-o','LineWidth',1,'MarkerSize',3);
    %         grid on
    %         legend([h1 h2],'���Խ�','��ͯ��','��----��','Location','best')
    %     title('����������')
    
    
    %     revenue1_m_max
    %     c(id2)
    %
    
end
%%ȫ������ֵ����
cr=10;   %ָ���ڼ���ʵ��
figure
plot(c,revenue1_m(cr,:),'k-','LineWidth',2);
[r_max,I] = max(revenue1_m(cr,:))
hold on
h1=  plot(c(I),r_max,'k*','LineWidth',2,'MarkerSize',8);
xlim([0  11E+4])
xlabel('������')
ylabel('ƽ������')
legend(h1,'����')
grid on

figure
plot(c,revenue1_m(cr,:),'k-','LineWidth',2);
hold on
h1=  plot(c(I),r_max,'k*','LineWidth',2,'MarkerSize',8);
xlim([2.7E+4  4.4E+4])
xlabel('������')
ylabel('ƽ������')
legend(h1,'����')
grid on

%% C-c���
figure
plot(DistanceGg,'k-o','LineWidth',1,'MarkerSize',8);
% title('f(c)')
grid on


%% �Ƚϵ������Ĳ�����������Ĳ���
figure
h1=plot(C_dd,'k-o','LineWidth',1,'MarkerSize',3);
hold on
h2=plot(c(id2),'k:o','LineWidth',1,'MarkerSize',3);
% h3=plot(mean(d(1:sj,:),2),'k:o','LineWidth',2,'MarkerSize',3);
h=plot(C_dd-c(id2),'k-o','LineWidth',2,'MarkerSize',3);
legend([h1 h2 h],'�㷨��','������','�㷨��-������','Location','best')
% legend([h1 h2 h3 h],'������','������','��������ֵ','������--������','Location','best')
% title('���������ϵı�����͵�����-������')
xlabel('��  ��ʵ��')
ylabel('������')
grid on

figure
h1=plot(revenue2_m_max,'k-o','LineWidth',1,'MarkerSize',3);
hold on
h2=plot(revenue1_m_max,'k:o','LineWidth',1,'MarkerSize',3);
h=plot(revenue2_m_max-revenue1_m_max','k-o','LineWidth',2,'MarkerSize',3);
legend([h1 h2 h],'�㷨��','������','�㷨��-������','Location','best')
xlabel('��  ��ʵ��')
ylabel('ƽ������')
% title('�����������������')
grid on
%
% figure
% title('���������ϵı�����͵�����')
% boxplot([(C_dd)', (c(id2))'])
cOpti1=mean(C_dd);    %������
cOpti2=mean(c(id2));   %������
%
% %% ��������ݱ仯����
%
% figure
% plot(c,revenue_m','LineWidth',1);
% grid on
% xlabel('������');
% ylabel('ƽ������');
% title('ƽ����������ݱ仯����')


%% ��֤��������������ԣ������µ�����
%��֤sj����
for rseed=1:sj   %��ͬ���������������
%     rng(rseed+sj);      % ����������ӣ���+sj����Ŀ��������һ���ֵ����ݲ�ͬ
%     d(rseed+sj,:)= lognrnd(m,v,1,n);      %���ɵ���Ļ���������
    B_1=0;                      %�������������ڶ���Ļ�����
    B_2=0;                      %�������������ڶ���Ļ�����
    revenueTest1=[];
    revenueTest2=[];
    for ts=1:n               %ts--�۲�ĵڼ�����   n-��������
        %������
        if (d(rseed+sj,ts)+sum(B_1))<cOpti1
            revenueTest1(ts)=(d(rseed+sj,ts)+sum(B_1))*rate-cOpti1*fare;
            B_1=0;    %�������ڶ���Ļ�����Ϊ0
        else
            B_1=syl([d(rseed+sj,ts) B_1],cOpti1);%%�������Ķ����������ã�����ڵ���ĵ�������ǰһ���������
            bc=sum(B_1.*(p*(1+r).^(0:(length(B_1)-1))));       %bc--����Ĳ����ܶ�
            revenueTest1(ts)=cOpti1*(rate-fare)-bc;
        end
        
        %������
        if (d(rseed+sj,ts)+sum(B_2))<cOpti2
            B_2=0;    %�������ڶ���Ļ�����Ϊ0
            revenueTest2(ts)=(d(rseed+sj,ts)+sum(B_2))*rate-cOpti2*fare;
        else
            B_2=syl([d(rseed,ts) B_2],cOpti2);%%�������Ķ����������ã�����ڵ���ĵ�������ǰһ���������
            bc=sum(B_2.*(p*(1+r).^(0:(length(B_2)-1))));       %bc--����Ĳ����ܶ�
            revenueTest2(ts)=cOpti2*(rate-fare)-bc;
        end

        
        if mod(ts,20)==0
            fprintf('���飬��ǰ�ǵ�%d-%d�Σ���%d-%d��\n', ts,rseed-sj,n,sj)
        end
    end
    revenueMtest1(rseed)=mean(revenueTest1);                   %������Ķ��������ƽ��ֵ
    revenueMtest2(rseed)=mean(revenueTest2);                   %������Ķ��������ƽ��ֵ
    
end
ratio=(revenueMtest1-revenueMtest2)./revenueMtest1;

figure
h1=plot(revenueMtest1,'k-o','LineWidth',1,'MarkerSize',3);
hold on
h2=plot(revenueMtest2,'k:o','LineWidth',1,'MarkerSize',3);
% h3=plot(revenue1_m_max,'r-o','LineWidth',1,'MarkerSize',3);
h=plot(revenueMtest1-revenueMtest2,'k-o','LineWidth',2,'MarkerSize',3);
legend([h1 h2 h],'�㷨��','������','�㷨��-������','Location','best')
xlabel('��  ��ʵ��')
ylabel('ƽ������')
% legend([h1 h2 h h3],'������','������','��---��','ѵ����','Location','best')
% title('������͵��������������ϵ�ƽ������')
grid on

figure
hist(revenueMtest1'- revenueMtest2')
title('��֤�����ϵĵ�����ͱ�����������')

figure
subplot(1,2,1);
plot(d(1,:),'k-','LineWidth',1,'MarkerSize',3);
grid on
% title('����������');
subplot(1,2,2);
[f,xi]=ksdensity(d(1,:));
plot(xi,f,'k-','LineWidth',1,'MarkerSize',3);
grid on
% title('�������ֲ��ܶ�');
hold on


fprintf('�����������Ų���ƽ��Ϊ��%.0f����׼��Ϊ��%.0f��%.0f�ڵ�ƽ������Ϊ��%.0f\n', cOpti1,std(C_dd),n,mean(revenueMtest1))
fprintf('�����������Ų���ƽ��Ϊ��%.0f����׼��Ϊ��%.0f��%.0f�ڵ�ƽ������Ϊ��%.0f\n', cOpti2,std(c(id2)),n,mean(revenueMtest2))
