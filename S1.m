%此程序研究加大补偿下的需求可延期满足的报童模型--不设置最后期限
%多期下只收单一价格货物（当天发送，推迟到第二天则提供赔偿）的最优仓容问题---总量
%遍历法确定最优仓容，迭代法的解与之比较，表明迭代法可行
%试验最优值的稳定性--不同的随机种子
%采用报童-迭代法(波尔查诺二分法)


clc
clear
close all

%% 参数设置
sj=10;     %产生随机数的次数--旧、新数据各sj次
n=5000;    %观测多少天
rate=1;       %运价
fare=0.6;   %仓容单位价格
p=0.2;      %每延后一天的基准补偿率
r=0.5;       %每延后一天的补偿累加率      补偿率=p*(1+r)^(k-1)     k为延后的天数
D1=0.2E4;     %试验仓容的最低值
D2=10.2E4;     %试验仓容的最高值,D1和D2的取值是根据到货的分布
Gap=1000;         %仓容每次的添加值
cs=(D2-D1)/Gap;    %试验的次数cs
c=D1+Gap:Gap:D2; %试验的仓容序列，与试验的次数对应
m=10;v=0.5;           %每天到来货物重量的lognrnd参数：均值和方差（其实是对于的正态分别的）
Co=fare;                      %供过于求的成本--针对于每单位的货物
Cu=rate-fare+p;          %供不应求的成本--针对于每单位的货物
SL=Cu/(Cu+Co);            %报童问题的服务水平

rng(1);      % 设置随机种子
d= lognrnd(m,v,sj*2,n);      %生成到货的货物总重量序列,共sj*2次(行），前sj次用于求解，后sj次用于验证

%% 遍历法和迭代法到最优仓容

for rseed=1:sj    %不同的随机种子下试验--数据不同
    %% 遍历法：所有的仓容（一定的歩距）去试验，得到最优仓容
    for cr=1:cs                %c(cr)-订舱量，cs-循环的次数
        B=0;                     %当天的遗留量
        for ts=1:n               %ts--观察的第几天数   n-共多少天
            if (d(rseed,ts)+sum(B))<c(cr)
                revenue1(cr,ts)=(d(rseed,ts)+sum(B))*rate-c(cr)*fare;    %以实际运送出去的货物量为准
                 revenue1_bb(cr,ts)=(d(rseed,ts)+sum(B))*rate-c(cr)*fare;  
                B=0;    %得到遗留给下一期的货物量，一定要放在计算revenue1语句的后边
            else
                %%遗留到第二天的货物量为：前一天的到货量+遗留到前一天的货物量-当前的仓容
                B=syl([d(rseed,ts) B],c(cr));%%如果当天的订舱量不够用（针对于当天的到货量和前一天的遗留）,syl()是函数-剩余量
                bc=B*(p*(1+r).^(0:(length(B)-1)))';       %bc--当天的补偿总额，遗留给下一期的货物量在本期补偿
                revenue1(cr,ts)=c(cr)*(rate-fare)-bc;
                revenue1_bb(cr,ts)=c(cr)*(rate-fare)-sum(B)*p;
            end
            if mod(ts,100)==0&mod(cr,100)==0
                fprintf('遍历求解，当前是第%d-%d-%d次，共%d-%d-%d次\n', ts,cr,rseed,n,cs,sj)
            end
        end
        revenue1_m(rseed,cr)=mean(revenue1(cr,:));                   %当前仓容量下的多期收益的平均值
        revenue1_m_bb(rseed,cr)=mean(revenue1_bb(cr,:));                   %当前仓容量下的多期收益的平均值        
    end
    [revenue1_m_max,id2]=max(revenue1_m,[],2);   %id2--遍历解，找到sj次随机产生数据的最大仓容值和位置
    [revenue1_m_bb_max,id2_bb]=max(revenue1_m_bb,[],2);   %id2--遍历解，找到sj次随机产生数据的最大仓容值和位置    
    
    
    
    
    
    %%  迭代法(波尔查诺二分法）----不是挨个去试
    
    ddcs=1;  %迭代次数
    a=D1;      %D1是明显的“不可能为最佳仓容”的较小值
    b=D2;      %D2是明显的“不可能为最佳仓容”的较大值
    c_dd(rseed,ddcs)=(a+b)/2;      %试验仓容
    Distance(rseed,ddcs)=Fc(c_dd(rseed,ddcs),d(rseed,:),SL,r);   %求取当前试验仓容c_dd(rseed,ddcs)、当前到货序列d(rseed,:)，服务水平SL下的C-c
    threshold=0.1;    %原始阈值
    times=1;          %设置阈值的次数
    while abs(Distance(rseed,ddcs))>threshold         %当当前的“C-c”还较大
        if  Fc(c_dd(rseed,ddcs),d(rseed,:),SL,r)*Fc(b,d(rseed,:),SL,r)<0
            a=c_dd(rseed,ddcs);
        else
            b=c_dd(rseed,ddcs);
        end
        c_dd(rseed,ddcs+1)=(a+b)/2;
        Distance(rseed,ddcs+1)=Fc(c_dd(rseed,ddcs+1),d(rseed,:),SL,r);
        ddcs=ddcs+1;
        if ddcs>20*times       %新的“阈值”下再次运行20次
            threshold=threshold+0.1;      %提升“阈值”5个单位
            times=times+1;
        end
    end
%     % '尝试解','报童解'接近情况图
%         figure
%         yx=(find(c_dd(rseed,:)~=0));       %c_dd中当前行中有效（yx）的列（有些列为0）
%         h1=plot(c_dd(rseed,yx),'k-o','LineWidth',1,'MarkerSize',3);
%         hold on
%         h2=plot(Distance(rseed,yx),'k:o','LineWidth',1,'MarkerSize',3);
%         grid on
%         legend([h1 h2],'订舱量','f(c)','Location','best')
%         title('迭代求解过程')
    
    C_dd(rseed)=c_dd(rseed,max(find(c_dd(rseed,:)>0))) ; %每次迭代实验的最优仓容值
    DistanceGg(rseed)=Distance(rseed,max(find(Distance(rseed,:)~=0))) ;     %得以“过关”的“C-c”
    
    
    
    %%%迭代解基于当前最优仓容下的收益
    B=0;
    for ts=1:n              %ts--观察的第几天数   n-共多少天
        if (d(rseed,ts)+sum(B))<C_dd(rseed)
            revenue2(ts)=(d(rseed,ts)+sum(B))*rate-C_dd(rseed)*fare;
            B=0;
        else
            B=syl([d(rseed,ts) B],C_dd(rseed));%%如果当天的订舱量不够用（针对于当天的到货量和前一天的遗留）
            bc=sum(B.*(p*(1+r).^(0:(length(B)-1))));       %bc--当天的补偿总额，遗留给下一期的货物量在本期补偿
            revenue2(ts)=C_dd(rseed)*(rate-fare)-bc;
        end
        
    end
    revenue2_m_max(rseed)=mean(revenue2);     %当前随机数据，迭代法最优仓容下的平均收益
    %         figure
    %         h1=plot(c_dd(rseed,find(c_dd(rseed,:)>0)),'k-o','LineWidth',1,'MarkerSize',3);
    %         hold on
    %         h2=plot(2:length(find(C(rseed,:)>0))+1,C(rseed,find(C(rseed,:)>0)),'r-o','LineWidth',1,'MarkerSize',3);
    %         grid on
    %         legend([h1 h2],'尝试解','报童解','迭----遍','Location','best')
    %     title('迭代求解过程')
    
    
    %     revenue1_m_max
    %     c(id2)
    %
    
end
%%全局最优值曲线
cr=10;   %指定第几次实验
figure
plot(c,revenue1_m(cr,:),'k-','LineWidth',2);
[r_max,I] = max(revenue1_m(cr,:))
hold on
h1=  plot(c(I),r_max,'k*','LineWidth',2,'MarkerSize',8);
xlim([0  11E+4])
xlabel('订舱量')
ylabel('平均收益')
legend(h1,'最大点')
grid on

figure
plot(c,revenue1_m(cr,:),'k-','LineWidth',2);
hold on
h1=  plot(c(I),r_max,'k*','LineWidth',2,'MarkerSize',8);
xlim([2.7E+4  4.4E+4])
xlabel('订舱量')
ylabel('平均收益')
legend(h1,'最大点')
grid on

%% C-c情况
figure
plot(DistanceGg,'k-o','LineWidth',1,'MarkerSize',8);
% title('f(c)')
grid on


%% 比较迭代法的仓容与遍历法的仓容
figure
h1=plot(C_dd,'k-o','LineWidth',1,'MarkerSize',3);
hold on
h2=plot(c(id2),'k:o','LineWidth',1,'MarkerSize',3);
% h3=plot(mean(d(1:sj,:),2),'k:o','LineWidth',2,'MarkerSize',3);
h=plot(C_dd-c(id2),'k-o','LineWidth',2,'MarkerSize',3);
legend([h1 h2 h],'算法解','遍历解','算法解-遍历解','Location','best')
% legend([h1 h2 h3 h],'迭代解','遍历解','到货量均值','迭代解--遍历解','Location','best')
% title('经验数据上的遍历解和迭代解-订舱量')
xlabel('第  次实验')
ylabel('订舱量')
grid on

figure
h1=plot(revenue2_m_max,'k-o','LineWidth',1,'MarkerSize',3);
hold on
h2=plot(revenue1_m_max,'k:o','LineWidth',1,'MarkerSize',3);
h=plot(revenue2_m_max-revenue1_m_max','k-o','LineWidth',2,'MarkerSize',3);
legend([h1 h2 h],'算法解','遍历解','算法解-遍历解','Location','best')
xlabel('第  次实验')
ylabel('平均收益')
% title('经验数据上最大收益')
grid on
%
% figure
% title('经验数据上的遍历解和迭代解')
% boxplot([(C_dd)', (c(id2))'])
cOpti1=mean(C_dd);    %迭代解
cOpti2=mean(c(id2));   %遍历解
%
% %% 收益随仓容变化曲线
%
% figure
% plot(c,revenue_m','LineWidth',1);
% grid on
% xlabel('订舱量');
% ylabel('平均收益');
% title('平均收益随仓容变化曲线')


%% 验证迭代法解的优良性，采用新的数据
%验证sj次数
for rseed=1:sj   %不同的随机种子下试验
%     rng(rseed+sj);      % 设置随机种子，“+sj”的目的是与上一部分的数据不同
%     d(rseed+sj,:)= lognrnd(m,v,1,n);      %生成当天的货物总重量
    B_1=0;                      %迭代解遗留到第二天的货物量
    B_2=0;                      %遍历解遗留到第二天的货物量
    revenueTest1=[];
    revenueTest2=[];
    for ts=1:n               %ts--观察的第几天数   n-共多少天
        %迭代解
        if (d(rseed+sj,ts)+sum(B_1))<cOpti1
            revenueTest1(ts)=(d(rseed+sj,ts)+sum(B_1))*rate-cOpti1*fare;
            B_1=0;    %遗留到第二天的货物量为0
        else
            B_1=syl([d(rseed+sj,ts) B_1],cOpti1);%%如果当天的订舱量不够用（针对于当天的到货量和前一天的遗留）
            bc=sum(B_1.*(p*(1+r).^(0:(length(B_1)-1))));       %bc--当天的补偿总额
            revenueTest1(ts)=cOpti1*(rate-fare)-bc;
        end
        
        %遍历解
        if (d(rseed+sj,ts)+sum(B_2))<cOpti2
            B_2=0;    %遗留到第二天的货物量为0
            revenueTest2(ts)=(d(rseed+sj,ts)+sum(B_2))*rate-cOpti2*fare;
        else
            B_2=syl([d(rseed,ts) B_2],cOpti2);%%如果当天的订舱量不够用（针对于当天的到货量和前一天的遗留）
            bc=sum(B_2.*(p*(1+r).^(0:(length(B_2)-1))));       %bc--当天的补偿总额
            revenueTest2(ts)=cOpti2*(rate-fare)-bc;
        end

        
        if mod(ts,20)==0
            fprintf('试验，当前是第%d-%d次，共%d-%d次\n', ts,rseed-sj,n,sj)
        end
    end
    revenueMtest1(rseed)=mean(revenueTest1);                   %迭代解的多期收益的平均值
    revenueMtest2(rseed)=mean(revenueTest2);                   %遍历解的多期收益的平均值
    
end
ratio=(revenueMtest1-revenueMtest2)./revenueMtest1;

figure
h1=plot(revenueMtest1,'k-o','LineWidth',1,'MarkerSize',3);
hold on
h2=plot(revenueMtest2,'k:o','LineWidth',1,'MarkerSize',3);
% h3=plot(revenue1_m_max,'r-o','LineWidth',1,'MarkerSize',3);
h=plot(revenueMtest1-revenueMtest2,'k-o','LineWidth',2,'MarkerSize',3);
legend([h1 h2 h],'算法解','遍历解','算法解-遍历解','Location','best')
xlabel('第  次实验')
ylabel('平均收益')
% legend([h1 h2 h h3],'迭代解','遍历解','迭---遍','训练解','Location','best')
% title('遍历解和迭代解在新数据上的平均收益')
grid on

figure
hist(revenueMtest1'- revenueMtest2')
title('验证数据上的迭代解和遍历解的收益差')

figure
subplot(1,2,1);
plot(d(1,:),'k-','LineWidth',1,'MarkerSize',3);
grid on
% title('到货量序列');
subplot(1,2,2);
[f,xi]=ksdensity(d(1,:));
plot(xi,f,'k-','LineWidth',1,'MarkerSize',3);
grid on
% title('到货量分布密度');
hold on


fprintf('迭代法的最优仓容平均为：%.0f，标准差为：%.0f，%.0f期的平均收益为：%.0f\n', cOpti1,std(C_dd),n,mean(revenueMtest1))
fprintf('遍历法的最优仓容平均为：%.0f，标准差为：%.0f，%.0f期的平均收益为：%.0f\n', cOpti2,std(c(id2)),n,mean(revenueMtest2))
