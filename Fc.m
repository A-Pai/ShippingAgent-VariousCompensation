function Answer= Fc(c,d,SL,r)
%求取在SL、r下最优订舱量与订舱量尝试值c的差值
%c:订舱量尝试值
%d：到货量序列
%SL :服务水平
%r :补偿累加率  
%Answer         最优订舱量与订舱量尝试值c的差值
n=length(d);
B=0;
D(1)=d(1);
for ts=1:n-1               %ts--观察的第几天数   n-共多少天
    if (d(ts)+sum(B))<c
        B=0;                                                    %遗留到第二天的货物量为0
    else
        B=syl([d(ts) B],c);%%如果当天的订舱量不够用（针对于当天的到货量和前一天的遗留）--剩余量
    end
    B_zs=B*((1+r).^(0:(length(B)-1)))';                      %根据补偿累加率折算的剩余量总值
    D(ts+1)=d(ts+1)+B_zs;
end
% [f,xi]=ksdensity(D);
% plot(xi,f,'k-','LineWidth',1,'MarkerSize',3);
% cdfplot(D);
% hold on
C=prctile(D,SL*100);  %由“第二天承运量”概率分布得到的最优仓容
Answer=C-c;
