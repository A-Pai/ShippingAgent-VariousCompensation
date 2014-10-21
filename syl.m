function Answer= syl(d,c)
%求取承运量序列经过“当天”后的剩余量序列（遗留到第二天）--按照先进先出的原则，即存的时间越长的货物越先发
%d       当天的承运量序列，下标“1-n”分别是“当天的到货量-前（n-1）天的遗留量 ”              
%c                  当天的订舱量
%Answer         剩余量序列
B_ls=fliplr(cumsum(fliplr(d)));
B_ls=B_ls-c;
k=max(find(B_ls>0));
B0=B_ls(1:k);
Answer=[fliplr(diff(fliplr(B0))) B_ls(k)];
