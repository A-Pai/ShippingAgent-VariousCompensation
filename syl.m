function Answer= syl(d,c)
%��ȡ���������о��������족���ʣ�������У��������ڶ��죩--�����Ƚ��ȳ���ԭ�򣬼����ʱ��Խ���Ļ���Խ�ȷ�
%d       ����ĳ��������У��±ꡰ1-n���ֱ��ǡ�����ĵ�����-ǰ��n-1����������� ��              
%c                  ����Ķ�����
%Answer         ʣ��������
B_ls=fliplr(cumsum(fliplr(d)));
B_ls=B_ls-c;
k=max(find(B_ls>0));
B0=B_ls(1:k);
Answer=[fliplr(diff(fliplr(B0))) B_ls(k)];
