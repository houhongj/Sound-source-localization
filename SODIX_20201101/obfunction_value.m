function [fun, df] = obfunction_value(x0, x, fun_syms, dfun_syms)
%Ŀ�꺯��ֵ���䵼��ֵ
%minimize.m  f�践��������
%x0, ��ʼֵ
%x, �Ա���
%fun_syms Ŀ�꺯��

% n = length(x0);
% fun_syms = (1-x(1))^2 + 2*(x(2)^2)+ (x(3)^2)+(2.5-x(4))^2;%Ŀ�꺯��
% fun = (1-x{1})^2 + 2*(x{2}^2)+ (x{3}^2)+(2.5-x{4})^2;

fun = double(subs(fun_syms, x, x0));

% for i = 1 : n
%     dfun_syms(i) = diff(fun_syms, x(i));%��ƫ����
% %     df(i) = diff(obf, x{i});
% end
df = double(subs(dfun_syms, x, x0));
df = real(df);
% aa=1;





end