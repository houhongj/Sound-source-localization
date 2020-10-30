function [fun, df] = obfunction_value(x0, x, fun_syms, dfun_syms)
%目标函数值及其导数值
%minimize.m  f需返回两参数
%x0, 初始值
%x, 自变量
%fun_syms 目标函数

% n = length(x0);
% fun_syms = (1-x(1))^2 + 2*(x(2)^2)+ (x(3)^2)+(2.5-x(4))^2;%目标函数
% fun = (1-x{1})^2 + 2*(x{2}^2)+ (x{3}^2)+(2.5-x{4})^2;

fun = double(subs(fun_syms, x, x0));

% for i = 1 : n
%     dfun_syms(i) = diff(fun_syms, x(i));%各偏导数
% %     df(i) = diff(obf, x{i});
% end
df = double(subs(dfun_syms, x, x0));
df = real(df);
% aa=1;





end