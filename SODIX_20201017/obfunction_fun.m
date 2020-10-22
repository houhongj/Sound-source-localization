function [fun, df] = obfunction_fun(x0,fun_syms, dfun_syms,x)
%目标函数值及其导数值



fun = double(subs(fun_syms, x, x0));

df = double(subs(dfun_syms, x, x0));
df = real(df);

aa=1;





end