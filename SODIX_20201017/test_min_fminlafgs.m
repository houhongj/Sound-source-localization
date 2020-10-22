% 测试Lbfgs hhj202010

close all;
clc
clear all
addpath('../rodyo-FEX-minimize-2f6ee87');
%% 
% x0 = [0];
% sym x
% options = setoptimoptions('Algorithm','fminlbfgs');
% options.PlotFcn = @optimplotfval;%fminsearch  fminlbfgs
% options.Display = 'plot';%  ,'Display','plot' 
% a = minimize(@obfunction_value1,x0,[],[],[],[],[],[],[],options)
% a = minimize(@obfunction_value1,x0)
% function [f] = obfunction_value1(x)
% f =  sin(x) + 3;
% end
%% 
% 'GradObj','on'与[f,g] = obfunction_value1(x)保持一致
% x0 = [0];
% sym x
% options = setoptimoptions('Algorithm','fminlbfgs','GradObj','on');
% a = minimize(@obfunction_value1,x0,[],[],[],[],[],[],[],options)
% % a = minimize(@obfunction_value1,x0)
% function [f,g] = obfunction_value1(x)
% f =  sin(x) + 3;
% g = cos(x);
% end
%%
% %上下约束、注意起始点在范围内
% %输出的参数
% x0 = [-0.5];
% sym x
% % ('AlwaysHonorConstraints', 'none');    %bounds
% options = setoptimoptions('Algorithm','fminlbfgs','GradObj','on','AlwaysHonorConstraints', 'bounds');%fminsearch  fminlbfgs
% [sol, fval, exitflag, output, grad] = minimize(@obfunction_value1,x0,[],[],[],[],[-1],[0],[],options)
% % a = minimize(@obfunction_value1,x0)
% function [f,g] = obfunction_value1(x)%会多次调用，还是另用函数
% f =  sin(x) + 3;
% g = cos(x);
% end
%% 
% 多元函数
x0 = [0 0 0 0];
x=sym('x',[1,4]);
y=sym('y',[1,4]);
% ('AlwaysHonorConstraints', 'none');    %bounds

options = setoptimoptions('Algorithm','fminlbfgs','GradObj','on','AlwaysHonorConstraints', 'bounds'...
    ,'Display','iter','MaxIter',5 );%fminsearch  fminlbfgs,    'PlotFcn',@optimplotfval
% options.Display = 'plot';%  
b=1;
c=[2,1;2,4];
[sol, fval, exitflag, output, grad] = minimize(@obfunction_value1,x0,[],[],[],[],[],[],[],options,b,c,y)
% [sol, fval, exitflag, output, grad] = minimize(@obfunction_value1,x0,[],[],[],[],[-10 -10 -10 -10],[2 2 2 2],[],options)
% a = minimize(@obfunction_value1,x0)
function [f,g] = obfunction_value1(x,b,c,y)%会多次调用，还是另用函数
% f =  (1-x(1))^2 + 2*(x(2)^2)+ (x(3)^2)+(2.5-x(4))^2;
% df(1) = -2*(1-x(1));
% df(2) = 4*x(2);
% df(3) = 2*x(3);
% df(4) = -2*(2.5-x(4));
% g = df;

if 1
    y=ones(1,4);
end
f =  (1-x(1))^2 + 2*(x(2)^2)+ (x(3)^2)+(2.5-x(4))^2+x(1)*x(2)+b+c(1)+c(2)+c(3)+c(4);
df(1) = -2*(1-x(1))+x(2);
df(2) = 4*x(2)+x(1);
df(3) = 2*x(3);
df(4) = -2*(2.5-x(4));
g = df;
mm=1;
end

%%


