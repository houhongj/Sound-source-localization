% 测试conjungate_gradient hhj20200710
% 无法解决复数问题
close all;
clc
clear all

% test conjungate gradient method
% 
addpath('../rodyo-FEX-minimize-2f6ee87');
%% define function and variable

%f = xs^2+2*ys^2-2*xs*ys + 2*ys + 2;
% f = (x1-1)^4 + (x1 - x2)^2;
% f = (1-x1)^2 + 2*(x2 - x1^2)^2;
% f = (1-x1)^2 + 2*(x2^2)+ (x3^2)+(2.5-x4)^2-2*x1*x3+4*x3*x4+3*x1*x2;
% f = (1-x1)^2 + 2*(x2^2)+ (x3^2)+(2.5-x4)^2+x1*x2;
% f = fa(x);
% x = {x1, x2 x3 x4};

% initial value
x0 = [0 0 0 0];
% tolerance
epsilon = 1e-1;

%% call conjungate gradient method
%for i=1:4
%     eval(['syms x',num2str(i)]);
%     eval(['x{i}=x',num2str(i)]);
%     xbb{i}={['syms x',num2str(i)]};
%     x(i)=sym (['x',num2str(i)]);
% end

% f1 = (1-x1)^2 + 2*(x2^2)+ (x3^2)+(2.5-x4)^2+x1*x2;
% show_detail = true;
% [bestf, bestx, count] = conjungate_gradient(f1, x, x0, epsilon, show_detail);
% % print result
% fprintf('bestx = %s, bestf = %f, count = %d\n', num2str(bestx), bestf, count);

%% 共轭梯度

% leng=10;
% [X, fX, i] = minimize(x0, @obfunction_value222, leng, x, f);


%% 低内存拟牛顿约束算法
x=sym('x',[1,4]);
% sym x
x1 =0;
% f = @(x)(1-x(1))^2 + 2*(x(2)^2)+ (x(3)^2)+(2.5-x(4))^2+x(1)*x(2);
% f = (1-x(1))^2 + 2*(x(2)^2)+ (x(3)^2)+(2.5-x(4))^2+x(1)*x(2);
% options = optimset('GradObj','on');
options = setoptimoptions('Algorithm','fminlbfgs','GradObj','off');

% options = optimset('TolFun', 1e-5, 'TolX', 1e-7);
a = minimize(@obfunction_value22,x0,[],[],[],[],[],[],[],options)
% a = minimize(@myfun,x0)



b=1;
%%
function [fun] = obfunction_value22(x0, x)
% function [fun, dfun] = obfunction_value222(x0, x, f)
f = (1-x(1))^2 + 2*(x(2)^2)+ (x(3)^2)+(2.5-x(4))^2+x(1)*x(2);
fun = f;
% fun = double(subs(f, x, x0));
n=length(x);


dfun = double(subs(f, x, x0));
a=1;
end
      function [f] = myfun(x)
      f = sin(x) + 3;
% 	    if ( nargout > 1 ), g = cos(x); end
      end