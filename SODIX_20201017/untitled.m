function [x,fval,exitflag,output,grad,hessian] = untitled(x0)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('fminunc');
%% Modify options setting
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'PlotFcn', {  @optimplotx @optimplotfval });
options = optimoptions(options,'Algorithm', 'quasi-newton');
options = optimoptions(options,'SpecifyObjectiveGradient', true);
options = optimoptions(options,'Hessian', 'off');
[x,fval,exitflag,output,grad,hessian] = ...
fminunc(@obfunction_value,x0,options);
