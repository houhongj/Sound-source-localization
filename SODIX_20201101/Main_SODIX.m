
%利用lbfgs计算，速度快
%正则化

%FIXME--exitflag为-2
%FIXME--传递函数？单频信号？白噪声？负号？利用单频点无法生成特定指向性线源
%FIXME--白噪声功率，wgn，特定频率下的功率
close all;
clc
clear all
tic
clearvars;
% addpath('../Acoustic-Beamforming-master');
% addpath('../GPExp-master\code\gpml3.5\util');
addpath('../rodyo-FEX-minimize-2f6ee87');
syms fun_syms


%% 初始化
x_range = 1*[-2 5];
bf_freq = 400;
c = 340;
mic_rang = [10 170];  %麦克风范围
R=10;                 %麦克风半径距离
mic_ang = 10;           %麦克风间隔
num_source_class = [1 1 1 1];%源类型：进气口1；机身2；喷嘴3，喷流4
source_info_loc = [-1.320 0.04 1.06 2];%源位置[-5.4 -2.2 0 2] [-1.320 0.04 1.06 2]
source_info_num = [3 3 3 30];%每组源个数
wavelength = c/bf_freq;
source_info_inte = [1/20*wavelength 1/20*wavelength 1/20*wavelength 2/20*wavelength];%每组源间隔
% Regularization_factor = [0 0];
% Regularization_factor_loop = [0 0];
Regularization_factor_loop = [0 0;0.00001 0;0.0001 0;0.001 0;0.01 0;0.1 0];

y_range = 1*[0 0];
z_range = 0;
dBrange = 30;
coherent = true;% coherent = false;
SNR =5;
% 1000hz---0.136   200hz---0.68    400--0.34   200hz用0.4以下100hz  -0.2lamuda 200----0.3
inte_1=0.4*wavelength;
x_t(:,1) = (x_range(1):inte_1:x_range(2));
x_t(:,2) = 0;
x_t(:,3) = 0;
num_syms_x = size(x_t,1);        %-未知点源个数
fs = 44.1e3;
t_sampling = 5;
%%
poldeg = (mic_rang(1):mic_ang:mic_rang(2));
mic_pos(:,1)=R*cos(poldeg/180*pi);
mic_pos(:,2) = 0;
mic_pos(:,3)=R*sin(poldeg/180*pi);

[source_info_1] = source_setup(num_source_class, source_info_loc,source_info_num, 120, x_range, source_info_inte);
[direct_amp, relative_ang] = source_direct(source_info_1, mic_pos, mic_rang, mic_ang, R);

%% 
source_info = source_info_1(:,1:5);
% [p, Fs] = simulateArraydata_2(source_info, mic_pos, coherent, c,44.1e3, 5, [0 0
% 0],direct_amp);  %单频信号
% [p, Fs] = simulateArraydata(source_info, mic_pos, c,44.1e3, 5, [0 0 0]); 

p = zeros(size(mic_pos,1),fs*t_sampling);
for i=1:size(source_info,1)
    [p1, Fs] = simulateArraydata_SODIX(source_info(i,:), mic_pos, coherent, c, fs, t_sampling, [0 0 0],direct_amp(i,:)); 
    p = p+p1;
end


[CSM, freqs] = developCSM(p.', bf_freq, bf_freq, Fs, 0.5, 0.5);

%% 目标函数
for loop=1:size(Regularization_factor_loop,1)
    Regularization_factor = Regularization_factor_loop(loop,:);
[fun_syms, x, dfun_syms,djm0] = obfunction_SODIX(num_syms_x, x_t, CSM, mic_pos, freqs, ...
    Regularization_factor, c);

%% fminunc

% options = optimoptions('fminunc');
% options = optimoptions(options,'Display', 'iter');
% options = optimoptions(options,'PlotFcn', {  @optimplotx @optimplotfval });
% options = optimoptions(options,'Algorithm', 'quasi-newton');
% options = optimoptions(options,'SpecifyObjectiveGradient', true);
% options = optimoptions(options,'Hessian', 'off');
% [X,fval,exitflag,output,grad,hessian] = ...
% fminunc(@obfunction_value,x0,options,x,fun_syms, dfun_syms);

 %% 计算慢的方法
%  fminsearch
%  CG
%% 
x0 = djm0;
% x0=1*ones(num_syms_x, size(mic_pos,1));
leng = min(num_syms_x*size(mic_pos,1),100);
options = setoptimoptions('Algorithm','fminlbfgs','AlwaysHonorConstraints', 'bounds','GradObj','on'...
    ,'Display','iter','MaxIter',leng ,'GoalsExactAchieve',1,'HessUpdate','lbfgs' );%fminsearch  fminlbfgs,    'PlotFcn',@optimplotfval ,'Display','iter'
 
lb = zeros(num_syms_x, size(mic_pos,1));%x0,[],[],[],[],[lb],[],[], options  ,'AlwaysHonorConstraints', 'bounds'

[X, fval, bestX,bestf, exitflag, output] = minimize( @obfunction_value, x0,[],[],[],[],lb,[],[], options, ...
    x,fun_syms, dfun_syms);
bestX = bestX(:,1:output.iterations-1);
bestf = bestf(:,1:output.iterations-1);
% [sol, fval, exitflag, output, grad] = minimize( @obfunction_fun, x0,[],[],[],[],[],[],[], options, ...
%     fun_syms, dfun_syms,x);
%% 画图

xt_relative_ang1 = xt_direct(x_t, mic_pos, mic_rang, mic_ang);
BB = 20*log10((X)/2e-5);

BB=rot90(BB);
xt_relative_ang=rot90(xt_relative_ang1);
maxSPL=max(max(BB));
% BB(BB<=maxSPL-dBrange)=NaN;
BB(BB<=maxSPL-dBrange)=maxSPL-dBrange;
figure();
x_loc1=x_t(:,1).';
x_loc = repmat(x_loc1,size(mic_pos,1),1);
y_loc1=(mic_rang(1):mic_ang:mic_rang(2)).';
y_loc = repmat(y_loc1,1,size(x_t,1));
contourf( x_loc, xt_relative_ang, BB,'LineStyle','none');
hold on

numsou = num_source_class.*source_info_num;
numsou1= [1,num_source_class(2)*sum(numsou(1))+1,num_source_class(3)*sum(numsou(1:2))+1,num_source_class(4)*sum(numsou(1:3))+1];
plot([source_info_1([numsou1(1),numsou1(2),numsou1(3),numsou1(4),end],1) source_info_1([numsou1(1),numsou1(2),numsou1(3),numsou1(4),end],1)], get(gca, 'YLim'), '-r')

hold on
plot(get(gca, 'XLim'), [source_info_1(:,7) source_info_1(:,7)],  '-r')
hold on;
axis([[min(min(x_loc)) max(max(x_loc))] [min(min(xt_relative_ang)) max(max(xt_relative_ang))]]);
load mycolor.mat;
colormap(mycolor);
colorbar;
caxis([maxSPL-dBrange maxSPL]);
axis xy;
title(['Regularization-factor: ' num2str(Regularization_factor)]);
name =  strrep(strrep(char(datetime),':','-'),' ','-');
save(name ,'BB','bestX','xt_relative_ang','dBrange','x_loc', 'source_info_1', 'R','X','Regularization_factor')
name_fig1 = strcat(name,'_fig1');
name_fig2 = strcat(name,'_fig2');
saveas(1,name_fig1,'jpg');
saveas(gcf,name_fig2,'jpg');
toc
end

