close all;
clc
clear all
tic

clearvars;
addpath('../Acoustic-Beamforming-master');
addpath('../GPExp-master\code\gpml3.5\util');
% addpath('../rodyo-FEX-minimize-2f6ee87');
syms fun_syms
%% 初始化
x_range = 1*[-2 0];
bf_freq = 400;
mic_rang = [10 160];  %麦克风范围
R=5;                 %麦克风半径距离
mic_ang = 20;           %麦克风间隔

y_range = 1*[0 0];
z_range = 0;
c = 340;
dBrange = 30;
coherent = true;% coherent = false;
SNR =5;
wavelength = c/bf_freq;
% 1000hz---0.136   200hz---0.68    400--0.34   200hz用0.4以下100hz  -0.2lamuda 200----0.3
inte_1=0.4*wavelength;
x_t(:,1) = [x_range(1):inte_1:x_range(2)];
x_t(:,2) = 0;
x_t(:,3) = 0;
num_syms_x = size(x_t,1);        %变量个数  -未知点源个数
%%
% angle1 = (mic_rang(1):mic_ang:20);
% angle2 = (25:mic_ang*0.5:75);
% angle3 = (80:mic_ang:mic_rang(2));
% poldeg = [angle1 angle2 angle3];
poldeg = (mic_rang(1):mic_ang:mic_rang(2));
mic_pos(:,1)=R*cos(poldeg/180*pi);
mic_pos(:,2) = 0;
mic_pos(:,3)=R*sin(poldeg/180*pi);

[source_info_1] = source_setup(bf_freq, 100, wavelength);

[direct_amp, relative_ang] = source_direct_2(source_info_1, mic_pos, mic_rang, mic_ang, R);
source_info = source_info_1(:,1:5);
[p, Fs] = simulateArraydata_2(source_info, mic_pos, coherent, c,48e3, 10, [0 0 0],direct_amp); 
[CSM, freqs] = developCSM(p.', bf_freq, bf_freq, Fs, 1, 0.5);
% CSM_Ref = zeros(size(mic_pos,1),size(mic_pos,1));

%% 目标函数
% fun_syms = 0;
[fun_syms, x, dfun_syms,djm0] = obfunction_2(num_syms_x, x_t, CSM, mic_pos, freqs, ...
    z_range, c);
% xt_relative_ang = xt_direct_2(x_t, mic_pos, mic_rang, mic_ang, R);
x0 = djm0;

%% 
leng = min(num_syms_x*size(mic_pos,1),50);
% leng = 5;        %最大搜索布数
% [X, fX, i, bestx] = minimize(x0, @obfunction_value, leng, x, fun_syms, dfun_syms);
[X, fX, i, bestX] = minimize(x0, @obfunction_value, leng, x, fun_syms, dfun_syms);
% s=1;
% a=1;
% while 1
%     [X, fX, i] = minimize(x0, @obfunction_value, leng, x, fun_syms, dfun_syms);
%     
%     xt_relative_ang = xt_direct_2(x_t, mic_pos, mic_rang, mic_ang, R);
%     BB = 20*log10((abs(X).^2)/2e-5);
%     % BB=flipud(BB);
%     BB=rot90(BB);
%     xt_relative_ang=rot90(xt_relative_ang);
%     maxSPL=max(max(BB));
%     % BB(BB<=maxSPL-dBrange)=NaN;
%     
% %     figure();
%     subplot(2,5,s)
%     x_loc1=x_t(:,1).';
%     x_loc = repmat(x_loc1,size(mic_pos,1),1);
%     y_loc1=[mic_rang(1):mic_ang:mic_rang(2)].';
%     y_loc = repmat(y_loc1,1,size(x_t,1));
%     contourf( x_loc, xt_relative_ang, BB,'LineStyle','none');
%     % imagesc( x_loc, xt_relative_ang.', BB.');
%     % contourf( x_loc, y_loc, djm0.');
%     hold on;
%     % axis([[-5 15] [mic_rang(1) mic_rang(2)]]);
%     axis([x_range [min(min(xt_relative_ang)) max(max(xt_relative_ang))]]);
%     colorbar;caxis([maxSPL-dBrange maxSPL]);
%     axis xy;
%     
% %     a=input('是否继续？y/n,\n','s');
%     while a~='y'&&a~='n'
%         disp('输入错误du\n');
%         a=input('是否继续？y/n,\n','s');
%     end
%     if a=='n'|| s==10
%         break;
%     end
%     x0 = X;
%     s = s+1;
% end
%% 
% xt_relative_ang = xt_direct_2(x_t, mic_pos, mic_rang, mic_ang, R);
%% 画图
%
xt_relative_ang1 = xt_direct_2(x_t, mic_pos, mic_rang, mic_ang, R);
BB = 20*log10((abs(X).^2)/2e-5);

BB=rot90(BB);
xt_relative_ang=rot90(xt_relative_ang1);
maxSPL=max(max(BB));
% BB(BB<=maxSPL-dBrange)=NaN;

figure();
x_loc1=x_t(:,1).';
x_loc = repmat(x_loc1,size(mic_pos,1),1);
y_loc1=[mic_rang(1):mic_ang:mic_rang(2)].';
y_loc = repmat(y_loc1,1,size(x_t,1));
contourf( x_loc, xt_relative_ang, BB,'LineStyle','none');
hold on
plot([source_info_1(:,1) source_info_1(:,1)], get(gca, 'YLim'), '-y')
% hold on
% plot([-source_info_1(:,1) -source_info_1(:,1)], get(gca, 'YLim'), '-y')
% imagesc( x_loc, xt_relative_ang.', BB.');
% contourf( x_loc, y_loc, djm0.');
hold on;
% axis([[-5 15] [mic_rang(1) mic_rang(2)]]);
axis([[min(min(x_loc)) max(max(x_loc))] [min(min(xt_relative_ang)) max(max(xt_relative_ang))]]);
colorbar;caxis([maxSPL-dBrange maxSPL]);
axis xy;

toc