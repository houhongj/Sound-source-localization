%利用保存的工作区文件，绘图
close all;
clc
clear all



load("2020-10-29 11-18-49.mat")
num_interval = 10;
num_fig_row = 2;
num_com = size(bestX,2);
num_fig = floor(num_com/num_interval);

num_fig_line = floor(num_fig/num_fig_row);
% x_loc1=x_t(:,1).';
% x_loc = repmat(x_loc1,size(mic_pos,1),1);
maxSPL=max(max(BB));
% BB(BB<=maxSPL-dBrange)=NaN;
BB(BB<=maxSPL-dBrange)=maxSPL-dBrange;
figure()
contourf( x_loc, xt_relative_ang, BB,'LineStyle','none');
hold on
plot([source_info_1(:,1) source_info_1(:,1)], get(gca, 'YLim'), '-y')
hold on;
axis([[min(min(x_loc)) max(max(x_loc))] [min(min(xt_relative_ang)) max(max(xt_relative_ang))]]);
load mycolor.mat;
colormap(mycolor);
colorbar;caxis([maxSPL-dBrange maxSPL]);
axis xy;
%% 
figure()
for i=1:num_fig
    
    X = reshape(bestX(:,i+10*(i-1)),size(x_loc,2),size(x_loc,1));
    BB = 20*log10((abs(X).^2)/2e-5);
    BB=rot90(BB);
    maxSPL=max(max(BB));
    % BB(BB<=maxSPL-dBrange)=NaN;
    subplot(num_fig_row,num_fig_line,i)
    contourf( x_loc, xt_relative_ang, BB,'LineStyle','none');
    hold on
    plot([source_info_1(:,1) source_info_1(:,1)], get(gca, 'YLim'), '-y')
    hold on;
    axis([[min(min(x_loc)) max(max(x_loc))] [min(min(xt_relative_ang)) max(max(xt_relative_ang))]]);
    load mycolor.mat;
    colormap(mycolor);
    colorbar;caxis([maxSPL-dBrange maxSPL]);
    axis xy;
    
end



