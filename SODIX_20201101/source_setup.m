function [source_info,num_source_info] = source_setup(num_source_class, source_info_loc,source_info_num, souSPL, x_range,source_info_inte)


%源设置
%源类型：进气口1；机身2；喷嘴3，喷流4
% source_info = zeros(1,11);
% num_source_class = [1 1 1 1];%源类型：进气口1；机身2；喷嘴3，喷流4
jet_gas = 2;%喷流噪声源曲线方差5
% source_info_set(:,1) = [-5.4 -2.2 0 2];%源位置x
source_info_set(:,1) = source_info_loc;%源位置[-1.320 0.04 1.06 2]
source_info_set(:,2) = [0 0 0 0];%源位置y
source_info_set(:,3) = [0 0 0 0];%源位置z
% source_info_set(:,4) = ones(4,1)*bf_freq;%源频率
source_info_set(:,4) = ones(4,1)*0;%源频率
source_info_set(:,5) = ones(4,1)*souSPL;%源声压
source_info_set(:,6) = [35 35 35 35];%源方向性指向形状——高斯方差
source_info_set(:,7) = [120 90 60 45];%源方向性指向最大值——高斯平均
source_info_set(:,8) = [1 1 1 1];%源方向性幅度
source_info_set(:,9) = [1 2 3 4];%源方向性组别
source_info_set(:,10) = source_info_num;%每组源个数
source_info_set(:,11) = source_info_inte;%每组源间隔
%%%%%%不懂？喷流源设置间隔大的时候，出现角度偏差？个数大的时候也
%% 点源
for I=1:3
    if num_source_class(I)==1
        num_source = source_info_set(I,10);
%         source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,1) = [-source_info_set(I,11)+source_info_set(I,1),source_info_set(I,1),source_info_set(I,11)+source_info_set(I,1)];
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,1) = (source_info_set(I,1):source_info_set(I,11):source_info_set(I,1)+source_info_set(I,11)*(source_info_set(I,10)-1));
%         source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,4) = ones(source_info_set(I,10),1)*bf_freq;
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,4) = ones(source_info_set(I,10),1)*source_info_set(I,4);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,5) = ones(source_info_set(I,10),1)*source_info_set(I,5);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,6) = ones(source_info_set(I,10),1)*source_info_set(I,6);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,7) = ones(source_info_set(I,10),1)*source_info_set(I,7);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,8) = ones(source_info_set(I,10),1)*source_info_set(I,8);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,9) = ones(source_info_set(I,10),1)*source_info_set(I,9);
           
    end
    
end
% num_source_info = size(source_info,1);
if sum(num_source_class(1:3))==0
    num_source_info = 0;
else
    num_source_info = size(source_info,1);
end
    
%% 喷流源
J=4;

jet_locx = zeros(1,source_info_set(4,10));
if num_source_class(J)==1
%     sou_rand =sort(2*(rand(source_info_set(4,10),1)-0.5));
%     source_info_1(:,1) = sou_rand;
    num_source = source_info_set(J,10);
    for I=1:num_source
        jet_locx(I) = (I-5)*source_info_set(J,11)+source_info_set(J,1);    %3为高斯均值前面的源个数
%         I=I+1;
    end
    jet_amp =  1*gaussmf(jet_locx,[jet_gas source_info_set(J,1)]);
    jet_amp1 =  1*gaussmf((x_range(1):0.1:x_range(2)),[jet_gas source_info_set(J,1)]);
    figure(99)
    plot(jet_locx,jet_amp, '*')
    hold on
    plot((x_range(1):0.1:x_range(2)),jet_amp1)
    source_info(num_source_info+1:num_source_info+num_source,1) = jet_locx.';
%     source_info(num_source_info+1:num_source_info+num_source,4) = ones(source_info_set(J,10),1)*bf_freq;
    source_info(num_source_info+1:num_source_info+num_source,4) = ones(source_info_set(J,10),1)*source_info_set(J,4);
    source_info(num_source_info+1:num_source_info+num_source,5) = ones(source_info_set(J,10),1)*source_info_set(J,5);
    source_info(num_source_info+1:num_source_info+num_source,6) = ones(source_info_set(J,10),1)*source_info_set(J,6);%源方向性指向形状——高斯方差
    source_info(num_source_info+1:num_source_info+num_source,7) = ones(source_info_set(J,10),1)*source_info_set(J,7);
    source_info(num_source_info+1:num_source_info+num_source,8) = jet_amp;
    source_info(num_source_info+1:num_source_info+num_source,9) = ones(source_info_set(J,10),1)*source_info_set(J,9);
end


num_source_info = size(source_info,1);
source_info_tem = zeros(1,size(source_info,2));

for s=1:num_source_info
    
    if isempty(find(source_info(s,:) == 0, 1))
        source_info_tem = [source_info_tem;source_info(s,:)];
    else
        if size(find(source_info(s,:) == 0),2) ~= size(source_info,2)
            source_info_tem = [source_info_tem;source_info(s,:)];
        end
    end
end
source_info = source_info_tem(2:end,:);
num_source_info = size(source_info,1);


end