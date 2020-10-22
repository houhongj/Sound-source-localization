function [source_info,num_source_info] = source_setup(bf_freq, souSPL, wavelength)



%Դ����
%Դ���ͣ�������1������2������3������4

num_source_class = [1 0 0 0];%Դ���ͣ�������1������2������3������4
jet_gas = 5;%��������Դ���߷���
% source_info_set(:,1) = [-5.4 -2.2 0 2];%Դλ��x
source_info_set(:,1) = [-1.320 0.04 1.06 2];%Դλ��[-1.320 0.04 1.06 2]
source_info_set(:,2) = [0 0 0 0];%Դλ��y
source_info_set(:,3) = [0 0 0 0];%Դλ��z
source_info_set(:,4) = ones(4,1)*bf_freq;%ԴƵ��
source_info_set(:,5) = ones(4,1)*souSPL;%Դ��ѹ
source_info_set(:,6) = [35 35 35 35];%Դ������ָ����״������˹����
source_info_set(:,7) = [120 90 60 45];%Դ������ָ�����ֵ������˹ƽ��
source_info_set(:,8) = [1 1 1 1];%Դ�����Է���
source_info_set(:,9) = [1 2 3 4];%Դ���������
source_info_set(:,10) = [3 3 3 25];%ÿ��Դ����
source_info_set(:,11) = [1/20*wavelength 1/20*wavelength 1/20*wavelength 10/20*wavelength];%ÿ��Դ���
%% ��Դ
for I=1:3
    if num_source_class(I)==1
        num_source = source_info_set(I,10);
%         source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,1) = [-source_info_set(I,11)+source_info_set(I,1),source_info_set(I,1),source_info_set(I,11)+source_info_set(I,1)];
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,1) = (source_info_set(I,1):source_info_set(I,11):source_info_set(I,1)+source_info_set(I,11)*(source_info_set(I,10)-1));
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,4) = ones(source_info_set(I,10),1)*bf_freq;
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,5) = ones(source_info_set(I,10),1)*source_info_set(I,5);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,6) = ones(source_info_set(I,10),1)*source_info_set(I,6);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,7) = ones(source_info_set(I,10),1)*source_info_set(I,7);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,8) = ones(source_info_set(I,10),1)*source_info_set(I,8);
        source_info(1+(I-1)*num_source:num_source+(I-1)*num_source,9) = ones(source_info_set(I,10),1)*source_info_set(I,9);
           
    end
    
end
num_source_info = size(source_info,1);
%% ����Դ
J=4;
jet_locx = zeros(1,source_info_set(4,10));
if num_source_class(J)==1
    
    num_source = source_info_set(J,10);
    for I=1:num_source
        jet_locx(I) = (I-3)*source_info_set(J,11)+source_info_set(J,1);    %3Ϊ��˹��ֵǰ���Դ����
%         I=I+1;
    end
    jet_amp =  1*gaussmf(jet_locx,[jet_gas source_info_set(J,1)]);
    jet_amp1 =  1*gaussmf((1:0.1:11),[jet_gas source_info_set(J,1)]);
    figure(99)
    plot(jet_locx,jet_amp, '*')
    hold on
    plot((1:0.1:11),jet_amp1)
    source_info(num_source_info+1:num_source_info+num_source,1) = jet_locx.';
    source_info(num_source_info+1:num_source_info+num_source,4) = ones(source_info_set(J,10),1)*bf_freq;
    source_info(num_source_info+1:num_source_info+num_source,5) = ones(source_info_set(J,10),1)*source_info_set(J,5);
    source_info(num_source_info+1:num_source_info+num_source,6) = ones(source_info_set(J,10),1)*source_info_set(J,6);
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