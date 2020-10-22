function [direct_amp, relative_ang] = source_direct_2(source_info_1, mic_info, mic_rang, mic_ang, R)
%Դ�����Խ�ģ
% source_info_1 ��6��7�п��Ʒ��������ߣ���8��Ϊ���ȣ���9��Ϊ���


num_source = size(source_info_1,1);


figure(num_source+1)
plot(mic_info(:,1),mic_info(:,3),'*') 
axis equal;axis([[-R R] [0 R]]);
% if size(source_info_1,2)<9
%     source_info_1(:,9)=ones(num_source,1);
% end
% Initialize the data 
relative_ang = zeros(size(source_info_1,1),size(mic_info,1));
length_sou_mic = zeros(size(source_info_1,1),size(mic_info,1));
cos_ang = zeros(size(source_info_1,1),size(mic_info,1));
direct_amp = zeros(size(source_info_1,1),size(mic_info,1));
direct_amp1 = zeros(size(source_info_1,1),180);
for s=1:num_source
    if source_info_1(s,1)==0
        relative_ang(s,:) = (mic_rang(1):mic_ang:mic_rang(2));
    else
        for m=1:size(mic_info,1)
            length_sou_mic(s,m) = sqrt( sum((mic_info(m, :) - source_info_1(s, 1:3)).^2) );
%             cos_ang(s,m) =acos( (source_info_1(s,1)^2+length_sou_mic(s,m)^2-R^2)/(2*(source_info_1(s,1))*length_sou_mic(s,m)));
%             relative_ang(s,m) = 180-cos_ang(s,m)*180/pi;
            % ���Ϊ���������нǼ��㣬������ά�ռ����
            loc_mic_tem = mic_info(m,1)-source_info_1(s,1);
            cos_ang(s,m) =acos( loc_mic_tem/length_sou_mic(s,m));
            relative_ang(s,m) = cos_ang(s,m)*180/pi;
        end
    end
    
    % ��ʱ����Ϊ��˹����  s==2 || s==3
    if source_info_1(s,9)==1|| source_info_1(s,9)==2 || source_info_1(s,9)==3|| source_info_1(s,9)==4
        direct_amp(s,:) =  source_info_1(s,8)*gaussmf(relative_ang(s,:),source_info_1(s,6:7));
        direct_amp1(s,:) =  source_info_1(s,8)*gaussmf((1:1:180),source_info_1(s,6:7));
        
%     elseif source_info_1(s,9)==2    %���2-�޸�����
%         direct_amp(s,:) =  source_info_1(s,8)*gaussmf(relative_ang(s,:),source_info_1(s,6:7));
%         direct_amp1(s,:) =  source_info_1(s,8)*gaussmf([1:1:180],source_info_1(s,6:7));
%         
%     else  %%���3-�޸�����
%         direct_amp(s,:) =  source_info_1(s,8)*gaussmf(relative_ang(s,:),source_info_1(s,6:7));
    end

X=[source_info_1(s,1) ,source_info_1(s,1)+0.5*R*source_info_1(s,8)*cos(source_info_1(s,7)/180*pi)];%0.5*R*Ϊ����ϵ��
Y=[source_info_1(s,3) ,source_info_1(s,3)+0.5*R*source_info_1(s,8)*sin(source_info_1(s,7)/180*pi)];
figure(num_source+1)
hold on
line(X,Y);

figure(s)
plot(relative_ang(s,:),direct_amp(s,:),'*')
hold on
plot((1:1:180),direct_amp1(s,:))
end




end