function [relative_ang] = xt_direct_2(x_t, mic_info, mic_rang, mic_ang, R)
%ªªÀ„∑ΩœÚ



num_source = size(x_t,1);
if size(x_t,2)<3
    x_t(:,3)=zeros(num_source,1);
end



for s=1:num_source
    if x_t(s,1)==0
        relative_ang(s,:) = [mic_rang(1):mic_ang:mic_rang(2)];
    else
        for m=1:size(mic_info,1)
            length_sou_mic(s,m) = sqrt( sum((mic_info(m, :) - x_t(s, 1:3)).^2) );
%             cos_ang(s,m) =acos( (x_t(s,1)^2+length_sou_mic(s,m)^2-R^2)/(2*(x_t(s,1))*length_sou_mic(s,m)));
%             relative_ang(s,m) = 180-cos_ang(s,m)*180/pi;
            
            loc_mic_tem = mic_info(m,1)-x_t(s,1);
            cos_ang(s,m) =acos( loc_mic_tem/length_sou_mic(s,m));
            relative_ang(s,m) = cos_ang(s,m)*180/pi;
        end
    end
    
end


%% 
a=1;
end