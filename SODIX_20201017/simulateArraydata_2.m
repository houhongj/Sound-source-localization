function [p, Fs] = simulateArraydata_2(source_info, mic_info, coherent, c, Fs, duration, U, direct_amp)
% hhj0920
%����Դ�ķ����ԣ��ڲ�ͬ��˷�ǿ��ǰ������㷽��ϵ��
% function [p, Fs] = simulateArraydata(source_info, mic_info, c, Fs, duration, U)
%simulateArraydata   Simulates pressure data for microphone array. 
%   simulateArraydata(source_info, mic_info, c, Fs, duration, U) generates
%   synthetic microphone array data. Option to include wind by using vector
%   U only for tonal source!
%   
%   simulateArraydata(source_info, mic_info) minimum requirement for call
%   
%   source_info(l, :) is [x y z freq SPL_at_array] for each row, 
%   i.e. several rows are more than one source. This generates tonal noise
%   source for frequency freq. If freq = 0 will output white noise.
%   Multiple white noise sources currently not working
%

%   Anwar Malgoezar, May 2018. 
%   Group ANCE

% Add coherent parameter for create the coherent or incoherent
% signals---nfl 2020.5.16

N_source = size(source_info, 1);
N_mic = size(mic_info, 1);
f_min = min(source_info(:, 4));% Ƶ��min
f_max = max(source_info(:, 4));

% Array center
x_ac = mean(mic_info,1);  %�޷���������������ѹ1��������2��ָ����

if nargin < 3
    % coherent signal
    coherent = true;
end
if nargin < 4
    % speed of sound
    c = 343;
end
if nargin < 5
    % sample freq to be 20 times the highest frequency, if only white noise
    % sources set to 50 kHz
    if f_max == 0
        Fs = 50e3;
    else
        Fs = 20*f_max;
    end
end
if nargin < 6
    if f_max == 0
        duration = 1;
    else
    % duration of signal to be 5 periods of lowest frequency
        duration = 5/f_min;
    end
end
if nargin < 7
    U = [0 0 0];
end
% See Pieter's derivation for steering vector with convection
% Ref.Phased array beamforming applied to wind tunnel and fly-over tests 
%  hhj--20200528
M = U/c;
beta_sq = 1 - dot(M, M);

% compute sample points
t = 0:1/Fs:(duration-1/Fs);
N_samples = length(t);

p = zeros(N_mic, N_samples);
p1 = zeros(N_mic, N_samples);
% r_ac = 1;
for I = 1:N_source
    r_ac = norm(x_ac-source_info(I, 1:3));
    
    % White noise or tonal
    if source_info(I, 4)==0
        % Before generating white noise, need to determine the shift in
        % samples for max and min distance compared to array center.
        r_max = sqrt(max(sum((mic_info-source_info(I, 1:3)).^2, 2)));
        

        max_shift = round(Fs*(r_max-r_ac)/c) * ((r_max-r_ac)>0);
        % Fs*(r_max-r_ac)/c �ӳٲ�����
%         bb=((r_max-r_ac)>0);???   ((r_ac-r_min>0)) hhj--20200528
        r_min = sqrt(min(sum((mic_info-source_info(I, 1:3)).^2, 2)));
        min_shift = round(Fs*(r_ac-r_min)/c) * ((r_ac-r_min>0));
        
        % Remember power in wgn is in dBW, i.e. 10*log10(P), for an 
        % impedance = rho*c = 1. For acoustics we have P = p_rms^2/(rho*c)
        % dBW - 20*log10(p_0) = dBSPL. At the center of the array Wc
        Wc = wgn(min_shift + Fs*duration + max_shift, 1, ...
                 source_info(I, 5)+10*log10(2e-5^2)).';
        % �����������Ĵ���ָ����ѹ���İ����� hhj--20200528
        for J = 1:N_mic

            r = sqrt( sum((mic_info(J, :) - source_info(I, 1:3)).^2) );
            ph = (r-r_ac)/c;
            N_shift = round(ph*Fs);
            p(J, :) = p(J, :) + Wc(1+min_shift+N_shift:Fs*duration+N_shift+min_shift)*r_ac/r;
        end
        
    else% ��Ƶ�ź�
        % From SPL = 20*log10(amp/2e-5), we scale by the distance from
        % source to array center to obtain given SPL at center.
        amp = r_ac*2e-5*10^(source_info(I, 5)/20);
%         amp = 2e-5*10^(source_info(I, 5)/20);
%         amp = 2e-5*10^(source_info(I, 5)/20);
        % in order to get the incoherent signals, add the random
        % phase to every source signal.---nfl 2020.5.16
        if coherent
            randphase = zeros([1,N_samples]);
        else
            randphase = randn(1,N_samples);
        end
        for J = 1:N_mic

%             r = sqrt( sum((mic_info(J, :) - source_info(I, 1:3)).^2) );
%             ph = r/c; %UΪ0ʱ������ʹ�õļ򻯹�ʽ
%             r1(J, :) = sqrt( dot(M, mic_info(J, :) - source_info(I, 1:3))^2 + ...
%                 beta_sq*dot(mic_info(J, :) - source_info(I, 1:3), mic_info(J, :) - source_info(I, 1:3)) );
            r = sqrt( dot(M, mic_info(J, :) - source_info(I, 1:3))^2 + ...
                beta_sq*dot(mic_info(J, :) - source_info(I, 1:3), mic_info(J, :) - source_info(I, 1:3)) );
            ph = ( -dot(M, mic_info(J, :) - source_info(I, 1:3)) + r ) / (c * beta_sq);
            
            p1(J, :) = p1(J, :) + sqrt(2)*amp*cos(2*pi*(source_info(I, 4)*(t))+randphase)/r;
            %hhj  ������ϵ��
            p(J, :) = p(J, :) + direct_amp(I,J)*sqrt(2)*amp*cos(2*pi*(source_info(I, 4)*(t-ph))+randphase)/r;
%             p(J, :) = p(J, :) + sqrt(2)*amp*cos(2*pi*(source_info(I, 4)*(t-ph))+randphase)/r;
            % Simulate dipole, still needs works
%             p(J, :) = p(J, :) + (0.0001*(2*pi*source_info(I, 4)/c)^2*sin(acos(source_info(I, 3)/r)))*sqrt(2)*amp*cos(2*pi*source_info(I, 4)*(t-ph))/r;
        end
        
    end
end    

end