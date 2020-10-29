% Modify the regularization conditions and derivatives，x
function [fun_syms, x, dfun_syms, djm0] = obfunction_4(num_syms, x_t, CSM, mic_positions, frequencies, ...
    plane_distance, c)

syms fun_syms 

N_mic = size(mic_positions, 1);
N_freqs = length(frequencies);
x_t(:,3)=plane_distance;


x=sym('x',[num_syms,N_mic]);
fun_syms = 0;
%目标函数
reverseStr = '';
G = zeros( size(x_t, 1),N_mic);
r_ti = zeros( size(x_t, 1),N_mic);
for K = 1:N_freqs
    msg = sprintf('\tBeamforming %d/%d frequency points...\n', K, N_freqs);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    k = 2*pi*frequencies(K)/c;
    
    for I = 1:N_mic
         r_ti(:, I) = sqrt( (x_t(:,1) - mic_positions(I,1)).^2 + ...
                     (x_t(:,2) - mic_positions(I,2)).^2 + ...
                     (x_t(:,3) - mic_positions(I,3)).^2 );%扫描点与每个麦克风的距离
                 
    end
    %发动机实际测试需修改导向矢量，包括U   ref.2008
    G = exp(-1i*k*(r_ti))./(r_ti);
    %% 构建最小化目标函数
%     ref.2012--An extended formulation of the SODIX method with application to aeroengine broadband noise


    CSM_mod = ((x.*G).')*(x.*conj(G));
%     CSM_mod = ((A.*conj(G)).')*(A.*(G));
    err_C = CSM-CSM_mod;
    
%     fun = norm(err_C,'fro')^2;
%     fun = sum(sum(err_C.*conj(err_C)));
fun = sum(sum( (real(CSM) - real(CSM_mod)).^2 )) + ...
    sum(sum( (imag(CSM) - imag(CSM)).^2 ));
% fun = sum(sum( (real(CSM) - real(CSM_mod)).^2 )) ;
    %% 目标函数导数
    %参见公式6，三部分
%     GG_mod = (G.*conj(G));%
%     err_C_0diag = err_C-diag(diag(err_C));
    DrmGmod = x.*G;
    DrlGmod = conj(G);

%     dfun1=-4*real(DrmGmod*real(err_C_0diag).*DrlGmod);
%     dfun2=-4*imag(DrmGmod*imag(err_C_0diag).*DrlGmod);
%     dfun3=-4*(x).*G.*conj(G)*diag(diag(err_C));
%     dfun = dfun1+dfun2+dfun3;
    
    
    %ref.2019，
    dfun=-4*real((DrmGmod*(err_C.')).*DrlGmod);





    %% 平滑函数/正则化条件
%     ref.Noise Source Analysis of an Aeroengine with a New Inverse Method SODIX
    %平滑函数1 ---更均匀的方向性/相邻麦克风变换不能太大
    
    A_1(:,3:N_mic+2) = x;
    A_1(:,N_mic+3:N_mic+4) =0;
    A_2 = A_1(:,2:end-1);
    Djm_1 = A_2(:,2:end-1);
    Djm_0 = A_2(:,1:end-2);
    Djm_2 = A_2(:,3:end);
    Djm = Djm_1-0.5*(Djm_0+Djm_2);
% 
    G1 = sum(sum( (Djm).^2 ));
% 
    dDjm_1 = A_1(:,1:end-4);
    dDjm_2 = A_1(:,2:end-3);
    dDjm_3 = A_1(:,3:end-2);
    dDjm_4 = A_1(:,4:end-1);
    dDjm_5 = A_1(:,5:end);
    dG1 = 0.5*dDjm_1-2*dDjm_2+3*dDjm_3-2*dDjm_4+0.5*dDjm_5;
%     
%     
%     %平滑函数2 ---平滑源强度沿发动机轴线的变化
    A_3(3:num_syms+2,:) = x;
    A_3(num_syms+3:num_syms+4,:) =0;
    A_4 = A_3(2:end-1,:);
    Dj_1m = A_4(2:end-1,:);
    Dj_0m = A_4(1:end-2,:);
    Dj_2m = A_4(3:end,:);
    Dj_m = Dj_1m-0.5*(Dj_0m+Dj_2m);

    G2 = sum(sum( (Dj_m).^2 ));
    
    dDj_1m = A_3(1:end-4,:);
    dDj_2m = A_3(2:end-3,:);
    dDj_3m = A_3(3:end-2,:);
    dDj_4m = A_3(4:end-1,:);
    dDj_5m = A_3(5:end,:);
    dG2    = 0.5*dDj_1m-2*dDj_2m+3*dDj_3m-2*dDj_4m+0.5*dDj_5m;

    %% 
%     fun_syms = fun+0.0001*size(x_t,1)*G1+0.05*G2;
%     dfun_syms = dfun+0.0001*size(x_t,1)*dG1+0.05*dG2;
%     fun_syms = fun+0.0001*size(x_t,1)*G1;
%     dfun_syms = dfun+0.0001*size(x_t,1)*dG1;
    fun_syms = fun;
    dfun_syms = dfun;
    %% 初始值
%     ref. Advancements in the source localization method SODIX and application to short cowl engine data
    rjm = (1./(4*pi*r_ti)).^2;
    sumrjm = sum(rjm,1);
    djm0_t = ((real(diag(CSM)).')./sumrjm).^0.5;%原文中为没有平方项，0.5，现在求解平方项。0.25
    djm0 = repmat(djm0_t,num_syms ,1);
    


end

end
