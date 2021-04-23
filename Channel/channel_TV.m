function [CH_Response, CH_Info] = channel_TV(p, K_abs)
% =========================================================================
% -- Function to generate a time-variant THz channel
% =========================================================================

% -- Function: [CH_Response, CH_Info] = channel_TV(p,K_abs)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       K_abs: Molecular absorption coefficient, matrix, size(p.Nsubc, p.Nsubb)
%              (number of subbands @ each subcarrier, number of subcarriers)

% -- Output Arguments:
%       CH_Response: A struct that contains two channel responses:
%           1) CH_Response.H: H(t,f) time-variant frequency domain response, a 5D-array of size (p.Qr, p.Qt, num_subcarries, time_realization, num_freqs_per_subcarr),
%               i.e. H(MIMO Rx SA, MIMO Tx SA, number of subcarries, number of sample per frame of the input signal, number of sub-bands inside each subcarrier)
%           2) CH_Response.h: h(t,tau) time-variant delay domain response, a 3D-Cell array of size (p.Qr, p.Qt, num_subcarries)
%               and each cell member is a matrix of size (time_realization, max_Delay_spread),
%               i.e. h{MIMO Rx SA, MIMO Tx SA, number of subcarries},
%               and max_Delay_spread refers to the time of arrival, delay, of the last ray in the last cluster
%       CH_Info: A struct that contains channel statistical parameters such as
%           1) CH_Info.Bc = Coherence bandwidth of the (time, delay) domain channel, a 4D-array of size (p.Qr, p.Qt, num_subcarries, time_realization)
%           2) CH_Info.Tau_rms = Root mean square of delay spread for the (time, delay) domain channel, a 4D-array of size (p.Qr, p.Qt, num_subcarries, time_realization)
%           3) CH_Info.Tc = Coherence time of the (time, delay) domain channel for each subcarrier, a vector of size(1, num_subcarries)
%           4) CH_Info.fd_max = Maximum Doppler shift due to mobility for the (time, delay) domain channel for each subcarrier, a vector of size(1, num_subcarries)

%=================================================

% -- (c) 2021 Simon Tarboush, Hadi Sarieddeen, Hui Chen, 
%             Mohamed Habib Loukil, Hakim Jemaa, 
%             Mohamed-Slim Alouini, Tareq Y. Al-Naffouri

% -- e-mail: simon.w.tarboush@gmail.com; hadi.sarieddeen@kaust.edu.sa; hui.chen@kaust.edu.sa;
%            mohamedhabib.loukil@kaust.edu.sa; hakim.jemaa@kaust.edu.sa;
%            slim.alouini@kaust.edu.sa; tareq.alnaffouri@kaust.edu.sa

% =========================================================================

% S. Tarboush, H. Sarieddeen, H. Chen, M.-H. Loukil, H. Jemaa, M.-S. Alouini, and T. Y. Al-Naffouri, 
%  "TeraMIMO:  A  channel  simulator for  wideband  ultra-massive  MIMO  terahertz  communications," 
%  arXivpreprint arXiv:2104.11054, 2021.

% =========================================================================

%% Channnel matrix and initialization parameters

time_realization = p.nSamplesperFrame;
num_freqs_per_subcarr = p.nFreq(1);
num_subcarries = p.nFreq(2);


T_samp = 1/p.BW_sub;

h_LoS_NLoS = cell(p.Qr, p.Qt, num_subcarries);

CH_Info.Bc = zeros(p.Qr, p.Qt, num_subcarries, time_realization);
CH_Info.Tau_rms = zeros(p.Qr, p.Qt, num_subcarries, time_realization);
CH_Info.Tc = zeros(1, num_subcarries);
fd_max = zeros(1, num_subcarries);
CH_Info.fd_max = zeros(1, num_subcarries);

%% Main loop

% Loop over subcarries
for indx_subc = 1:num_subcarries
    
    [CH_Info.Tc(1, indx_subc), fd_max(1, indx_subc)] = ...
        get_Tc_FdMax(p.Fc_sub(1,indx_subc), p.Vel, p.c, p.TCApproxMode);
    CH_Info.fd_max(1, indx_subc) = fd_max(1, indx_subc);
    
    if strcmp(p.channelType,'Multipath') || strcmp(p.channelType,'Multipath+LoS')
        % MP generation
        [ToAClusters, ToARays, ToAMPtot, AoD_NLOS, AoA_NLOS] = ...
            get_ClustersRays(p, indx_subc);
    end
    
    % Loops over SA MIMO arrays
    for mr = 1:p.Mr
        for nr = 1:p.Nr
            for mt = 1:p.Mt
                for nt = 1:p.Nt
                    
                    [AlphaLoS,  AlphaLoSAmp] = get_PathLoss(p, indx_subc, mr, nr, mt, nt, K_abs);
                    [AoD_SV_LOS, AoA_SV_LOS, AoD_BF, AoA_BF] = get_AoD_AoA_LoS(p, mr, nr, mt, nt);
                    if strcmp(p.channelType,'LoS') || strcmp(p.channelType,'Multipath+LoS')
                        Gt_LOS = get_SectorGain(p, AoD_SV_LOS, 'LoS');
                        Gr_LOS = get_SectorGain(p, AoA_SV_LOS, 'LoS');
                        
                        if  strcmp(p.WaveModelAE,'Plane')
                            
                            a_sv_t_LoS = get_ArrayResponse(p, indx_subc, p.Mat, p.Nat, AoD_SV_LOS, 'SV','T');
                            a_bf_t_LoS = get_ArrayResponse(p, indx_subc, p.Mat, p.Nat, AoD_BF, 'BF','T');
                            At_LoS = a_sv_t_LoS.*a_bf_t_LoS;
                            
                            a_sv_r_LoS = get_ArrayResponse(p, indx_subc, p.Mar, p.Nar, AoA_SV_LOS, 'SV','R');
                            a_bf_r_LoS = get_ArrayResponse(p, indx_subc, p.Mar, p.Nar, AoA_BF, 'BF','R');
                            Ar_LoS = a_sv_r_LoS.*a_bf_r_LoS;
                            
                            At_eq_LoS = sum(sum(At_LoS,2),1);
                            Ar_eq_LoS = sum(sum(Ar_LoS,2),1);
                            
                            H_LoS_tmp = Gr_LOS*Gt_LOS*AlphaLoS.*Ar_eq_LoS.*At_eq_LoS;
                        else
                            % AE SWM need to add SV and/or BF
                            error('This wavemodel for AE level isn''t implemented. Supported Options are: SWM/PWM');
                        end
                        H_LoS_tmp1 = zeros(1,1,1,1,num_freqs_per_subcarr);
                        H_LoS_tmp1(1,1,1,1,:) = H_LoS_tmp;
                        
                    end
                    
                    if strcmp(p.channelType,'Multipath') || strcmp(p.channelType,'Multipath+LoS')
                        AlphaNLoS = get_NLoS(p, indx_subc, AlphaLoSAmp, ToAClusters, ToARays, ToAMPtot);
                        
                        Gt_NLOS = get_SectorGain(p, AoD_NLOS, 'NLoS');
                        Gr_NLOS = get_SectorGain(p, AoA_NLOS, 'NLoS');
                        
                        At_eq_NLoS_tot = zeros(p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr);
                        Ar_eq_NLoS_tot = zeros(p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr);
                        
                        for indx_clusters = 1: p.nClustersVec(indx_subc)
                            for indx_rays = 1: p.nRaysVec(indx_subc)
                                
                                a_sv_t_NLoS = ...
                                    get_ArrayResponse(p, indx_subc, p.Mat, p.Nat, AoD_NLOS(:,indx_clusters,indx_rays,:),'SV','T');
                                a_bf_t_NLoS = ...
                                    get_ArrayResponse(p, indx_subc, p.Mat, p.Nat, AoD_BF, 'BF','T');
                                At_NLoS = a_sv_t_NLoS.*a_bf_t_NLoS;
                                
                                a_sv_r_NLoS = ...
                                    get_ArrayResponse(p, indx_subc, p.Mar, p.Nar, AoA_NLOS(:,indx_clusters,indx_rays,:), 'SV','R');
                                a_bf_r_NLoS = ...
                                    get_ArrayResponse(p, indx_subc, p.Mar, p.Nar, AoA_BF, 'BF','R');
                                Ar_NLoS = a_sv_r_NLoS.*a_bf_r_NLoS;
                                
                                At_eq_NLoS = sum(sum(At_NLoS,2),1);
                                Ar_eq_NLoS = sum(sum(Ar_NLoS,2),1);
                                
                                At_eq_NLoS_tot(indx_clusters,indx_rays,:) = At_eq_NLoS;
                                Ar_eq_NLoS_tot(indx_clusters,indx_rays,:) = Ar_eq_NLoS;
                            end
                        end
                        
                        H_NLoS_tmp = sum(sum(Gr_NLOS.*Gt_NLOS.*AlphaNLoS.*Ar_eq_NLoS_tot.*At_eq_NLoS_tot,2),1);
                        H_NLoS_tmp1 = zeros(1,1,1,1,num_freqs_per_subcarr);
                        H_NLoS_tmp1(1,1,1,1,:) = H_NLoS_tmp;
                        
                    end
                    
                    switch p.channelType
                        case 'LoS'
                            
                            error('This model isn''t implemented');
                            
                        case 'Multipath'
                            
                            h_trun1 = FD2IR(p,H_NLoS_tmp1);
                            h_nlos = squeeze(cell2mat(h_trun1));
                            h_nlos_len = length(h_nlos);
                            if h_nlos_len == 0
                                error('Invalid channel generation: There are no MP components!!');
                            end
                            DEL = find(abs(h_nlos) > 0)*T_samp;
                            h_val = h_nlos(abs(h_nlos) > 0);
                            
                        case 'Multipath+LoS'
                            
                            error('This model isn''t implemented');
                            
                        otherwise
                            error('This type of channel isn''t supported!!');
                    end
                    % Get (time, delay) domain channel
                    h_nlos_los_tmp = get_TV_DelayDomain(p, h_val, DEL, ...
                        T_samp, fd_max(1, indx_subc));
                    
                    % Get channel delay-domain statistical parameters: Bc, Tau_rms
                    for indx_time = 1:time_realization
                        PWR = abs(h_nlos_los_tmp(indx_time,abs(h_nlos_los_tmp(indx_time,:)) > 0)).^2;
                        [CH_Info.Bc((mr-1)*p.Nr+nr,(mt-1)*p.Nt+nt, indx_subc, indx_time),...
                            CH_Info.Tau_rms((mr-1)*p.Nr+nr,(mt-1)*p.Nt+nt, indx_subc, indx_time)] ...
                            = get_Bc_TauRms(PWR, DEL, p.BCApproxMode);
                    end
                    
                    % Get final channel response in the (time, delay) domain
                    % Misalignment coefficient is changed from time realization to another
                    h_LoS_NLoS{(mr-1)*p.Nr+nr,(mt-1)*p.Nt+nt,indx_subc,:} = ...
                        repmat(p.h_ma,1,size(h_nlos_los_tmp,2)).*h_nlos_los_tmp;
                end
            end
        end
    end
end

%% Channel response

% h(time, delay) domain implementation
CH_Response.h = h_LoS_NLoS;

% H(time, frequency) domain implementation
CH_Response.H = IR2FR(p, CH_Response.h);

end