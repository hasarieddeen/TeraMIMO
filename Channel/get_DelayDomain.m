function h = get_DelayDomain(p, h_taps_vec, delay_vec, delay_LoS, T_sampling)
% =========================================================================
% -- Function to generate time-invariant delay domain THz channel response
% =========================================================================

% -- Function: h = get_delaydomain(p, h_taps_vec, delay_vec, delay_LoS, T_sampling)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       h_taps_vec: Vector contains the path gain of each tap of the MP channel
%       delay_vec: Vector contains the delay (sec) of each tap of the MP channel
%       delay_LoS: ToA of LoS path
%       T_sampling: Sampling frequency

% -- Output Arguments:
%       h: h(tau) time-invariant delay domain response, a 3D-Cell Array of size (p.Qr, p.Qt, num_subcarries);
%          each cell member is a vector of size (1, max_Delay_spread),
%          i.e. h{MIMO Rx SA, MIMO Tx SA, number of subcarries},
%          and max_Delay_spread refers to the time of arrival, delay, of the last ray in the last cluster

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


Ch_Delay = round(delay_vec/T_sampling);
Ch_Delay_LoS = round(delay_LoS/T_sampling)+1;
Lh = max(Ch_Delay);
num_ch_taps = length(Ch_Delay);

Path_amp = zeros(1,Lh);
for indx_ch_tap = 1:num_ch_taps
    indx_del = Ch_Delay(indx_ch_tap);
    Path_amp(1,indx_del)= Path_amp(1,indx_del)+h_taps_vec(1,indx_ch_tap);
end

h = zeros(1,Lh);
for indx_ch_len = 1:Lh
    if strcmp(p.addrandomness,'On')
        h(1, indx_ch_len) = Path_amp(1, indx_ch_len)*abs(randn+1i*randn)/sqrt(2);
        
        if (strcmp(p.channelType,'LoS')||strcmp(p.channelType,'Multipath+LoS'))
            if indx_ch_len == Ch_Delay_LoS
                % No randomness for LoS component
                h(1,indx_ch_len) = Path_amp(1, indx_ch_len);
            end
        end
    elseif strcmp(p.addrandomness,'Off')
        h(1,indx_ch_len) = Path_amp(1, indx_ch_len);
        % If we remove the randomness option then no need for the for loop
        % Direct implementation: h = Path_amp;
    else
        error('Options are: On/Off');
    end
end

end