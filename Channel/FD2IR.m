function h = FD2IR(p, H, trun_threshold)
% =========================================================================
% -- Function to convert frequency domain channel into delay domain channel response
%     by adding a threshold to eleminate non-resolvable paths
%     Used method for truncation: All paths below the peak of the channel power by a threshold are eliminated
% =========================================================================

% -- Function: h = FD2IR(p, H, trun_threshold)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       H: Frequency domain channel
%       trun_threshold: A threshold for truncating the NLoS compared to LoS in the case of "Multipath+LoS" channel,
%        where the default truncation threshold is defined in generate_channel_param _TIV or _TV "p.threshold_dB"

% -- Output Arguments:
%       h: Delay domain response

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

num_rx = size(H,1);
num_tx = size(H,2);
num_subcarries = size(H,3);
h = cell(num_rx,num_tx,num_subcarries);

for rx_ind = 1:num_rx
    for tx_ind = 1:num_tx
        for subb_ind = 1:num_subcarries
            
            h_temp = ifft(H(rx_ind,tx_ind,subb_ind,:));
            h_temp_abs = abs(h_temp);
            
            if exist('trun_threshold','var') == 0
                thresh_temp = 10^(p.threshold_dB/20) * max(h_temp_abs);
            else
                thresh_temp = trun_threshold;
            end
            
            h_trueest = [];
            h_trueest(h_temp_abs > thresh_temp) = h_temp(h_temp_abs > thresh_temp);
            h{rx_ind,tx_ind,subb_ind,:} = h_trueest;
            
        end
    end
end

end

