function H = IR2FR(p, h)
% =========================================================================
% -- Function to convert h(t,tau), the time-variant delay domain channel
%    response, to H(t,f), the time-variant frequency domain channel response,
%    using Fourier transform
% =========================================================================

% -- Function: H = IR2FR(p, h)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       h: h(t,tau) time-variant delay domain response

% -- Output Arguments:
%       H: H(t,f) time-variant frequency domain response, a 5D-Array of size (p.Qr, p.Qt, num_subcarries, time_realization, num_freqs_per_subcarr),
%        i.e. H(MIMO Rx SA, MIMO Tx SA, number of subcarries, number of sample per frame of the input signal, number of sub-bands in each subcarrier).

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


num_rx = size(h,1);
num_tx = size(h,2);
num_subcarries = size(h,3);
time_realization = p.nSamplesperFrame;
num_freqs_per_subcarr = p.Nsub_b;

H = zeros(num_rx, num_tx, num_subcarries, time_realization, num_freqs_per_subcarr);

for idx_Rx  = 1:num_rx
    for idx_Tx = 1:num_tx
        for idx_Sub = 1:num_subcarries
            
            h_temp = h{idx_Rx,idx_Tx,idx_Sub};
            
            for idx_Time = 1:time_realization
                H(idx_Rx,idx_Tx,idx_Sub,idx_Time,:) = fft(h_temp(idx_Time,:),num_freqs_per_subcarr);
            end
            
        end
    end
end
