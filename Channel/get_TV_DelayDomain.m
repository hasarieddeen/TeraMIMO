function h = get_TV_DelayDomain(p, h_taps_vec, delay_vec, T_sampling, fd_max)
% =========================================================================
% -- Function to generate a time-variant delay domain THz channel response for
%    multipath components following a specific Doppler spectrum
% =========================================================================

% -- Function: h = get_tv_delaydomain(p, h_taps_vec, delay_vec, T_sampling, fd_max)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       h_taps_vec: Vector containing the path gain of each tap of the MP channel
%       delay_vec: Vector containing the delay (sec) of each tap of the MP channel
%       T_sampling: Sampling frequency.
%       fd_max: Maximum Doppler shift.

% -- Output Arguments:
%       h: h(t,tau) time-variant delay domain response, a 3D-Cell Array of size (p.Qr, p.Qt, num_subcarries)
%          and each cell member is a matrix of size (time_realization, max_Delay_spread),
%          i.e. h{MIMO Rx SA, MIMO Tx SA, number of subcarries};
%          max_Delay_spread refers to the time of arrival, delay, of the last ray in the last cluster

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
Lh = max(Ch_Delay);
num_ch_taps = length(Ch_Delay);

FrameLen = p.nSamplesperFrame;
% FrameLen*T_sampling = time duration of the signal
numBLKperTC = FrameLen*T_sampling*fd_max;
% Number of elements in the frame over which the channel will give significant variance based on fd_max
num_changes = floor(numBLKperTC);
% floor(numBLKperTC) to ovoid f = fmax

Path_amp = zeros(1,Lh);
for indx_ch_tap = 1:num_ch_taps
    indx_del = Ch_Delay(indx_ch_tap);
    Path_amp(1,indx_del)= Path_amp(1,indx_del)+h_taps_vec(1,indx_ch_tap);
end

% Initialize channel h(t,tau) (time,delay) domain
h = zeros(FrameLen,Lh);

nu_freq = -num_changes:num_changes;
nu_len = length(nu_freq);

if num_changes == 0
    S_nu = 1;% PSD of Doppler spectrum
else
    switch p.DopplerSpecShape
        case 'Jakes'
            
            S_nu = 1./(pi*fd_max*sqrt(1-(nu_freq/numBLKperTC).^2));
        case 'Flat'
            
            S_nu=ones(1,nu_len)/(2*fd_max);
        otherwise
            
            error('The selected Doppler spectrum shape isn''t supported. Options: Jakes, Flat, ...')
    end
    S_nu = S_nu/sum(S_nu);
    
    figure(),plot(nu_freq./numBLKperTC*fd_max,S_nu,'r--x');xlabel('Frequency (Hz)'); ylabel('Amplitude');
    title(strcat('PSD of Doppler Spectrum, fd-max = ',num2str(fd_max/1e3),' (KHz)'));
    
end

% Generate h(nu,tau) (doppler, delay) domain
h_del_doppler = zeros(nu_len,Lh);

for indx_ch_len = 1:Lh
    for indx_nu_doppler = 1:nu_len
        h_del_doppler(indx_nu_doppler,indx_ch_len) = Path_amp(1,indx_ch_len)*sqrt(S_nu(nu_len));
    end
end

% Go back to h(t,tau) (time, delay) domain by using inverse Fourier transform
for indx_ch_len = 1:Lh
    h(:,indx_ch_len) = exp(1j*2*pi/FrameLen*(0:FrameLen-1).'*(-num_changes:num_changes)) ...
        *h_del_doppler(:,indx_ch_len);
end

end