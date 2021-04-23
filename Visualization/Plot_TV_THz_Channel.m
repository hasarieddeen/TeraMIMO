function Plot_TV_THz_Channel(p, H, h, fd_max_subcarriers ,num_of_plots)

% =========================================================================
% -- Function to generate plots for the frequency and delay domain time-variant THz channel,
%   and to check the time correlation between channel taps
% =========================================================================

% -- Function: Plot_TV_THz_Channel(p, H, h, fd_max_subcarriers, num_of_plots)

% -- Input Arguments:
%       p: channel struct that contains the channel parameters
%       H: H(t,f) TV frequency domain response
%       h: h(t,tau) TV delay domain response
%       fd_max_subcarriers: maximum Doppler shift for a subcarrier
%       num_of_plots: number of plots for time observations

% -- Output Arguments:

% =========================================================================

% -- (c) 2021 Simon Tarboush, Hadi Sarieddeen, Hui Chen, 
%             Mohamed Habib Loukil, Hakim Jemaa, 
%             Mohamed-Slim Alouini, Tareq Y. Al-Naffouri

% -- e-mail: simon.w.tarboush@gmail.com; hadi.sarieddeen@kaust.edu.sa; hui.chen@kaust.edu.sa;
%            mohamedhabib.loukil@kaust.edu.sa; hakim.jemaa@kaust.edu.sa;
%            slim.alouini@kaust.edu.sa; tareq.alnaffouri@kaust.edu.sa

% =========================================================================

% S. Tarboush, H. Sarieddeen, H. Chen, M.-H. Loukil, H. Jemaa, M.-S. Alouini, and T. Y. Al-Naffouri, 
%  "TeraMIMO: A channel simulator for wideband ultra-massive MIMO terahertz communications,"
%  arXivpreprint arXiv:2104.11054, 2021.

% =========================================================================

Ts = 1/p.BW_sub;
num_freq_subband = p.Nsubc;
vec_time = randi(p.nSamplesperFrame,1,num_of_plots);

for rx_ind = 1 : size(h,1)
    for tx_ind = 1: size(h,2)
        for subb_ind = 1: size(h,3)
            
            h_temp = h{rx_ind,tx_ind,subb_ind};
            
            for time_ind = 1: length(vec_time)
                
                f = p.freq(:,subb_ind)/1e12;% Frequency vector (THz)
                
                %%%%%% Magnitude of Frequency Domain %%%%%
                HH = H(rx_ind,tx_ind,subb_ind,vec_time(time_ind),:);
                H_temp1 = fft(h_temp(vec_time(time_ind),:),num_freq_subband);
                figure();
                plot(f,mag2db(squeeze(abs(HH))),...
                    'Color','b','Marker','o');hold on;
                plot(f,mag2db(abs(H_temp1)),...
                    'Color','r','Marker','x');
                title('Channel Frequency Response');
                xlabel('Frequency (THz)');
                ylabel('Path Gain (dB)');
                legend('Freq. Domain','From Delay Domain');
                hold off;
                
                %%%%%% Magnitude of Delay Domain %%%%%
                h1 = ifft(H(rx_ind,tx_ind,subb_ind,vec_time(time_ind),:));
                t_vec1 = (1:length(h1))*Ts/1e-9;
                t_vec = (1:length(h_temp(vec_time(time_ind),:)))*Ts/1e-9;
                figure();
                stem(t_vec1, squeeze(abs(h1)), ...
                    'Color','b',...
                    'LineStyle','-.','Marker','o');
                hold on;
                stem(t_vec, squeeze(abs(h_temp(vec_time(time_ind),:))),...
                    'Color','r', ...
                    'LineStyle','--','Marker','*');
                title('Delay-Domain Channel'); % Impulse response
                xlabel('$\tau (nsec) $','Interpreter','latex');
                ylabel('Magnitude');
                legend('From Freq. Domain','Delay Domain');
                hold off;
                
            end
        end
    end
end
%%%%%% Correlation between paths %%%%%
% Select maximum path gain only
t = 0:Ts:Ts*(p.nSamplesperFrame-1);
for rx_ind = 1 : size(h,1)
    for tx_ind = 1: size(h,2)
        for subb_ind = 1: size(h,3)
            
            h_temp2 = h{rx_ind,tx_ind,subb_ind};
            % find max path gain
            [~, indx_h_max_g ] = max(abs(h_temp2(1,:)));
            h_max_g = squeeze(h_temp2(:,indx_h_max_g));
            
            Acn = xcorr(h_max_g,'biased');
            Acn_pos = Acn(p.nSamplesperFrame:end).';
            Acn_pos = Acn_pos/max(Acn_pos);
                       
            switch p.DopplerSpecShape
                case 'Jakes'
                    Ac_th = besselj(0,2*pi*fd_max_subcarriers(1,subb_ind)*t);
                case 'Flat'
                    Ac_th = sinc(2*fd_max_subcarriers(1,subb_ind)*t);
                otherwise
                    error('The Selected Doppler Spectrum Shape isn''t supported, Options: Jakes, Flat, ...')
            end
            
            figure();
            plot(t,real(Acn_pos),'b');hold on;
            plot(t,real(Ac_th),'r');
            title('Autocorrelation Function for Time-Variant Channel');
            xlabel('time (msec)');ylabel('Amplitude');
            legend('Simulation','Theory');
            grid on;hold off;
        end
    end
end

