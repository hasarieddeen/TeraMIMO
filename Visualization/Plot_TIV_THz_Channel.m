function Plot_TIV_THz_Channel(p, H, h)
% =========================================================================
% -- Function to generate plots for the frequency and delay domain time-invariant THz channel
% =========================================================================

% -- Function: Plot_TIV_THz_Channel(p, H, h)

% -- Input Arguments:
%       p: channel struct that contains the channel parameters
%       H: H(f) TIV frequency domain response
%       h: h(tau) TIV delay domain response

% -- Output Arguments:

%=================================================

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
rx_ind = 1;
tx_ind = 1;
if size(H,3) <= 3
    slec_ind = 1:size(H,3);
else
    slec_ind = [1, ceil(size(H,3)/2) size(H,3)];
end
for indx = 1: length(slec_ind)
    subc_ind = slec_ind(indx);
    f = p.freq(:,subc_ind)/1e9;% Frequency vector (THz)
    
    %%%%%% Magnitude of Frequency Domain %%%%%
    HH = abs(H(rx_ind,tx_ind,subc_ind,:));
    H_temp1 = abs(fft(h{rx_ind,tx_ind,subc_ind},size(H,4)));
    figure();
    plot(f,mag2db(squeeze(HH)),...
        'Color','b');hold on;
    plot(f,mag2db(squeeze(H_temp1)),...
        'Color','r');
    title(strcat('Channel frequency response, subcarrier ',num2str(subc_ind), ...
        ', (Rx,Tx) pair (',num2str(rx_ind),',',num2str(tx_ind),')'));
    xlabel('Frequency (GHz)');
    ylabel('Path Gain (dB)');
    legend('Frequency Domain','From Delay Domain');
    hold off;
    
    %%%%%% Phase of Frequency Domain %%%%%
    HH_ph = angle(H(rx_ind,tx_ind,subc_ind,:));
    H_temp1_ph = angle(fft(h{rx_ind,tx_ind,subc_ind},size(H,4)));
    figure();
    plot(f,rad2deg(squeeze(HH_ph)),...
        'Color','b','LineStyle','-.','Marker','o');hold on;
    plot(f,rad2deg(squeeze(H_temp1_ph)),...
        'Color','r','LineStyle','--','Marker','*');
    title(strcat('Channel phase response, subcarrier ',num2str(subc_ind), ...
        ', (Rx,Tx) pair (',num2str(rx_ind),',',num2str(tx_ind),')'));
    xlabel('Frequency (GHz)');
    ylabel('Phase Response (Degree)');
    legend('Frequency Domain','From Delay Domain');
    hold off;
    
    %%%%%% PDP of Delay Domain %%%%%
    h1 = ifft(H(rx_ind,tx_ind,subc_ind,:));
    t_vec = (1:length(h{rx_ind,tx_ind,subc_ind}))*Ts/1e-9;
    t_vec1 = (1:length(h1))*Ts/1e-9;
    figure();
    stem(t_vec1, squeeze(mag2db(abs(h1))), ...
        'Color','b',...
        'LineStyle','-.','Marker','o');
    hold on;
    stem(t_vec, squeeze(mag2db(abs(h{rx_ind,tx_ind,subc_ind}))),...
        'Color','r', ...
        'LineStyle','--','Marker','*');
    title(strcat('Channel power delay profile, subcarrier ',num2str(subc_ind), ...
        ', (Rx,Tx) pair (',num2str(rx_ind),',',num2str(tx_ind),')'));
    xlabel('$\tau (nsec) $','Interpreter','latex');
    ylabel('PDP (dB)');
    set(gca, 'Ydir', 'reverse');
    legend('From Freq. Domain','Delay Domain');
    hold off;
    
    %%%%%% Magnitude of Delay Domain %%%%%
    h1 = ifft(H(rx_ind,tx_ind,subc_ind,:));
    t_vec = (1:length(h{rx_ind,tx_ind,subc_ind}))*Ts/1e-9;
    t_vec1 = (1:length(h1))*Ts/1e-9;
    figure();
    stem(t_vec1, squeeze(abs(h1)), ...
        'Color','b',...
        'LineStyle','-.','Marker','o');
    hold on;
    stem(t_vec, squeeze(abs(h{rx_ind,tx_ind,subc_ind})),...
        'Color','r', ...
        'LineStyle','--','Marker','*');
    title(strcat('Channel delay-domain, subcarrier ',num2str(subc_ind), ...
        ', (Rx,Tx) pair (',num2str(rx_ind),',',num2str(tx_ind),')'));
    xlabel('$\tau (nsec) $','Interpreter','latex');
    ylabel('Magnitude');
    legend('From Freq. Domain','Delay Domain');
    hold off;
end
