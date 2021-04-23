% =========================================================================
% -- Script to compare absorption coefficient methods
% =========================================================================

% S. Tarboush, H. Sarieddeen, H. Chen, M.-H. Loukil, H. Jemaa, M.-S. Alouini, and T. Y. Al-Naffouri, 
%  "TeraMIMO: A channel simulator for wideband ultra-massive MIMO terahertz communications,"
%  arXivpreprint arXiv:2104.11054, 2021.

% =========================================================================

%% 
clear; clc;
close all;


%% Add Paths
path(pathdef); addpath(pwd); cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd Molecular_Absorption;addpath(genpath(pwd)); cd ..; 
cd Visualization;addpath(genpath(pwd)); cd ..;
cd Examples;

%% Absorption coefficient calculations
% Options: /'Hitran' /'Approx1' /'Approx2'

% To simulate Approx1 over a wider range remove the following condition
% if (p.Fc - p.BW/2) < 0.275e12 || (p.Fc + p.BW/2) > 0.4e12
% the same applies for Approx2

p_upd.absorptionType = 'Hitran';
p_ch = update_channel_param_TIV(p_upd);
p_ch.d_tx_rx;
K_abs_hit = compute_Abs_Coef(p_ch);
save('K_abs_hit.mat','K_abs_hit');

p_upd.absorptionType = 'Approx1';
p_ch = update_channel_param_TIV(p_upd);
p_ch.d_tx_rx;
K_abs_app1 = compute_Abs_Coef(p_ch);
save('K_abs_app1.mat','K_abs_app1');

p_upd.absorptionType = 'Approx2';
p_ch = update_channel_param_TIV(p_upd);
p_ch.d_tx_rx;
K_abs_app2 = compute_Abs_Coef(p_ch);
save('K_abs_app2.mat','K_abs_app2');

%% Visualize results

x1 = 0.2e12/1e9;
x2 = 0.4e12/1e9;
figure();
semilogy(p_ch.freq/1e9,K_abs_hit,'Color','r');grid on; hold on;
semilogy(p_ch.freq/1e9,K_abs_app1,'Color','b');
semilogy(p_ch.freq/1e9,K_abs_app2,'Color','k');
axis([p_ch.freq(1)/1e9 p_ch.freq(end)/1e9 1e-5 0.2]);
xlabel('Frequecny (GHz)');
ylabel('K_{abs} Absorption Coefficient (m^{-1})');
title('Absorption Coefficient Calculation Methods');
legend('HITRAN','Approx1','Approx2','Location','southeast');
plot([x1 x1],get(gca,'ylim'),'Color','k','LineStyle','-.','LineWidth',1.5);
text(x1+5,1e-4,'\leftarrow 275 (GHz)');
plot([x2 x2],get(gca,'ylim'),'Color','k','LineStyle','-.','LineWidth',1.5);
text(x2+5,1e-4,'\leftarrow 400 (GHz)');
dim2 = [.53 .5 .07 .2];
annotation('rectangle',dim2);
annotation('line',[.4 .53],[.7 .5])
annotation('line',[.4 .53],[.9 .7])
ax=axes;
set(ax,'position',[0.2,0.7,0.2,0.2]);
box(ax,'on')
semilogy(p_ch.freq/1e9,K_abs_hit,'r','parent',ax);grid on;hold on;
semilogy(p_ch.freq/1e9,K_abs_app1,'b','parent',ax);
semilogy(p_ch.freq/1e9,K_abs_app2,'k','parent',ax);
set(ax,'xlim',[310,340],'ylim',[1e-3,2e-2]);