% =========================================================================
% -- Script to generate a time-variant (TV) THz channel
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

%%
clear; clc;
close all;


%% Add Paths
path(pathdef); addpath(pwd);
cd Channel; addpath(genpath(pwd)); cd ..;
cd Molecular_Absorption;addpath(genpath(pwd)); cd ..; 
cd Visualization;addpath(genpath(pwd)); cd ..; 


%% Initialize Parameters
p_ch = generate_channel_param_TV();


%% Calculation of Absorption Coefficient 
K_abs = compute_Abs_Coef(p_ch);


%% Call Channel
[CH_Response, CH_Info] = channel_TV(p_ch, K_abs); 


%% Visualize Channel
num_ch_obsrv = 3;
Plot_TV_THz_Channel(p_ch, CH_Response.H, CH_Response.h, CH_Info.fd_max, num_ch_obsrv);

