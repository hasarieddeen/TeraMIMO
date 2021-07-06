% =========================================================================
% -- Script to visualize beam split (uniform linear array: ULA)
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

%% Initialize parameters
p_ch = generate_channel_param_TIV();

% p.Fc = 0.3e12;
% p.BW = 0.05e12;
% p.Nsub_c = 1; 
% p.Nsub_b = 5; 
% p.Mat = 8;
% p.deltaMt = p.lambdac(1)/2;
% p.BeamSplitEffect = 'On';

% Beamforming angles
AoD_BF1 = [0 ; deg2rad(45)];
AoD_BF2 = [0 ; deg2rad(105)];
AoD_BF3 = [0 ; deg2rad(135)];

theta = 0:0.2:180;
phi = zeros(1,length(theta));

theta_all = [phi; theta];
len_phi_vec = size(theta_all,2);

%%

G1 = zeros(p_ch.nFreq(1), len_phi_vec);
G2 = zeros(p_ch.nFreq(1), len_phi_vec);
G3 = zeros(p_ch.nFreq(1), len_phi_vec);
for indx_phi = 1:len_phi_vec

    a_sv1 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, deg2rad(theta_all(:,indx_phi)), 'SV','T');
    a_bf1 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, AoD_BF1, 'BF','T');
    
    a_sv2 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, deg2rad(theta_all(:,indx_phi)), 'SV','T');
    a_bf2 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, AoD_BF2, 'BF','T');
    
    a_sv3 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, deg2rad(theta_all(:,indx_phi)), 'SV','T');
    a_bf3 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, AoD_BF3, 'BF','T');
    
    G1(:, indx_phi) = abs(squeeze(sum(sum(a_sv1.*a_bf1,2),1)).')/sqrt(p_ch.Mat*p_ch.Nat);
    G2(:, indx_phi) = abs(squeeze(sum(sum(a_sv2.*a_bf2,2),1)).')/sqrt(p_ch.Mat*p_ch.Nat);
    G3(:, indx_phi) = abs(squeeze(sum(sum(a_sv3.*a_bf3,2),1)).')/sqrt(p_ch.Mat*p_ch.Nat);
    
end
figure();plot(theta_all(2,:), G1);
figure();plot(theta_all(2,:), G2);
figure();plot(theta_all(2,:), G3);

%% Polar Plot

figure(); polarplot(deg2rad(theta_all(2,:)),G1(3,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G2(3,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G3(3,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G1(1,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G2(2,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G3(2,:));

thetalim([0 180]);
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'left';
%legend('SA_1/f_c (\theta = -45^\circ)', 'SA_2/f_c (\theta = 15^\circ)','SA_3/f_c (\theta = 45^\circ)', ...
%    'SA_1/f_1 (\theta = -45^\circ)',  'SA_2/f_2 (\theta = 15^\circ)', 'SA_3/f_2 (\theta = 45^\circ)', ...
%    'Location', 'southoutside');

%% Change number of antennas by updating parameters

p_ch.Mat = 32;
p_ch = update_channel_param_TIV(p_ch);

%%

G1 = zeros(p_ch.nFreq(1), len_phi_vec);
G2 = zeros(p_ch.nFreq(1), len_phi_vec);
G3 = zeros(p_ch.nFreq(1), len_phi_vec);
for indx_phi = 1:len_phi_vec

    a_sv1 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, deg2rad(theta_all(:,indx_phi)), 'SV','T');
    a_bf1 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, AoD_BF1, 'BF','T');
    
    a_sv2 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, deg2rad(theta_all(:,indx_phi)), 'SV','T');
    a_bf2 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, AoD_BF2, 'BF','T');
    
    a_sv3 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, deg2rad(theta_all(:,indx_phi)), 'SV','T');
    a_bf3 = get_ArrayResponse(p_ch, 1, p_ch.Mat, p_ch.Nat, AoD_BF3, 'BF','T');
    
    G1(:, indx_phi) = abs(squeeze(sum(sum(a_sv1.*a_bf1,2),1)).')/sqrt(p_ch.Mat*p_ch.Nat);
    G2(:, indx_phi) = abs(squeeze(sum(sum(a_sv2.*a_bf2,2),1)).')/sqrt(p_ch.Mat*p_ch.Nat);
    G3(:, indx_phi) = abs(squeeze(sum(sum(a_sv3.*a_bf3,2),1)).')/sqrt(p_ch.Mat*p_ch.Nat);
    
end
figure();plot(theta_all(2,:), G1);
figure();plot(theta_all(2,:), G2);
figure();plot(theta_all(2,:), G3);

%% Polar plot

figure();polarplot(deg2rad(theta_all(2,:)),G1(3,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G2(3,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G3(3,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G1(1,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G2(2,:));
hold on;polarplot(deg2rad(theta_all(2,:)),G3(2,:));

thetalim([0 180]);
pax = gca;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'left';
%legend('SA_1/f_c (\theta = -45^\circ)', 'SA_2/f_c (\theta = 15^\circ)','SA_3/f_c (\theta = 45^\circ)', ...
%    'SA_1/f_1 (\theta = -45^\circ)',  'SA_2/f_2 (\theta = 15^\circ)', 'SA_3/f_2 (\theta = 45^\circ)', ...
%    'Location', 'southoutside');