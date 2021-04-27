function p  = generate_channel_param_TV()
% =========================================================================
% -- Function to generate the required time-variant THz channel parameters
% =========================================================================

% -- Function: p  = generate_channel_param_TV()

% -- Input Arguments:
%       NULL (no input is required)

% -- Output Arguments:
%       p: Channel struct that contains the channel parameters

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

%% Constants

p.c = 2.9979e8;             % Speed of light in vacuum
p.T0 = 296;                 % Reference temperature (Kelvin)
p.T = 298.15;               % System temperature (Kelvin) 
p.p0 = 1;                   % Standard pressure (atm)
p.p = p.p0;                 % System pressure (atm)
p.PLE = 2;                  % Path loss exponent

%% Define channel type

p.channelType = 'Multipath';
% Options: /'LoS' /'Multipath' /'Multipath+LoS'

%% Transmission parameters (Fc and BW)

p.Fc = 0.3e12;                  % Center frequency of the transmission bandwidth (Hz)
p.BW = 0.01e12;                 % Total channel bandwidth (Hz)
p.Nsubb = 1;                    % Number of subcarriers to divide the total Bandwidth (K-subcarriers)
p.BW_sub = p.BW/p.Nsubb;        % Subcarrier bandwidth (Hz)
p.Nsubc = 2^9;                  % Number of sub-bands in each subcarrier

p.freq = zeros(p.Nsubc, p.Nsubb);
p.nFreq = size(p.freq);

p.fstart_sub = p.Fc-p.BW/2+p.BW_sub/2;
p.fstop_sub = p.fstart_sub+p.BW-p.BW_sub/2;
p.Fc_sub = p.fstart_sub:p.BW_sub:p.fstop_sub;  % Center frequency of each subcarrier (Hz)

for indx_sub = 1:p.Nsubb
    p.fstep = p.BW_sub/p.Nsubc;
    p.fstart = p.Fc_sub(indx_sub)-p.BW_sub/2+p.fstep/2;
    p.fstop = p.fstart+p.BW_sub-p.fstep/2;
    p.freq (:,indx_sub)= (p.fstart:p.fstep:p.fstop).';
end

p.lambdac = p.c./p.Fc_sub;     % Wavelength at center frequency of each subcarrier (m)
p.lambdak = p.c./p.freq;       % Wavelength for each sub-band center frequency inside every subcarrier (m)

%% Absorption Coefficient Calculation

p.absorptionType = 'Hitran';
% Options are: /'Hitran' /'Approx1' /'Approx2'
%  Hitran: Exact absorption coefficient, valid in the frequency band: [0.1-10] THz
% Approx1: First approximation of absorption coefficient, valid in the frequency band: [275-400] GHz
% Approx2: Second approximation of absorption coefficient, valid in the frequency band: [100-450] GHz

if strcmp(p.absorptionType,'Hitran')
    
    p.h = 6.6262e-34;           % Planck constant
    p.Kb = 1.3806e-23;          % Boltzmann constant
    p.R = 8.3144;               % Gas constant
    p.Na = 6.0221e23;           % Avogadro constant
    p.Tstp = 273.15;            % Temperature at standard pressure
    
    % Molecules of transmission medium
    p.molecules = {'N2','O2','H2O','CO2','CH4'};
    % Molecules ratio
	p.moleculesRatio = [76.5450 20.946 1.57 0.033 0.906]/100;
    
	if (sum(p.moleculesRatio) ~= 1)
		error('Sum of molecules ratio should be 100%');
	end
    % Note: for more molecules, check "Molecular_Absorption" folder ---> "Data" subfolder.
elseif strcmp(p.absorptionType,'Approx1')
    
    p.phi = 50;      % Relative humidity
    if (p.Fc - p.BW/2) < 0.275e12 || (p.Fc + p.BW/2) > 0.4e12
        error('Error: approx1, approximation of abs_coef is valid in: [275-400] GHz');
    end
elseif strcmp(p.absorptionType,'Approx2')
        
    p.phi = 50;      % Relative humidity
    if (p.Fc - p.BW/2) < 0.1e12 || (p.Fc + p.BW/2) > 0.45e12
        error('Error: approx2, approximation of abs_coef is valid in: [100-450] GHz');
    end   
else
    error('This method for computing absorption coefficient isn''t implemented, options are: Hitran, Approx1, Approx2');
end

%% UM-MIMO transceiver design (AoSA structure)

% Number of subarrays (SAs)

p.Mt = 1; % Number of transmitter SAs (row)
p.Nt = 1; % Number of transmitter SAs (column)

p.Mr = 1; % Number of Receiver SAs (row)
p.Nr = 1; % Number of Receiver SAs (column)

p.Qt = p.Mt*p.Nt; % Total number of transmitter SAs
p.Qr = p.Mr*p.Nr; % Total number of receiver SAs

% Number of antenna elements (AEs) inside each SA

p.Mat = 1; % Number of transmitter AEs (row) inside each SA
p.Nat = 1; % Number of transmitter AEs (column) inside each SA

p.Mar = 1; % Number of receiver AEs (row) inside each SA
p.Nar = 1; % Number of receiver AEs (column) inside each SA

p.Qat = p.Mat*p.Nat; % Total number of transmitter AEs inside each SA
p.Qar = p.Mar*p.Nar; % Total number of receiver AEs inside each SA

% Define SA Spacing

% Deltaxx is defined from center of SA(i) to center of SA(i+1)
% default 1e-2

p.DeltaMt = 1e-2;    % Spacing between rows of SAs @Tx
p.DeltaNt = 1e-2;    % Spacing between columns of SAs @Tx

p.DeltaMr = 1e-2;    % Spacing between rows of SAs @Rx
p.DeltaNr = 1e-2;    % Spacing between columns of SAs @Rx

% Define AE Spacing

% deltaxx is defined from center of AE(i) to center of AE(i+1)
% default 1e-8

p.deltaMt = 1e-8;    % Spacing between rows of AEs @Tx
p.deltaNt = 1e-8;    % Spacing between columns of AEs @Tx

p.deltaMr = 1e-8;    % Spacing between rows of AEs @Rx
p.deltaNr = 1e-8;    % Spacing between columns of AEs @Rx

%% Design geometry

% Define local/global position and Euler angles
p.positionTx = [0; 0; 0];     % Tx center 3D positions (global coordinates)
p.eulerTx = [0; 0; 0];     % Tx Euler rotation angles, following ZYX intrinsic rotation

p.positionRx = [1; 0; 0];     % Rx center 3D positions (global coordinates)
p.eulerRx = [pi; 0; 0];       % Rx Euler rotation angles, following ZYX intrinsic rotation

% Define rotaion matrix

% Update Tx geometry (global --> local)
% Rotation matrix is real orthogonal (inv = transpose)

% global = Rotm*local;
% local = (Rotm^-1)*global;

p.RotmTx = eul2rotm(p.eulerTx', 'ZYX');
p.unitDirTx = (p.positionRx-p.positionTx)/norm(p.positionRx-p.positionTx);
% unit direction vector (global) from Tx to Rx
p.unitDirLocalTx = p.RotmTx'*p.unitDirTx;
% unit direction vector (local) from Tx to Rx

% Update Rx geometry (global --> local)

p.RotmRx = eul2rotm(p.eulerRx', 'ZYX');
p.unitDirRx = (p.positionTx-p.positionRx)/norm(p.positionTx-p.positionRx);
% unit direction vector (global) from Rx to Tx
p.unitDirLocalRx = p.RotmRx'*p.unitDirRx;
% unit direction vector (local) from Rx to Tx

% Transmission distance (m)
p.d_tx_rx = norm(p.positionTx-p.positionRx);

% Define angle-of-departure/arrival (AoD/AoA) and Analog beamforming (BF) angles
% AoD/AoA (based on local coordinates)

p.AoD = [atan2(p.unitDirLocalTx(2), p.unitDirLocalTx(1)); ...
    acos(p.unitDirLocalTx(3))];
% Physical angle-of-departure AoD [Azimuth; Elevation] (Rad)
p.AoA = [atan2(p.unitDirLocalRx(2), p.unitDirLocalRx(1)); ...
    acos(p.unitDirLocalRx(3))];
% Physical angle-of-arrival AoA [Azimuth; Elevation] (Rad)

% Analog BF @ (Tx, Rx) only on SA level

p.IdealABF = 'On';
if strcmp(p.IdealABF,'On')
    % Ideal analog BF
    p.AoD_BF = p.AoD;
    p.AoA_BF = p.AoA;
else
    % This error leads to add phase uncertainty in phase shifters
    p.err1 = 0.1;
    p.err2 = 0.1;
    p.AoD_BF = p.AoD + p.err1;
    p.AoA_BF = p.AoA + p.err2;
end

%% Parameters of Saleh-Valenzuela (S-V) model for multipath (MP)

% This channel model is a statistical model (many clusters, with many rays whitin each cluster),
% NO RAY-TRACING in this version

% For MP generation, we need to compute: 1) path gain, 2) AoD/AoA, 3) time-of-arrival (ToA)
% Note: Path Gain (dB) = PathLoss_LoS (dB) + Additional Components (dB)

% Define additional path gain components for MP (clusters and rays)
% Double exponential decay profile
% These default parameters, @ p.Fc=0.3 THz, are based on measurements from Ref. [1] and [2]

% Ref [1]: S. Priebe and T. Kurner, "Stochastic modeling of THz indoor radio channels",
% IEEE Trans. Wireless Commun., vol. 12, no. 9, pp. 4445-4455, Sep. 2013.
% Ref [2]: S. Priebe, M. Kannicht, M. Jacob, and T. Kurner, "Ultra broadband indoor channel measurements and calibrated ray tracing propagation modeling at THz frequencies",
% IEEE J. Commun. and Networks, vol.15, no. 6, pp. 547-558, Dec. 2013.
% Ref [3]: C. Lin and G. Y. Li, "Indoor Terahertz Communications: How Many Antenna Arrays Are Needed?"
% in IEEE Transactions on Wireless Commun., vol. 14, no. 6, pp. 3097-3107, June 2015.

% Important note: These MP parameters are only valid to 0.55 THz
if (p.Fc + p.BW/2) < 0.55e12
    p.clusterDecayFactor = 3.12*1e-9;       % Cluster arrival decay factor (nsec)
    p.rayDecayFactor = 0.91*1e-9;           % Ray arrival decay factor (nsec)
else
    p.clusterDecayFactor = 0;       % Cluster arrival decay factor (nsec)
    p.rayDecayFactor = 0;           % Ray arrival decay factor (nsec)
end

% In general, these parameters are frequency-dependent and wall-material dependent
% due to lack of measurements we fix these values in this code

p.clusterDecayFactorVec = repmat(p.clusterDecayFactor,p.nFreq(2),1);    % Vector form for every subcarrier
p.rayDecayFactorVec = repmat(p.rayDecayFactor,p.nFreq(2),1);            % Vector form for every subcarrier

% Define AoD/AoA for clusters and rays
% Model: Total azimuth/elevation = 1) mean azimuth/elevation (clusters) + 2) azimuth/elevation for each ray within a cluster
% 1) Mean azimuth/elevation (clusters) ---> uniform distribution
% 2) Azimuth/elevation for each ray within a cluster ---> Gaussian Mixture Model (GMM)

% These default parameters, @ p.Fc=0.3 THz, are based on measurements from Ref. [1] and [2]
% Gaussian Mixture Model (GMM)
% Three parameters are required: a) means, b) sigmas, c) weights

% a) Means = zero
p.AoDAzimuthMean = [0;0];
p.AoDElevationMean = [0;0];
p.AoAAzimuthMean = [0;0];
p.AoAElevationMean = [0;0];

% b) Sigmas
if (p.Fc + p.BW/2) < 0.55e12
    p.AoDAzimuthSigma = reshape([1.776;20.225],1,1,[]);
    p.AoDElevationSigma = reshape([1.811;12.201],1,1,[]);
    p.AoAAzimuthSigma = reshape([1.895;5.487],1,1,[]);
    p.AoAElevationSigma = reshape([2.198;7.297],1,1,[]);
else
    p.AoDAzimuthSigma = reshape([0;0],1,1,[]);
    p.AoDElevationSigma = reshape([0;0],1,1,[]);
    p.AoAAzimuthSigma = reshape([0;0],1,1,[]);
    p.AoAElevationSigma = reshape([0;0],1,1,[]);
end

p.AoDAzimuthSigmaVec = repmat(p.AoDAzimuthSigma,p.nFreq(2),1);
p.AoDElevationSigmaVec = repmat(p.AoDElevationSigma,p.nFreq(2),1);
p.AoAAzimuthSigmaVec = repmat(p.AoAAzimuthSigma,p.nFreq(2),1);
p.AoAElevationSigmaVec = repmat(p.AoAElevationSigma,p.nFreq(2),1);

% c) Weights
if (p.Fc + p.BW/2) < 0.55e12
    p.AoDAzimuthWeight = [0.475 0.525];
    p.AoDElevationWeight = [0.429 0.571];
    p.AoAAzimuthWeight = [0.583 0.417];
    p.AoAElevationWeight = [0.724 0.276];
else
    % Sum of weights should be equal to 1 (default)
    p.AoDAzimuthWeight = [0.5 0.5];
    p.AoDElevationWeight = [0.5 0.5];
    p.AoAAzimuthWeight = [0.5 0.5];
    p.AoAElevationWeight = [0.5 0.5];  
end

p.AoDAzimuthWeightVec = repmat(p.AoDAzimuthWeight,p.nFreq(2),1);
p.AoDElevationWeightVec = repmat(p.AoDElevationWeight,p.nFreq(2),1);
p.AoAAzimuthWeightVec = repmat(p.AoAAzimuthWeight,p.nFreq(2),1);
p.AoAElevationWeightVec = repmat(p.AoAElevationWeight,p.nFreq(2),1);

% Generate realizations from GMM distribution
% gmdistribution(mean, sigma, weight)
for indx_subb = 1:p.nFreq(2)
    p.gmAoDAzimuth{indx_subb} = gmdistribution(p.AoDAzimuthMean, p.AoDAzimuthSigmaVec(indx_subb,:,:), p.AoDAzimuthWeightVec(indx_subb,:));
    p.gmAoDElevation{indx_subb} = gmdistribution(p.AoDElevationMean, p.AoDElevationSigmaVec(indx_subb,:,:), p.AoDElevationWeightVec(indx_subb,:));
    p.gmAoAAzimuth{indx_subb} = gmdistribution(p.AoAAzimuthMean, p.AoAAzimuthSigmaVec(indx_subb,:,:), p.AoAAzimuthWeightVec(indx_subb,:));
    p.gmAoAElevation{indx_subb} = gmdistribution(p.AoAElevationMean, p.AoAElevationSigmaVec(indx_subb,:,:), p.AoAElevationWeightVec(indx_subb,:));
end

% Define Time-of-Arrival (ToA) for clusters and rays
p.ToAType = 'Exponential';  % Two options: 'Exponential' (default)/'Paraboloid' 
if strcmp(p.ToAType,'Exponential')
    % ToA is exponentially distributed conditioned on the ToA of the previous clusters/rays.
    
    % Inter-cluster and intra-cluster arrival time
    if (p.Fc + p.BW/2) < 0.55e12
        % These default parameters, @ p.Fc=0.3 THz, are based on measurements from Ref. [1] and [2]
        % In general, these parameters are frequency-dependent and wall-material dependent
        % due to lack of measurements we fix these values in this code
        p.clusterArrivalRate = 0.13/1e-9 ;      % Cluster arrival rate (nsec^-1)
        p.rayArrivalRate = 0.37/1e-9;           % Ray arrival rate (nsec^-1)
    else
        p.clusterArrivalRate = Inf;      % Cluster arrival rate (nsec^-1)
        p.rayArrivalRate = Inf;          % Ray arrival rate (nsec^-1)
    end
    p.clusterArrivalRateVec = repmat(p.clusterArrivalRate,p.nFreq(2),1);
    p.rayArrivalRateVec = repmat(p.rayArrivalRate,p.nFreq(2),1);
    
    % Check over realizations, LoS always before MP
    p.tau_LoS = p.d_tx_rx/p.c; % From (center of array @Tx) to (center of array @Rx)
    if (strcmp(p.channelType,'Multipath+LoS'))
        for indx_subb = 1:p.nFreq(2)
            if (p.Fc + p.BW/2) < 0.55e12
                if p.tau_LoS > 1/p.clusterArrivalRateVec(indx_subb)
                    error('Paramters aren''t consistant: delay of LoS should be less than the delay of first cluster');
                end
            end
        end
    end
    
    if (p.Fc + p.BW/2) < 0.55e12
        
        p.nClustersVec = max(1, poissrnd(1./p.clusterArrivalRateVec*1e9));     % Number of clusters in each subcarrier
        p.nRaysVec = max(1, poissrnd(1./p.rayArrivalRateVec*1e9));             % Number of rays in each subcarrier
    else
        
        p.nClustersVec = zeros(p.nFreq(2),1);
        p.nRaysVec = zeros(p.nFreq(2),1);
    end
    
elseif strcmp(p.ToAType,'Paraboloid')
    % Ref [1]: S. Priebe, M. Jacob, and T. Kurner, "AoA, AoD and ToA characteristics of scattered multipath clusters for THz indoor channel modeling",
    % in European Conf. Antennas Propagation, Rome, Italy, Apr. 2011.
    error('Still not implemented in this version');
    %     ToAAzimuthMean = [0.570;-12.525];
    %     ToAAzimuthSigma = reshape([1.196;14.430],1,1,[]);
    %     ToAAzimuthWeight = [0.661;1.108];
    
    %     ToAElevationMean = [0928;-3.386];
    %     ToAElevationSigma = reshape([1.018;11.594],1,1,[]);
    %     ToAElevationWeight = [0.705;0.804];

    %     p.gmToAAzimuth = gmdistribution(ToAAzimuthMean,ToAAzimuthSigma,ToAAzimuthWeight);
    %     p.gmToAElevation = gmdistribution(ToAElevationMean,ToAElevationSigma,ToAElevationWeight);
else
    error('Supported options of ToA are: (Exponential/Paraboloid)');
end

%% Tx, Rx antenna gain (frequency-independent)

% Ref [1]: Constantine A. Balanis, "Antenna Theory: Analysis and Design" 4th edition Eq.(2-26)
% assuming Directivity = Gain i.e. e_cd = e_c * e_d = 1
% No losses, i.e. perfect conduction efficiency & perfect dielectric efficiency
% Also, perfect reflection (mismatch) efficiency e_r (i.e. e_r=1)
% Reflection coefficient equal to zero (matching between source & antenna)

p.psi_azi = deg2rad(rad2deg(2*sqrt(pi)));           % half-power beamwidth in azimuth-plane, examples: rad2deg(2*sqrt(pi)) 27.7 60 120
p.psi_elev = deg2rad(rad2deg(2*sqrt(pi)));          % half-power beamwidth in elevation-plane, examples: rad2deg(2*sqrt(pi)) 27.7 30 60

% input value in (degree), output in (Rad)
% default rad2deg(2*sqrt(pi))

p.GaindBi = pow2db(4*pi/(p.psi_azi*p.psi_elev));    % Gain (dBi)
p.Gain = 10.^(p.GaindBi/10);                        % Gain (scalar value)

%% Spherical and Planar Wave Model

p.WaveModelSA = 'Plane'; %'Sphere'/'Plane'
p.WaveModelAE = 'Plane';  %'Sphere'/'Plane'
% Supported Combinations for this version are SA/AE: Plane/Plane, Sphere/Plane, Sphere/Sphere (NO steering vector)

%% Beam split parameters

p.BeamSplitEffect = 'Off'; % 'Off' 'On'

%% Compute Bc, Tc, and TV delay domain response parameters

p.BCApproxMode = 'Most Popular';
% Options are: /'Most Popular'/'Pessimistic'/'Optimistic'/

switch p.channelType
    case 'LoS'
        p.threshold_dB = -40;
    case 'Multipath'
        p.threshold_dB = -40;
    case 'Multipath+LoS'
        p.threshold_dB = -40;
        p.gap_los_nlos = -40;
        % Ref [1]: H. Yuan, N. Yang, K. Yang, C. Han and J. An, "Hybrid Beamforming for Terahertz Multi-Carrier Systems Over Frequency Selective Fading,"
        % in IEEE Transactions on Communications, vol. 68, no. 10, pp. 6186-6199, Oct. 2020.
        % Sec. II.B. Frequency Selective THz Channel Model
        % The gap between the LOS and NLOS path gains in THz channels (e.g., more than 15 dB on average) is more significant than that in mmWave channels.
    otherwise
        error('This channel type isn''t implemented !!');
end

p.TCApproxMode = 'Geometric Mean';
% Options are: /'Geometric Mean'/'Pessimistic'/'Optimistic'

p.Vel = 40;                     % Velocity of the UE/BS (km/hr)
p.DopplerSpecShape = 'Jakes';   % Options are: 'Jakes', 'Flat'
p.nSamplesperFrame = 60000;     % Length of input signal in samples

%% Misalignment parameters

p.Misalignment = 'Off';  % 'Off' 'On'

% Input parameters for misalignment
% This misalignment model will be updated in the second version of TeraMIMO

p.areaRx = (p.Mr*p.DeltaMr)*(p.Nr*p.DeltaNr);
p.radius_Rx = sqrt(p.areaRx/pi);

p.theta_M = 2*asind(1/pi*(2.782/p.Mt));
p.theta_N = 2*asind(1/pi*(2.782/p.Nt));
p.areaTx = 2*sind(p.theta_M/2)*p.d_tx_rx*2*sind(p.theta_N/2)*p.d_tx_rx;
p.radius_Tx = sqrt(p.areaTx/pi);

p.sigma_s = 0.1;     % Variance of misalignment
p.realizations = p.nSamplesperFrame;

if(strcmp(p.Misalignment, 'On'))
    p.h_ma = get_Misalignment_Fading(p.radius_Rx, p.radius_Tx, p.sigma_s, p.realizations); 
    % Misalignment coefficient
else
    p.h_ma = ones(p.realizations,1); % No misalignment coefficient, default as 1.
end

end
