function p  = update_channel_param_TIV(q)
% =========================================================================
% -- Function to update the required time-invariant THz channel parameters
% =========================================================================

% -- Function: p = update_channel_param_TIV()

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

if isfield(q,'T')
    p.T = q.T;              % System temperature (Kelvin) 
else 
    p.T = p.T0;
end
p.p0 = 1;                   % Standard pressure (atm)
if isfield(q,'p')
    p.p = q.p;              % System pressure (atm)
else
    p.p = p.p0;
end
if isfield(q,'PLE')
    p.PLE = q.PLE;
else
    p.PLE = 2;              % Path loss exponent
end

%% Define channel type

if isfield(q,'channelType')
    p.channelType = q.channelType;
else
    p.channelType = 'LoS';
    % Options: /'LoS' /'Multipath' /'Multipath+LoS'
end
if isfield(q,'addrandomness')
    p.addrandomness = q.addrandomness;
else
    p.addrandomness = 'Off';    % 'On', default is 'Off'
end

%% Transmission parameters (Fc & BW)

if isfield(q,'Fc')
    p.Fc = q.Fc;         % Center frequency of the transmission bandwidth (Hz)
else
    p.Fc = 0.275e12;
end
if isfield(q,'BW')
    p.BW = q.BW;         % Total channel bandwidth (Hz)
else
    p.BW = 0.01e12;
end
if isfield(q,'Nsubb')
    p.Nsub_c = q.Nsub_c;   % Number of subcarriers to divide the total Bandwidth (K-subcarriers)  
else
    p.Nsub_c = 1;
end
p.BW_sub = p.BW/p.Nsub_c; % Subcarrier bandwidth (Hz)
if isfield(q,'Nsubc')
    p.Nsub_b = q.Nsub_b;   % Number of sub-bands in each subcarrier
else
    p.Nsub_b = 2^9;
end
p.freq = zeros(p.Nsub_b, p.Nsub_c);
p.nFreq = size(p.freq);

p.fstart_sub = p.Fc-p.BW/2+p.BW_sub/2;
p.fstop_sub = p.fstart_sub+p.BW-p.BW_sub/2;
p.Fc_sub = p.fstart_sub:p.BW_sub:p.fstop_sub;  % Center frequency of each subcarrier (Hz)

for indx_subc = 1:p.nFreq(2)
    p.fstep = p.BW_sub/p.Nsub_b;
    p.fstart = p.Fc_sub(indx_subc)-p.BW_sub/2+p.fstep/2;
    p.fstop = p.fstart+p.BW_sub-p.fstep/2;
    p.freq (:,indx_subc)= (p.fstart:p.fstep:p.fstop).';
end

p.lambdac = p.c./p.Fc_sub;     % Wavelength at center frequency of each subcarrier (m)
p.lambdak = p.c./p.freq;       % Wavelength for each sub-band center frequency inside every subcarrier (m)

%% Absorption Coefficient Calculation

if isfield(q,'absorptionType')
    p.absorptionType = q.absorptionType;
else
    p.absorptionType = 'Hitran';
    % Options are: /'Hitran' /'Approx1' /'Approx2'
    %  Hitran: Exact absorption coefficient, valid in the frequency band: [0.1-10] THz
    % Approx1: First approximation of absorption coefficient, valid in the frequency band: [275-400] GHz
    % Approx2: Second approximation of absorption coefficient, valid in the frequency band: [100-450] GHz
end

if strcmp(p.absorptionType,'Hitran')
    
    p.h = 6.6262e-34;           % Planck constant
    p.Kb = 1.3806e-23;          % Boltzmann constant
    p.R = 8.3144;               % Gas constant
    p.Na = 6.0221e23;           % Avogadro constant
    p.Tstp = 273.15;            % Temperature at standard Pressure
    
    if (isfield(q,'molecules') && isfield(q,'moleculesRatio'))
        
        p.molecules =  q.molecules;
        p.moleculesRatio = q.moleculesRatio;
		if (sum(p.moleculesRatio) ~= 1)
			error('Sum of molecules ratio should be 100%');
		end
	else
        
        % Molecules of transmission medium
        p.molecules = {'N2','O2','H2O','CO2','CH4'};
        % Molecules ratio
		p.moleculesRatio = [76.5450 20.946 1.57 0.033 0.906]/100;
    
		if (sum(p.moleculesRatio) ~= 1)
			error('Sum of molecules ratio should be 100%');
		end
        % Note: for more molecules, check "Molecular_Absorption" folder ---> "Data" subfolder.
    end
    
elseif strcmp(p.absorptionType,'Approx1')
    if isfield(q,'phi')
        p.phi = q.phi;            % Relative humidity
    else
        p.phi = 50;
    end
    if (p.Fc - p.BW/2) < 0.275e12 || (p.Fc + p.BW/2) > 0.4e12
        error('Error: approx1, approximation of abs_coef is valid in: [275-400] GHz');
    end
elseif strcmp(p.absorptionType,'Approx2')
    if isfield(q,'phi')
        p.phi = q.phi;            % Relative humidity
    else
        p.phi = 50;
    end
    if (p.Fc - p.BW/2) < 0.1e12 || (p.Fc + p.BW/2) > 0.45e12
        error('Error: approx2, approximation of abs_coef is valid in: [100-450] GHz');
    end
else
    error('This method for computing absorption coefficient isn''t implemented, options are: Hitran, Approx1, Approx2');
end

%% UM-MIMO transceiver design (AoSA structure)

% Number of subarrays (SAs)

if isfield(q,'Mt')
    p.Mt = q.Mt; % Number of transmitter SAs (row)
else
    p.Mt = 1;
end
if isfield(q,'Nt')
    p.Nt = q.Nt; % Number of transmitter SAs (column)
else
    p.Nt = 1;
end

if isfield(q,'Mr')
    p.Mr = q.Mr; % Number of Receiver SAs (row)
else
    p.Mr = 1;
end
if isfield(q,'Nr')
    p.Nr = q.Nr; % Number of Receiver SAs (column)
else
    p.Nr = 1; 
end

p.Qt = p.Mt*p.Nt; % Total number of transmitter SAs
p.Qr = p.Mr*p.Nr; % Total number of receiver SAs

% Number of antenna elements (AEs) inside each SA

if isfield(q,'Mat')
    p.Mat = q.Mat; % Number of transmitter AEs (row) inside each SA
else
    p.Mat = 1;
end
if isfield(q,'Nat')
    p.Nat = q.Nat; % Number of transmitter AEs (column) inside each SA
else
    p.Nat = 1;
end

if isfield(q,'Mar')
    p.Mar = q.Mar; % Number of receiver AEs (row) inside each SA
else
    p.Mar = 1;
end
if isfield(q,'Nar')
    p.Nar = q.Nar; % Number of receiver AEs (column) inside each SA
else
    p.Nar = 1;
end

p.Qat = p.Mat*p.Nat; % Total number of transmitter AEs inside each SA
p.Qar = p.Mar*p.Nar; % Total number of receiver AEs inside each SA

% Define SA Spacing

% Deltaxx is defined from center of SA(i) to center of SA(i+1)
% default 1e-2

if isfield(q,'DeltaMt')
    p.DeltaMt = q.DeltaMt;    % Spacing between rows of SAs @Tx
else
    p.DeltaMt = 1e-2;
end
if isfield(q,'DeltaNt')
    p.DeltaNt = q.DeltaNt;    % Spacing between columns of SAs @Tx
else
    p.DeltaNt = 1e-2;
end

if isfield(q,'DeltaMr')
    p.DeltaMr = q.DeltaMr;    % Spacing between rows of SAs @Rx
else
    p.DeltaMr = 1e-2;
end
if isfield(q,'DeltaNr')
    p.DeltaNr = q.DeltaNr;    % Spacing between columns of SAs @Rx
else
    p.DeltaNr = 1e-2;
end

% Define AE Spacing

% deltaxx is defined from center of AE(i) to center of AE(i+1)
% default 1e-8

if isfield(q,'deltaMt') 
    p.deltaMt = q.deltaMt;    % Spacing between rows of AEs @Tx
else
    p.deltaMt = 1e-8;
end
if isfield(q,'deltaNt') 
    p.deltaNt = q.deltaNt;    % Spacing between columns of AEs @Tx
else
    p.deltaNt = 1e-8;
end

if isfield(q,'deltaMr') 
    p.deltaMr = q.deltaMr;    % Spacing between rows of AEs @Rx
else
    p.deltaMr = 1e-8;
end
if isfield(q,'deltaNr') 
    p.deltaNr = q.deltaNr;    % Spacing between columns of AEs @Rx
else
    p.deltaNr = 1e-8;
end

%% Geometry design

% Define local/global position and Euler angles
if isfield(q,'positionTx') 
    p.positionTx = q.positionTx;     % Tx center 3D positions (global coordinates)
else
    p.positionTx = [0; 0; 0];
end
if isfield(q,'eulerTx')
    p.eulerTx = q.eulerTx;           % Tx Euler rotation angles, following ZYX intrinsic rotation
else
    p.eulerTx = [0; 0; 0];
end

if isfield(q,'positionRx')
    p.positionRx = q.positionRx;     % Rx center 3D positions (global coordinates)
else
    p.positionRx = [1; 0; 0];
end
if isfield(q,'eulerRx')
    p.eulerRx = q.eulerRx;           % Rx Euler rotation angles, following ZYX intrinsic rotation
else
    p.eulerRx = [pi; 0; 0];
end

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

% Define angle-of-departure/arrival (AoD/AoA) && Analog beamforming (BF) angles
% AoD/AoA (based on local coordinates)

p.AoD = [atan2(p.unitDirLocalTx(2), p.unitDirLocalTx(1)); ...
    acos(p.unitDirLocalTx(3))];
% Physical angle-of-departure AoD [Azimuth; Elevation] (Rad)
p.AoA = [atan2(p.unitDirLocalRx(2), p.unitDirLocalRx(1)); ...
    acos(p.unitDirLocalRx(3))];
% Physical angle-of-arrival AoA [Azimuth; Elevation] (Rad)

% Analog BF @ (Tx, Rx) only on SA level

if isfield(q,'IdealABF')
    p.IdealABF = q.IdealABF;
else
    p.IdealABF = 'On';
end

if strcmp(p.IdealABF,'On')
    % Ideal analog BF
    p.AoD_BF = p.AoD;
    p.AoA_BF = p.AoA;
else
    if isfield(q,'err1')
        % This error adds phase uncertainty in phase shifters
        p.err1 = q.err1;
    else
        p.err1 = 0.1;
    end
    if isfield(q,'err2')
        % This error adds phase uncertainty in phase shifters
        p.err2 = q.err2;
    else
        p.err2 = 0.1;
    end
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

% Important note: These MP parameters are only valid for 0.55 THz
if (p.Fc + p.BW/2) < 0.55e12
    if isfield(q,'clusterDecayFactor')
        p.clusterDecayFactor = q.clusterDecayFactor;       % Cluster decay factor (nsec)
    else
        p.clusterDecayFactor = 3.12*1e-9;
    end
    if isfield(q,'rayDecayFactor') 
        p.rayDecayFactor = q.rayDecayFactor;               % Ray decay factor (nsec)
    else
        p.rayDecayFactor = 0.91*1e-9;
    end
else
    p.clusterDecayFactor = 0;       % Cluster arrival decay factor (nsec)
    p.rayDecayFactor = 0;           % Ray arrival decay factor (nsec)
end

% In general, these parameters are frequency-dependent and wall-material dependent
% due to lack of measurements we fix these values in this code

p.clusterDecayFactorVec = repmat(p.clusterDecayFactor,p.nFreq(2),1);    % Vector form for every subcarrier
p.rayDecayFactorVec = repmat(p.rayDecayFactor,p.nFreq(2),1);            % Vector form for every subcarrier

% Define AoD/AoA for clusters and rays
% Model: total Azi/Elev = 1) mean Azi/Elev (clusters) + 2) Azi/Elev for each ray within a cluster
% 1) Mean azimuth/elevation (clusters) ---> dniform distribution
% 2) Azimuth/elevation for each ray within a cluster ---> Gaussian Mixture Model (GMM)

% These default parameters, @ p.Fc=0.3 THz, are based on measurements from Ref. [1] and [2]
% Gaussian Mixture Model (GMM)
% Three parameters are required: a) means, b) sigmas, c) weights

% a) Means = zero
if isfield(q,'AoDAzimuthMean')
    p.AoDAzimuthMean = q.AoDAzimuthMean;
else
    p.AoDAzimuthMean = [0;0];
end
if isfield(q,'AoDElevationMean') 
    p.AoDElevationMean = q.AoDElevationMean;
else
    p.AoDElevationMean = [0;0];
end
if isfield(q,'AoAAzimuthMean') 
    p.AoAAzimuthMean = q.AoAAzimuthMean;
else
    p.AoAAzimuthMean = [0;0];
end
if isfield(q,'AoAElevationMean') 
    p.AoAElevationMean = q.AoAElevationMean;
else
    p.AoAElevationMean = [0;0];
end

% b) Sigmas
if (p.Fc + p.BW/2) < 0.55e12
    if isfield(q,'AoDAzimuthSigma')
        p.AoDAzimuthSigma = q.AoDAzimuthSigma;
    else
        p.AoDAzimuthSigma = reshape([1.776;20.225],1,1,[]);
    end
    if isfield(q,'AoDElevationSigma')
        p.AoDElevationSigma = q.AoDElevationSigma;
    else
        p.AoDElevationSigma = reshape([1.811;12.201],1,1,[]);
    end
    if isfield(q,'AoAAzimuthSigma')
        p.AoAAzimuthSigma = q.AoAAzimuthSigma;
    else
        p.AoAAzimuthSigma = reshape([1.895;5.487],1,1,[]);
    end
    if isfield(q,'AoAElevationSigma')
        p.AoAElevationSigma = q.AoAElevationSigma;
    else
        p.AoAElevationSigma = reshape([2.198;7.297],1,1,[]);
    end
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
    if isfield(q,'AoDAzimuthWeight') 
        p.AoDAzimuthWeight = q.AoDAzimuthWeight;
    else
        p.AoDAzimuthWeight = [0.475 0.525];
    end
    if isfield(q,'AoDElevationWeight')
        p.AoDElevationWeight = q.AoDElevationWeight;
    else
        p.AoDElevationWeight = [0.429 0.571];
    end
    if isfield(q,'AoAAzimuthWeight')
        p.AoAAzimuthWeight = q.AoAAzimuthWeight;
    else
        p.AoAAzimuthWeight = [0.583 0.417];
    end
    if isfield(q,'AoAElevationWeight')
        p.AoAElevationWeight = q.AoAElevationWeight;
    else
        p.AoAElevationWeight = [0.724 0.276];
    end
else
    % sum of weights should be 1 (default)
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
for indx_subc = 1:p.nFreq(2)
    p.gmAoDAzimuth{indx_subc} = gmdistribution(p.AoDAzimuthMean, p.AoDAzimuthSigmaVec(indx_subc,:,:), p.AoDAzimuthWeightVec(indx_subc,:));
    p.gmAoDElevation{indx_subc} = gmdistribution(p.AoDElevationMean, p.AoDElevationSigmaVec(indx_subc,:,:), p.AoDElevationWeightVec(indx_subc,:));
    p.gmAoAAzimuth{indx_subc} = gmdistribution(p.AoAAzimuthMean, p.AoAAzimuthSigmaVec(indx_subc,:,:), p.AoAAzimuthWeightVec(indx_subc,:));
    p.gmAoAElevation{indx_subc} = gmdistribution(p.AoAElevationMean, p.AoAElevationSigmaVec(indx_subc,:,:), p.AoAElevationWeightVec(indx_subc,:));
end

% Define time of arrival (ToA) for clusters and rays
if isfield(q,'ToAType') 
    p.ToAType = q.ToAType;  
else
    p.ToAType = 'Exponential';  % Two options: 'Exponential' (default)/'Paraboloid' 
end
if strcmp(p.ToAType,'Exponential')
    % ToA is exponentially distributed conditioned on the ToA of the previous clusters/rays.
    
    % Inter-Cluster and Intra-Cluster Arrival Time
    if (p.Fc + p.BW/2) < 0.55e12
        % These default parameters, @ p.Fc=0.3 THz, are based on measurements from Ref. [1] and [2]
        % In general, these parameters are frequency-dependent and wall-material dependent
        % due to lack of measurements we repeat these values
        if isfield(q,'clusterArrivalRate')
            p.clusterArrivalRate = q.clusterArrivalRate;      % Cluster arrival rate (nsec^-1)
        else
            p.clusterArrivalRate = 0.13/1e-9;
        end
        if p.d_tx_rx > p.c/p.clusterArrivalRate
            % Update the cluster arrival rate so that we don't show an
            % error but we generate channel realizations where their
            % parameters aren't based on measuremnts
            % here we say that always the first cluster comes after 3 nsec
            p.clusterArrivalRate = 1/(p.d_tx_rx/p.c+3e-9) ;
        end
        if isfield(q,'rayArrivalRate')
            p.rayArrivalRate = q.rayArrivalRate;              % Ray arrival rate (nsec^-1)
        else
            p.rayArrivalRate = 0.37/1e-9;
        end
    else
        p.clusterArrivalRate = Inf;      % Cluster arrival rate (nsec^-1)
        p.rayArrivalRate = Inf;          % Ray arrival rate (nsec^-1)
    end
    p.clusterArrivalRateVec = repmat(p.clusterArrivalRate,p.nFreq(2),1);
    p.rayArrivalRateVec = repmat(p.rayArrivalRate,p.nFreq(2),1);
    
    % Check over realizations, LoS always before MP
    p.tau_LoS = p.d_tx_rx/p.c; % From (center of array @Tx) to (center of array @Rx)
    if (strcmp(p.channelType,'Multipath+LoS'))
        for indx_subc = 1:p.nFreq(2)
            if (p.Fc + p.BW/2) < 0.55e12
                if p.tau_LoS > 1/p.clusterArrivalRateVec(indx_subc)
                    error('Paramters aren''t consistant: delay of LoS should be less than the delay of First Cluster');
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

if isfield(q,'psi_azi')
    p.psi_azi = q.psi_azi;           % half-power beamwidth in azimuth-plane, examples: rad2deg(2*sqrt(pi)) 27.7 60 120
else
    p.psi_azi = deg2rad(rad2deg(2*sqrt(pi)));
end
if isfield(q,'psi_elev')
    p.psi_elev = q.psi_elev;         % half-power beamwidth in elevation-plane, examples: rad2deg(2*sqrt(pi)) 27.7 30 60
else
    p.psi_elev = deg2rad(rad2deg(2*sqrt(pi)));
end

% input value in (degree), output in (Rad)
% default rad2deg(2*sqrt(pi))

if isfield(q,'GaindBi')
    p.GaindBi = q.GaindBi;    % Gain (dBi)
else
    p.GaindBi = pow2db(4*pi/(p.psi_azi*p.psi_elev));
end
p.Gain = 10.^(p.GaindBi/10);  % Gain (scalar value)

%% Spherical and planar wave model

% Supported combinations for this version are (SA/AE): Plane/Plane, Sphere/Plane, Sphere/Sphere (NO steering vector)
if isfield(q,'WaveModelSA') 
    p.WaveModelSA = q.WaveModelSA; 
else
    p.WaveModelSA = 'Plane';    %'Sphere'/'Plane'
end
if isfield(q,'WaveModelAE') 
    p.WaveModelAE = q.WaveModelAE;
else
    p.WaveModelAE = 'Plane';  %'Sphere'/'Plane'
end

%% Beam split parameters

if isfield(q,'BeamSplitEffect')
    p.BeamSplitEffect = q.BeamSplitEffect;    
else
    p.BeamSplitEffect = 'Off'; % 'Off' 'On'
end

%% Compute Bc and TIV delay domain response parameters

if isfield(q,'BCApproxMode')
    p.BCApproxMode = q.BCApproxMode;
else
    p.BCApproxMode = 'Most Popular';
    % Options are: /'Most Popular'/'Pessimistic'/'Optimistic'/
end
switch p.channelType
    case 'LoS'
        if isfield(q,'threshold_dB')
            p.threshold_dB = q.threshold_dB;
        else
            p.threshold_dB = -40;
        end
    case 'Multipath'
        
        if isfield(q,'threshold_dB')
            p.threshold_dB = q.threshold_dB;
        else
            p.threshold_dB = -40;
        end
    case 'Multipath+LoS'
        if isfield(q,'threshold_dB')
            p.threshold_dB = q.threshold_dB;
        else
            p.threshold_dB = -40;
        end
        if isfield(q,'gap_los_nlos')
            p.gap_los_nlos = q.gap_los_nlos;
        else
            p.gap_los_nlos = -40;
        end
        
        % Ref [1]: H. Yuan, N. Yang, K. Yang, C. Han and J. An, "Hybrid Beamforming for Terahertz Multi-Carrier Systems Over Frequency Selective Fading,"
        % in IEEE Transactions on Communications, vol. 68, no. 10, pp. 6186-6199, Oct. 2020.
        % Sec. II.B. Frequency Selective THz Channel Model
        % The gap between the LOS and NLOS path gains in THz channels (e.g., more than 15 dB on average) is more significant than that in mmWave channels.
    otherwise
        error('This channel type isn''t implemented !!');
end

%% Misalignment parameters

if isfield(q,'Misalignment')
    p.Misalignment = q.Misalignment;
else
    p.Misalignment = 'Off';  % 'Off'
end

% Input parameters for misalignment
% This misalignment model will be updated in the second version of TeraMIMO

p.areaRx = (p.Mr*p.DeltaMr)*(p.Nr*p.DeltaNr);
p.radius_Rx = sqrt(p.areaRx/pi);

p.theta_M = 2*asind(1/pi*(2.782/p.Mt));
p.theta_N = 2*asind(1/pi*(2.782/p.Nt));
p.areaTx = 2*sind(p.theta_M/2)*p.d_tx_rx*2*sind(p.theta_N/2)*p.d_tx_rx;
p.radius_Tx = sqrt(p.areaTx/pi);

if isfield(q,'sigma_s')
    p.sigma_s = q.sigma_s;
else
    p.sigma_s =  0.1;   % variance of Misalignment
end
p.realizations = 1;     % 1 realization for 1 channel generation (TIV)

if(strcmp(p.Misalignment, 'On'))
    p.h_ma = get_Misalignment_Fading(p.radius_Rx, p.radius_Tx, p.sigma_s, p.realizations);
    % misalignment coefficient
else
    p.h_ma = 1; % no misalignment coefficient, default is 1.
end

end