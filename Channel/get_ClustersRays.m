function [ToAClusters, ToARays, ToAMPtot, AoD_NLOS, AoA_NLOS] = get_ClustersRays(p, indx_subc)
% =========================================================================
% -- Function to generate the required multipath (MP) parameters following the modified S-V model to fit the THz band peculiarities,
%    where MP components are in the form of multiple clusters, and multiple rays inside each cluster. 
% =========================================================================

% -- Function: [ToAClusters, ToARays, ToAMPtot, AoD_NLOS, AoA_NLOS] = get_ClustersRays(p, indx_subc)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       indx_subc: Index of the subcarrier

% -- Output Arguments:
%       ToAClusters: Time-of-arrival (ToA) of the clusters following "Exponential Distribution", a 3D-Array of size (p.nClustersVec(indx_subc), 1, num_freqs_per_subcarr),
%        i.e. size (Number of clusters in a specific subcarrier (based on indx_subc), 1, Number of sub-bands inside each subcarrier)
%       ToARays: Time-of-arrival (ToA) of the rays inside a Cluster following "Exponential Distribution", a 3D-Array of size (p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr),
%        i.e. size (Number of clusters in a specific subcarrier (based on indx_subc), Number of rays in a Specific subcarrier (based on indx_subc), Number of sub-bands inside each subcarrier)
%       ToAMPtot: Time-of-arrival (ToA) of each rays inside each cluster with reference to the ToA of first ray in the first cluster, 
%        a 3D-Array of size (p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr)
%       AoD_NLOS: Angle-of-departure, a 4D-Array of size (2, p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr),
%        where the first dimension refers to [Azimuth; Elevation] angles
%       AoA_NLOS: Angle-of-arrival, a 4D-Array of size (2, p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr),
%        where the first dimension refers to [Azimuth; Elevation] angles 

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

% -- References: 
%       Ref [1]: C. Lin and G. Y. Li, "Indoor Terahertz Communications: How Many Antenna Arrays Are Needed?,"
%                 in IEEE Transactions on Wireless Communications, vol. 14, no. 6, pp. 3097-3107, June 2015, doi: 10.1109/TWC.2015.2401560.

% =========================================================================
%% Initialize Outputs

num_freqs_per_subcarr = p.nFreq(1);
ToAClusters = zeros(p.nClustersVec(indx_subc), 1, num_freqs_per_subcarr);
ToARays =  zeros(p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr);

AoD_NLOS = zeros(2, p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr); % 2 in first dimension referes to [Azimuth; Elevation]
AoA_NLOS = zeros(2, p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr); % 2 in first dimension referes to [Azimuth; Elevation]

%% Generate ToA of clusters and rays

if strcmp(p.ToAType,'Exponential')
    % ToA Clusters Generation
    ToACluster1freq = cumsum( ...
        random('Exponential',1/p.clusterArrivalRateVec(indx_subc),p.nClustersVec(indx_subc),1),1);
    if ~strcmp(p.channelType,'Multipath') 
        if (p.Fc + p.BW/2) < 0.55e12
            while (ToACluster1freq(1) <= p.tau_LoS)
                ToACluster1freq = cumsum( ...
                    random('Exponential',1/p.clusterArrivalRateVec(indx_subc),p.nClustersVec(indx_subc),1),1);
            end
        end        
    end
    tmp_ToACluster = repmat(ToACluster1freq,1,num_freqs_per_subcarr);
    
    % Rays ToA generation
    if p.nRaysVec(indx_subc) ~= 1
        ToARays1freq = zeros(p.nClustersVec(indx_subc),p.nRaysVec(indx_subc),1);
        ToARays1freq(:,2:end,:) = cumsum( ...
            random('Exponential',1/p.rayArrivalRateVec(indx_subc),p.nClustersVec(indx_subc),p.nRaysVec(indx_subc)-1,1),2);
        ToARays = repmat(ToARays1freq,1,1,num_freqs_per_subcarr);
    end
    
    % Clusters/rays ToA generation
    ToAClusters(:,1,:) = tmp_ToACluster;
    ToAMPtot = repmat(ToAClusters,1,p.nRaysVec(indx_subc),1) + ToARays;

elseif strcmp(p.ToAType,'Paraboloid') 
    error('Still Not Implementad in this version');
else
    error('Supported options of ToA are: (Exponential/Paraboloid)'); 
end

%% Generate AoD/AoA of clusters and rays
% Angle means for each cluster

% AoD mean in azimuth/elevation
AoDAzimuthMean = unifrnd(-pi,pi,p.nClustersVec(indx_subc),1);
AoDElevationMean = unifrnd(-pi/2,pi/2,p.nClustersVec(indx_subc),1);

% AoA mean in azimuth/elevation
AoAAzimuthMean = unifrnd(-pi,pi,p.nClustersVec(indx_subc),1);
AoAElevationMean = unifrnd(-pi/2,pi/2,p.nClustersVec(indx_subc),1);

% Angles for each ray within a cluster
if (p.Fc + p.BW/2) < 0.55e12
    % AoD ray in azimuth
    tmp_AoDA1 = reshape(random(p.gmAoDAzimuth{indx_subc},p.nClustersVec(indx_subc)*p.nRaysVec(indx_subc)), ...
        p.nClustersVec(indx_subc),p.nRaysVec(indx_subc));
    tmp_AoDA2 = zeros(p.nClustersVec(indx_subc),p.nRaysVec(indx_subc),1);
    tmp_AoDA2(:,:,1) = tmp_AoDA1;
    tmp_AoDA = repmat(tmp_AoDA2,1,1,num_freqs_per_subcarr);
    AoDAzimuth = repmat(AoDAzimuthMean,1,p.nRaysVec(indx_subc),num_freqs_per_subcarr)+ tmp_AoDA;
    
    % AoD ray in elevation
    tmp_AoDE1 = reshape(random(p.gmAoDElevation{indx_subc},p.nClustersVec(indx_subc)*p.nRaysVec(indx_subc)), ...
        p.nClustersVec(indx_subc),p.nRaysVec(indx_subc));
    tmp_AoDE2 = zeros(p.nClustersVec(indx_subc),p.nRaysVec(indx_subc),1);
    tmp_AoDE2(:,:,1) = tmp_AoDE1;
    tmp_AoDE = repmat(tmp_AoDE2,1,1,num_freqs_per_subcarr);
    AoDElevation = repmat(AoDElevationMean,1,p.nRaysVec(indx_subc),num_freqs_per_subcarr)+ tmp_AoDE;
    
    % AoA ray in azimuth
    tmp_AoAA1 = reshape(random(p.gmAoAAzimuth{indx_subc},p.nClustersVec(indx_subc)*p.nRaysVec(indx_subc)), ...
        p.nClustersVec(indx_subc),p.nRaysVec(indx_subc));
    tmp_AoAA2 = zeros(p.nClustersVec(indx_subc),p.nRaysVec(indx_subc),1);
    tmp_AoAA2(:,:,1) = tmp_AoAA1;
    tmp_AoAA = repmat(tmp_AoAA2,1,1,num_freqs_per_subcarr);
    AoAAzimuth = repmat(AoAAzimuthMean,1,p.nRaysVec(indx_subc),num_freqs_per_subcarr) + tmp_AoAA;
    
    % AoA ray in elevation
    tmp_AoAE1 = reshape(random(p.gmAoAElevation{indx_subc},p.nClustersVec(indx_subc)*p.nRaysVec(indx_subc)), ...
        p.nClustersVec(indx_subc),p.nRaysVec(indx_subc));
    tmp_AoAE2 = zeros(p.nClustersVec(indx_subc),p.nRaysVec(indx_subc),1);
    tmp_AoAE2(:,:,1) = tmp_AoAE1;
    tmp_AoAE = repmat(tmp_AoAE2,1,1,num_freqs_per_subcarr);
    AoAElevation = repmat(AoAElevationMean,1,p.nRaysVec(indx_subc),num_freqs_per_subcarr) + tmp_AoAE;
    
    % Final azimuth/elevation AoD and truncate angles
    AoD_NLOS(1,:,:,:) = truncate_Angle(AoDAzimuth,-pi,pi,'open','close');
    AoD_NLOS(2,:,:,:) = truncate_Angle(AoDElevation,-pi/2,pi/2,'close','close');
    
    % Final azimuth/elevation AoA and truncate angles
    AoA_NLOS(1,:,:,:) = truncate_Angle(AoAAzimuth,-pi,pi,'open','close');
    AoA_NLOS(2,:,:,:) = truncate_Angle(AoAElevation,-pi/2,pi/2,'close','close');
end

end