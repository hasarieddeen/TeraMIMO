function AlphaNLoS = get_NLoS(p, indx_subc, AlphaLoSAmp, ToAClusters, ToARays, ToAMPtot)
% =========================================================================
% -- Function to compute the NLoS total complex gain using frequency domain implementation
% =========================================================================

% -- Function: AlphaNLoS = get_NLoS(p, indx_subc, AlphaLoSAmp, ToAClusters, ToARays, ToAMPtot)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       indx_subc: Index of the subcarrier
%       AlphaLoSAmp: Magnitude of path loss for the LoS component
%       ToAClusters: Time-of-arrival (ToA) of the clusters
%       ToARays: Time-of-arrival (ToA) of the rays inside a cluster
%       ToAMPtot: Time-of-arrival (ToA) of each rays inside each cluster with a reference to the ToA of first ray in the first cluster

% -- Output Arguments:
%       AlphaNLoS: NLoS complex gain, a 3D-Array of size (p.nClustersVec(indx_subc), p.nRaysVec(indx_subc), num_freqs_per_subcarr)

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

num_freqs_per_subcarr = p.nFreq(1);

PathGain_Spread_Absorp = repmat(AlphaLoSAmp,p.nClustersVec(indx_subc),p.nRaysVec(indx_subc),1);
PathGain_DoubleDecay = repmat(exp(-ToAClusters/2/p.clusterDecayFactorVec(indx_subc)),...
    1,p.nRaysVec(indx_subc),1).*exp(-ToARays/2/p.rayDecayFactorVec(indx_subc));

AlphaNLoS_Phase_singlefreq = truncate_Angle(unifrnd(0,2*pi,p.nClustersVec(indx_subc),...
    p.nRaysVec(indx_subc),1),0,2*pi,'close','open');
AlphaNLoS_Phase = repmat(AlphaNLoS_Phase_singlefreq,1,1,num_freqs_per_subcarr);

freq_new = zeros(1,1,num_freqs_per_subcarr);
freq_new(1,1,:) = p.freq(:,indx_subc);
NLOS_PH_Delay = exp(-1j.*2.*pi.*repmat(freq_new,p.nClustersVec(indx_subc), p.nRaysVec(indx_subc),1).*ToAMPtot);

AlphaNLoS = PathGain_Spread_Absorp.*PathGain_DoubleDecay.*exp(1j*AlphaNLoS_Phase).*NLOS_PH_Delay;
end