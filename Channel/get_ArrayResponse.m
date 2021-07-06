function a = get_ArrayResponse(p, indx_subc, M, N, AngleIn, Psi_type, TRx)
% =========================================================================
% -- Function to compute array response (AE level) for antenna structure using local coordinates
% =========================================================================

% -- Function: a = get_ArrayResponse(p, indx_subc, M, N, AngleIn, Psi_type, TRx)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       indx_subc: Index of the subcarrier
%       M: Number of receiver/transmitter SAs (rows) based on TRx
%       N: Number of receiver/transmitter SAs (columns) based on TRx
%       AngleIn: AoD/AoA for LoS/NLoS; in the case of NLoS, only for one ray in a cluster
%       Psi_type: Define whether to calculate a beamsteering or beamforming vector 'SV', 'BF', as follows:
%        1) @ Tx SV ---> AoD, @ Rx SV ---> AoA
%        2) @ Tx BF ---> AoD_BF, @ Rx BF ---> AoA_BF
%       TRx: Defines the direction of the link, 'T' at Tx side, 'R' at Rx side

% -- Output Arguments:
%       a: Array response, a 3D-Array of size(M, N, num_freqs_per_subcarr)

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

%% Initialize Output
num_freqs_per_subcarr = p.nFreq(1);
a = zeros(M, N, num_freqs_per_subcarr);

%% 
% p.deltaNt @Tx, p.deltaNr @Rx
% p.deltaMt @Tx, p.deltaMr @Rx
if strcmp(TRx, 'T')
    
    delta1 = p.deltaNt; 
    delta2 = p.deltaMt; 
elseif strcmp(TRx, 'R')
    
    delta1 = p.deltaNr;
    delta2 = p.deltaMr;
else
    
    error('TRx has only two options: T/R');
end

if strcmp(Psi_type,'SV')
    lambda = p.lambdak(:,indx_subc);
elseif strcmp(Psi_type,'BF')
    % Add Beam Split Effect
    % Beam Split is caused by fixed beamformer vector (frequency-independent)
    
    if strcmp(p.BeamSplitEffect,'On')
        
        lambda = repmat(p.lambdac(indx_subc), num_freqs_per_subcarr, 1);
    elseif strcmp(p.BeamSplitEffect,'Off')
        
        lambda = p.lambdak(:,indx_subc);
    else
        
        error('BeamSplitEffect has only two options: On/Off');
    end
else
    
    error('Options are only beamsteering vector (SV) or analog beamforming (BF)')
end

for m = 1:M
    for n = 1:N
        
        Psi = delta1*(n-1-(N-1)/2)*sin(squeeze(AngleIn(1,:,:,:))).*sin(squeeze(AngleIn(2,:,:,:))) + ...
            delta2*(m-1-(M-1)/2)*cos(squeeze(AngleIn(2,:,:,:)));
        
        if strcmp(Psi_type,'SV')
            
            a(m,n,:) = exp(1j*2*pi*Psi./lambda);
        elseif strcmp(Psi_type,'BF')
            
            a(m,n,:) = exp(-1j*2*pi*Psi./lambda);
        else
            
    error('Options are only beamsteering vector (SV) or analog beamforming (BF)')
        end
    end
end

% Normalization
if strcmp(Psi_type,'SV')
    
    SV_norm = sqrt(M*N);
    a = a/SV_norm;
elseif strcmp(Psi_type,'BF')
    
    BF_norm = 1;
    a = a/BF_norm;
else
    
    error('Options are only beamsteering vector (SV) or analog beamforming (BF)')
end

end
