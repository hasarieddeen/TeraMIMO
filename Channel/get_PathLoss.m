function [AlphaLoS, AlphaLoSAmp, Tau_LoS] = get_PathLoss(p, indx_subc, mr, nr, mt, nt, K_abs)
% =========================================================================
% -- Function to compute the LoS path loss. This function supports multiple options based on the used model: spherical wave model (SWM) or plane wave model (PWM)
%    These computations are on the level of SA/AE. Options are: (PWM/PWM, SWM/PWM, SWM/SWM)
% =========================================================================

% -- Function: [AlphaLoS, AlphaLoSAmp, Tau_LoS] = get_PathLoss(p, indx_subc, mr, nr, mt, nt, K_abs)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       indx_subc: Index of the subcarrier
%       mr: Index of the number of the receiver SAs (rows)
%       nr: Index of the number of the receiver SAs (columns)
%       mt: Index of the number of the transmitter SAs (rows)
%       nt: Index of the number of the transmitter SAs (columns)
%       K_abs: Molecular absorption coefficient, matrix, size(p.Nsub_b, p.Nsub_c)

% -- Output Arguments:
%       AlphaLoS: Path loss for the LoS component; the size is based on the used model for the SA/AE
%        1) PWM/PWM: 3D-Array of size(1,1, num_freqs_per_subcarr)
%        2) SWM/PWM: 3D-Array of size(1,1, num_freqs_per_subcarr)
%        3) SWM/SWM: 3D-Array of size(p.Qar,p.Qat, num_freqs_per_subcarr),
%           where p.Qar/p.Qat is the total number of receiver/transmitter AEs inside each SA
%       AlphaLoSAmp: Magnitude of path loss for LoS component; the size is based on the used model for the SA/AE
%        1) PWM/PWM: 3D-Array of size(1,1, num_freqs_per_subcarr)
%        2) SWM/PWM: 3D-Array of size(1,1, num_freqs_per_subcarr)
%        3) SWM/SWM: 3D-Array of size(p.Qar,p.Qat, num_freqs_per_subcarr),
%            where p.Qar/p.Qat is the total number of receiver/transmitter AEs inside each SA
%       Tau_LoS: Time-of-arrival (ToA) of the LoS component; the size is based on the used model for the SA/AE
%        1) PWM/PWM: Scalar value (1x1)
%        2) SWM/PWM: Scalar value (1x1)
%        3) SWM/SWM: Matrix of size(p.Qar,p.Qat)

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

if((strcmp(p.WaveModelSA,'Sphere')) && strcmp(p.WaveModelAE,'Plane'))

    positionLocalSARx = [0; (nr-1-(p.Nr-1)/2)*p.DeltaNr; (mr-1-(p.Mr-1)/2)*p.DeltaMr];
    positionGlobalSARx = p.RotmRx*positionLocalSARx + p.positionRx;

    positionLocalSATx = [0; (nt-1-(p.Nt-1)/2)*p.DeltaNt; (mt-1-(p.Mt-1)/2)*p.DeltaMt];
    positionGlobalSATx = p.RotmTx*positionLocalSATx + p.positionTx;

    d_SAAEs = norm((positionGlobalSARx-positionGlobalSATx));
    
elseif ((strcmp(p.WaveModelSA,'Plane')) && strcmp(p.WaveModelAE,'Plane'))

    d_AoA = p.DeltaNr*(nr-1-(p.Nr-1)/2)*sin(p.AoA(1))*sin(p.AoA(2)) + ...
        p.DeltaMr*(mr-1-(p.Mr-1)/2)*cos(p.AoA(2));
    d_AoD = p.DeltaNt*(nt-1-(p.Nt-1)/2)*sin(p.AoD(1))*sin(p.AoD(2)) + ...
        p.DeltaMt*(mt-1-(p.Mt-1)/2)*cos(p.AoD(2));
    
    d_SAAEs = -(d_AoA+d_AoD)+ p.d_tx_rx;

elseif ((strcmp(p.WaveModelSA,'Sphere')) && strcmp(p.WaveModelAE,'Sphere'))
    
    d_SAAEs = zeros(p.Qar,p.Qat);
    
    positionLocalSARx = [0; (nr-1-(p.Nr-1)/2)*p.DeltaNr; (mr-1-(p.Mr-1)/2)*p.DeltaMr];
    positionLocalSATx = [0; (nt-1-(p.Nt-1)/2)*p.DeltaNt; (mt-1-(p.Mt-1)/2)*p.DeltaMt];

    for aemr = 1:p.Mar
        for aenr = 1:p.Nar
            for aemt = 1:p.Mat
                for aent = 1:p.Nat
                    
                    positionLocalAERx = [0; (aenr-1-(p.Nar-1)/2)*p.deltaNr; (aemr-1-(p.Mar-1)/2)*p.deltaMr] + positionLocalSARx;
                    positionGlobalAERx = p.RotmRx*positionLocalAERx + p.positionRx;

                    positionLocalAETx = [0; (aent-1-(p.Nat-1)/2)*p.deltaNt; (aemt-1-(p.Mat-1)/2)*p.deltaMt] + positionLocalSATx;
                    positionGlobalAETx = p.RotmTx*positionLocalAETx + p.positionTx;
                    
                    d_SAAEs((aemr-1)*p.Nar+aenr,(aemt-1)*p.Nat+aent) = norm((positionGlobalAERx-positionGlobalAETx))- p.d_tx_rx;
                end
            end
        end
    end
else
    error('This wavemodel on SA/AE level isn''t implemented. Supported options are: SWM/SWM, SWM/PWM, PWM/PWM');
end

D_eff = d_SAAEs; % Distance of SAs/AEs (based on model) between Rx & Tx
Tau_LoS = D_eff./p.c; % Deterministic delay of the LoS path

if strcmp(p.WaveModelAE,'Plane')
    
    AlphaLoSAmp = zeros(1, 1, num_freqs_per_subcarr);
    AlphaLoSPhase = zeros(1, 1, num_freqs_per_subcarr);
    
    AlphaLoSAmp(1,1,:) = (p.c./(4*pi*p.freq(:,indx_subc).*D_eff)).^(p.PLE/2)...
        .*exp(-K_abs(:,indx_subc).*D_eff/2);
    AlphaLoSPhase(1,1,:) = exp(-1j*2*pi*p.freq(:,indx_subc).*Tau_LoS);
    
elseif strcmp(p.WaveModelAE,'Sphere')
    
    AlphaLoSAmp = zeros(p.Qar, p.Qat, num_freqs_per_subcarr);
    AlphaLoSPhase = zeros(p.Qar, p.Qat, num_freqs_per_subcarr);
    
    for idx_freq = 1:num_freqs_per_subcarr   
        % Absorption and Spreading Loss
        AlphaLoSAmp(:, :, idx_freq) = (p.c./(4*pi*p.freq(idx_freq,indx_subc).*D_eff)).^(p.PLE/2) ...
            .*exp(-K_abs(idx_freq,indx_subc).*D_eff/2);
        AlphaLoSPhase(:, :, idx_freq) = exp(-1j*2*pi*p.freq(idx_freq,indx_subc).*Tau_LoS);
        
    end
else
    error('This wavemodel for AE level isn''t Implemented. Supported options are: SWM/PWM');
end

% Get path loss for LOS
% AlphaLoS = SpreadingLoss.*AbsorptionLoss .* AlphaLoSPhase
% AlphaLoS = AlphaLoSAmp .* AlphaLoSPhase

AlphaLoS = AlphaLoSAmp.*AlphaLoSPhase;

end