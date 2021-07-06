function [AoD_SV, AoA_SV, AoD_BF, AoA_BF] = get_AoD_AoA_LoS(p, mr, nr, mt, nt)
% =========================================================================
% -- Function to compute local AoD/AoA, where this function supports multiple options based on the model, Spherical Wave Model (SWM) or Plane Wave Model (PWM).
%    these computations are on the level of SA/AE. Options are (PWM/PWM, SWM/PWM, SWM/SWM)
% =========================================================================

% -- Function: [AoD_SV, AoA_SV, AoD_BF, AoA_BF] = get_AoD_AoA_LoS(p, mr, nr, mt, nt)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       mr: Index of the number of the receiver SAs (rows)
%       nr: Index of the number of the receiver SAs (columns)
%       mt: Index of the number of the transmitter SAs (rows)
%       nt: Index of the number of the transmitter SAs (columns)

% -- Output Arguments:
%       AoD_SV: Angle-of-departure for beamsteering vector (SV); size is based on the used model for the SA/AE
%        1) PWM/PWM: Vector of size(2, 1), (Azimuth; Elevation)
%        2) SWM/PWM: Vector of size(2, 1), (Azimuth; Elevation)
%        3) SWM/SWM: 3D-Array of size(2, p.Qar, p.Qat), (Azimuth; Elevation),
%            where p.Qar/p.Qat is the total number of receiver/transmitter AEs inside each SA
%       AoA_SV: Angle-of-arrival for beamsteering vector (SV); size is based on the used model for the SA/AE
%        1) PWM/PWM: Vector of size(2, 1), (Azimuth; Elevation)
%        2) SWM/PWM: Vector of size(2, 1), (Azimuth; Elevation)
%        3) SWM/SWM: 3D-Array of size(2, p.Qar, p.Qat), (Azimuth; Elevation),
%            where p.Qar/p.Qat is the total number of receiver/transmitter AEs inside each SA
%       AoD_BF: Angle-of-departure for analog beamforming (BF, SA level); size is based on the used model for the SA/AE
%        1) PWM/PWM: Vector of size(2, 1), (Azimuth; Elevation)
%        2) SWM/PWM: Vector of size(2, 1), (Azimuth; Elevation)
%        3) SWM/SWM: 3D-Array of size(2, p.Qar, p.Qat), (Azimuth; Elevation),
%            where p.Qar/p.Qat is the total number of receiver/transmitter AEs inside each SA
%       AoA_BF: Angle-of-arrival for analog beamforming (BF, SA level); size is based on the used model for the SA/AE
%        1) PWM/PWM: Vector of size(2, 1), (Azimuth; Elevation)
%        2) SWM/PWM: Vector of size(2, 1), (Azimuth; Elevation)
%        3) SWM/SWM: 3D-Array of size(2, p.Qar, p.Qat), (Azimuth; Elevation),
%            where p.Qar/p.Qat is the total number of receiver/transmitter AEs inside each SA

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

if((strcmp(p.WaveModelSA,'Sphere')) && strcmp(p.WaveModelAE,'Plane'))
    % SWM for SA, compute different AoA/AoD for each SA pair
    positionLocalSARx = [0; (nr-1-(p.Nr-1)/2)*p.DeltaNr; (mr-1-(p.Mr-1)/2)*p.DeltaMr];
    positionGlobalSARx = p.RotmRx*positionLocalSARx + p.positionRx;
    
    positionLocalSATx = [0; (nt-1-(p.Nt-1)/2)*p.DeltaNt; (mt-1-(p.Mt-1)/2)*p.DeltaMt];
    positionGlobalSATx = p.RotmTx*positionLocalSATx + p.positionTx;
    
    unitDirTx = (positionGlobalSARx-positionGlobalSATx)/norm(positionGlobalSARx-positionGlobalSATx);
    unitDirLocalTx = p.RotmTx'*unitDirTx;
    
    unitDirRx = (positionGlobalSATx-positionGlobalSARx)/norm(positionGlobalSATx-positionGlobalSARx);
    unitDirLocalRx = p.RotmRx'*unitDirRx;
    
    % SV (physical) AoD/AoA
    AoD_SV = [atan2(unitDirLocalTx(2), unitDirLocalTx(1)); acos(unitDirLocalTx(3))];
    AoA_SV = [atan2(unitDirLocalRx(2), unitDirLocalRx(1)); acos(unitDirLocalRx(3))];
    
    % Analog BF AoD/AoA
    AoD_BF = p.AoD_BF;
    AoA_BF = p.AoA_BF;
    
elseif ((strcmp(p.WaveModelSA,'Plane')) && strcmp(p.WaveModelAE,'Plane'))
    % PWM for SA, compute same DoA/DoD for each SA pair
    % SV AoD/AoA
    AoD_SV = p.AoD;
    AoA_SV = p.AoA;
    
    % Analog BF AoD/AoA
    AoD_BF = p.AoD_BF; 
    AoA_BF = p.AoA_BF; 
    
elseif ((strcmp(p.WaveModelSA,'Sphere')) && strcmp(p.WaveModelAE,'Sphere'))

    AoD_SV = zeros(2, p.Qar, p.Qat);
    AoA_SV = zeros(2, p.Qar, p.Qat);
    
    % SV AoD/AoA
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
                    
                    unitDirTx = (positionGlobalAERx-positionGlobalAETx)/norm(positionGlobalAERx-positionGlobalAETx);
                    unitDirLocalTx = p.RotmTx'*unitDirTx;   
                    
                    unitDirRx = (positionGlobalAETx-positionGlobalAERx)/norm(positionGlobalAETx-positionGlobalAERx);
                    unitDirLocalRx = p.RotmRx'*unitDirRx;
                    
                    % SV AoD/AoA
                    AoD_SV(:,(aemr-1)*p.Nar+aenr,(aemt-1)*p.Nat+aent) = [atan2(unitDirLocalTx(2), unitDirLocalTx(1)); acos(unitDirLocalTx(3))];
                    AoA_SV(:,(aemr-1)*p.Nar+aenr,(aemt-1)*p.Nat+aent) = [atan2(unitDirLocalRx(2), unitDirLocalRx(1)); acos(unitDirLocalRx(3))];
                end
            end
        end
    end
    
    % Analog BF AoD/AoA
    AoD_BF = repmat(p.AoD_BF, 1, p.Qar, p.Qat);
    AoA_BF = repmat(p.AoA_BF, 1, p.Qar, p.Qat);
    
else
    error('This wavemodel on SA/AE level isn''t Implemented. Supported options are: SWM/SWM, SWM/PWM, PWM/PWM');
end

end