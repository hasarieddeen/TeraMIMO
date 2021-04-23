function G = get_SectorGain(p, angle_in, Flag)
% =========================================================================
% -- Function to compute antenna gains based on the ideal sector model (ISM)
% =========================================================================

% -- Function: G = get_SectorGain(p, angle_in, Flag)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       angle_in: Input angles, maybe AoA/AoD for LoS/NLoS
%       Flag: Options are: 'LoS'/'NLoS'

% -- Output Arguments:
%       G: Amplitude antenna gain, size based on "Flag"
%        1) LoS, (SA/AE): (SWM/PWM, PWM/PWM) scalar size(1,1), (SWM/SWM) matrix, size(p.Qar, p.Qat) for Rx/Tx.
%        2) NLoS (PWM only for both SA/AE): 3D-Array, size(p.nClusters_per_subcarr, p.nRays_per_subcarr, num_freqs_per_subcarr).
 
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

if strcmp(Flag,'LoS')
    % Here, we support spherical and planar AE wave models
    G = zeros(size(angle_in,2),size(angle_in,3));
    
    for indx_aer = 1:size(angle_in,2)
        for indx_aet = 1:size(angle_in,3)
            
            if ((-p.psi_azi/2) <= angle_in(1,indx_aer,indx_aet) && ...
                    angle_in(1,indx_aer,indx_aet) <= (p.psi_azi/2)) && ...
                    ((pi/2 - p.psi_elev/2) <= angle_in(2,indx_aer,indx_aet) && ...
                    angle_in(2,indx_aer,indx_aet) <= (pi/2 + p.psi_elev/2))
                G(indx_aer,indx_aet) = sqrt(p.Gain);
            else
                G(indx_aer,indx_aet) = 0;
            end
            
        end
    end
elseif strcmp(Flag,'NLoS')
    G = zeros(size(angle_in,2),size(angle_in,3),size(angle_in,4));

     for indx_clusters = 1:size(angle_in,2)
         for indx_rays = 1:size(angle_in,3)
            for indx_freq = 1:size(angle_in,4)
                
                 if ((-p.psi_azi/2) <= angle_in(1,indx_clusters,indx_rays,indx_freq)&& ...
                         angle_in(1,indx_clusters,indx_rays,indx_freq) <= (p.psi_azi/2)) && ...
                         ((pi/2 - p.psi_elev/2) <= (pi/2 + angle_in(2,indx_clusters,indx_rays,indx_freq)) && ...
                         (pi/2 + angle_in(2,indx_clusters,indx_rays,indx_freq)) <= (pi/2 + p.psi_elev/2))
                    G(indx_clusters,indx_rays,indx_freq) = sqrt(p.Gain);
                 else
                     G(indx_clusters,indx_rays,indx_freq) = 0;
                 end
             end
         end
     end
else
    error('Options are: "LoS" or "NLoS"')
end
end