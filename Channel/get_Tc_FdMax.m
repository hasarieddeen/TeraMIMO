function [Tc, fd_max] = get_Tc_FdMax(fc,V,c,ApproxMode)
% =========================================================================
% -- Function to compute coherence time and maximum Doppler shift of the time-variant delay domain channel response
% =========================================================================

% -- Function: [Tc, fd_max] = get_Tc_FdMax(fc,V,c,ApproxMode)

% -- Input Arguments:
%       fc: Center frequency of Tx signal at each subcarrier (Hz)
%       V: Velocity of the BS/UE (km/hr)
%       c: Speed of light in vacuum
%       ApproxMode: Options are: 'Optimistic'/'Pessimistic'/'Geometric Mean'

% -- Output Arguments:
%       Tc: Coherence time, scalar value
%       fd_max: Maximum Doppler shift, scalar value

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
%       Ref [1]: B. Sklar, "Rayleigh fading channels in mobile digital communication systems. I. Characterization,"
%        IEEE Communications magazine. 1997 Jul;35(7):90-100.

% =========================================================================

fd_max=V/c*fc*(1000/3600);% Maximum Doppler shift

switch ApproxMode
    
    case 'Optimistic'
        
        Tc=1/fd_max;
    case 'Pessimistic'
        
        Tc=9/(16*pi*fd_max);
    case 'Geometric Mean'
        
        Tc=sqrt(9/(16*pi))/fd_max;
    otherwise
        
        error('Your choice for model isn''t true. Options are: "Geometric Mean", "Pessimistic", "Optimistic"');
end

end