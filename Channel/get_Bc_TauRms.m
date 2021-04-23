function [Bc, Tau_rms] = get_Bc_TauRms(Pscal, Tau ,ApproxMode)
% =========================================================================
% -- Function to compute coherence bandwidth and root mean square of delay
%     spread of the delay domain channel response
% =========================================================================

% -- Function: [Bc, Tau_rms] = get_Bc_TauRms(Pscal, Tau ,ApproxMode)

% -- Input Arguments:
%       Pscal: Vector containing the power of each tap in the MP channel response
%       Tau: Vector containing the associated delay of each tap in the MP channel response
%       ApproxMode: Options are: 'Optimistic'/'Pessimistic'/'Most Popular'

% -- Output Arguments:
%       Bc: Coherence bandwidth, scalar value
%       Tau_rms: Root mean square of delay spread, scalar value

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

mean_excess_delay=Pscal*Tau.'/sum(Pscal);
mean_squared_excess_delay=Pscal*(Tau.^2).'/sum(Pscal);
Tau_rms=sqrt(mean_squared_excess_delay-mean_excess_delay.^2);

switch ApproxMode
    case 'Optimistic'
        
        Bc=1/(2*pi*Tau_rms);
    case 'Pessimistic'
        
        Bc=1/(50*Tau_rms);
    case 'Most Popular'
        
        Bc=1/(5*Tau_rms);
    otherwise
        
        error('Your choice for model isn''t true. Options are: "Most Popular", "Pessimistic", "Optimistic"');
end

end

