function K_abs = compute_Abs_Coef(p)
% =========================================================================
% -- Function to select the calculation method of molecular absorption coefficient
% =========================================================================

% -- Function: K_abs = compute_Abs_Coef(p)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%          "p.absorptionType" defines the method of computation
%          Options are: Hitran, Approx1, and Approx2
%          1) Hitran: Exact absorption coefficient, valid in the frequency band: [0.1-10] THz
%          2) Approx1: First approximation of absorption coefficient, valid in the frequency band: [275-400] GHz
%          3) Approx2: Second approximation of absorption coefficient, valid in the frequency band: [100-450] GHz

% -- Output Arguments:
%       K_abs: Molecular absorption coefficient, a matrix of size(p.Nsubc, p.Nsubb),
%              (number of subbands @ each subcarrier, number of subcarriers)

%=================================================

% -- (c) 2021 Simon Tarboush, Hadi Sarieddeen, Hui Chen, 
%             Mohamed Habib Loukil, Hakim Jemaa, 
%             Mohamed-Slim Alouini, Tareq Y. Al-Naffouri

% -- e-mail: simon.w.tarboush@gmail.com; hadi.sarieddeen@kaust.edu.sa; hui.chen@kaust.edu.sa;
%            mohamedhabib.loukil@kaust.edu.sa; hakim.jemaa@kaust.edu.sa;
%            slim.alouini@kaust.edu.sa; tareq.alnaffouri@kaust.edu.sa

% =========================================================================

% S. Tarboush, H. Sarieddeen, H. Chen, M.-H. Loukil, H. Jemaa, M.-S. Alouini, and T. Y. Al-Naffouri, 
%  "TeraMIMO: A channel simulator for wideband ultra-massive MIMO terahertz communications,"
%  arXivpreprint arXiv:2104.11054, 2021.

% =========================================================================

%% Select Computation Method

switch p.absorptionType
    case 'Hitran'
        % 1) Hitran: Exact absorption coefficient is valid in the frequency band: [0.1-10] THz
        K_abs = Abs_Coef_Hitran(p);
    case 'Approx1'
        % 2) Approx1: Approximation of absorption coefficient is valid in the frequency band: [275-400] GHz
        K_abs = Abs_Coef_Approx1(p);
    case 'Approx2'
        % 3) Approx2: Approximation of absorption coefficient is valid in the frequency band: [100-450] GHz
        K_abs = Abs_Coef_Approx2(p);
    otherwise
         error('This method for computing absorption coefficient isn''t implemented, options are: Hitran, Approx1, Approx2');
end

end