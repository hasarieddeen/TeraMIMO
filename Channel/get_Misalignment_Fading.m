function MisCoef = get_Misalignment_Fading(a, w_d, sigma_s, N)
% =========================================================================
% -- Function to generate misalignment coefficients for the THz channel;
%    this function supports both time-invariant and time-variant channels
% =========================================================================

% -- Function MisCoef:  = get_Misalignment_Fading(a, w_d, sigma_s, N)

% -- Input Arguments:
%       a: Radius of the RX effective area
%       w_d: Maximum radius of the beam @ distance d
%       sigma_s: Standard deviation of the pointing error displacement at the RX
%       N: number of realizations

% -- Output Arguments:
%       MisCoef: Misalignment coefficient

% -- Example:
%       hp = get_Misalignment_Fading(0.01, 0.05, 0.01, 5);

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
%       Ref [1]: A. A. Boulogeorgos, E. N. Papasotiriou and A. Alexiou, "Analytical Performance Assessment of THz Wireless Systems,"
%                 in IEEE Access, vol. 7, pp. 11436-11453, 2019, doi: 10.1109/ACCESS.2019.2892198
%       Ref [2]: A. A. Farid and S. Hranilovic, "Outage Capacity Optimization for Free-Space Optical Links With Pointing Errors,"
%                 in Journal of Lightwave Technology, vol. 25, no. 7, pp. 1702-1710, July 2007, doi: 10.1109/JLT.2007.899174.

% =========================================================================
%% Initialize Parameters

% u: Intermediate variable
u = sqrt(pi)*a/(sqrt(2)*w_d);
% A_0: Fraction of the collected power at r = 0
A_0 = erf(u)^2;
% w_eq2: w_eq^2, related to w_d^2
w_eq2 = w_d^2*sqrt(pi)*erf(u)/(2*u*exp(-u^2));
% gamma: Ratio between the equivalent beam width radius and the std of r at RX
gamma = sqrt(w_eq2)/2/sigma_s;

if u > 1
    disp('Beamwidth too Small!');
    hp = zeros(N,1);
elseif gamma > 5        % set a limit for gamma
    hp = ones(N,1)*A_0;
else
    % resolution for generating discrete cdf
    resolution = A_0/1000;
    hp_sample = resolution:resolution:A_0;
    f_hp = gamma^2/A_0^(gamma^2)*hp_sample.^(gamma^2-1); % vector
    % generate CDF
    xCDF=hp_sample;
    yCDF=cumsum(f_hp)/sum(f_hp);
    x1=rand(N,1);
    x1(x1<yCDF(1)) = yCDF(1); % interpolation resolution
    hp = interp1(yCDF,xCDF,x1);
end

MisCoef = hp/A_0;
end