function K_Gas = Gas_Abs_Coef(p, data, indx_sb)
% =========================================================================
% -- Function to compute exact molecular absorption coefficients in the THz band for a specific gas
%     using the "HITRAN" database, valid in the frequency range: [0.1-10] THz. 
% =========================================================================

% -- Function: K_Gas = Gas_Abs_Coef(p, data, indx_sb)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters
%       data: A struct that contains the information to compute molecular absorption coefficient for a specific gas
%       indx_sb: Index of a subcarrier

% -- Output Arguments:
%       K_Gas: Molecular absorption coefficient for a specific gas, a vector of size(p.Nsubc, 1),
%              (number of subbands @ each subcarrier, 1)

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

% -- References: 
%       Ref [1]: J. M. Jornet and I. F. Akyildiz, "Channel Modeling and Capacity Analysis for Electromagnetic Wireless Nanonetworks in the Terahertz Band,"
%                in IEEE Transactions on Wireless Communications, vol. 10, no. 10, pp. 3211-3221, October 2011, doi: 10.1109/TWC.2011.081011.100545.

% =========================================================================

R = p.R/101325;                                     % Gas constant (m^3 atm/K/mol)
f = repmat(p.freq(:,indx_sb),1,length(data.fc0));   % Frequency

% Total # of molecules (Eq.4 Ref. [1])
% The abundance is not accounted here
% It is accounted in line intensity data.S (Eq.4 Ref. [1])
Q = (p.p/(R*p.T))*data.Q*p.Na;

% Position of resonant f = zero_pos + linear_p shift (Eq.6 Ref. [1])
fc = repmat((data.fc0+data.delta*(p.p/p.p0))',p.nFreq(1),1);

% Lorentz half-width (Eq.7 Ref. [1])
% The first alpha is computed using q^{i,j},
% and the second one assumes that the abundance is already multiplied with
alpha = repmat((((1-data.q).*(data.alphaAir)+data.q.*(data.alphaGas)).*((p.p/p.p0)*((p.T0/p.T).^data.gamma)))',p.nFreq(1),1);

% Van Vleck-Weisskopf asymmetric line shape (NOTE THAT WE REMOVED THE TERM 100*p.c) (Eq.8 Ref. [1])
F = (1/pi) * alpha.*f./fc.* (1./((f-fc).^2+(alpha).^2) + 1./((f+fc).^2+(alpha).^2));

% Line shape (Eq.9 Ref. [1])
G = f./fc.* tanh((p.h*p.c*f)/(2*p.Kb*p.T)) ./ tanh((p.h*p.c*fc)/(2*p.Kb*p.T)) .*F;

% Absorption cross section = line intensity * line shape (Eq.5 Ref. [1])
sigma = repmat((data.S)',p.nFreq(1),1) .* G; 
% Note that in this version of the code, the Line Intensity is computed only @ T0 = 296 K

% Absorption coefficient (Eq.3 Ref. [1])
K_Gas = sum(p.p/p.p0*p.Tstp/p.T*Q*sigma, 2);

end