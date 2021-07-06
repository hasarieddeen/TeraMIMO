function K_abs = Abs_Coef_Approx1(p)
% =========================================================================
% -- Function to compute approximate molecular absorption coefficients in the sub-THz band
%     using an approximation that is valid in the frequency range: [275-400] GHz 
% =========================================================================

% -- Function: K_abs = Abs_Coef_Approx1(p)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters

% -- Output Arguments:
%       K_abs: Molecular absorption coefficient, a matrix of size(p.Nsub_b, p.Nsub_c),
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

% -- References:
%       Ref [1]: J. Kokkoniemi, J. Lehtomaki and M. Juntti, "Simplified molecular absorption loss model for 275–400 gigahertz frequency band,"
%                12th European Conference on Antennas and Propagation (EuCAP 2018), London, UK, 2018, pp. 1-5, doi: 10.1049/cp.2018.0446.

% =========================================================================

%% Constants
c1 = 10.835;c2 = 12.664;
p1 = 5.54e-37;p2 = -3.94e-25;
p3 = 9.06e-14;p4 = -6.36e-3;
g1 = 0.2205;g2 = 0.1303;g3 = 0.0294;g4 = 0.4093;g5 = 0.0925;
g6 = 2.014;g7 = 0.1702;g8 = 0.0303;g9 = 0.537;g10 = 0.0956;
q1 = 6.1121;q2 = 1.0007;q3 = 3.46e-6;
q4 = 17.502;q5 = 273.15;q6 = 32.18;

pressurehPa = (p.p)*1013.25;     % Pressure (hPa)

% The saturated water vapor partial pressure at temperature (T)
% (Buck's equation: "Improved magnus form approximation of saturation vapor pressure") (Eq. 15 Ref. [1])
Pw = q1*(q2+q3*pressurehPa) * exp(q4 * (p.T-q5) / (p.T-q6) );

% Volume mixing ratio of the water vapor (Eq.14 Ref. [1])
v = (p.phi/100)*(Pw/pressurehPa);

Av = g1*v*(g2*v+g3);
Bv = (g4*v+g5)^2;
Cv = g6*v*(g7*v+g8);
Dv = (g9*v+g10)^2;

y1 = Av./ (Bv + ((p.freq)/p.c/100 - c1).^2);             % (Eq.11 Ref. [1])
y2 = Cv./ (Dv + ((p.freq)/p.c/100 - c2).^2);             % (Eq.12 Ref. [1])

g = p1*((p.freq).^3) + p2*((p.freq).^2) + p3*(p.freq) + p4;     % (Eq.16 Ref. [1])

K_abs = (y1+y2)+g;                                              % (Eq.13 Ref. [1])

end