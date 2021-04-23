function K_abs = Abs_Coef_Approx2(p)
% =========================================================================
% -- Function to compute approximate molecular absorption coefficients in the sub-THz band
%     using an approximation that is valid in the frequency range: [100-450] GHz
% =========================================================================

% -- Function: K_abs = Abs_Coef_Approx2(p)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters

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

% -- References: 
%       Ref [1]: J. Kokkoniemi, J. Lehtomaki and M. Juntti, "A line-of-sight channel model for the 100-450 gigahertz frequency band,"
%                 arXiv preprint arXiv:2002.04918. Feb, 2020.

% =========================================================================

%% Constants
q1 = 6.1121;q2 = 1.0007;q3 = 3.46e-6;
q4 = 17.502;q5 = 273.15;q6 = 32.18;
a1 = 5.159e-5;a2 = -6.65e-5;a3 = 0.0159;
b1 = -2.09e-4;b2 = 0.05;
c1 = 0.1925;c2 = 0.135;c3 = 0.0318;
d1 = 0.4241;d2 = 0.0998;
e1 = 0.2251;e2 = 0.1314;e3 = 0.0297;
f1 = 0.4127;f2 = 0.0932;
g1 = 2.053;g2 = 0.1717;g3 = 0.0306;
h1 = 0.5394;h2 = 0.0961;
i1 = 0.177;i2 = 0.0832;i3 = 0.0213;
j1 = 0.2615;j2 = 0.0668;
k1 = 2.146;k2 = 0.1206;k3 = 0.0277;
l1 = 0.3789;l2 = 0.0871;
p1 = 3.96;p2 = 6.11;p3 = 10.84;
p4 = 12.68;p5 = 14.65;p6 = 14.94;
ga = 0.0157;gb = 2e-4;gc = 0.915e-112;gd = 9.42;

pressurehPa = (p.p)*1013.25;     % Pressure (hPa)

% The saturated water vapor partial pressure at temperature (T)
% (Buck's equation:"Improved magnus form approximation of saturation vapor pressure") (Eq. 15 Ref. [1])
Pw = q1*(q2+q3*pressurehPa) * exp(q4 * (p.T-q5) / (p.T-q6) );

% Volume mixing ratio of the water vapor (Eq.14 Ref. [1])
mu = (p.phi/100)*(Pw/pressurehPa);

A = a1*(1-mu)*(a2*(1-mu)+a3);
B = (b1*(1-mu)+b2)^2;
C = c1*mu*(c2*mu+c3);
D = (d1*mu+d2)^2;
E = e1*mu*(e2*mu+e3);
F = (f1*mu+f2)^2;
G = g1*mu*(g2*mu+g3);
H = (h1*mu+h2)^2;
I = i1*mu*(i2*mu+i3);
J = (j1*mu+j2)^2;
K = k1*mu*(k2*mu+k3);
L = (l1*mu+l2)^2;

y1 = A./(B + ((p.freq)/p.c/100 - p1 ).^2);
y2 = C./(D + ((p.freq)/p.c/100 - p2 ).^2);
y3 = E./(F + ((p.freq)/p.c/100 - p3 ).^2);
y4 = G./(H + ((p.freq)/p.c/100 - p4 ).^2);
y5 = I./(J + ((p.freq)/p.c/100 - p5 ).^2);
y6 = K./(L + ((p.freq)/p.c/100 - p6 ).^2);
g = mu/ga*(gb+gc*(p.freq).^gd);

K_abs = (y1+y2+y3+y4+y5+y6)+g;

end