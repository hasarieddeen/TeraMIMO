function GasData = load_Gas_Data(moleculeName, moleculeRatio, c)
% =========================================================================
% -- Function to load data of "HITRAN" database from the attached .csv file for a specific gas
%     For more information about the "data" struture, please see "header.txt" file in "Data" subfolder
% =========================================================================

% -- Function: GasData = load_Gas_Data(moleculeName, moleculeRatio, c)

% -- Input Arguments:
%       moleculeName: Gas name
%       moleculeRatio: Ratio of a gas in the transmission medium
%       c: Speed of light in vaccum

% -- Output Arguments:
%       GasData: A struct that contains the information to compute the molecular absorption coefficient for a specific gas

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

% Read .csv data file
M = csvread(sprintf('Data/%s.csv',moleculeName));

% Molecule percentage
GasData.Q = moleculeRatio;         %[ ]

% Zero-pressure position of the resonance frequency
GasData.fc0 = c*100*M(:,2);        %[Hz]

% Absorption strength using line intensity
GasData.S = c/100*M(:,3);          %[Hz.m^2/molecule]

% Linear pressure shift
GasData.delta = c*100*M(:,4);      %[Hz]

% Temperature coefficient
GasData.gamma = M(:,5);            %[ ]

% Air half-widths
GasData.alphaAir = c*100*M(:,6);   %[Hz]

% Self-broadened half-widths
GasData.alphaGas = c*100*M(:,7);   %[Hz]

% Mixing ratio for the isotopologue i of gas g
% gas ratio * isotopologue abundance
GasData.q = moleculeRatio*M(:,8);  %[ ]

end