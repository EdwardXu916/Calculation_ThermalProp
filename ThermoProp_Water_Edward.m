% Title: Calculate the other three variables based on two knowing
%        thermodynamic variables. And plot the T(l)-s Diagram of Water.
% Based on: MATLAB program from << Chemical, Biochemical, and Engineering
%                                  Thermodynamics >>
% Version: 3.0, Edward Xu, 18.4.24
% SubTitle: Calculate the other three variables based on two knowing
%           thermodynamic variables.

% Input
% T(l)
% p(m)

% Constants ------------------------------------------------------------------
R = 8.314472;                % J/(mol*K), Universial Gas Constant
M = 18 / 1000;               % kg / mol , Molar Mass
R_G = R / M;                 % J/(K*kg) , Gas Constant - Water
AAA = 32.24;                 % C_p1, heat capacity calculation parameter of Water
BBB = 0.1923E-2;             % C_p2, heat capacity calculation parameter of Water
CCC = 1.055E-5;              % C_p3, heat capacity calculation parameter of Water
DDD = -3.595E-9;             % C_p4, heat capacity calculation parameter of Water
T_c = 374 + 273.15;          % K   647.15
p_c = 22064 * 1000;          % Pa 
% h_c = 2085.9 * 1000;       % J/kg
% s_c = 4.4092 * 1000;       % J/(kg*K)
% v_c = 0.003106             % m^3/kg
% OMEGA = 0.176;
% T_boil = 243.4;            % ???
T_ref = 0 + 273.15;         % K
p_ref = 101.325 * 1000;      % Pa

% Part1: Peng-Robinson Constant Calculation,
% Reference: << Chemical, Biochemical, and Engineering Thermodynamics, 5e >>
%            E6.4-2, P221, Ch6
a_Tc = 0.457235529 * (R_G * T_c)^2 ./ p_c;       % Critical Point Restriction "a(T_c)"
KAPPA = 0.4069;                                  % Dependent on OMEGA(working substance), Temperature-independent parameter in PR-EOS
%(KAPPA = 0.37464 + (1.54226 - 0.26992 * OMEGA) * OMEGA;)
T_r = T(l) ./ T_c;                               % Reduced Temerature
ALPHASqrt = 1 + KAPPA * (1 - sqrt(T_r));
ALPHA = ALPHASqrt^2;                             % Temperature-dependent parameter in PR-EOS
a_T = a_Tc * ALPHA;                              % Temperature-dependent parameter in PR-EOS (dimensions depend on equation)
b = 0.077796074 * R_G * T_c ./ p_c;              % m^3/mol, Critical Point Restriction "b", Temperature-independent parameter in PR-EOS
DADT = - a_Tc * KAPPA * ALPHASqrt ./ sqrt(T_c*T(l));
A = a_T * p(m) ./ ((R_G * T(l))^2);              % Parameters for Cubic Form of Equation of State, dimensionless form of EOS parameter a_T
B = b * p(m) ./ (R_G * T(l));                    % Parameters for Cubic Form of Equation of State, dimensionless form of EOS parameter b
[Z_g,Z_l] = ZZroot2(A,B);

% Part2: Solve Peng-Robinson EOS to get compressibility fa_Tctor
% Reference: << Chemical, Biochemical, and Engineering Thermodynamics, 5e >>
%            E6.4-4, P222-223, Ch6
% Root = SolveCubic(Para_TcF(1),Para_TcF(2),Para_TcF(3));
% Root,
% Z(1) = max(ZZ); % Vapor Phase, most compressible;
% Z(2) = min(ZZ); % Liquid Phase, least compressible;
Z(1) = Z_g;
Z(2) = Z_l;

% Part3: Solve for Peng-Robinson compressibility factor
TT = T(l);
pp = p(m);
DepartHS = SolveDepartHS(a_T,b,B,Z,DADT,TT,R,R_G);
DH(1) = DepartHS(1);
DH(2) = DepartHS(2);
DS(1) = DepartHS(3);
DS(2) = DepartHS(4);

% Part4: Enthalpy and Entropy of Ideal Gas Change from Reference State to (T(l),p(m))
H_IG = AAA * (T(l) - T_ref) + BBB * (T(l)^2 - T_ref^2)./2 + ...
       CCC * (T(l)^3 - T_ref^3)./3 + DDD * (T(l)^4 - T_ref^4)./4;
S_IG = AAA * log(T(l) ./ T_ref) + BBB * (T(l) - T_ref) + ...
       CCC * (T(l)^2 - T_ref^2)./2 + DDD * (T(l)^3 - T_ref^3)./3;
S_IG = S_IG - R * log(p(m)/p_ref);
H = DH + H_IG;
S = DS + S_IG;

%% Optput the result.
Compressibility = Z;
%{
fprintf('Temperature in this ondition is %4.1f K.\n', T(l));                          % K
fprintf('Pressure in this condition is %4.1f kPa.\n', p(m));                          % kPa
fprintf('Enthalpy of saturated vapor in this condition is %f J/mol.\n', H(1));     % J/mol
fprintf('Enthalpy of saturated liquid in this condition is %f J/mol.\n', H(2));    % J/mol
fprintf('Entropy of saturated vapor in this condition is %f J/mol*K.\n', S(1));    % J/mol*K
fprintf('Entropy of saturated liquid in this condition is %f J/mol*K.\n', S(2));   % J/mol*K

SpecifyVolume = Z * 1e3 * R_G * T(l) ./ p(m);
fprintf('Specify Volume of saturated vapor in this condition is %f.\n',SpecifyVolume(1)); % m^3/kmol
fprintf('Specify Volume of saturated liquid in this condition is %f.\n',SpecifyVolume(2)); % m^3/kmol

fprintf('Compressibility Factor of saturated vapor in this condition is %f.\n',Z(1));
fprintf('Compressibility Factor of saturated liquid in this condition is %f.\n',Z(2));
%}

% Define SubFunction area -------------------------------------------------

% SubFunction1 SolveFugacity:
function Fugacity = SolveFugacity(A,B,Z,pp)
ParaFuga2 = A/B/sqrt(8);
for i=1:2
    ParaFuga1(i) = log((Z(i)+(1+sqrt(2))*B)/(Z(i)+(1-sqrt(2))*B));
    LogFugacityCoeff(i) = (Z(i)-1)-log(Z(i)-B)-ParaFuga1(i)*ParaFuga2;
    FugacityCoeff(i) = exp(LogFugacityCoeff(i));
    Fugacity(i) = FugacityCoeff(i) * pp;
end
clear i;
end

% SubFunction2 SolveDepartHS:
function DepartHS = SolveDepartHS(a_T,b,B,Z,DADT,TT,R,R_G)
for j=1:2,
    ParaFuga1(j) = log((Z(j) + (1+sqrt(2))*B) ./ (Z(j) + (1-sqrt(2))*B));
    DH(j) = (TT * DADT - a_T) * ParaFuga1(j)/b/sqrt(8);
    DH(j) = R * TT * (Z(j)-1) + DH(j) * R ./ R_G;
    % Gas Enthalpy Departure , EQN 6.4-29
    DS(j) = DADT * ParaFuga1(j) ./ b ./ sqrt(8);
    DS(j) = R * log((Z(j)-B)) + DS(j) * R ./ R_G;
    % Gas Entropy Departure , EQN 6.4-30
end
DepartHS = [DH DS];
clear j;
end

%{
% SubFunction3 ZZroot: Solve the Equation of State.
function [ZZ,D] = ZZroot(A,B)
V(1) = 1;
V(2) = - 1 + B;
V(3) = A - B * (3 * B + 2);
V(4) = B * (B * B + B - A);
ZZ = roots(V);

v(1) = 2/27 * V(2)^3 - V(2)*V(3)/3 + V(4);
v(2) = V(3) - V(2)^2/3;
D = v(1)^2/4 + v(2)^3/27;

for m = 1:3 % Get rid off the imag root -----------------------------------
    if imag(ZZ(m)) ~= 0
        ZZ(m) = 0;
    end
end

ZZ = sort(ZZ);

if abs(ZZ(1)) < 1e-8
    ZZ(1) = ZZ(3);
end
if abs(ZZ(3)) < 1e-8
    ZZ(3) = ZZ(1);
end

end
%}

% SubFunction4 ZZroot2: Solve the Equation of State using Determining Equation.
function [Z_g,Z_l] = ZZroot2(A,B)
c1 = B - 1;
c2 = A - 3 * B^2 - 2 * B;
c3 = B^3 + B^2 - A*B;
q = 2/27 * c1^3 - c1*c2/3 +c3;
r = c2 - c1^2 / 3;
D = q^2 / 4 + r^3 / 27;                                %
OMEGA1 = (-1 + sqrt(3) * 1i) / 2;
OMEGA2 = (-1 - sqrt(3) * 1i) / 2;
Z1 = (-q/2 + sqrt((q/2)^2 + (r/3)^3))^(1/3) + ...
     (-q/2 - sqrt((q/2)^2 + (r/3)^3))^(1/3);
Z2 = OMEGA1 * (-q/2 + sqrt((q/2)^2 + (r/3)^3))^(1/3) + ...
     OMEGA2 * (-q/2 - sqrt((q/2)^2 + (r/3)^3))^(1/3);
Z3 = OMEGA2 * (-q/2 + sqrt((q/2)^2 + (r/3)^3))^(1/3) + ...
     OMEGA1 * (-q/2 - sqrt((q/2)^2 + (r/3)^3))^(1/3);
ZZZ = [Z1 Z2 Z3]                                      % 
ZZZ = ZZZ(imag(ZZZ)==0);                               %
if D > 0
    Z_g = ZZZ(1);
    Z_l = 0;
else
    Z_g = max(ZZZ);
    Z_l = min(ZZZ);
end
end





