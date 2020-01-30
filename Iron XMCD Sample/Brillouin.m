function [ Bj ] = Brillouin(B,T,J,g)
%Brillouin Function
%   Calculates the vlaues of the Brillouin function for different B
%   Parameters:
%   1) B - magnetic field (array of values)
%   2) J - total orbital angular momentum
%   3) T - temperature

kB = 8.6173324e-5; %Boltzmann constant eV/K
mu_B = 5.7883818066e-5; %Bohr magneton eV/T

x=g*J*(mu_B*B)./(kB*T);

%Makes two constants to make the equation below simpler
%tJ is the constant used in the second half equal to 1/2J
tJ = 1/(2*J);
%tJp is the constant used in teh first half equal to (2J+1)/2J
tJp = (2*J+1)*tJ;

%Calculates Bj
Bj = g*J*(tJp*coth(tJp*x)-tJ*coth(x*tJ));

end

