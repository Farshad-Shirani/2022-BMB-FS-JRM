function F = construct_F1(solution, variableLocation, Dx)
%--- This function returnes the value of the F1 component of ODE vector field. The ordering is in x1 direction.

global D1 D2
sigma1_squared = D1;
sigma2_squared = D2;

N1 = solution(:,1);
N2 = solution(:,4);

solution_minusTwo = reflect(solution,[2 0]); % solution(i-2)  reflect also applies the reflecting boundary condition
solution_minusOne = reflect(solution,[1 0]); % solution(i-1)  
solution_plusOne = reflect(solution,[-1 0]); % solution(i+1)  
solution_plusTwo = reflect(solution,[-2 0]); % solution(i+2)

dS = ( -solution_plusTwo + 8*solution_plusOne - 8*solution_minusOne + solution_minusTwo ) / (12 * Dx);
dN1 = dS(:,1);
dQ1 = dS(:,2);
dV1 = dS(:,3);
dN2 = dS(:,4);
dQ2 = dS(:,5);
dV2 = dS(:,6);

d2S = ( -solution_plusTwo + 16*solution_plusOne - 30*solution + 16*solution_minusOne - solution_minusTwo ) / (12 * Dx^2);
d2S = d2S'; % reorders 2-dim array d2S so that linear indexing d2S(:) corresponds to current ordeing in U

F = d2S(:);
F(variableLocation(1,:)) = F(variableLocation(1,:)) * sigma1_squared;
F(variableLocation(2,:)) = ( F(variableLocation(2,:)) + 2 * dN1 .* dQ1 ./ N1 ) * sigma1_squared;
F(variableLocation(3,:)) = ( F(variableLocation(3,:)) + 2 * dN1 .* dV1 ./ N1 + 2 * dQ1.^2 ) * sigma1_squared;
F(variableLocation(4,:)) = F(variableLocation(4,:)) * sigma2_squared;
F(variableLocation(5,:)) = ( F(variableLocation(5,:)) + 2 * dN2 .* dQ2 ./ N2 ) * sigma2_squared;
F(variableLocation(6,:)) = ( F(variableLocation(6,:)) + 2 * dN2 .* dV2 ./ N2 + 2 * dQ2.^2 ) * sigma2_squared;
end