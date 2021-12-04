function F = construct_F1(solution, variableLocation, Dx)
%--- This function returnes the value of the F1 component of ODE vector field. The ordering is in x1 direction.

global D
sigma_squared = D;

N = solution(:,1);

solution_minusTwo = reflect(solution,[2 0]); % solution(i-2)  reflect also applies the reflecting boundary condition
solution_minusOne = reflect(solution,[1 0]); % solution(i-1)  
solution_plusOne = reflect(solution,[-1 0]); % solution(i+1)  
solution_plusTwo = reflect(solution,[-2 0]); % solution(i+2)

dS = ( -solution_plusTwo + 8*solution_plusOne - 8*solution_minusOne + solution_minusTwo ) / (12 * Dx);
dN = dS(:,1);
dQ = dS(:,2);
dV = dS(:,3);

d2S = ( -solution_plusTwo + 16*solution_plusOne - 30*solution + 16*solution_minusOne - solution_minusTwo ) / (12 * Dx^2);
d2S = d2S'; % reorders 2-dim array d2S so that linear indexing d2S(:) corresponds to current ordeing in U

F = d2S(:);
F(variableLocation(1,:)) = F(variableLocation(1,:)) * sigma_squared;
F(variableLocation(2,:)) = ( F(variableLocation(2,:)) + 2 * dN .* dQ ./ N ) * sigma_squared;
F(variableLocation(3,:)) = ( F(variableLocation(3,:)) + 2 * dN .* dV ./ N + 2 * dQ.^2 ) * sigma_squared;

end