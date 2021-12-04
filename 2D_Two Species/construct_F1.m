function F = construct_F1(solution, variableLocation, Dx)
%--- This function returnes the value of the F1 component of ODE vector field. The ordering is in x1 direction.

global D1 D2
sigma1_squared = D1(1,1);
sigma2_squared = D2(1,1);

N1 = reshape(solution(:,:,1), [], 1);
N2 = reshape(solution(:,:,4), [], 1);

solution_minusTwo = reflect(solution,[2 0 0]); % solution(i-2,j)  reflect also applies the reflecting boundary condition
solution_minusOne = reflect(solution,[1 0 0]); % solution(i-1,j)  
solution_plusOne = reflect(solution,[-1 0 0]); % solution(i+1,j)  
solution_plusTwo = reflect(solution,[-2 0 0]); % solution(i+2,j)

dS = ( -solution_plusTwo + 8*solution_plusOne - 8*solution_minusOne + solution_minusTwo ) / (12 * Dx);
dN1 = reshape(dS(:,:,1), [], 1);
dQ1 = reshape(dS(:,:,2), [], 1);
dV1 = reshape(dS(:,:,3), [], 1);
dN2 = reshape(dS(:,:,4), [], 1);
dQ2 = reshape(dS(:,:,5), [], 1);
dV2 = reshape(dS(:,:,6), [], 1);

d2S = ( -solution_plusTwo + 16*solution_plusOne - 30*solution + 16*solution_minusOne - solution_minusTwo ) / (12 * Dx^2);
d2S = permute(d2S,[3,1,2]); % reorders 3-dim array d2S so that linear indexing d2S(:) corresponds to current ordeing in U

F = d2S(:);
F(variableLocation(1,:)) = F(variableLocation(1,:)) * sigma1_squared;
F(variableLocation(2,:)) = ( F(variableLocation(2,:)) + 2 * dN1 .* dQ1 ./ N1 ) * sigma1_squared;
F(variableLocation(3,:)) = ( F(variableLocation(3,:)) + 2 * dN1 .* dV1 ./ N1 + 2 * dQ1.^2 ) * sigma1_squared;
F(variableLocation(4,:)) = F(variableLocation(4,:)) * sigma2_squared;
F(variableLocation(5,:)) = ( F(variableLocation(5,:)) + 2 * dN2 .* dQ2 ./ N2 ) * sigma2_squared;
F(variableLocation(6,:)) = ( F(variableLocation(6,:)) + 2 * dN2 .* dV2 ./ N2 + 2 * dQ2.^2 ) * sigma2_squared;
end