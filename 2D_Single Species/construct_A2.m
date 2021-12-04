function A = construct_A2(solution, variableLocation, Dx)
%--- This function gives the sparse matrix resulted from discritization in x2 direction & linearization
%--- solution is already permutred according to x2-direction ordering. Therefore, the code is identical to "construct_A1" if the equations
%--- are  symmetric in each direction, with the only difference being the use of "circshift" instead of "reflect" to apply periodic BC
%--- in x2-direction instead of reflecting BC. Here we should use D(2,2) instead of D(1,1) as in construct_A1.

global D 
sigma_squared = D(2,2);

segmentSize = numel(solution(:,:,1)); % number of (A)ij segments
numNonZero = 1 * (5 + 10 + 14) * segmentSize; %number of non-zero elements in A.
odeSystemSize = numel(solution);

solution_minusTwo = circshift(solution,[2 0 0]); % solution(i-2,j)  circshift also applies the periodic boundary condition
solution_minusOne = circshift(solution,[1 0 0]); % solution(i-1,j)  
solution_plusOne = circshift(solution,[-1 0 0]); % solution(i+1,j)  
solution_plusTwo = circshift(solution,[-2 0 0]); % solution(i+2,j)  


dS = ( -solution_plusTwo + 8*solution_plusOne - 8*solution_minusOne + solution_minusTwo ) / (12 * Dx);

%---creating linear arrays
N1 = reshape(solution(:,:,1), [], 1);
dQ = reshape(dS(:,:,2), [], 1);
dN_N = reshape(dS(:,:,1), [], 1) ./ N1;
dQ_N = dQ ./ N1;
dV_N = reshape(dS(:,:,3), [], 1) ./ N1;

rows = zeros(numNonZero, 1); % row indices of nonzero elements of A
columns = zeros(numNonZero, 1); % column indicis of nonzero elements of A
elements = zeros(numNonZero, 1); % nonzero elements of A

segmentCounter = 0;
%% diagonal blocks of A (center block of (A)ij)===========================================================
firstColumn = variableLocation(1,:); % indices of the first (left) columns of the center block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -2 * dN_N .* dQ_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -2 * dN_N .* dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma_squared;
segmentCounter = segmentCounter + 1;


%% 1st lower diagonal of A (first left block of (A)ij)====================================================
firstColumn = circshift(variableLocation(1,:), 1); % indices of the first (left) columns of the first left block of (A)ij. circshift applies the periodic boundary condition

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dQ_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN_N ) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * (-8) / (12*Dx) * dQ * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN_N ) * sigma_squared;
segmentCounter = segmentCounter + 1;


%% 2nd lower diagonal of A (secend left block of (A)ij)======================================================
firstColumn = circshift(variableLocation(1,:), 2); % indices of the first (left) columns of the second left block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dQ_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN_N ) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 / (12*Dx) * dQ * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN_N ) * sigma_squared;
segmentCounter = segmentCounter + 1;


%% 1st upper diagonal of A (first right block of (A)ij)======================================================
firstColumn = circshift(variableLocation(1,:), -1); % indices of the first (left) columns of the first right block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dQ_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN_N ) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * 8 / (12*Dx) * dQ * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN_N ) * sigma_squared;
segmentCounter = segmentCounter + 1;


%% 2nd upper diagonal of A (secend right block of (A)ij)=====================================================
firstColumn = circshift(variableLocation(1,:), -2); % indices of the first (left) columns of the second right block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dQ_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN_N ) * sigma_squared;
segmentCounter = segmentCounter + 1;


segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * (-1) / (12*Dx) * dQ * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN_N ) * sigma_squared;
segmentCounter = segmentCounter + 1;


%% Constructing sparse A ==================================================================================

A = sparse(rows, columns, elements, odeSystemSize, odeSystemSize);

end

