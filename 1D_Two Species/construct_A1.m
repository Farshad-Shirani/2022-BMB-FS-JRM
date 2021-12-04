function A = construct_A1(solution, variableLocation, Dx)
%--- This function gives the sparse matrix resulted from discritization in x1 direction & linearization

global D1 D2
sigma1_squared = D1;
sigma2_squared = D2;

segmentSize = numel(solution(:,1)); % number of (A)ij segments
numNonZero = 2 * (5 + 10 + 14) * segmentSize; %number of non-zero elements in A.
odeSystemSize = numel(solution);

solution_minusTwo = reflect(solution,[2 0]); % solution(i-2)  reflect also applies the reflecting boundary condition
solution_minusOne = reflect(solution,[1 0]); % solution(i-1)  
solution_plusOne = reflect(solution,[-1 0]); % solution(i+1)  
solution_plusTwo = reflect(solution,[-2 0]); % solution(i+2)  


dS = ( -solution_plusTwo + 8*solution_plusOne - 8*solution_minusOne + solution_minusTwo ) / (12 * Dx);

%---creating linear arrays
N1 = solution(:,1);
N2 = solution(:,4);
dQ1 = dS(:,2);
dQ2 = dS(:,5);
dN1_N1 = dS(:,1)./ N1;
dQ1_N1 = dQ1 ./ N1;
dV1_N1 = dS(:,3) ./ N1;
dN2_N2 = dS(:,4) ./ N2;
dQ2_N2 = dQ2 ./ N2;
dV2_N2 = dS(:,6)./ N2;

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
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -2 * dN1_N1 .* dQ1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -2 * dN1_N1 .* dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) =  -2 * dN2_N2 .* dQ2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = -2 * dN2_N2 .* dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma2_squared;
segmentCounter = segmentCounter + 1;


%% 1st lower diagonal of A (first left block of (A)ij)====================================================
firstColumn = reflect(variableLocation(1,:), 1); % indices of the first (left) columns of the first left block of (A)ij. reflect applies the reflecting boundary condition

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dQ1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN1_N1 ) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * (-8) / (12*Dx) * dQ1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN1_N1 ) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dQ2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN2_N2 ) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 4 * (-8) / (12*Dx) * dQ2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN2_N2 ) * sigma2_squared;
segmentCounter = segmentCounter + 1;


%% 2nd lower diagonal of A (secend left block of (A)ij)======================================================
firstColumn = reflect(variableLocation(1,:), 2); % indices of the first (left) columns of the second left block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dQ1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN1_N1 ) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 / (12*Dx) * dQ1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN1_N1 ) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dQ2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN2_N2 ) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 4 / (12*Dx) * dQ2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN2_N2 ) * sigma2_squared;
segmentCounter = segmentCounter + 1;


%% 1st upper diagonal of A (first right block of (A)ij)======================================================
firstColumn = reflect(variableLocation(1,:), -1); % indices of the first (left) columns of the first right block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dQ1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN1_N1 ) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * 8 / (12*Dx) * dQ1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN1_N1 ) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dQ2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN2_N2 ) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 4 * 8 / (12*Dx) * dQ2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN2_N2 ) * sigma2_squared;
segmentCounter = segmentCounter + 1;


%% 2nd upper diagonal of A (secend right block of (A)ij)=====================================================
firstColumn = reflect(variableLocation(1,:), -2); % indices of the first (left) columns of the second right block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dQ1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN1_N1 ) * sigma1_squared;
segmentCounter = segmentCounter + 1;


segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dV1_N1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * (-1) / (12*Dx) * dQ1 * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN1_N1 ) * sigma1_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(4,:); % 4th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dQ2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;


segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(5,:); % 5th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN2_N2 ) * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 3; % 4th column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dV2_N2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 4; % 5th column in each block
elements(segmentBegining : segmentEnding) = 4 * (-1) / (12*Dx) * dQ2 * sigma2_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(6,:); % 6th row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 5; % 6th column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2  * (-1) / (12*Dx) * dN2_N2 ) * sigma2_squared;
segmentCounter = segmentCounter + 1;

%% Constructing sparse A ==================================================================================

A = sparse(rows, columns, elements, odeSystemSize, odeSystemSize);

end

