function [xVertex, alph] = geoDist(x1,x2,deltaX,N)

%% Calculate the total distance between the points
L = hypot(x1(1)-x2(1),x1(2)-x2(2));

%% Set the initial bounds for the alpha value

alpha_low = 0;
alpha_high = 10;

%% Loop through the binary search 40 times.
% You could alternatively use a while condition
for i = 1:40
    
    % Determine the mid point of the interval
    alph = (alpha_low+alpha_high)/2;
    
    % Calculate the vertex positions
    dx = deltaX*(alph.^(0:N-1));
    
    % Set the final total length of the boundary
    L_calc = sum(dx);         
    
    % Perform the binary criterion
    if L_calc > L
        alpha_high = alph;
    else
        alpha_low = alph;
    end
    
end

%% Perform a check to confirm the cells aren't growing/shrinking too fast
if alph>2 || alph<1/2
    error('Bad Cells');
end

%% With the distribution of the cell sizes determine vertex positions
x = x1(1) + (x2(1)-x1(1))*(cumsum([0 dx])/L);
y = x1(2) + (x2(2)-x1(2))*(cumsum([0 dx])/L);

xVertex = [x' y'];

%% Set the final position to be x2 so that the vertex span precisely the range provided
xVertex(end,:) = x2;