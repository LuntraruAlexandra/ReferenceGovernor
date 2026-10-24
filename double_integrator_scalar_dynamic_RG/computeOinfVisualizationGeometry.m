function [Oinf,vMin,vMax,numSlices,referenceValues,slices,colorMap] = ...
    computeOinfVisualizationGeometry(Oinf,H,h,n,m,lpOptions)

%% 11. ADMISSIBLE REFERENCE INTERVAL
% This part is not actually related to RG
% Used this for the geometric representation of Oinf
% Minimize the first reference component, then maximize it.
% The plotting and RG example below assume a scalar reference (m = 1).
referenceObjective = [1; zeros(n,1)];
[~,minimumObjective,minExitFlag] = linprog( ...
    referenceObjective,H,h,[],[],[],[],lpOptions);
[~,maximumNegativeObjective,maxExitFlag] = linprog( ...
    -referenceObjective,H,h,[],[],[],[],lpOptions);

if minExitFlag ~= 1 || maxExitFlag ~= 1
    error('Unable to compute the admissible reference bounds.');
end

vMin = minimumObjective;
vMax = -maximumNegativeObjective;
Oinf.vMin = vMin;
Oinf.vMax = vMax;
fprintf('\nAdmissible reference interval:\n');
fprintf('[%.6f, %.6f]\n',vMin,vMax);

%% 12. COMPUTE CROSS-SECTIONS OF O_INF*

% The following geometry is intentionally specific to n=2, m=1.
assert(n == 2 && m == 1, 'Cross-section plots require n=2, m=1.');

% For fixed scalar v: H(:,1:m)*v + H(:,m+1:m+n)*x <= h.
% Therefore, Hx*x <= hx defines a polygon in the (x1,x2) plane.
numSlices = 20;
endpointMargin = 1e-4;

if vMax-vMin <= 2*endpointMargin
    error('The admissible reference interval is too narrow.');
end

referenceValues = linspace(vMin+endpointMargin,vMax-endpointMargin,numSlices);
slices = cell(numSlices,1);
Hx = H(:,m+1:m+n);

for sliceIndex = 1:numSlices
    fixedReference = referenceValues(sliceIndex);
    hx = h-H(:,1:m)*fixedReference;
    vertices = [];

    % Intersect every pair of boundary lines.
    for rowA = 1:size(Hx,1)-1
        for rowB = rowA+1:size(Hx,1)
            lineMatrix = Hx([rowA rowB],:);
            lineBounds = hx([rowA rowB]);

            % Skip parallel or nearly parallel lines.
            if rcond(lineMatrix) < 1e-11
                continue;
            end

            candidateVertex = lineMatrix\lineBounds;

            % Keep only intersections satisfying every inequality.
            if all(Hx*candidateVertex <= hx+1e-7)
                vertices = [vertices;candidateVertex']; %#ok<AGROW>
            end
        end
    end

    if size(vertices,1) < 3
        continue;
    end

    % Remove duplicate vertices and order the convex boundary.
    vertices = uniquetol(vertices,1e-7,'ByRows',true);
    if size(vertices,1) < 3
        continue;
    end

    hullIndices = convhull(vertices(:,1),vertices(:,2));
    slices{sliceIndex} = vertices(hullIndices,:);
end

Oinf.referenceValues = referenceValues;
Oinf.slices = slices;

%% 13-14. PLOTTING IS COMBINED INTO ONE FIGURE AT THE END
% The slice geometry is preserved; no intermediate figures are opened.
colorMap = parula(numSlices);

%% 15. NUMERICAL VERIFICATION OF SLICE VERTICES

maximumViolation = -inf;
for sliceIndex = 1:numSlices
    sliceVertices = slices{sliceIndex};
    if isempty(sliceVertices)
        continue;
    end

    for vertexIndex = 1:size(sliceVertices,1)
        augmentedState = [referenceValues(sliceIndex); ...
                          sliceVertices(vertexIndex,:)'];
        constraintViolation = max(H*augmentedState-h);
        maximumViolation = max(maximumViolation,constraintViolation);
    end
end

if isfinite(maximumViolation)
    fprintf('\nMaximum constraint violation: %.3e\n',maximumViolation);
else
    warning('No valid cross-sections were generated.');
end

disp('O_inf* computation and slice visualization completed.');

%% 16. CHECK GEOMETRIC BOUNDS OF THE SLICES

allVertices = [];
for sliceIndex = 1:numSlices
    if ~isempty(slices{sliceIndex})
        allVertices = [allVertices;slices{sliceIndex}]; %#ok<AGROW>
    end
end

if ~isempty(allVertices)
    fprintf('\nGeometric bounds of the slices:\n');
    fprintf('x1: [%.8f, %.8f]\n', ...
            min(allVertices(:,1)),max(allVertices(:,1)));
    fprintf('x2: [%.8f, %.8f]\n', ...
            min(allVertices(:,2)),max(allVertices(:,2)));
    fprintf('Theoretical x1 bounds: [-1, 1]\n');
    fprintf('Theoretical x2 bounds: [-0.1, 0.1]\n');
else
    warning('Geometric bounds cannot be checked: no slice vertices exist.');
end
end
