clear;
clc;
close all;

%% =========================================================
% REFERENCES
%
% [1] E. G. Gilbert, K. T. Tan (1991)
% "Linear Systems with State and Control Constraints:
%  The Theory and Application of Maximal Output Admissible Sets"
%
% Algorithms 3.1 / 3.2:
% Iterative construction of O_inf and LP-based redundancy tests.
%
% [2] E. G. Gilbert, I. Kolmanovsky, K. T. Tan (1995)
% "Discrete-Time Reference Governors and the Nonlinear
%  Control of Systems with State and Control Constraints"
%
% Theorem 3.1, equations (17)-(22):
% Finitely determined inner approximation O_inf* using
% the tightened steady-state constraint set Y(epsilon).
%
% [3] Garone, Di Cairano, Kolmanovsky (2017)
% "Reference and Command Governors for Systems
%  with Constraints: A Survey on Theory and Applications"
%
% Example: double integrator with controller
% u = 0.9170*(v-x1) - 1.6821*x2
%
% OBJECTIVES:
% Compute O_inf* for z = [w; x1; x2], with constant w in predictions.
% Remove redundant inequalities.
% Visualize 3D slices and their 2D contours.
%
% A dynamic scalar Reference Governor is implemented after the O_inf* computation.
%% =========================================================

%% 1. DISCRETE-TIME DOUBLE-INTEGRATOR MODEL

Ts = 0.1;

% Model used in the paper's example:
% x1(k+1) = x1(k) + 0.1*x2(k)
% x2(k+1) = x2(k) + 0.1*u(k)
Ad = [1 Ts;
      0 1];
Bd = [0;
      Ts];

%% 2. CONTROLLER GIVEN IN THE PAPER

% u(k) = 0.9170*(v(k)-x1(k)) - 1.6821*x2(k)
% Equivalently: u(k) = -K*x(k) + K(1)*v(k)
K = [0.9170 1.6821];

% Closed-loop dynamics: x(k+1) = Ah*x(k) + Bh*v(k)
Ah = Ad - Bd*K;
Bh = Bd*K(1);

% Dimensions inferred from the closed-loop model.
% n: plant-state dimension; m: admissible-reference dimension.
n = size(Ah,1);
m = size(Bh,2);

assert(size(Ah,2) == n, 'Ah must be square.');
assert(size(Bh,1) == n, 'Bh must have n rows.');

closedLoopEigenvalues = eig(Ah);

disp('Controller gain K from the paper:'); disp(K);
disp('Matrix Ah:'); disp(Ah);
disp('Matrix Bh:'); disp(Bh);
disp('Closed-loop eigenvalues:'); disp(closedLoopEigenvalues);

if any(abs(closedLoopEigenvalues) >= 1)
    error('The closed-loop system is not Schur stable.');
end

%% 3. SYSTEM CONSTRAINTS

% |x1| <= 1, |x2| <= 0.1, |u| <= 0.1
stateMin = [-1; -0.1];
stateMax = [ 1;  0.1];
inputMin = -0.1;
inputMax =  0.1;

% Constrained output: y = [x1; x2; u] = Ch*x + Dh*v
% y was chosen to contain all constraints
Ch = [eye(n);
      -K];
Dh = [zeros(n,m);
      K(1)];

% p: number of constrained outputs (position, velocity, input).
p = size(Ch,1);
assert(all(size(Dh) == [p m]), 'Dh must be p-by-m.');

outputMin = [stateMin; inputMin];
outputMax = [stateMax; inputMax];

% Y = {y : F*y <= g}
F = [eye(p);
     -eye(p)];

% q: number of scalar output inequalities F*y <= g.
q = size(F,1);
g = [outputMax;
     -outputMin];

disp('Matrix Ch:'); disp(Ch);
disp('Matrix Dh:'); disp(Dh);

%% 4. AUGMENTED SYSTEM - REFERENCE [2]

% z = [v; x1; x2] for this example (dimension m+n)
% Assume that v remains constant when computing O_inf:
% v(k+1) = v(k)
% x(k+1) = Ah*x(k) + Bh*v(k)
% z(k+1) = AG*z(k), y(k) = CG*z(k)
AG = [eye(m) zeros(m,n);
      Bh     Ah];
CG = [Dh Ch];

disp('Matrix AG:'); disp(AG);
disp('Matrix CG:'); disp(CG);
disp('Eigenvalues of AG:'); disp(eig(AG));

%% 5. STEADY-STATE OUTPUT - REFERENCE [2]

% For constant v:
% x_eq = (I-Ah)^(-1)*Bh*v
% y_eq = H0*v
% H0 = Dh + Ch*(I-Ah)^(-1)*Bh
H0 = Dh + Ch*((eye(n)-Ah)\Bh);
disp('Steady-state output matrix H0:'); disp(H0);

%% 6. STEADY-STATE CONSTRAINT TIGHTENING - THEOREM 3.1

% Considered F instead of H and g instead of y:
% [2], equation (17): Y(epsilon) = {y : F*y <= g-epsilon}
% Only the steady-state output is tightened.
% Predicted transient outputs must satisfy the original set Y.
epsilon = 0.01;

if epsilon <= 0 || epsilon >= min(g)
    error('The selected epsilon is not admissible.');
end

% H0*v belongs to Y(epsilon): F*H0*v <= g-epsilon
% See [2], equations (21)-(22).
Hsteady = [F*H0 zeros(q,n)];
hSteady = g - epsilon;

%% 7. INITIALIZE O_0* - REFERENCES [1], [2]

% O_0* = {z : H0*v belongs to Y(epsilon), CG*z belongs to Y}
H = [Hsteady;
     F*CG];
h = [hSteady;
     g];
% These matrices for writing the two inequalities in 21/22 from [2]
% in a single inequality

% Metadata for each inequality:
% constraintType = 0: steady-state inequality
% constraintType = 1: predicted-output inequality
constraintType = [zeros(q,1);
                  ones(q,1)];

% predictionStep = -1: steady-state inequality
% predictionStep =  0: output constraint at k = 0
predictionStep = [-ones(q,1);
                   zeros(q,1)];

% Index of the corresponding original inequality in F
constraintIndex = [(1:q)';
                   (1:q)'];

%% 8. LINEAR PROGRAMMING OPTIONS AND INITIAL REDUCTION

maxIterations = 500;
redundancyTolerance = 1e-8;
zeroTolerance = 1e-10;
lpOptions = optimoptions('linprog','Display','none');

initialInequalityCount = size(H,1);
[H,h,constraintType,predictionStep,constraintIndex] = ...
    removeRedundantInequalities(H,h,constraintType,predictionStep, ...
                                constraintIndex,lpOptions, ...
                                redundancyTolerance,zeroTolerance);

fprintf('\nInitial inequalities: %d\n',initialInequalityCount);
fprintf('Inequalities after reducing O_0*: %d\n',size(H,1));

%% 9. ITERATIVE CONSTRUCTION OF O_INF*

% Based on [1], Algorithms 3.1 / 3.2, and [2], Theorem 3.1.
% At each prediction step k, test F*CG*AG^k*z <= g.
% Stop once every new inequality is implied by the current set.
AGpower = AG;  % Constraints at k=0 are already included.
tauStar = NaN;

for k = 1:maxIterations
    Hnew = F*CG*AGpower; % again Hnew * y - hNew is f
    hNew = g;
    allRedundant = true;

    % Test each new inequality against the current polytope.
    for j = 1:size(Hnew,1)
        % Extract each new inequality one by one
        objectiveRow = Hnew(j,:);
        upperBound = hNew(j);

        % Maximize objectiveRow*z subject to H*z <= h.
        [~,objectiveMin,exitFlag] = linprog( ...
            -objectiveRow',H,h,[],[],[],[],lpOptions);

        if exitFlag ~= 1
            error('LP failed at k=%d, inequality %d (flag=%d).', ...
                  k,j,exitFlag);
        end

        maximumValue = -objectiveMin;
        % now check if the new maximum is greater than before
        % if it does go over the constraint then consider it non redundant
        if maximumValue > upperBound + redundancyTolerance
            allRedundant = false;
        end
    end
    % If all new inequalities are redundant then finite determination is reached.
    if allRedundant
        tauStar = k-1;
        fprintf('\nFinite determination achieved.\n');
        fprintf('tau* = %d\n',tauStar);
        break;
    end

    % Append constraints for prediction step k.
    % This appends only the non redundant inequalities
    H = [H;Hnew];
    h = [h;hNew];
    % Update the information about new constraints
    constraintType = [constraintType;ones(q,1)];
    predictionStep = [predictionStep;k*ones(q,1)];
    constraintIndex = [constraintIndex;(1:q)'];
    countBeforeReduction = size(H,1);

    % The above code still contains redundant inequalities
    % because only if a single one was foun to be non redundant then it 
    % would append all of the inequalities at that step
    % so an additional function to remove all redundant ineqs
    [H,h,constraintType,predictionStep,constraintIndex] = ...
        removeRedundantInequalities(H,h,constraintType,predictionStep, ...
                                    constraintIndex,lpOptions, ...
                                    redundancyTolerance,zeroTolerance);

    fprintf('k = %3d | before reduction = %3d | after reduction = %3d\n', ...
            k,countBeforeReduction,size(H,1));
    AGpower = AG*AGpower;  % AG^k -> AG^(k+1)
end

if isnan(tauStar)
    error('Finite determination was not established within maxIterations.');
end

%% 10. FINAL POLYTOPE REPRESENTATION

% O_inf* = {z : H*z <= h}, where z = [v; x1; x2] (m+n components).
Oinf.H = H;
Oinf.h = h;
Oinf.AG = AG;
Oinf.CG = CG;
Oinf.H0 = H0;
Oinf.epsilon = epsilon;
Oinf.tauStar = tauStar;
Oinf.constraintType = constraintType;
Oinf.predictionStep = predictionStep;
Oinf.constraintIndex = constraintIndex;

fprintf('\n====================================\n');
fprintf('FINAL O_INF* RESULT\n');
fprintf('====================================\n');
fprintf('Final number of inequalities: %d\n',size(H,1));
fprintf('tau* = %d\n',tauStar);
disp('Matrix H:'); disp(H);
disp('Vector h:'); disp(h);

inequalityTable = table( ...
    (1:size(H,1))', ...
    constraintType, ...
    predictionStep, ...
    constraintIndex, ...
    'VariableNames', ...
    {'InequalityID', 'ConstraintType', ...
     'PredictionStep', 'ConstraintIndex'});

disp('Non-redundant inequalities:');
disp(inequalityTable);
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

%% 17. DYNAMIC SCALAR REFERENCE GOVERNOR - [2], EQUATIONS (10)-(11)

% Augmented state xG(k) = [w(k); x1(k); x2(k)].
% w(k) is the previously accepted reference; r(k) is the desired reference.
% xG(k+1) = AG*xG(k) + BG*K(k)*(r(k)-w(k)).
% y(k) = CG*xG(k) must satisfy F*y(k) <= g.
% Choose the largest kappa in [0,1] such that xG(k+1) belongs to O_inf*.
% The plant evolves using w(k), not the newly computed w(k+1).

BG = [eye(m); zeros(n,m)];
Oinf.BG = BG;

simulationDuration = 40;                 % seconds
simulationSteps = round(simulationDuration/Ts);
time = (0:simulationSteps)*Ts;

% Desired position reference: zero initially, then a step to 0.8.
% Change 0.8 to 1.2 to illustrate an inadmissible steady-state target.
desiredReference = zeros(1,simulationSteps+1);
desiredReference(time >= 1) = 0.8;

% Initial augmented state must belong to O_inf*.
initialReference = 0;
initialPlantState = zeros(n,1);
xG = [initialReference; initialPlantState];

% This simulation and the 2D slice plots are for m = 1, n = 2.
assert(m == 1 && n == 2, ...
    'The plotting and scalar RG example require m=1 and n=2.');
feasibilityTolerance = 1e-7;

if any(H*xG > h+feasibilityTolerance)
    error('Initial augmented state is outside O_inf*.');
end

augmentedHistory = zeros(n+m,simulationSteps+1);
outputHistory = zeros(p,simulationSteps+1);
kappaHistory = zeros(1,simulationSteps);
nextStateViolation = zeros(1,simulationSteps);
augmentedHistory(:,1) = xG;

for step = 1:simulationSteps
    r = desiredReference(step);
    w = xG(1);
    referenceError = r-w;

    % In the article's timing convention, y(k) uses the current w(k).
    outputHistory(:,step) = CG*xG;

    % xG_next = baseNext + direction*kappa.
    baseNext = AG*xG;
    direction = BG*referenceError;

    % Enforce H*(baseNext + direction*kappa) <= h.
    % Each row yields: coefficient_i*kappa <= margin_i.
    coefficients = H*direction;
    margins = h-H*baseNext;
    kappaLower = 0;
    kappaUpper = 1;

    for rowIndex = 1:size(H,1)
        coefficient = coefficients(rowIndex);
        margin = margins(rowIndex);

        if coefficient > 1e-12
            kappaUpper = min(kappaUpper,margin/coefficient);
        elseif coefficient < -1e-12
            kappaLower = max(kappaLower,margin/coefficient);
        elseif margin < -feasibilityTolerance
            error('RG infeasible at step %d: zero coefficient.',step);
        end
    end

    if kappaLower > kappaUpper+feasibilityTolerance
        error('No feasible kappa in [0,1] at step %d.',step);
    end

    % Largest feasible kappa; numerical clipping to [0,1].
    kappa = max(0,min(1,kappaUpper));
    if kappa < kappaLower-feasibilityTolerance
        error('Numerically infeasible kappa at step %d.',step);
    end

    xGnext = baseNext+direction*kappa;
    nextStateViolation(step) = max(H*xGnext-h);
    if nextStateViolation(step) > feasibilityTolerance
        error('Next augmented state violates O_inf* at step %d.',step);
    end

    kappaHistory(step) = kappa;
    xG = xGnext;
    augmentedHistory(:,step+1) = xG;
end
outputHistory(:,end) = CG*xG;

% Physical constraints: x1, x2 and u are all included in y.
physicalResidual = F*outputHistory-g;
physicalConstraintViolation = max(physicalResidual(:));

fprintf('\n===== DYNAMIC REFERENCE GOVERNOR =====\n');
fprintf('Maximum O_inf* violation: %.3e\n',max(nextStateViolation));
fprintf('Maximum physical constraint violation: %.3e\n', ...
        physicalConstraintViolation);
fprintf('Final desired reference: %.4f\n',desiredReference(end));
fprintf('Final accepted reference: %.4f\n',augmentedHistory(1,end));
fprintf('Final position: %.4f\n',augmentedHistory(2,end));

Oinf.rg.time = time;
Oinf.rg.desiredReference = desiredReference;
Oinf.rg.augmentedHistory = augmentedHistory;
Oinf.rg.outputHistory = outputHistory;
Oinf.rg.kappaHistory = kappaHistory;

%% 18. ONLY FIGURE: RG TRAJECTORY AND O_INF* CROSS-SECTIONS
% Plot the RG trajectory together with the constant-reference cross-sections.
% Mark with X the exact intersections with each plotted plane w = constant.
% Intersections are linearly interpolated between consecutive RG samples.

fig = figure('Color',[1 1 1], ...
       'Name','Oinf slices and RG trajectory', ...
       'NumberTitle','off','WindowStyle','normal', ...
       'MenuBar','figure','ToolBar','figure', ...
       'Position',[100 100 1150 750]);
ax = gca;
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
% Explicit light theme, independent of MATLAB dark-mode defaults.
ax.Color = [1 1 1];
ax.XColor = [0 0 0];
ax.YColor = [0 0 0];
ax.ZColor = [0 0 0];
ax.GridColor = [0.65 0.65 0.65];
ax.GridAlpha = 0.25;
ax.Projection = 'perspective';

% Draw O_inf* cross-sections first, with low opacity.
for sliceIndex = 1:numSlices
    sliceVertices = slices{sliceIndex};
    if isempty(sliceVertices)
        continue;
    end

    fixedReference = referenceValues(sliceIndex);
    fill3(sliceVertices(:,1),sliceVertices(:,2), ...
          fixedReference*ones(size(sliceVertices,1),1), ...
          colorMap(sliceIndex,:), ...
          'FaceAlpha',0.18, ...
          'EdgeColor',colorMap(sliceIndex,:), ...
          'EdgeAlpha',0.45, ...
          'LineWidth',0.8, ...
          'HandleVisibility','off');
end

% Extract the augmented-state trajectory in the plot coordinate order.
positionTrajectory = augmentedHistory(2,:);
velocityTrajectory = augmentedHistory(3,:);
referenceTrajectory = augmentedHistory(1,:);

% Plot a single, clearly visible red trajectory.
trajectoryHandle = plot3(positionTrajectory,velocityTrajectory, ...
      referenceTrajectory,'-','Color',[0.85 0.08 0.08],'LineWidth',3, ...
      'DisplayName','RG trajectory');

% Compute intersections of each trajectory segment with each slice plane.
intersectionPoints = zeros(0,3);
intersectionTimes = zeros(0,1);
intersectionPlane = zeros(0,1);
intersectionTolerance = 1e-10;

for sliceIndex = 1:numSlices
    if isempty(slices{sliceIndex})
        continue;
    end
    planeReference = referenceValues(sliceIndex);

    for step = 1:simulationSteps
        wA = referenceTrajectory(step);
        wB = referenceTrajectory(step+1);

        % Skip segments parallel to the plane, including segments in it.
        if abs(wB-wA) <= intersectionTolerance
            continue;
        end

        alpha = (planeReference-wA)/(wB-wA);
        if alpha < -intersectionTolerance || alpha > 1+intersectionTolerance
            continue;
        end
        alpha = max(0,min(1,alpha));

        xCross = positionTrajectory(step) + alpha* ...
                 (positionTrajectory(step+1)-positionTrajectory(step));
        vCross = velocityTrajectory(step) + alpha* ...
                 (velocityTrajectory(step+1)-velocityTrajectory(step));
        tCross = time(step) + alpha*(time(step+1)-time(step));

        candidate = [xCross vCross planeReference];
        % Avoid duplicate intersections at trajectory vertices.
        if ~isempty(intersectionPoints) && ...
           any(vecnorm(intersectionPoints-candidate,2,2) < 1e-8)
            continue;
        end

        intersectionPoints(end+1,:) = candidate; %#ok<AGROW>
        intersectionTimes(end+1,1) = tCross; %#ok<AGROW>
        intersectionPlane(end+1,1) = sliceIndex; %#ok<AGROW>
    end
end

if ~isempty(intersectionPoints)
    crossHandle = plot3(intersectionPoints(:,1),intersectionPoints(:,2), ...
          intersectionPoints(:,3),'x', ...
          'Color',[0 0 0], ...
          'MarkerSize',11,'LineWidth',2.5, ...
          'DisplayName','Trajectory / slice intersections');
else
    crossHandle = plot3(nan,nan,nan,'x', ...
          'Color',[0 0 0], ...
          'MarkerSize',11,'LineWidth',2.5, ...
          'DisplayName','Trajectory / slice intersections');
end

startHandle = plot3(positionTrajectory(1),velocityTrajectory(1), ...
      referenceTrajectory(1),'o','MarkerSize',9, ...
      'MarkerFaceColor',[0.1 0.65 0.25], ...
      'MarkerEdgeColor','k','DisplayName','Start');
endHandle = plot3(positionTrajectory(end),velocityTrajectory(end), ...
      referenceTrajectory(end),'o','MarkerSize',9, ...
      'MarkerFaceColor',[0.95 0.6 0.05], ...
      'MarkerEdgeColor','k','DisplayName','End');

xlabel('x_1 - Position');
ylabel('x_2 - Velocity');
zlabel('w - Accepted reference');
title('O_\infty^* cross-sections and RG trajectory');
legend([trajectoryHandle crossHandle startHandle endHandle], ...
       'Location','northeastoutside');
view(38,23);
xlim([-1.1 1.1]);
ylim([-0.12 0.12]);
zlim([vMin-0.05 vMax+0.05]);
pbaspect([2 1 2]);
% Use the standard MATLAB figure toolbar to rotate the 3D view.

% Store intersection data for later analysis.
Oinf.rg.intersectionPoints = intersectionPoints;
Oinf.rg.intersectionTimes = intersectionTimes;
Oinf.rg.intersectionPlaneIndex = intersectionPlane;

fprintf('\nRG trajectory intersections with slice planes: %d\n', ...
        size(intersectionPoints,1));
if ~isempty(intersectionPoints)
    intersectionTable = table(intersectionTimes,intersectionPoints(:,3), ...
        intersectionPoints(:,1),intersectionPoints(:,2), ...
        'VariableNames',{'Time','Reference','Position','Velocity'});
    disp(intersectionTable);
end
%% LABEL START, END AND DESIRED REFERENCE

% Start coordinates: (x1, x2, w)
startPoint = [augmentedHistory(2,1), ...
              augmentedHistory(3,1), ...
              augmentedHistory(1,1)];

% End coordinates: (x1, x2, w)
endPoint = [augmentedHistory(2,end), ...
            augmentedHistory(3,end), ...
            augmentedHistory(1,end)];

% Final desired reference
rDesired = desiredReference(end);

% Start label
text(startPoint(1), startPoint(2), startPoint(3), ...
    sprintf('(%.3f, %.3f, %.3f)', startPoint), ...
    'FontSize', 10, ...
    'FontWeight', 'bold', ...
    'Color', [0 0.5 0], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

% End label
text(endPoint(1), endPoint(2), endPoint(3), ...
    sprintf('(%.3f, %.3f, %.3f)', endPoint), ...
    'FontSize', 10, ...
    'FontWeight', 'bold', ...
    'Color', [0.8 0.2 0], ...
    'BackgroundColor', 'w', ...
    'Margin', 3);

% Desired reference label
text(endPoint(1), endPoint(2), endPoint(3) + 0.08, ...
    sprintf('r = %.3f', rDesired), ...
    'FontSize', 10, ...
    'FontWeight', 'bold', ...
    'Color', 'b', ...
    'BackgroundColor', 'w', ...
    'Margin', 3);
hold off


%% =========================================================
% 18. TIME EVOLUTION AND CONSTRAINT VERIFICATION
%% =========================================================

% Extract simulation data
t = time;

r = desiredReference;
w = augmentedHistory(1,:);
x1 = augmentedHistory(2,:);
x2 = augmentedHistory(3,:);
u = outputHistory(3,:);

% Create one figure with four subplots
figure('Color','w', ...
       'Name','RG - Time Evolution and Constraints', ...
       'NumberTitle','off');

tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

%% 1. Position and references
nexttile;
hold on; grid on; box on;

plot(t,r,'g--','LineWidth',1.5);
plot(t,w,'b','LineWidth',2);
plot(t,x1,'r','LineWidth',2);

yline(stateMax(1),'m--','LineWidth',1.2);
yline(stateMin(1),'m--','LineWidth',1.2);

ylabel('Position / Reference');
title('Position Tracking');

legend('Desired reference r', ...
       'Accepted reference w', ...
       'Position x_1', ...
       'Position constraints', ...
       'Location','best');

xlim([t(1) t(end)]);

%% 2. Velocity
nexttile;
hold on; grid on; box on;

plot(t,x2,'b','LineWidth',2);

yline(stateMax(2),'r--','LineWidth',1.3);
yline(stateMin(2),'r--','LineWidth',1.3);

ylabel('Velocity x_2');
title('Velocity Constraint');

legend('Velocity x_2', ...
       'Maximum velocity', ...
       'Minimum velocity', ...
       'Location','best');

xlim([t(1) t(end)]);

%% 3. Command
nexttile;
hold on; grid on; box on;

stairs(t,u,'b','LineWidth',1.8);

yline(inputMax,'r--','LineWidth',1.3);
yline(inputMin,'r--','LineWidth',1.3);

ylabel('Control input u');
title('Control Input Constraint');

legend('Command u', ...
       'Maximum input', ...
       'Minimum input', ...
       'Location','best');

xlim([t(1) t(end)]);

%% 4. Reference Governor parameter
nexttile;
hold on; grid on; box on;

stairs(t(1:end-1),kappaHistory,'b','LineWidth',1.8);

yline(1,'r--','LineWidth',1.2);
yline(0,'r--','LineWidth',1.2);

ylabel('\kappa(k)');
xlabel('Time [s]');
title('Reference Governor Parameter K');

legend('\K(k)', ...
       'Upper bound', ...
       'Lower bound', ...
       'Location','best');

ylim([-0.05 1.05]);
xlim([t(1) t(end)]);

sgtitle('Dynamic Reference Governor - Time Evolution');

%% =========================================================
% NUMERICAL CONSTRAINT VERIFICATION
%% =========================================================

tol = 1e-7;

% Maximum violations (positive values mean violation)
positionViolation = max([ ...
    max(x1-stateMax(1)), ...
    max(stateMin(1)-x1), 0]);

velocityViolation = max([ ...
    max(x2-stateMax(2)), ...
    max(stateMin(2)-x2), 0]);

inputViolation = max([ ...
    max(u-inputMax), ...
    max(inputMin-u), 0]);

fprintf('\n====================================\n');
fprintf('CONSTRAINT CHECK\n');
fprintf('====================================\n');

fprintf('Position constraint violation: %.3e\n', ...
        positionViolation);

fprintf('Velocity constraint violation: %.3e\n', ...
        velocityViolation);

fprintf('Input constraint violation: %.3e\n', ...
        inputViolation);

if max([positionViolation,velocityViolation,inputViolation]) <= tol
    fprintf('All physical constraints are satisfied.\n');
else
    fprintf('WARNING: Physical constraint violation detected.\n');
end