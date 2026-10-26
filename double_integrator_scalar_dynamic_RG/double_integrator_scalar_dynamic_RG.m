clear;
clc;
close all;

%%
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
%% 

%% 1. DISCRETE-TIME DOUBLE-INTEGRATOR MODEL 


% Model used in the paper's example [3]:
% x1(k+1) = x1(k) + 0.1*x2(k)
% x2(k+1) = x2(k) + 0.1*u(k)
Ts = 0.1;


Ad = [1 0.1;
      0 1];
Bd = [0;
      0.1];
[n,m]= size(Bd);


%% 2. CONTROLLER 

% u(k) = 0.9170*(v(k)-x1(k)) - 1.6821*x2(k)--> this is from the paper
% Equivalently: u(k) = -K*x(k) + K(1)*v(k)

% Tests were also done with this controller:
% Kd = [0.9170 1.6821]; --> from the paper

% Closed-loop dynamics: x(k+1) = Ah*x(k) + Bh*v(k)

Q = 10 * eye(n)
R = 2 * eye(m);


Kd = dlqr(Ad,Bd,Q,R)
Ah = Ad - Bd*Kd;
Bh = Bd*Kd(1);

% Dimensions inferred from the closed-loop model.
% n: plant-state dimension; m: admissible-reference dimension.
n = size(Ah,1);
m = size(Bh,2);

eigenvalues = eig(Ah);
disp('LQR controller gain:'); disp(Kd);
disp('Matrix Ah:'); disp(Ah);
disp('Matrix Bh:'); disp(Bh);
disp('Closed-loop eigenvalues:'); disp(eigenvalues);

if any(abs(eigenvalues) >= 1)
    error('The closed-loop system is not Schur stable.');
end

%% 3. SYSTEM CONSTRAINTS

% |x1| <= 1, |x2| <= 0.1, |u| <= 0.1
stateMin = [-1; -0.1];
stateMax = [ 1;  0.1];
inputMin = -0.1;
inputMax =  0.1;

% Constrained output: y = [x1; x2; u] = Ch*x + Dh*v

% u = -Kd * x + Kd(1)*v

% y was chosen to contain all constraints
Ch = [eye(n);
      -Kd];
Dh = [zeros(n,m);
      Kd(1)];

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



%% 7. INITIALIZE O_0* - REFERENCES [1], [2]

% H0*v belongs to Y(epsilon): F*H0*v <= g-epsilon
% See [2], equations (21)-(22).
% O_0* = {z : H0*v belongs to Y(epsilon), CG*z belongs to Y}

% f(CG*AGt *xG) = F*CG *(AGt* xG) - g <= 0
% f([H0 0]xG) --> F*H0 *xG  - g <= -eps
H = [F*H0 zeros(q,n);
     F*CG];
h = [g - epsilon;
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
[Oinf,vMin,vMax,numSlices,referenceValues,slices,colorMap] = ...
    computeOinfVisualizationGeometry(Oinf,H,h,n,m,lpOptions);

%% 17. DYNAMIC SCALAR REFERENCE GOVERNOR - [2], EQUATIONS (10)-(11)

% Augmented state xG(k) = [w(k); x1(k); x2(k)].
% w(k) is the previously accepted reference; r(k) is the desired reference.
% xG(k+1) = AG*xG(k) + BG*K(k)*(r(k)-w(k)).
% y(k) = CG*xG(k) must satisfy F*y(k) <= g.
% Choose the largest K in [0,1] such that xG(k+1) belongs to O_inf*.
% The plant evolves using w(k), not the newly computed w(k+1).

BG = [eye(m); zeros(n,m)];
Oinf.BG = BG;

simulationDuration = 40;                 % seconds
simulationSteps = round(simulationDuration/Ts);
time = (0:simulationSteps)*Ts;

% Desired position reference: zero initially, then a step to 0.5
% Change 0.5 to 1.2 to illustrate an inadmissible steady-state target
desiredReference = zeros(1,simulationSteps+1);
desiredReference(time >= 1) = 0.5;

% Initial augmented state must belong to O_inf*.
xG = zeros(n+1,1);

feasibilityTolerance = 1e-7;

if any(H*xG > h+feasibilityTolerance)
    error('Initial augmented state is outside O_inf*.');
end

augmentedHistory = zeros(n+m,simulationSteps+1);
outputHistory = zeros(p,simulationSteps+1);
KHistory = zeros(1,simulationSteps);
nextStateCheck = zeros(1,simulationSteps);
augmentedHistory(:,1) = xG;

for step = 1:simulationSteps
    r = desiredReference(step);
    w = [eye(m) zeros(m,n)] * xG;

    % This is [I 0]*xG actually for this case

    % In the article's timing convention, y(k) uses the current w(k).
    outputHistory(:,step) = CG*xG;

    % xG_next = AG*xG + BG*K*(r-w)

    % Enforce H*(AG*xG + BG*K*(r-w)) <= h.
    % Each row yields: coefficient_i*K <= margin_i.
    coefficients = H*BG*(r-w);
    margins = h-H*AG*xG;
    KLower = 0;
    KUpper = 1; 

    for rowIndex = 1:size(H,1)
        coefficient = coefficients(rowIndex);
        margin = margins(rowIndex);
        
        if coefficient > 1e-12
            KUpper = min(KUpper,margin/coefficient);
        elseif coefficient < -1e-12
            KLower = max(KLower,margin/coefficient);
        elseif margin < -feasibilityTolerance
            error('RG infeasible at step %d: zero coefficient.',step);
        end
    end

    if KLower > KUpper+feasibilityTolerance
        error('No feasible K in [0,1] at step %d.',step);
    end

    % Largest feasible K; keeping it to interval --> [0,1].
    K = max(0,min(1,KUpper));
    if K < KLower-feasibilityTolerance
        error('Numerically infeasible K at step %d.',step);
    end

    xGnext = AG*xG+BG*K*(r-w);
    nextStateCheck(step) = max(H*xGnext-h);
    if nextStateCheck(step) > feasibilityTolerance
        error('Next augmented state violates O_inf* at step %d.',step);
    end

    KHistory(step) = K;
    xG = xGnext;
    augmentedHistory(:,step+1) = xG;
end
outputHistory(:,end) = CG*xG;



Oinf.rg.time = time;
Oinf.rg.desiredReference = desiredReference;
Oinf.rg.augmentedHistory = augmentedHistory;
Oinf.rg.outputHistory = outputHistory;
Oinf.rg.KHistory = KHistory;

%% 18. ONLY FIGURE: RG TRAJECTORY AND O_INF* CROSS-SECTIONS
% Plot the RG trajectory together with the constant-reference cross-sections.
% Mark with X the exact intersections with each plotted plane w = constant.
% Intersections are linearly interpolated between consecutive RG samples.

Oinf = plotRGTrajectoryAndSlices(Oinf,slices,numSlices,referenceValues, ...
    colorMap,augmentedHistory,simulationSteps,time,vMin,vMax,desiredReference);

%%
% 18. TIME EVOLUTION AND CONSTRAINT VERIFICATION

% Extract simulation data
t = time;

r = desiredReference;
w = augmentedHistory(1,:);
x1 = augmentedHistory(2,:);
x2 = augmentedHistory(3,:);
u = outputHistory(3,:);

plotRGTimeEvolution(t,r,w,x1,x2,u,stateMax,stateMin, ...
    inputMax,inputMin,KHistory);
