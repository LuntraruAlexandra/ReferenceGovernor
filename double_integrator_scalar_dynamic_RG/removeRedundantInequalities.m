%% =========================================================
% Inspired by [1], Algorithm 3.2.
% For each inequality a_i*z <= b_i, test whether the other
% inequalities imply it by solving:
% maximize a_i*z, subject to a_j*z <= b_j for j ~= i.
% If the maximum is <= b_i + tolerance, remove the inequality.
% Constraints are removed sequentially so that equivalent
% inequalities are not all deleted at the same time.
%% =========================================================

function [H,h,constraintType,predictionStep,constraintIndex] = ...
    removeRedundantInequalities(H,h,constraintType,predictionStep, ...
                                constraintIndex,lpOptions, ...
                                redundancyTolerance,zeroTolerance)

    rowIndex = 1;
    while rowIndex <= size(H,1)
        objectiveRow = H(rowIndex,:);
        upperBound = h(rowIndex);

        % Temporarily remove the inequality being tested.
        % Why? 
        % To test if it is redundant we come with the question:
        % Are constraints still guaranteed if it is no longer part of the
        % set?
        keepRows = true(size(H,1),1);
        keepRows(rowIndex) = false;
        Hremaining = H(keepRows,:);
        hRemaining = h(keepRows);

        % Handle inequalities with nearly zero coefficients.
        if norm(objectiveRow,inf) < zeroTolerance
            if upperBound >= -redundancyTolerance
                isRedundant = true;
            else
                error('An infeasible zero-coefficient inequality was found.');
            end
        else
            % max a*z = -min(-a*z)
            % here looking for the worst case possible with linprog
            
            % If the worst case is not achieved then adding back the
            % inequality would be useless
            [~,objectiveMin,exitFlag] = linprog( ...
                -objectiveRow',Hremaining,hRemaining, ...
                [],[],[],[],lpOptions);

            if exitFlag == 1
                maximumValue = -objectiveMin;
                isRedundant = maximumValue <= ...
                              upperBound+redundancyTolerance;
            elseif exitFlag == -3
                % Without this inequality, the LP is unbounded.
                % Therefore, the inequality must be retained.
                isRedundant = false;
            else
                error('LP failed while reducing inequality %d (flag=%d).', ...
                      rowIndex,exitFlag);
            end
        end

        if isRedundant
            H(rowIndex,:) = [];
            h(rowIndex) = [];
            constraintType(rowIndex) = [];
            predictionStep(rowIndex) = [];
            constraintIndex(rowIndex) = [];
        else
            rowIndex = rowIndex+1;
        end
    end
end

