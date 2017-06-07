% -----------------------------------------
% Graph-Cuts for F-Formation (GCFF)
% 2015 - University of Verona
% Written by Francesco Setti
% -----------------------------------------
%
% EVALGROUPS computes precision and recall scores, as well as the number of
% true positives, false positives and false negatives, given a set of
% detected groups and a set of ground truth groups.
%
% INPUT:
% ======
%  group      := <G-elements> cell array, each element contains a group
%                detected by tested algorithm, each one defined by the
%                array of subjects' ID of individuals belonging to it.
%  GT         := <T-elements> cell array, each element contains a group
%                provided by ground truth, each one defined by the array of
%                subjects' ID of individuals belonging to it.
%  TH         := scalar, it defines the threshold used to determine wheter
%                a detected group matches the GT group, or not. Older
%                versions of ff_evalgroups are contained in this function
%                and the relation is:
%                  - TH=2/3  <=>  cardinality criteria
%                  - TH=1    <=>  complete overlap criteria
%  cardmode   := this parameter can either be 0/1 or ~/'cardmode'. If it is
%                0 or does not exist, the computation will not take into
%                account the different group cardinalities, otherwise it
%                will be done for each cardinality separately.
%
% OUTPUT:
% =======
%  precision  := the ability to detect 'only' the ground truth (or the
%                ability not to generate false positives).
%                   Pr = TP / (TP+FP)
%  recall     := the ability to detect 'all' the ground truth (or the
%                ability not to generate false negatives).
%                   Re = TP / (TP+FN)
%  TP         := number of True Positives
%  FP         := number of False Positives
%  FN         := number of False Negatives
%
% NB: if cardmode is set OFF, all the output are scalars; otherwise they
% are all vectors with different sizes depending on the maximum cardinality
% of group and GT inputs.
%

function [precision,recall,TP,FP,FN] = ff_evalgroups(group,GT,TH,cardmode)


%% INITIALIZATION
% Default values definition and input data manipulation

% DEFAULT: TH = 2/3
if ~exist('TH','var')
    TH = 2/3 ;
elseif ischar(TH)
    switch TH
        case 'card'
            TH = 2/3 ;
        case 'all'
            TH = 1 ;
        otherwise
            error('TH input parameter cannot be set. It can be either a double between 0 and 1 or a string ''card'' or ''all''.')
    end
elseif isnumeric(TH) && (TH<0 || TH>1)
    error('TH value out of bound. It can be either a double between 0 and 1 or a string ''card'' or ''all''.')
end

% DEFAULT: cardmode = 0
if ~exist('cardmode','var') || strcmp(cardmode,'') || cardmode==0
    cardmode = 0 ;
elseif cardmode==1 || strcmp(cardmode,'cardmode')
    cardmode = 1 ;
else
    error('CARDMODE input parameter cannot be set. It can be either a boolean (0/1) or a string ''cardmode''.')
end


% % CHECK: group
% if isempty(group)
%     warning('group_empty','Input variable GROUPS is empty! Precision should be NaN!')
% end
%
% % CHECK: GT
% if isempty(GT)
%     warning('GT_empty','Input variable GT is empty! Recall should be NaN!')
% end


%% PROCESSING

switch cardmode

    case 0  % CARDMODE: OFF

        % Degenerate cases:
        if isempty(group) && isempty(GT)
            TP = 0 ;
            FP = 0 ;
            FN = 0 ;
            precision = 1 ;
            recall = 1 ;
            return ;
        elseif isempty(group) && ~isempty(GT)
            TP = 0 ;
            FP = 0 ;
            FN = numel(GT) ;
            precision = 1 ;
            recall = 0 ;
            return ;
        elseif ~isempty(group) && isempty(GT)
            TP = 0 ;
            FP = numel(group) ;
            FN = 0 ;
            precision = 0 ;
            recall = 1 ;
            return ;
        end

        % Initialize the 'True Positive'
        TP = 0 ;

        % For each GT group
        for ii = 1:length(GT)

            % Select the GT group and its cardinality (= number of elements)
            GTtmp  = GT{ii} ;
            GTcard = length(GTtmp) ;

            % for each detected group
            for jj = 1:size(group,2)

                % Select the detected group and its cardinality (= number of elements)
                grouptmp  = group{jj} ;
                groupcard = length(grouptmp) ;

                % Find the intersection between the GT and detected groups
                % and its cardinality.
                inters     = intersect(GTtmp,grouptmp) ;
                interscard = length(inters) ;

                if ge(interscard/max(GTcard,groupcard),TH)
                    TP = TP + 1 ;
                    break
                end
            end
        end

        % Define the number of false positives and negatives.
        FP = size(group,2) - TP ;
        FN = size(GT,2) - TP ;

        % Compute precision and recall.
        precision = TP ./ (TP+FP) ;
        recall = TP ./ (TP+FN) ;

    case 1  % CARDMODE: ON

        % Detect the number of groups in detected and ground truth.
        for ii = 1:length(GT)
            c_GT(ii) = numel(GT{ii}) ;
        end
        for jj = 1:length(group)
            c_group(jj) = numel(group{jj}) ;
        end

        if ~isempty(GT)
            GTcardmax    = max(c_GT) ;
            h_GT = hist(c_GT,1:GTcardmax) ;
        end
        if ~isempty(group)
            groupcardmax = max(c_group) ;
            h_group = hist(c_group,1:groupcardmax) ;
        end

        if isempty(group) && isempty(GT)        % GT = G = []
            precision = NaN ;
            recall    = NaN ;
            TP        = NaN ;
            FP        = NaN ;
            FN        = NaN ;

        elseif isempty(group) && ~isempty(GT)   % G = []
            TP        = zeros(GTcardmax,1) ;
            FP        = 0 ;
            FN        = h_GT ;

        elseif ~isempty(group) && isempty(GT)   % GT = []
            TP        = zeros(1,groupcardmax) ;
            FP        = h_group ;
            FN        = 0 ;

        else                                    % no []

            % Initialize the 'True Positives'
            TP = zeros(GTcardmax,groupcardmax) ;

            % For each GT group
            for ii = 1:length(GT)

                % Select the GT group and its cardinality (= number of elements)
                GTtmp  = GT{ii} ;
                GTcard = numel(GTtmp) ;
                if GTcard > GTcardmax
                    GTcardmax = GTcard ;
                end

                % for each detected group
                for jj = 1:length(group)

                    % Select the detected group and its cardinality (= number of elements)
                    grouptmp  = group{jj} ;
                    groupcard = numel(grouptmp) ;
                    if groupcard > groupcardmax
                        groupcardmax = groupcard ;
                    end

                    % Find the intersection between the GT and detected groups
                    % and its cardinality.
                    inters     = intersect(GTtmp,grouptmp) ;
                    interscard = numel(inters) ;
                    if ge(interscard/max(GTcard,groupcard),TH)
                        TP(GTcard,groupcard) = TP(GTcard,groupcard) + 1 ;
                        break
                    end

                end

            end

            % Define the number of false positives and negatives.
            FP = h_group - sum(TP,1) ;
            FN = h_GT - sum(TP,2)' ;

        end

        % Define precision.
        precision = sum(TP,1) ./ (sum(TP,1)+FP) ;

        % Define recall.
        recall = sum(TP,2)' ./ (sum(TP,2)'+FN) ;

end
