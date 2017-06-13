% -----------------------------------------
% Graph-Cuts for F-Formation (GCFF)
% 2015 - University of Verona
% Written by Francesco Setti
% -----------------------------------------
%
% This script is just an example of how to use the GCFF code to run
% experiments on provided data, in case you do not have groundtruths.
%


%% INITIALIZATION

% Clean the workspace
clear all,
close all,
clc,

% Set data folder
seqpath = '../data/CoffeeBreak/Seq1' ;
% NB: edit here your own path to data!!!


% Load features
load(fullfile(seqpath,'features.mat')) ;

%Load settings
load(fullfile(seqpath,'settings_gc.mat') ) ;


% % Load groundtruth
% load(fullfile(seqpath,'groundtruth.mat')) ;
%
%
% % If only some frames are annotated, delete all the others from features.
% [~,indFeat] = intersect(timestamp,GTtimestamp) ;
% timestamp = timestamp(indFeat) ;
% features  = features(indFeat) ;
%
% % Initialize evaluation variables
% TP = zeros(1,length(timestamp)) ;
% FP = zeros(1,length(timestamp)) ;
% FN = zeros(1,length(timestamp)) ;
% precision = zeros(1,length(timestamp)) ;
% recall = zeros(1,length(timestamp)) ;


for idxFrame = 1:length(timestamp)
    gg = gc( features{idxFrame}, stride, mdl ) ;
    groups{idxFrame} = [] ;
    for ii = 1:max(gg)+1
        groups{idxFrame}{ii} = features{idxFrame}(gg==ii-1,1) ;
    end

    if ~isempty(groups{idxFrame})
        groups{idxFrame} = ff_deletesingletons(groups{idxFrame}) ;
    end
%     if ~isempty(GTgroups{idxFrame})
%         GTgroups{idxFrame} = ff_deletesingletons(GTgroups{idxFrame}) ;
%     end
%     [precision(idxFrame),recall(idxFrame),TP(idxFrame),FP(idxFrame),FN(idxFrame)] = ff_evalgroups(groups{idxFrame},GTgroups{idxFrame},'card') ;

    % DISPLAY RESULTS
    % Frame:
    fprintf('Frame: %d/%d\n', idxFrame, length(timestamp))
    % Found:
    fprintf('   FOUND:-- ')
    if ~isempty(groups{idxFrame})
        for ii=1:size(groups{idxFrame},2)
            fprintf(' %i',groups{idxFrame}{ii})
            fprintf(' |')
        end
    else
        fprintf(' No Groups ')
    end
    fprintf('\n')
%     % GT:
%     fprintf('   GT   :-- ')
%     if ~isempty(GTgroups{idxFrame})
%         for ii=1: size(GTgroups{idxFrame},2)
%             fprintf(' %i',GTgroups{idxFrame}{ii})
%             fprintf(' |')
%         end
%     else
%         fprintf(' No Groups ')
%     end
%     fprintf('\n');

    figure(idxFrame),
    ff_plot_person_tv(features{idxFrame},'b',50), axis equal

    if ~isempty(groups{idxFrame})
        for ii=1:size(groups{idxFrame},2)
            for jj=1:size(groups{idxFrame}{ii},1)-1
                for kk=jj+1:size(groups{idxFrame}{ii},1)
                    xx = [features{idxFrame}(features{idxFrame}(:,1)==groups{idxFrame}{ii}(jj),2), features{idxFrame}(features{idxFrame}(:,1)==groups{idxFrame}{ii}(kk),2)] ;
                    yy = [features{idxFrame}(features{idxFrame}(:,1)==groups{idxFrame}{ii}(jj),3), features{idxFrame}(features{idxFrame}(:,1)==groups{idxFrame}{ii}(kk),3)] ;
                    line(xx,yy,'Color','g','LineWidth',2)
                end
            end
        end
    end
    pause

end

% pr = mean(precision) ;
% re = mean(recall) ;
%
% F1 = 2 * pr .* re ./ (pr+re) ;
%
%
% [~,indFeat] = intersect(timestamp,GTtimestamp) ;
% pr_avg = mean(precision(indFeat)) ;
% re_avg = mean(recall(indFeat)) ;
% F1_avg = 2 * pr_avg * re_avg / ( pr_avg + re_avg ) ;
% fprintf('Average Precision: -- %d\n',pr_avg)
% fprintf('Average Recall: -- %d\n',re_avg)
% fprintf('Average F1 score: -- %d\n',F1_avg)
