function [translationTable] = id2clust(BPGAtable,targetName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
hlp = BPGAtable(:, {'cluster', targetName});
hlp = hlp(~cellfun(@isempty,hlp{:,targetName}), :);

splitTarget = cellfun(@(x) strsplit(x, ';'), hlp{:,targetName}, 'UniformOutput', false);
clusterRepeats = repelem(hlp.cluster, cellfun(@numel, splitTarget));

targetList = [];
clusterList = [];
for i = 1:height(hlp)
    currentClusters = repmat(hlp.cluster(i), numel(splitTarget{i}), 1);
    targetList = [targetList; splitTarget{i}(:)];
    clusterList = [clusterList; currentClusters];
end

translationTable = table(targetList, clusterList, 'VariableNames', {targetName, 'cluster'});
translationTable = table2cell(translationTable)
translationTable;
end