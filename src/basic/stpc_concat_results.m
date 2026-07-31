function results = stpc_concat_results(primaryResults, comparisonResults)
%STPC_CONCAT_RESULTS Concatenate result structs with compatible top fields.
%
% STPC products contain selection/timing metadata that smoothing comparison
% products do not. MATLAB requires identical top-level field sets for
% struct concatenation, so add absent fields as [] while preserving the
% historical flat struct-array return interface.

if isempty(primaryResults)
    results = comparisonResults;
    return
end
if isempty(comparisonResults)
    results = primaryResults;
    return
end

allFields = unique([fieldnames(primaryResults); ...
    fieldnames(comparisonResults)], 'stable');
primaryResults = local_add_fields(primaryResults, allFields);
comparisonResults = local_add_fields(comparisonResults, allFields);
primaryResults = orderfields(primaryResults, allFields);
comparisonResults = orderfields(comparisonResults, allFields);
results = [primaryResults, comparisonResults];
end

function values = local_add_fields(values, allFields)
existing = fieldnames(values);
missing = setdiff(allFields, existing, 'stable');
for i = 1:numel(missing)
    [values.(missing{i})] = deal([]);
end
end
