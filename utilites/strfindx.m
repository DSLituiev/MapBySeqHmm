function out = strfindx(str, pattern)

cout = strfind(str, pattern);

out = ~all(cellfun(@isempty,cout));
