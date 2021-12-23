function assignhere(varname, vardata)
%ASSIGNHERE Assign variable in the caller workspace ("here").
assignin('caller', varname, vardata);
end

