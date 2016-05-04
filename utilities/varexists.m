function [var_exists, isfotf, thisclass]=varexists(varName)
% VAREXISTS Determine the type of variable in the workspace.

    if nargin == 0
        var_exists = false;
    else
        allVars  = evalin('base', 'who');
    
        var_exists = max(strcmp(allVars, varName));

        % Check output
        if isempty(var_exists)
            var_exists=false;
        end

        % Check class
        isfotf = false;

        if var_exists

            % Get this class name
            thisclass = evalin('base',strcat('class(',varName,')'));

            % Is it a FOTF object?
            if strcmp(thisclass, 'fotf')
                isfotf = true;
            end

        else

            thisclass='';

        end
    end
end