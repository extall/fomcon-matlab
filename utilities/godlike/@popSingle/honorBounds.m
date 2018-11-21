% FIXME: (Rody Oldenhuis) once we've switched to using objectiveFunctions(), we 
% can simply remove this whole method

function [newpop, newfit] = honorBounds(pop, newpop, newfit)

    % find violation sites
    outsiders1 = false; outsiders2 = false;
    if ~isempty(newpop)
        outsiders1 = newpop < pop.lb;
        outsiders2 = newpop > pop.ub;
    end

    % PSO requires more elaborate check
    if strcmpi(pop.algorithm, 'PSO')

        % rename for clarity
        velocity = pop.pop_data.velocities; % extract velocities
        velUb = (pop.ub - pop.lb)/5;        % upper bounds on velocity
        velLb = (pop.ub - pop.lb)/1e50;     % lower bounds on velocity

        % bounce against bounds
        if any(outsiders1(:) | outsiders2(:))
            newpop(outsiders1)   = newpop(outsiders1) - velocity(outsiders1);
            newpop(outsiders2)   = newpop(outsiders2) - velocity(outsiders2);
            velocity(outsiders1) = -velocity(outsiders1);
            velocity(outsiders2) = -velocity(outsiders2);
        end

        % limit velocity
        Velsign    = sign(velocity);
        outsiders1 = abs(velocity) > abs(velUb);
        outsiders2 = abs(velocity) < abs(velLb);
        if any(outsiders1(:)) || any(outsiders2(:))
            velocity(outsiders1) = Velsign(outsiders1).*velUb(outsiders1);
            velocity(outsiders2) = Velsign(outsiders2).*velLb(outsiders2);
        end

        % re-insert velocity
        pop.pop_data.velocities = velocity;

    % boundary violations in all other algorithms
    % are simply reinitialized
    else
        reinit = pop.lb + rand(pop.size, pop.dimensions).*(pop.ub-pop.lb);
        if any(outsiders1(:) | outsiders2(:))
            newpop(outsiders1) = reinit(outsiders1);
            newpop(outsiders2) = reinit(outsiders2);
            % also remove any function values
            newfit(any(outsiders1,2), :) = NaN;
            newfit(any(outsiders2,2), :) = NaN;
        end
    end

end
