classdef pbar < handle
%PBAR Progress bar with remaining time estimation%
%
% Usage: H = PBAR(MSG) creates a progress bar with a specified message,
%        UPDATE(H, N, NTOTAL, MSG) will update the progress bar, where N is
%        the current iteration and NTOTAL is the total number of iterations
%        to be performed.
%        DELETE(H) will remove the progress bar.
    
    properties (GetAccess='private', SetAccess='private')
        
        % Waitbar data
        msg                 % Message
        wbar                % Waitbar handle
        
        % Timers/counters
        iterCnt 		    % Counter
        beginTimer          % Start time (tic)
        
        updatesAtIter       % Array with iterations such that if
                            % iter > updatesAtIter(N) then the local
                            % iteration number will be increased
                            
        nextUpdate          % Next iteration number to update at
        
        % Iteration time calculation
        iterMemory          % Array in form [iterTime; iterNum]
        
		% Time data
		elapsedTime			% Elapsed time in seconds
		remainingTime		% Remaining time in seconds
        
        estimationMessage   % Estimated time message
                
    end
    
    methods
        
        % Initializer
        function pBar = pbar(msg)
           
            % Default message
            if nargin == 0
                msg = 'Please wait...';              
            end
            
            % Create a waitbar
            pBar.wbar = waitbar(0, msg, 'Name', 'Progress');
            
            % Initial parameters
			pBar.msg = msg;
			
			pBar.iterCnt 		= 0;
            pBar.iterMemory     = [];
			pBar.beginTimer 	= [];
			
			% Time data
			pBar.elapsedTime	= 0;
			pBar.remainingTime  = [];
            
            pBar.estimationMessage = 'estimating...';
            
        end
        
        % Destructor
        function delete(h)
            
            % Close the progress bar figure
            close(h.wbar);
                
        end
		
		function update(pBar, iter, totalIter, msg)
		%UPDATE(PBAR, ITER, TOTALITER, MSG) Updates PROGBAR display
		%	PBAR		- progress bar object
		%	ITER		- current iteration count
		%	TOTALITER 	- total iterations
		%	MSG			- message
    
            % On first call, create all necessary data structures
            if pBar.iterCnt == 0
               
                pBar.beginTimer = tic;
                pBar.nextUpdate = totalIter / 20; % 20 Updates, i.e. progress
                                                  % bar resolution
                                             
                pBar.updatesAtIter   = pBar.nextUpdate : ...     % Updates array
                                       pBar.nextUpdate : ...
                                       totalIter;
                
                pBar.iterCnt = pBar.iterCnt + 1;
                
                pBar.iterMemory(1, 1) = ...
                            toc(pBar.beginTimer);
                        
                        pBar.iterMemory(1, 2) = ...
                            iter;
       
            else
            
                % Otherwise "do the math"
                if iter > pBar.nextUpdate
                    
                    % Set new update margin
                    pBar.nextUpdate = pBar.updatesAtIter(pBar.iterCnt+1);
                    
                    % Given 6 points, approximate the curve
                    if pBar.iterCnt == 6
                        
                        % Second degree polynomial approximation
                        [poly2,S,mu] =  polyfit(pBar.iterMemory(:,2), ...
                                        pBar.iterMemory(:,1), ...
                                        2);

                        pBar.remainingTime = polyval(poly2, totalIter,S,mu);
                            
                        % Remaining time estimation message
                        pBar.estimationMessage = stohms(pBar.remainingTime);
                        
                    else

                        % Add new point: (time, iteration)
                        pBar.iterMemory(size(pBar.iterMemory,1)+1,:) = ...
                            [toc(pBar.beginTimer), iter];
                        
                    end
                    
                    % Set elapsed time
					pBar.elapsedTime = toc(pBar.beginTimer);
                    
                    % Check elapsed vs. estimated
                    if ~isempty(pBar.remainingTime) && ...
                        pBar.elapsedTime >= pBar.remainingTime
                        
                        % Second degree polynomial approximation
                        [poly2,S,mu] =  polyfit(pBar.iterMemory(:,2), ...
                                        pBar.iterMemory(:,1), ...
                                        2);

                        pBar.remainingTime = polyval(poly2, totalIter,S,mu);
                            
                        % Remaining time estimation message
                        pBar.estimationMessage = stohms(pBar.remainingTime);
                        
                    end
                    
                    % Change message if required
                    if nargin == 4
                        pBar.msg = msg;
                    end
      
                    % Get text
					estimatedLeft = [' (' stohms(pBar.elapsedTime) ' / ' pBar.estimationMessage ')'];
                
                    % Update the waitbar
                    waitbar(iter/totalIter, pBar.wbar, [pBar.msg estimatedLeft]);

                    pBar.iterCnt = pBar.iterCnt + 1;
                    
                end
            end
        end
    end
end

