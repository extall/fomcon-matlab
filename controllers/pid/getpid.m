function [Kp, Ki, lam, Kd, mu] = getpid(pid)
%GETPID Recover fractional or integer-order PID gains from a LTI model
%   
%   Usage: [Kp, Ki, lam, Kd, mu] = GETPID(MODEL), where MODEL may be an LTI 
%                                  system given either by TF, ZPK or FOTF.
%
%   See also: fracpid
	
	% Set every parameter to zero initially
	Kp = 0; Ki=0; Kd=0; lam=0; mu=0;
	
	% Determine controller type and return the corresponding parameters
	[a na b nb] = fotfparam(fotf(pid));
	
	if 		length(a) == 1 && length(na) == 1 && ...			% P controller
			length(b) == 1 && length(nb) == 1 && ...
			fleq(a,1) && fleq(na,0) && fleq(nb, 0)
				Kp = b;								
	elseif	length(a) == 1 && length(na) == 1 && ...			% PI controller
			length(b) == 2 && length(nb) == 2 && ...
			fleq(a,1) && ~fleq(na,0) && fleq(nb(2),0)
				Kp = b(1); Ki = b(2); lam = na;		
	elseif 	length(a) == 1 && length(na) == 1 && ...			% PD controller
			length(b) == 2 && length(nb) == 2 && ...
			fleq(a,1) && fleq(na,0) && fleq(nb(2),0)
				Kp = b(2); Kd = b(1); mu = nb(1);	
	elseif	length(a) == 1 && length(na) == 1 && ...			% PID controller
			length(b) == 3 && length(nb) == 3 && ...
			fleq(a,1) && fleq(nb(3),0)
				Kp = b(2); Ki = b(3); Kd = b(1);				
				lam = na; mu = nb(1)-lam;
	end
	
end

