function [array_tf] = oustapp(S, varargin)
%OUSTAPP Obtain an integer-order approximation of a fractional-order system.
%
% Usage: ARRAY_TF = OUSTAPP(G, wb, wh, N, method) 
%        
%        where G - fotf object,
%              wb, wh - lower and higher frequency bounds
%              N - approximation order
%              method - 'oust' for Oustaloup's method (default) or 
%                       'ref' for the refined Oustaloup method
%
%              (defaults: wb = 0.001, wh = 1000, N=5, method='oust')
%
%        See also: fsparam

	% Create the array
    array_tf = zpk;
	
    noOuts = size(S.fotf_array,1);
    noIns  = size(S.fotf_array,2);
    
	% Go through all of the FO TFs in the FOSS array
	for k=1:noOuts
		for l=1:noIns
			array_tf(k,l) = oustapp(getfotf(S,k,l),varargin{:});
		end
	end
	
end

