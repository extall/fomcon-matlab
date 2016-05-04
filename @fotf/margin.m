function [Gm,Pm,Wg,Wp] = margin(varargin)
%MARGIN Finds gain and phase margins and crossover frequencies.
%
%   Usage: [GM, PM, WG, WP] = MARGIN(G)
%          [GM, PM, WG, WP] = MARGIN(MAG, PH, W)
%
%          where GM, PM   - gain and phase margins,
%                WG, WP   - corresponding crossover frequencies,
%                G        - fractional-order transfer function,
%                MAG,PH,W - magnitude, phase and frequency vectors.


    if nargin == 1
        
        % Get margin for the system
        [mag, ph, w] = bode(varargin{1});
        [Gm,Pm,Wg,Wp] = margin(mag,ph,w);
        
    elseif nargin == 3
        
        % Get margin for magnitude, phase and w
        [Gm,Pm,Wg,Wp] = margin(varargin{1},varargin{2},varargin{3});
        
    else
        
        error('Wrong number of arguments.');
        
    end
        
end

