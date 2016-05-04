function varargout = bode(fdata)
%BODE Plot the experimentally observed frequency-domain data.
    
    varargout = {};

    % Convert to FRD
    mag = db2mag(fdata.mag);
    ph = deg2rad(fdata.phase);
    h = frd(mag.*exp(sqrt(-1)*ph), fdata.w);
    
    % Compute response
    if nargout == 0
        bode(h);
    else
        [varargout{1}, varargout{2}] = bode(h);
    end

end

