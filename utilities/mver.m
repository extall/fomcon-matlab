function mnum = mver()
%MVER MATLAB version number

    mversion = ver('matlab');
    mnum = str2num(mversion.Version);

end

