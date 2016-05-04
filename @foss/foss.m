classdef(CaseInsensitiveProperties = true) foss
    %FOSS Fractional-order state-space system
    %
    % Usage: S = FOSS(A,B,C,D,q) - creates a fractional-order state-space
    % system with matrices A, B, C, and D such that
    %
    %    q
    %   d  x(t) = A(t)x(t) + B(t)u(t)
    %      y(t) = C(t)x(t) + D(t)u(t).
    
    properties
        A       % Matrix A
        B       % Matrix B
        C       % Matrix C
        D       % Matrix D
        q       % Commensurate order
    end
    
    properties(SetAccess = private, Hidden = true)
        
        fotf_array   % Corresponding array of FO transfer functions
        
    end
    
    methods
        
        % Constructor
        function fss = foss(varargin)
            
            if nargin < 1
                error('FOSS:NotEnoughInputParameters', ...
                    'Not enough input arguments.');
            end
            
            if isa(varargin{1},'ss') && nargin == 1
                mySs  = varargin{1};
                fss.A = mySs.a;
                fss.B = mySs.b;
                fss.C = mySs.c;
                fss.D = mySs.d;
                fss.q = 1;
            elseif nargin == 5
                fss.A = varargin{1};
                fss.B = varargin{2};
                fss.C = varargin{3};
                fss.D = varargin{4};
                fss.q = varargin{5};
            else
                error('FOSS:WrongInputParameters', ...
                    'Wrong type or number of input parameters.');
            end
            
            % Create and populate the tf_array
            
            % Correct size
            fss.fotf_array = ...
                foss.build_fotf_array(fss.A,fss.B,fss.C,fss.D,fss.q);
            
        end
        
        % Redefine subsasgn
        function fss = subsasgn(fss, s, val)
            
            % Use the built-in function for assignment
            fss = builtin('subsasgn', fss, s, val);
            
            % Regenerate the internal fotf array on variable reassignment
            fss.fotf_array = ...
                foss.build_fotf_array(fss.A,fss.B,fss.C,fss.D,fss.q);
            
        end
        
        % Display
        function display(fss)
            
            disp('Fractional-order state-space system');
            disp('A =');
            disp(fss.A);
            disp('B =');
            disp(fss.B);
            disp('C =');
            disp(fss.C);
            disp('D =');
            disp(fss.D);
            disp('(comm. order) q=');
            disp(fss.q);
            
        end
        
    end
    
    methods (Static)
        
        % Generate the FOTF array from state space description in terms
        % of matrices A, B, C, and D, and commensurate order q
        function f_arr = build_fotf_array(A,B,C,D,q)
            f_arr = cell(size(C,1),size(B,2));
            
            % Use Control Systems toolbox conversion from
            % state space to transfer function form
            myFotf = tf(ss(A,B,C,D));
            for k=1:size(myFotf,1)
                for l=1:size(myFotf,2)
                    nowFotf = fotf(myFotf(k,l));
                    [a,na,b,nb,L] = fotfparam(nowFotf);
                    na=na*q; nb=nb*q;
                    myNewFotf = fotf(a,na,b,nb,L);
                    % Normalize the fotf
                    f_arr{k,l} = normalize(myNewFotf);
                end
            end
        end
        
    end
    
end

