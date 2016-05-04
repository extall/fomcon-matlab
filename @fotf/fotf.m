classdef fotf
    %FOTF Creates a new fractional-order transfer function object
    %
    %Usage G = FOTF(A, NA, B, NB, T)
    %      G = FOTF('s')
    %      G = FOTF(BPOLY, APOLY, T)
    %      G = FOTF(H)
    %
    %where G,H          - FOTF objects,
    %      A, NA, B, NB - vectors containing pole polynomical coefficients and
    %                     exponents and zero polynomial coefficients and
    %                     exponents respectively.
    %      BPOLY, APOLY - polynomial strings, e.g. 's^0.5+1e3', '1-5s^0.75'.
    %      T            - input/output delay [sec]   

    properties
       
        a           % Pole polynomial coefficients
        na          % Pole polynomial exponents
        b           % Zero polynomial coefficients
        nb          % Zero polynomial exponents
        ioDelay     % Input-output delay [sec]
        
    end
    
    methods
    
        % Constructor
        function G=fotf(a,na,b,nb,T)

            % Load configuration
            config = fomcon('config');
            epsi   = config.Core.General.Internal_computation_accuracy;
            
            % Empty object
            if nargin==0,

                G.a=[];
                G.na=[];
                G.b=[];
                G.nb=[];
                G.ioDelay=0;

            % FOTF object given
            elseif isa(a,'fotf')

                G=a;

            % Create FOTF from number
            elseif nargin==1 && isa(a,'double')

                G=fotf(1,0,a,0,0);

            % S term
            elseif nargin==1 && isa(a, 'char') && a == 's'

                G = fotf(1,0,1,1,0);
                
            % TF or ZPK object
            elseif nargin==1 && (isa(a,'tf') || isa(a,'zpk'))
                
                % Convert an integer-order transfer function to a fractional-order
                % transfer function
                G1 = a;
                
                [n,d]=tfdata(G1,'v');
                i1=find(abs(n)<epsi);
                i2=find(abs(d)<epsi);

                if length(i1)>0 && i1(1)==1
                    n=n(i1(1)+1:end);
                end

                if length(i2)>0 && i2(1)==1
                    d=d(i2(1)+1:end);
                end
                
                G=fotf(d,length(d)-1:-1:0,n,length(n)-1:-1:0);

                ioDelay = get(G1, 'ioDelay');

                if ~isempty(ioDelay) && ioDelay > 0
                    G.ioDelay = ioDelay;
                end
                
            % FOTF given by string
            elseif nargin>=2 && isa(a,'char') && isa(na,'char')
                
                G = newfotf(a,na);
                
                % Delay term
                if nargin==3 && isa(b,'double')
                    G.ioDelay = abs(b);
                end

            % FOTF given by coefficients
            else

                ii=find(abs(a)<epsi);
                a(ii)=[];
                na(ii)=[];
                ii=find(abs(b)<epsi);
                b(ii)=[];
                nb(ii)=[];

                % ioDelay
                if nargin == 5
                    G.ioDelay = abs(T);
                else
                    G.ioDelay = 0;
                end

                G.a=a;
                G.na=na;
                G.b=b;
                G.nb=nb;

            end

        end
    
    end

end