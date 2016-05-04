function display(G)
% DISPLAY  Display the fractional-order transfer function in the command line.

	% Load configuration parameters
    config = fomcon('config');
    numSigDig = config.Core.General.Model_significant_digits;

    strN=polydisp(G.b,G.nb);
    strD=polydisp(G.a,G.na);
    nn=length(strN);
    nd=length(strD);
    nm=max([nn,nd]);
    
    n1=(nm-nn)/2;
    n2=(nm-nd)/2;
    
    T = G.ioDelay;
    
    if T>0
        delay_s=[' exp(-' num2str(T,numSigDig) 's)'];
    else
        delay_s='';
    end
    
    num = [char(' '*ones(1,floor(n1))) strN];
    lin = [char('-'*ones(1,nm)) delay_s];
    den = [char(' '*ones(1,floor(n2))) strD];
    
    if isempty(G.a) && isempty(G.na) && isempty(G.b) && isempty(G.nb)
        
        disp('Empty fractional-order transfer function.');
        
    else
    
        disp('Fractional-order transfer function:');

        if ~(length(G.a)==1 && fleq(G.a,1) && G.na==0)
            disp(num);
            disp(lin);
            disp(den);
        else
            disp([num delay_s]);
        end
        
    end
    
end

function strP=polydisp(p,np)

    % Load configuration parameters
    config = fomcon('config');
    numSigDig = config.Core.General.Model_significant_digits;
    
    P='';
    [np,ii]=sort(np,'descend');
    p=p(ii);

    for i=1:length(p),
        P=[P,'+',num2str(p(i),numSigDig), 's^{',num2str(np(i),numSigDig),'}'];
    end

    P=P(2:end);
    P=strrep(P,'s^{0}','');
    P=strrep(P,'+-','-');
    P=strrep(P,'^{1}','');
    P=strrep(P,'+1s','+s');
    strP=strrep(P,'-1s','-s');

    if length(strP)>=2 && strcmp(strP(1:2),'1s')
        strP=strP(2:end);
    end

end