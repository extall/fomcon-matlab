function varargout = testsim(id, G, method, perc, N)
%TESTSIM Do OAT or Monte-Carlo time domain simulations for identified FOTF
%   The function allows to do either One-at-a-Time (OAT) parameter
%   simulation or Monte-Carlo simulation (when all parameters are
%   perturbed) in the time domain. The parameter uncertainty is supplied by
%   the user.
%
% Usage: [D, Y] = TESTSIM(ID, G, PERC, N), where
%   Input arguments: ID - fidata identification dataset.
%                    G  - either a FOTF model or FSPARAM structure with the
%                         identified model.
%                    METHOD - the method to be used in simulations. May
%                             either be 'oat' for the One-at-a-Time method,
%                             or 'monte-carlo'
%                    PERC - uncertainty of the parameter in percent, i.e.
%                           0.0 means zero uncertainty, and 1.0 an
%                           uncertainty within the 100% range of the
%                           parameter. When PERC is a vector, it must
%                           contain either 2 elements in the following
%                           form: [C_perc, E_perc], where K_perc is the
%                           global uncertainty for all model coefficients
%                           including the delay, and E_perc is that for all
%                           orders of the 's' operator. When PERC is a cell
%                           array, it must contain five entries in the
%                           following way: {Ca_p, Cna_p, Cb_p, Cnb_p, L_p},
%                           where the size of the vectors is determined
%                           by the size of the corresponding vectors of the
%                           FOTF object. See FOTFPARAM for details.
%                    N - (optional) the number of simulations to perform in
%                        case of the Monte-Carlo method. When the method is
%                        'oat', the number of simulations performed
%                        corresponds to the number of model parameters,
%                        which is also true for the 'monte-carlo' method,
%                        if this argument is omitted.
%
%   Output arguments: D - a cell array, each entry of which is a structure
%                         TSIM of the form
%                           TSIM.model - model with perturbed parameter(s),
%                           TSIM.resnorm - residual norm,
%                           TSIM.percfit - percent fit,
%                           TSIM.maxresid - maximum residual,
%                           TSIM.isstable - 0 or 1 - system stability,
%                           TSIM.resp - time-domain response to ID.u
%
%                     If output arguments are omitted, a plot will be drawn
%                     showing the response of the system to parameter
%                     variations.
%
%                     NB! This function can consume a significant amount of
%                     memory depending on the size of the identification
%                     dataset and desired amount of simulations!
%
% See also: fidata, fotf, fsparam, fotfparam

% Check input arguments
if nargin < 4
    error('TESTSIM:NotEnoughInputArguments', ...
        'Not enough input arguments.');
end

% Get the plant and its parameters
if ~isa(G, 'fsparam'), Gp = G; else Gp = G.plant; end
[a,na,b,nb,L] = fotfparam(Gp);

% Simulations amount
if nargin < 5 || strcmpi(method, 'oat')
    N = numel(a) + numel(na) + numel(b) + numel(nb) + numel(L);
end

% Create or retrieve the perturbation percentage
if isa(perc, 'cell')
    ap = perc{1}; nap = perc{2};
    bp = perc{3}; nbp = perc{4};
    Lp = perc{5};
else
    k1 = perc(1); k2 = perc(2);
    ap = k1*ones(1,numel(a)); nap = k2*ones(1,numel(na));
    bp = k1*ones(1,numel(b)); nbp = k2*ones(1,numel(nb));
    Lp = k1;
end

% Make sure every vector is a row vector
if size(ap,1)>size(ap,2), ap = ap'; end
if size(nap,1)>size(nap,2), nap = nap'; end
if size(bp,1)>size(bp,2), bp = bp'; end
if size(nbp,1)>size(nbp,2), nbp = nbp'; end

% Make sure we have valid entries
ap(ap<0)=0; ap(ap>1)=1;
nap(nap<0)=0; nap(nap>1)=1;
bp(bp<0)=0; bp(bp>1)=1;
nbp(nbp<0)=0; nbp(nbp>1)=1;
Lp(Lp<0)=0; Lp(Lp>1)=1;

% Construct a cell array containing parameter names (for OAT method)
paramNames = {};

% Current model index
kur = 1;
D = {};
switch(lower(method))
   
    case 'oat'
 
        len = init_display_progress(0,N);
        
        % Do 4 loops and a single simulation for the delay change
        for k=1:numel(a)
            
            pp = a;
            pp(k) = pp(k)+randsign(1).*rand(1).*ap(k).*pp(k);
            G_new = fotf(pp,na,b,nb,L);
            [y,res,pf,mr,ist] = do_simulation(G, G_new, id);
            tsim = struct;
            tsim.model = G_new; tsim.resnorm = res;
            tsim.percfit = pf; tsim.maxresid = mr;
            tsim.isstable = ist; tsim.resp = y;
            D{kur} = tsim;
            
            % Store parameter name
            paramNames{end+1} = ['a_{' num2str(numel(a)-k) '}'];
            
            % Update iteration count
            len = upd_display_progress(kur,N,len);
            kur = kur + 1;
        end
        
        for k=1:numel(na)
            
            pp = na;
            pp(k) = pp(k)+randsign(1).*rand(1).*nap(k).*pp(k);
            G_new = fotf(a,pp,b,nb,L);
            [y,res,pf,mr,ist] = do_simulation(G, G_new, id);
            tsim = struct;
            tsim.model = G_new; tsim.resnorm = res;
            tsim.percfit = pf; tsim.maxresid = mr;
            tsim.isstable = ist; tsim.resp = y;
            D{kur} = tsim;
            
            paramNames{end+1} = ['na_{' num2str(numel(na)-k) '}'];
            
            len = upd_display_progress(kur,N,len);
            kur = kur + 1;
        end
        
        for k=1:numel(b)
            
            pp = b;
            pp(k) = pp(k)+randsign(1).*rand(1).*bp(k).*pp(k);
            G_new = fotf(a,na,pp,nb,L);
            [y,res,pf,mr,ist] = do_simulation(G, G_new, id);
            tsim = struct;
            tsim.model = G_new; tsim.resnorm = res;
            tsim.percfit = pf; tsim.maxresid = mr;
            tsim.isstable = ist; tsim.resp = y;
            D{kur} = tsim;
            
            paramNames{end+1} = ['b_{' num2str(numel(b)-k) '}'];
            
            len = upd_display_progress(kur,N,len);
            kur = kur + 1;
        end
        
        for k=1:numel(nb)
            
            pp = nb;
            pp(k) = pp(k)+randsign(1).*rand(1).*nbp(k).*pp(k);
            G_new = fotf(a,na,b,pp,L);
            [y,res,pf,mr,ist] = do_simulation(G, G_new, id);
            tsim = struct;
            tsim.model = G_new; tsim.resnorm = res;
            tsim.percfit = pf; tsim.maxresid = mr;
            tsim.isstable = ist; tsim.resp = y;
            D{kur} = tsim;
            
            paramNames{end+1} = ['nb_{' num2str(numel(nb)-k) '}'];
            
            len = upd_display_progress(kur,N,len);
            kur = kur + 1;
        end
        
        pp = L;
        pp = pp+randsign(1).*rand(1).*Lp.*pp;
        G_new = fotf(a,na,b,nb,pp);
        [y,res,pf,mr,ist] = do_simulation(G, G_new, id);
        tsim = struct;
        tsim.model = G_new; tsim.resnorm = res;
        tsim.percfit = pf; tsim.maxresid = mr;
        tsim.isstable = ist; tsim.resp = y;
        D{kur} = tsim;
        
        paramNames{end+1} = 'L';
        
        upd_display_progress(kur,N,len);
        
    case 'monte-carlo'
        
        len = init_display_progress(0,N);
        
        % Create the new parameter matrices
        ppa = repmat(a,N,1);
        ppna = repmat(na,N,1);
        ppb = repmat(b,N,1);
        ppnb = repmat(nb,N,1);
        ppL  = repmat(L,N,1);
        
        % Perturb the parameters
        ppa = ppa + randsign(size(ppa)) .* ...
            rand(size(ppa)) .* repmat(ap,N,1) .* ppa;
        ppna = ppna + randsign(size(ppna)) .* ...
            rand(size(ppna)) .* repmat(nap,N,1) .* ppna;
        ppb = ppb + randsign(size(ppb)) .* ...
            rand(size(ppb)) .* repmat(bp,N,1) .* ppb;
        ppnb = ppnb + randsign(size(ppnb)) .* ...
            rand(size(ppnb)) .* repmat(nbp,N,1) .* ppnb;
        ppL = ppL + randsign(size(ppL)) .* ...
            rand(size(ppL)) .* repmat(Lp,N,1) .* ppL;
        
        % Do the simulation
        for k=1:N    
            
            G_new = fotf(ppa(k,:), ppna(k,:), ppb(k,:), ppnb(k,:), ppL(k));
            [y,res,pf,mr,ist] = do_simulation(G, G_new, id);
            tsim = struct;
            tsim.model = G_new; tsim.resnorm = res;
            tsim.percfit = pf; tsim.maxresid = mr;
            tsim.isstable = ist; tsim.resp = y;
            D{kur} = tsim;
            len = upd_display_progress(kur,N,len);
            kur = kur + 1;
            
        end
    
    otherwise
        error('TESTSIM:UnknownMethod', 'Unknown simulation method.');
end

% Move a few lines
fprintf(1,'\n\n');

% Just do the plot
if nargout < 1
    % Simulate original system to get the percent fit
    [e1, e2, origPercFit] = do_simulation(G,Gp,id);
    
    h = figure;
    subplot(2,1,1);
    for k=1:numel(D)
       st = D{k};
       if ~isempty(st.isstable) && st.isstable
        plot(id.t,st.resp);
        hold on;
       else
           warning('TESTSIM:UnstableSystem', ...
               ['Parameter variation may have resulted in an unstable ' ...
               'system - not showing.']);
       end
    end
    title('Test simulation runs of identified model');
    xlabel('Time [s]');
    ylabel('y(t) (under parameter uncertainties)');
    grid;
    
    myPercFit = [];
    for k=1:numel(D)
       st = D{k};
       if ~isempty(st.isstable) && st.isstable
            myPercFit(k) = st.percfit;
       else
            myPercFit(k) = 0;
       end
    end
    subplot(2,1,2);
    hh = bar(1:numel(D), myPercFit, 'r');
    hold on; plot([0 numel(D)+1], [origPercFit origPercFit], '--r');
    myYLim = ylim;
    myYLimA = myYLim(1); myYLimB = myYLim(2);
    myYLimA = myYLimA + sign(myYLimA)*25;
    myYLimB = myYLimB + sign(myYLimB)*25;
    ylim([myYLimA myYLimB]);
    xlim([0 numel(D)+1]);
    myPerfLabel = 'Test run #';
    if strcmpi(method, 'oat'), myPerfLabel = [myPerfLabel ...
            ' / Perturbed parameters']; end
    xlabel(myPerfLabel);
    ylabel('Percent fit [%]');
    
    % Set parameter names (OAT method only)
    if strcmpi(method, 'oat')
        
        for j=1:numel(D)
            x = j;
            y = myPercFit(j)+2;
            text(x,y,paramNames{j}, ...
                'Color','k',...
                'HorizontalAlignment','left',...
                'Rotation',90);
        end
    end
    
    set(h, 'NumberTitle', 'off');
    set(h, 'Name', 'System response due to parameter variations');
    
    varargout = {};
else
    varargout{1} = D;
end

end

function [y, resnorm, percfit, maxresid, isstab] = ...
    do_simulation(G, G_new, id)

% Check the sign of frac. powers before simulation
[a,na,b,nb,L] = fotfparam(G_new);
na(na<0) = 0; nb(nb<0) = 0;
G_new = fotf(a,na,b,nb,L);

% Determine the type of simulation (GL or approximated classic LTI)
if ~isa(G, 'fsparam')
    % Simulate system using GL based solver
    y = lsim(G_new, id.u, id.t);
else
    % Use the approximation
    Z = oustapp(G_new, G.w(1), G.w(2), G.N, G.approx);
    y = lsim(Z, id.u, id.t);
end

% Compute the desired characteristics
err = y - id.y;
resnorm = sum(err.^2);
percfit = 100*(1-(norm(err))/(norm(id.y-mean(id.y))));
maxresid = max(err);
isstab = isstable(G_new);

end

function len = init_display_progress(k,N)
    
    str = [num2str(k) ' of ' num2str(N) ' completed.'];
    len = length(str);
    fprintf(1, ['Simulation progress: ' str]);
    
end

function len1 = upd_display_progress(k,N,len)

    str = [num2str(k) ' of ' num2str(N) ' completed.'];
    len1 = length(str);
    for p=1:len, fprintf(1,'\b'); end
    fprintf(1, str);

end

