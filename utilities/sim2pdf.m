function sim2pdf(diagram)
%SIM2PDF Create a PDF file from a Simulink diagram
if nargin < 1
    error('SIM2PDF:NotEnoughInputArguments',...
          'Not enough input arguments.');
end

% Check whether the diagram is already opened
openMdl     = find_system('SearchDepth', 0);
openThisMdl = max(strcmpi(openMdl,diagram));

% If it is not opened, attempt to open it
if ~openThisMdl
    open_system(which(diagram));
end

% Save the system as PDF
saveas(get_param(diagram,'Handle'), [diagram '.pdf']);