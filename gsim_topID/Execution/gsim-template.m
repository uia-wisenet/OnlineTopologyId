function gsim(onlyplot, experimentIndex, niter, experimentClassName)
%
% This is the main file and entry point for all executions. It is
% recommended that Git is configured to ignore this file, so it will not be
% synced through github. This means that each researcher can freely set
% parameters here without affecting other researchers.
%
% Input:
%        onlyplot   0: Execute and plot. 
%                   1: Do not execute, just display the results of the last
%                   execution of experiment indicated by experimentIndex
%        experimentIndex  [optional] index of the experiment within the 
%                   file indicated by the variable experimentFile below. If
%                   experimentIndex not given, its default value specified
%                   below is used. If experimentIndex is a vector, then all
%                   the experiments whose indices are in experimentIndex
%                   are executed. 
%        niter      [optional] parameter passed to the experiment function.
%                   It can be used by the experiment function e.g. as the
%                   number of iterations of some optimization or Monte
%                   Carlo algorithm. However, gsim does not use this
%                   parameter at all. 
%        experimentClassName [optional] string containing the name of the
%                   file of experiments to be used. This option is
%                   useful if you do not want to edit gsim.m every time you
%                   want to execute experiments from different files.
%
% Examples:
% >> gsim(0)          % Executes the default experiment and plots the results
% >> gsim(1)          % Plots the results of the last execution of the default
%                       experiment
% >> gsim(0,1102,400) % Executes the experiment 1102 with 400 iterations. 
% >> gsim(0, 2005, 250, 'TemplateExperiments') % Executes the experiment 2005 in
%                                   % TemplateExperiments.m with 250 iterations.
%
% Please report any bug or inquiry to daniel.romero@uia.no.
%

% INITIALIZATIONS - DO NOT EDIT ===========================================
addpath(genpath('./Execution/'));
initializeGsim.initializePath;
global displaySettings

% EXECUTION PARAMETERS - EDIT =============================================
defaultExperimentClassName   = 'TemplateExperiments';
% defaultExperimentClassName = 'TutorialGfigureExperiments';

defaultExperimentIndex = 1001; 
defaultNiter = 200 ;

runningOnServer = 0;       % optionally set this variable to 1 to prevent a 
                           % server without graphical interface from
                           % attempting to plot the results. 

% GRAPHICAL REPRESENTATION PARAMETERS - EDIT ==============================

displaySettings.b_savePlots=0;
                           % 1: save plotted figures to .fig and
						   %    .pdf/.eps.
						   % 0: do not save plotted figures to
						   %    files
displaySettings.b_outputInPdf = 1;    
                           % 1: write figures to pdf files (recommended)
						   % 0: write figures to eps files
displaySettings.b_titleToCaption = 1; 
                           % Set to 1 when exporting a figure for a paper.  									  
displaySettings.v_windowPosition = [];% 1 x 4 vector indicating the screen 
                           % position and size of the window containing the
                           % figure. Use the same vector for all figures in
                           % your paper. For example, 
		                   %   [100   402   444   273];
						   % You can also adjust this parameter graphically
						   % as described in the manual.
displaySettings.b_dockedOn = 0;
                           % Set to 1 to display multiple figures on the
                           % same window. 					 
displaySettings.c_globalTranslationTable = {}; 
                           % similar to c_localTranslationTable described 
						   % in GFigure.m,	but this one is applied to all
						   % figures after applying c_localTranslationTable. 	   
displaySettings.b_inhibitPause = 0; 
                           % do not pause between figures when displaying a
                           % figure sequence.  
displaySettings.titleFontSize = []; 
                           % set to an integer to use a non-default font
                           % size for the title of figures.  

% DO NOT EDIT AFTER THIS POINT ============================================

% Input check
assert(exist('./gsim.m', 'file')>0, 'The gsim.m function should be called from the gsim folder.')
assert(nargin>=1, 'Not enough arguments. Usage:  >> gsim(onlyplot[,experimentIndex[,niter[,experimentClassName]]])');
assert(nargin<=4, 'Too many arguments.');

if nargin<4, experimentClassName  = defaultExperimentClassName;   end
if nargin<3 || isempty(niter),      niter = defaultNiter;         end 
if nargin<2, experimentIndex      = defaultExperimentIndex;       end

if onlyplot>1, 	warning('Onlyplot is greater than 1.'); end
assert(isnumeric(niter) && niter == round(niter) && niter >= 0, ...
		'niter must be a nonnegative integer');

% Create object of experiments
try
	experimentObj = feval(experimentClassName);
catch
	error ([experimentClassName '.m not found.'])
end

% Run the simulation
experimentObj.gsimExecute(experimentIndex, niter, onlyplot, runningOnServer)

end

