function varargout = gsinfo(varargin)
% GSINFO returns a list with the available experiment collections in the Experiments
% folder.
% GSINFO(ExperimentName) prints a list with the available experiments in
% the specified file.
% c_list = GSINFO(...) returns a cell structure with function names,
% commands, and HTML links 

% Developer: Luismi

try
	gsim
catch ME
	if contains(ME.message,'Not enough arguments.')
		% OK
	elseif contains(ME.message, 'from the gsim folder')
		error 'GSINFO should be called from the gsim folder.'
	else
		rethrow(ME)
	end
end


switch nargin
	case 0 %If GSINFO is called without arguments, 
		% list available experiment function sets
		my_list = dir('Experiments');
		cell_out = cell(0);
		for i = 1:length(my_list)
			entry = my_list(i).name;
			if contains(entry, 'Experiments') && ~contains(entry, '_data')
				aux = strrep(entry, 'Experiments', '');
				expName = strrep(aux, '.m', '');
				if not(isempty(expName))
					commandStr = sprintf('gsinfo %s', expName);
					listStr = sprintf('<a href="matlab:%s">%s</a>', commandStr, expName);
					cell_out = cat(1, cell_out, {expName, commandStr, listStr});
				end
			end
			% strjust(pad, strings, 'center')) % use to center
		end
		if nargout == 0
			fprintf('I found %d experiment function sets:\n', length(cell_out))
			for i = 1:size(cell_out,1)
				disp(strjust(pad(cell_out{i, end}), 'center'));
			end
		end
		
	case 1 %list available experiments in the specified file
		className = fixClassName(varargin{1});
		try
			methodList = methods(className);
			cell_out = cell(0);
			for i = 1:length(methodList)
				methodName = methodList{i};
				if contains(methodName, 'experiment_')
					numStr = strrep(methodName, 'experiment_', '');
					cell_out = cat(1, cell_out, getCommands(numStr, className, varargin{1}));
				end
			end
			nMethods = size(cell_out, 1);
		catch ME
			if isequal(ME.identifier, 'MATLAB:class:InvalidSuperClass')
				wraptext(ME.message);
				nMethods = -1;
			else
				rethrow(ME);
			end
		end
		
		if nargout == 0
			disp(helpfunc(className));
			if nMethods > 0
				fprintf('<a href="matlab:gsinfo %s">%s</a> contains %d experiments:\n', varargin{1}, className, nMethods)
				disp(cell_out)
			elseif nMethods == 0
				fprintf('<a href="matlab:gsinfo %s">%s</a> contains no experiments.\n', varargin{1}, className);
			else
				fprintf('<a href="matlab:gsinfo %s">%s</a> has some problem (look above).\n', varargin{1}, className);
				fprintf('<a href="matlab:open(''%s'')">Open %s</a>\n', className, varargin{1});
			end
% 			for i = 1:size(cell_out, 1)
% 				disp(cell_out(i, :));
% 			end
		end
		
	case 2 % show description for the specified experiment function
		className = fixClassName(varargin{1});
		numStr = varargin{2};
		methodName = ['experiment_' numStr];
		addpath Experiments
		assert(exist(className, 'class')>0, '%s is not a valid Experiment function set or is not available in the path', className);
		assert(ismethod(className, methodName), '%s does not contain such experiment %s', className, numStr);
		methodsList = methods(className);
		for i = 1:length(methodsList)
			if strcmp(methodsList{i}, methodName)
				cell_out = getCommands(numStr, className, varargin{1});
			end
		end
		if nargout == 0
			help([className '.' methodName]);
			disp(cell_out)
		end
		
	otherwise
		error 'usage: gsinfo [experimentName [experimentNumber]]'
end

if nargout == 1
	varargout{1} = cell_out;
end

end

function className_out = fixClassName(c_in)
    className_out = [c_in 'Experiments'];
	%addpath Experiments
	if exist(className_out, 'class')
		%OK
	else
		dir_list = dir('Experiments');
		for i_file = 1:length(dir_list)
			fileName = dir_list(i_file).name;
			if strcmpi(fileName, [className_out '.m'])
				warning 'found a case-insensitive match.'
				className_out = strrep(fileName, '.m', '');
				break
			end
		end
	end
	assert(exist(className_out, 'class')>0, ['%s is not a valid '...
		'ExperimentFunctionSet or is not available in the path'], ...
		className_out);
end

function c_out = getCommands(numStr, className, setName)
    methodName = ['experiment_' numStr];
    dataFileName = [className '_' numStr '.mat'];
	if exist(dataFileName, 'file')
		runStr = sprintf('<a href="matlab:fprintf(''gsim(0, %s, [], ''''%s'''')\\n'')">cmd:rerun exp_%s</a>', numStr, className, numStr);
		figStr = sprintf('<a href="matlab: gsim(1, %s, [], ''%s'') %%onlyplot">Show Fig. %s</a>', numStr, className, numStr);
	else
		runStr = sprintf('<a href="matlab: gsim(0, %s, [], ''%s'')">Run exp_%s</a>', numStr, className, numStr);
		figStr = '(no Fig)';
	end
	editStr = sprintf('<a href="matlab: matlab.desktop.editor.openAndGoToFunction(which(''%s''), ''%s'');">Edit exp_%s</a>', className, methodName, numStr);
	commandStr = ['gsinfo ' setName ' ' numStr];
	infoStr = sprintf('<a href="matlab: %s ">Info on exp_%s</a>', commandStr, numStr);
	c_out = {numStr, figStr, runStr, editStr, infoStr};	
end