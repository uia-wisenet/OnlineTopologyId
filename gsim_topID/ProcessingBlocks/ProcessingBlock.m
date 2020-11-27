classdef ProcessingBlock < matlab.mixin.Heterogeneous & handle
	% This class allows one to obtain information about parameters of
	% descendants (e.g. generators, samplers, estimators) in text form.
	% This is useful to make figures. 
	%
	% all classes extending ProcessingBlock shall allow their constructors to be
	% called without any parameter.
    %
	% EXTENDING THE CLASS PARAMETER
	%
	% - It is recommendable to write a constructor that invokes the 
	% constructor of ProcessingBlock (of the constructor of the parent class,
	% which is a sublcass of parameter) with all arguments. Example:
	%
	%	function obj = GraphGenerator(varargin)
	%		obj@ProcessingBlock(varargin{:});
	%	end
	%
	% - It is also recommendable to extend the methods printProperty and
	% (possibly) propertiesToPrint.
	%
	% DEVELOPER: Daniel Romero
	%

	
	
	properties
		c_replicatedVerticallyAlong = {};   % property names in this list 
		% are used to make legends In an array of objects of class 
		% ProcessingBlock, only the values of this property for the objects in the first
		% column are considered. If this property is empty for the entry
		% (n,1), then the one corresponding to the entry (m,1) is used,
		% where m is the largest index less than n such that the element
		% (m,1) has a non-empty c_replicatedVerticallyAlong field.
		%
		% If all legend entries have to have the same structure, just set
		% this field for the (1,1) element.
		
		c_replicatedHorizontallyAlong = {}; % property name in this list is 
		% used to make the label of the x-axis. For arrays of objects of
		% class ProcessingBlock, the only considered is the (1,1) entry.

		def_chars_per_line = 200;
	end
	
	methods % to be extended by subclasses

		function ch_out = printProperty(obj,ch_propertyName,b_justName)
			% ch_out is a string displaying the value of the property whose
			% name is ch_propertyName. 
			%
			% b_justName (optional; default = 0)
			% if b_justName = 1, just the name of the property is printed.
			% Otherwise, both the name and value must be given.
			%
			% Example: obj has a property with name numberOfNodes with
			% value 4. 
			%
			%    >> obj.printProperty('numberOfNodes')
			%
			%        'N = 4'
			%
			%    >> obj.printProperty('numberOfNodes',1)
			%
			%        'N'
			%
			%
			
			assert(numel(obj)==1,'OBJ must be scalar');
			
			if nargin<3
				b_justName = 0;
			end
			
			ch_out = '';
			
			
			% backwards compatibility
			if isprop(obj,'c_parsToPrint')&&isprop(obj,'c_patternToPrint')&&isprop(obj,'c_stringToPrint')
				c_properties = obj.propertiesToPrint;
				for k=1:length(c_properties)
					if strcmp(c_properties{k},ch_propertyName)
						if b_justName
							ch_out = obj.c_stringToPrint{k} ;
						else
							if ismethod(obj,[ch_propertyName '_print'])
								ch_out = eval(['obj.' ch_propertyName '_print']);
							else
								ch_out = my_sprintf( obj.c_patternToPrint{k} , obj.c_stringToPrint{k} , my_getfield(obj,obj.c_parsToPrint{k}) );
							end
							break;
						end
					end
					
				end
			else
				warning(sprintf('method printProperty needs to be extended by class %s in order to automatically create legends, titles, and axis labels. See class ProcessingBlock for more information',class(obj)));
			end
			
		end
				
		function c_out = propertiesToPrint(obj)		
			% c_out is a cell array of strings containing the name of the
			% properties to print. For example, if obj has properties
			% prop1, prop2, and prop3, then c_out can be {'prop1','prop3'}
			
			% default behaviour
			c_out = 'default';
			
			% backward compatibility
			if isprop(obj,'c_parsToPrint')
				warning(sprintf('Property c_parsToPrint is set by class %s, but this behaviour is outdated',class(obj)));			
				c_out = obj.c_parsToPrint;
			else
				warning(sprintf('method printProperty needs to be extended by class %s in order to automatically create legends, titles, and axis labels. See class ProcessingBlock for more information',class(obj)));			
			end
			
		end
		
	end
	
	methods % constructor
		
		function obj = ProcessingBlock(varargin)
			% varargin is a sequence of (property_name,property_value)
			% pairs
			% Example:
			%     car = ProcessingBlock('numberOfWheels',4,'year',2009,'color','blue')
			obj =  assignParametersByName(obj,varargin{:});
			
		end
	
	end
	
	methods(Static)
			
		function [xlab,xvalues] = getXLabel(varargin)
			%
			%
			%
			for k = 1:length(varargin)
				obj_array_now = varargin{k};				
				obj1=obj_array_now(1,1);
				if isempty(obj1.c_replicatedHorizontallyAlong) || ...
                        isempty(obj1.c_replicatedHorizontallyAlong{1})
					continue;					
				end
				xlab = obj1.printProperty(obj1.c_replicatedHorizontallyAlong{1},1);
				xvalues =[ obj_array_now(1,:).(obj1.c_replicatedHorizontallyAlong{1})];
				return
			end
			error('Object array not replicated horizontally');
		end
			
		function [xvalues] = getXAxis(varargin)	
			[~,xvalues] = ProcessingBlock.getXLabel(varargin{:});
		end
			
		function leg = getLegend(varargin)
			%
			%
			%
			assert(nargin>0,'ProcessingBlock.getLegend is a static method that takes at least one ProcessingBlock argument' );
			
			list = ProcessingBlock.getLegendList(varargin{:});
			
			leg ={};
			for k = 1:size(list,1)
				leg{k} = ProcessingBlock.strListToText( {list{k,:}},1000 );
			end
			
		end
				
		function tit = getTitle(varargin)
			global chars_per_line
			if isempty(chars_per_line)
				chars_per_line = ProcessingBlock.def_chars_per_line;
			end
			
			list = ProcessingBlock.getTitleList(varargin{:});
			tit = ProcessingBlock.strListToText(list,chars_per_line);
		end				
		
	end
	
	methods(Static,Access=private)
		
		function leg = getLegendList(varargin)
			% it makes a 2D cell array with the entries for the legend
			%
			% leg = getLegendList(OBJ_1,OBJ_2,...,OBJ_N)
			%
			%   OBJ_n   : array of objects of class ProcessingBlock
			%
			%   leg     : cell array where the m-th row is a cell array
			%             with strings printing the values of the
			%             parameters whose names are in the field
			%             c_replicatedVerticallyAlong of OBJ_n(m,1). Those
			%             parameters are supposed to be common to all
			%             objects in the m-th row, but this is not checked.
			%
			%
			last_col = 0;
			leg = {};
			for k = 1:nargin
				obj_array_now = varargin{k};
				obj_now = obj_array_now(1,1);
				leg_pars = obj_now.c_replicatedVerticallyAlong;
				if (size(obj_array_now,1)<2)&&(~isempty(leg_pars)&&(~isempty(leg_pars{1})))
					str = sprintf('Property c_replicatedVerticallyAlong is non-empty for an array of class %s not truly replicated vertically',class(obj_array_now));
					warning(str);
				end
				
				for row = 1:size(obj_array_now,1)
					obj_now = obj_array_now(row,1);
					if ~isempty(obj_now.c_replicatedVerticallyAlong)
						leg_pars = obj_now.c_replicatedVerticallyAlong;
					end
					
					for par = 1:length(leg_pars)
						par_now = leg_pars{par};
						if isempty(par_now)
							continue
						end
						str = obj_now.printProperty(par_now);
						leg{row,last_col + par} = str;
					end
					
				end
				
				last_col = size(leg,2);
				
			end
		end
			
		function tit = strListToText(list,chars_per_line)
			% this function takes a list of strings and concatenates them
			% with commas in between and inserting an EOL when the length
			% of the current line exceeds CHARS_PER_LINE
			
			tit = [];
			tit_len = 0;
			for k=1:length(list)
				if isempty(list{k})
					continue
				end
				% end of line
				if tit_len + length(list{k}) > chars_per_line
					tit = [tit sprintf('\n')];
					tit_len = length(list{k});										
				else					
					tit_len = tit_len + length(list{k});					
				end
				tit = [tit list{k}];
				% comma
				if k~=length(list)
					tit = [tit ', '];
				end
			end

		end
					
		function tit = getTitleFilter(no_list,varargin)
			% tit is a string, with ends of line if needed  
			% 
			% tit is a concatenation of parameter values separated by
			% commas.
			% 
			% the parameter values with indices in no_list are excluded
			%
			
			global chars_per_line
			if isempty(chars_per_line)
				chars_per_line = ProcessingBlock.def_chars_per_line;
			end
			
			list = ProcessingBlock.getTitleList(varargin{:});
			assert(max(no_list)<=length(list));
			inds = setxor(1:length(list),no_list);
			list = list(inds);
			tit = [];
			tit_len = 0;
			for k=1:length(list)
				% end of line
				if tit_len + length(list{k}) > chars_per_line
					tit = [tit sprintf('\n')];
					tit_len = length(list{k});										
				else					
					tit_len = tit_len + length(list{k});					
				end
				tit = [tit list{k}];
				% comma
				if k~=length(list)
					tit = [tit ', '];
				end
			end

		end
				
		function list = getTitleList(varargin)
			% list is a cell array with the strings to put in the title.
			% arguments are objects of the class ProcessingBlock

			list = {};
			for nn = 1:nargin
				obj1 = varargin{nn};				
				obj1 = obj1(1,1);
				props = obj1.propertiesToPrint;
				for k=1:length(props)
					par_name = props{k};
					if (~is_c_replicatedVerticallyAlong(obj1,par_name) )&&(~is_c_replicatedHorizontallyAlong(obj1,par_name) )
						str = obj1.printProperty(par_name);
						if ~isempty(str)
							list = {list{:},str};
						end
					end
				end
			end
		end
		
	end
	
	methods
				
		function obj_mat = replicate(obj,fieldname_1,fieldvalues_1,fieldname_2,fieldvalues_2)
			% input 
			%    FIELDVALUES_1: cell array with M values for the file in
			%                   ProcessingBlock called FIELDNAME_1
			%    FIELDNAME_1  : string with the name of the field
			%    
			%    FIELDVALUES_2: cell array with N values for the file in
			%                   ProcessingBlock called FIELDNAME_1
			%    FIELDNAME_2  : string with the name of the field
			%
			% output
			%    obj_mat      : MxN matrix where the (m,n)-th element has
			%                   all the values equal to those of the
			%                   current object except from the field
			%                   FIELDNAME_1, which has value
			%                   FIELDVALUES_1{m} and the field
			%                   FIELDNAME_2, which has value
			%                   FIELDVALUES_2{n}, that is, FIELDNAME_1
			%                   replicates the object vertically, and
			%                   FIELDNAME_2 does it horizontally
			%
			
			assert(numel(obj)==1,'REPLICATE not implemented for input array objects');
			% this  function is implemented for scalar objects. see
			% replicate_horizontally above
			
			% check that fieldvalues are cells. You can use num2cell
			assert( (isempty(fieldvalues_1)||iscell(fieldvalues_1))&&(isempty(fieldvalues_2)||iscell(fieldvalues_2) ));
			
			obj.c_replicatedVerticallyAlong=[obj.c_replicatedVerticallyAlong {fieldname_1}];
			obj.c_replicatedHorizontallyAlong=[obj.c_replicatedHorizontallyAlong {fieldname_2}];
			M = max(length(fieldvalues_1),1);
			N = max(length(fieldvalues_2),1);
			obj_mat = repmat(obj,M,N);
			for m=1:M
				for n=1:N
					obj_mat(m,n)=clone(obj);
					if ~isempty(fieldname_1)
						obj_mat(m,n).(fieldname_1) = fieldvalues_1{m};
					end
					if ~isempty(fieldname_2)
						obj_mat(m,n).(fieldname_2) = fieldvalues_2{n};
					end
				end
			end
			
		end
				
		function d = is_c_replicatedVerticallyAlong(obj,str)
			
			obj = obj(1,1);
			
			if isempty( obj.c_replicatedVerticallyAlong )
				d = 0;
			else
				d = 0;
				for k=1:length( obj.c_replicatedVerticallyAlong )
					if strcmp( str , obj.c_replicatedVerticallyAlong{k} )
						d=1;
						return
					end
				end
			end
			
		end
				
		function d = is_c_replicatedHorizontallyAlong(obj,str)
			
			obj = obj(1,1);
			
			if isempty( obj.c_replicatedHorizontallyAlong )
				d = 0;
			else
				d = 0;
				for k=1:length( obj.c_replicatedHorizontallyAlong )
					if strcmp( str , obj.c_replicatedHorizontallyAlong{k} )
						d=1;
						return
					end
				end
			end
			
		end
				
	end
	
	% Methods to replicate the object (form a matrix)
	methods(Sealed) % sealed allows these methods to be executed on heterogeneous
		% arrays. I do not exactly understand how this attribute works anyway.
		
		function copiedObject = clone(obj)
			% copiedObject is a replica of obj. 
			%
			% If obj is an M x N array, then copiedObject is also M x N
		    % where the (m,n)-th entry of copiedObject is a clone of
		    % the (m,n)-th  entry of obj.
			%
			% Motivation: Since ProcessingBlock inherits
			% from handle, setting obj2 = obj1 does not clone the ProcessingBlock
			% object, just a pointer. This method creates a copy of the
			% handle object, not of the pointer.
			%
			% This method creates a new ProcessingBlock object  copiedObject with
			% the same property values as obj.
			% 
			% If a property of obj is itself a ProcessingBlock, it is also
			% copied in a recursive fashion.
			%
		    
			for k1 = size(obj,1):-1:1
				for k2 = size(obj,2):-1:1
					copiedObject(k1,k2) = cloneSingleObject(obj(k1,k2));
				end
			end
			
			
		end
		
		function copiedObject = cloneSingleObject(obj)
			% copiedObject is a replica of obj. 
			%
			% Motivation: Since ProcessingBlock inherits
			% from handle, setting obj2 = obj1 does not clone the ProcessingBlock
			% object, just a pointer. This method creates a copy of the
			% handle object, not of the pointer.
			%
			% This method creates a new ProcessingBlock object  copiedObject with
			% the same property values as obj.
			% 
			% If a property of obj is itself a ProcessingBlock, it is also
			% copied in a recursive fashion.
		
		
			copiedObject = eval(class(obj));
			c_properties = properties(obj);
			for k=1:length(c_properties)
				%eval(['value = obj.' c_properties{k}]);
				value = obj.(c_properties{k});
				
				% see if the property is a ProcessingBlock or handle
				if isa(value,'ProcessingBlock')
					value = clone(value); % value was a pointer to the ProcessingBlock object in the k-th property of obj. Now it points to a new object. 
				elseif isa(value,'handle')
					warning('A property of this object contains a handle object which will not be copied. Just a pointer to it will be copied');
				end
				
				% check if the property is read only
				mp = findprop(obj,c_properties{k});
				if strcmp(mp.SetAccess,'none')
					continue
				end
				
				% copy property
				%eval(['copiedObject.' c_properties{k} ' = value;']);			
				copiedObject.(c_properties{k}) = value;
			end
		end
		
		% We overwrite repmat to clone ProcessingBlock objects
		function m_out = repmat(obj,nRows,nColumns)
			[S1_in,S2_in] = size(obj);		
			for k1 = nRows:-1:1
				for k2 = nColumns:-1:1
					m_out((k1-1)*S1_in+1:k1*S1_in,(k2-1)*S2_in+1:k2*S2_in) = clone(obj);
				end
			end		
		end

	end
	
end


