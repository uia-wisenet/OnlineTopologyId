classdef GFigure
	%
	% Objects of this class contain all information needed to plot a figure. Set
	% its properties appropriately and then plot invoking the method plot.
	% 
	% Examples: see function "GFigure.test" below.
	%
	% Vectors of objects GFigure are used to represent multiple figures.
	%
	% Each entry of this vector corresponds to a figure, with one or more
	% axes (subplots). It may directly contain the figure data in the properties 
	% 'm_X' and 'm_Y', or the subplots may be contained as GFigure objects in the
	% property 'm_multiplot'. By default, a subplot is created with the same
	% dimensions as the matrix m_multiplot. 
	%
	
	
	
	properties
		
		% General properties ----------------------------------------------
		axes;              % axes object if not null then plot it. Use 
		                   % F = GFigure.captureCurrentFigure()   
		                   % if you want to create an object of class
		                   % GFigure that contains the contents of the
		                   % current axes. 		   

		figureNumber = 1;  % number of the first figure		

		% Multiplot properties --------------------------------------------
		ch_multiplotType = 'subplot';
		                   % this can be:
						   %     'subplot' : a subplot is made with the
						   %                 same dimensions as the matrix
						   %                 of objects of class GFigure 
						   %                 stored in m_multiplot
						   %     'sequence': each element of m_multiplot
						   %                 is represented sequentially 
						   %                 and a pause command is
						   %                 executed between the display
						   %                 of two consecutive elements.
						   %                 Each element of the sequence 
						   %                 can have subplots.
		
        m_multiplot;       % Matrix/vector of objects of class GFigure
						   % If m_multiplot is not empyt, then m_X, m_Y and
						   % m_Z should be empty
		
		
		% Figure properties (disregarded in multiplot mode) ---------------
		
		% 2D and 3D plots
		m_X = [];          % In 2D plots:
		                   % Ncurves x Npoints matrix or 1 x Npoints vector
		                   % with  the X-coordinates. In the latter case,
		                   % it is understood that all curves share the
		                   % X-axis. Use method formMatrixByPaddingRows if
		                   % needed 
						   % In 3D plots: see m_Z
	    m_Y = [];          % In 2D plots:
		                   % Ncurves x Npoints matrix with the
		                   % Y-coordinates 
						   % In 3D plots: see m_Z
		
		v_xlim = [];       % 1 x 2 vectors with the limits for the x-axis
		v_ylim = [];       % 1 x 2 vectors with the limits for the y-axis
		
		b_logx = 0;        % set to 1 to have a logarithmic x-axis
        b_logy = 0;        % set to 1 to have a logarithmic y-axis
		
		m_errorBars = [];  % Ncurves x Npoints matrix with the semi-length 
		                   % of the error bars.If an element is NaN, no
		                   % error bar is represented for the corresponding
		                   % point. 
		
		% 3D plots only
		m_Z = [];          % in 3D plots, the displayed surface joins the 
		                   % points m_X(ind_x,ind_y), m_Y(ind_x,ind_y), and
						   % m_Z(ind_x,ind_y) for all ind_x and ind_y. The
						   % dimensions of m_X, m_Y, and m_Z is therefore
						   % numberOfPointsXAxis x numberOfPointsYAxis.		
		v_caxis =[];       % Color axis. This is the argument to the caxis 
		                   % command. If it is empty, the default [min
		                   % max] is applied. 
		b_viewFromTop = 0; % if ~= 0, the plot is viewed from the top -> 
		                   % only colors, but the Z axis becomes visible
		                   % by manually rotating the axes in the MATLAB
		                   % figure. 
							
		% text
		ch_title = '';     % the property name is self-explanatory
		ch_xlabel = '';    % the property name is self-explanatory
		ch_ylabel = '';    % the property name is self-explanatory
		ch_zlabel = '';    % the property name is self-explanatory
		c_legend = {};     % 1 x Ncurves cell array, where c_legend{n} is the  
		                   % char vector with the legend string of the n-th
		                   % curve
		ch_legendPosition = '';      
		                   % it can be 'northwest', 'south', etc
		v_legendPosition = [];  
		                   % this must be a vector with four entries 
		                   % which indicate the legend position within the 
						   % figure. 
		                   % For a figure that is being displayed, you can
		                   % move the legend with the graphical tool to the
						   % position you want and then use the function
						   % get_ch_legendPosition(fig) to get the vector
						   % corresponding to that legend position. 
		
		% colors and c_styles
		m_colorset =[0 0 0;1 0 0;0 0 .9 ;0 .65 0; .9 0 .9 ;.5 .5 0;0 .7 .7;...
			.5 0 1; 1 .5 0;  1 0 .5; 0 1 .5;0 1 0];
		                   % The n-th row is the RGB color of the n-th
		                   % curve. 
		colorPeriod = 12;  % The color of the curves repeats every 
		                   % colorPeriod curves. <colorPeriod> cannot be
		                   % larger than size(m_colorset,1)
			               % NOTE: color can also be set in the style
			               % strings, e.g. '-w' makes a solid white line.  
						   % For plot3 this is the only possibility with
						   % the current implementation. 
						 
		c_styles = {'-','-','-','-','-','-','-','-','-','-','-','-',...
			'--','--','--','--','--','--','--',...
			'-.','-.','-.','-.','-.','-.','-.',...
			':',':',':',':',':',':',':'};
		                   % c_styles{n} is the style of the n-th curve.
		                   % The format is the same as for MATLAB plot
		                   % function; type <help plot> in the command line
		                   % to see the description.
		           
		ch_plotType2D = 'plot'; 
		                   % other possibilities: 'stem' and 'bar'
		ch_plotType3D = 'imagesc'; 
		                   % other possibilities: 'plot3','surf' 
		
		ch_gridStyle = '--';
		                   % style for the grid. Other values are '' and ':'
		
		                   
		
		ch_caption = '';   % if it is not empty, this string is printed
		                   % and saved to a txt file. If it is empty
		                   % and the value of the global variable
		                   % <b_titleToCaption> is different from 0, then
		                   % <ch_caption> is set equal to <ch_title>.
		
		c_localTranslationTable = {}; 
		                   % If this is an N x 2 cell array of strings, the
						   % occurences of c_localTranslationTable(n,1) in
						   % either the title, X-label, Y-label, Z-label,
						   % or legend with c_localTranslationTable(n,2) for all
						   % n
						   %
						   % If this is an N x 3 cell aray of strings, the
						   % ocurrences of
						   % c_localTranslationTable(i,1) are replaced with
						   % c_localTranslationTable(i,3) in the title and
						   % caption but they are replaced with
						   % c_localTranslationTable(i,2) in the 
						   % legend, X-label and Y-label
						   %
		ch_interpreter = ''; 
		                   % can be tex, none
		                 
	end
	
	
	methods
		
		% Constructor
		function obj = GFigure(varargin)
			% 
			% Syntax: 
			% >>  F = GFigure(fieldName,value,fieldName,value...)
			%
	
			if nargin>0
								
				if mod(length(varargin),2)~=0
					error('the number of arguments must be even');
				end
				
				for k=1:2:length(varargin)-1
					obj.(varargin{k})=varargin{k+1};
				end
								
			end
		end
		
        % Main method of the class. 
		function plot(v_GFigure,figureNumber)		
			% 
			% This function displays and optionally saves the figures
			% corresponding to the GFigure objects in the vector
			% v_GFigure. 
			%
			% v_GFigure  : vector of objects of class GFigure
			% figureNumber: [optional] number of the window to contain the
			%               figure corresponding to v_GFigure(1)
			global displaySettings
			
			if nargin>1
				experimentIndex = figureNumber;
			else
				experimentIndex = v_GFigure(1,1,1).figureNumber;
			end			
			
			v_GFigure = v_GFigure(:);
				
			for k=1:length(v_GFigure)
				
				F_now = v_GFigure(k);				
				F_now.check_parameters;
			
				% Prepare a new window
				figure(experimentIndex+k-1)			
				if (~isempty(displaySettings))
					if isfield(displaySettings,'b_dockedOn') && displaySettings.b_dockedOn
						set(gcf,'WindowStyle','docked');
					else
						if isfield(displaySettings,'v_windowPosition')&&~isempty(displaySettings.v_windowPosition)
							set(gcf,'position',displaySettings.v_windowPosition);
						end
					end
				end
				
				% Plot the figure and save it in the file system
				if ~isempty(F_now.axes)
					try
						copyobj(F_now.axes,gcf);
						F_now.save_figure(experimentIndex,k);
					catch ME						
						switch class(ancestor(F_now.axes, 'figure'))
							case 'matlab.ui.Figure'
								display (ME.message);
								disp 'This error will appear when the GFigure object was created'
								disp 'from a figure and that figure was deleted. This will not happen if '
								disp 'the GFigure object has been loaded from a .MAT file.'
								
								error 'Axis handle was deleted (probably because its parent figure was closed).'
							otherwise
								error(ME.message);
						end
					end
				else
					plot_scalar_F(F_now,experimentIndex,k); 
				end
															
			end
			
			
		end	
		
	end
	
	methods(Access=private)
		
		function plot_scalar_F(F,experimentIndex,fig_index)
			% 
			% F corresponds to a figure.
			%
			% The name of the file where the figure is stored is
			%    figbasename_experimentIndex-fig_index.pdf
			% where figbasename is given by the globals
			%
			global displaySettings
			
			F = F.moveTitleToCaptionIfNeeded();
								
			% Handle multiplots if needed
			if isempty(F.m_multiplot) % then we plot a single axis
				subplot(1,1,1)
				F.plotAxis;
				F.print_caption;				
				F.save_figure(experimentIndex,fig_index);
			else % then we need to see which multiplot type we have
				switch(F.ch_multiplotType)
					case 'subplot'						
						nrows = size(F.m_multiplot,1);
						ncols = size(F.m_multiplot,2);
						for row = 1:nrows
							for col = 1:ncols
								subplot(nrows,ncols,(row-1)*ncols + col );
								plotAxis(F.m_multiplot(row,col,:));
							end
						end
						F(1,1,1).print_caption;
						F(1,1,1).save_figure(experimentIndex,fig_index);
					case 'sequence'
						Fseq = F.m_multiplot(:);
						for k=1:length(Fseq)
							plot_scalar_F(Fseq(k),experimentIndex,fig_index+k-1);
							if (~isempty(displaySettings))&&(isfield(displaySettings,'b_inhibitPause'))&&(displaySettings.b_inhibitPause)
								% not pausing
							else
								disp('Paused. Press a button to continue.')
								pause();
							end							
						end
								
					otherwise
						error('F.ch_multiplotType does not contain a valid string');
				end
				
			end
		end
		
		function F = moveTitleToCaptionIfNeeded(F)

			global displaySettings

			% parsing globals
			if ~isempty(displaySettings)
				if isfield(displaySettings,'b_titleToCaption')&&...
						displaySettings.b_titleToCaption&&...
						isfield(displaySettings,'b_savePlots')&&...
						displaySettings.b_savePlots
					if F.ch_title
						F.ch_caption = F.ch_title;
						F.ch_title = '';
					end
				end
			end
		end
			
		function plotAxis(F)
			% representation of the curves
			assert(numel(F)==1)
			F = F.translateAxis;

			% Actual representation
			if isempty(F(1,1,1).m_Z)
				plotAxis_2D(F);
			else
				plotAxis_3D(F);
			end
			
			% axis limits
			if ~isempty(F(1,1,1).v_xlim)
				xlim(F(1,1,1).v_xlim);
			end
			if ~isempty(F(1,1,1).v_ylim)
				ylim(F(1,1,1).v_ylim);
			end
			
			% text	
			title(F(1,1,1).ch_title);			
			global displaySettings					
			if (~isempty(displaySettings))&&(isfield(displaySettings,'titleFontSize'))&&(~isempty(displaySettings.titleFontSize))
				set(get(gca,'Title'),'FontSize',displaySettings.titleFontSize);
			end
			if ~isempty(F(1,1,1).ch_interpreter)
				set(get(gca,'Title'),'Interpreter',F(1,1,1).ch_interpreter);
			end
			
			% Axis labels
			xlabel(F(1,1,1).ch_xlabel);
			ylabel(F(1,1,1).ch_ylabel);
			zlabel(F(1,1,1).ch_zlabel);
			
			% Legend
			if ~isempty(F(1,1,1).c_legend)
				if ~isempty(F(1,1,1).ch_legendPosition)
					legend(F(1,1,1).c_legend,'Location',F(1,1,1).ch_legendPosition);
				else
					legend(F(1,1,1).c_legend);
				end
				if ~isempty(F(1,1,1).v_legendPosition)
					if ~isempty(F(1,1,1).ch_legendPosition)
						warning('ignoring F.ch_legendPosition');
					end
					cv = get(gcf,'Children');
					set(cv(1),'Position',F(1,1,1).v_legendPosition);
				end
				if ~isempty(F(1,1,1).ch_interpreter)
					set(get(gca,'Legend'),'Interpreter',F(1,1,1).ch_interpreter);
				end
			end
				
		end
			
		function plotAxis_2D(F)
			ncurves = size(F.m_Y,1);
			
			% complete F.c_styles if more curves to plot than
			% length(F.c_styles)
			if (length(F.c_styles)<ncurves)
				for k = length(F.c_styles)+1:ncurves
					F.c_styles{k}='-';
				end
			end
				
			% default 
			if isempty(F.m_X)
				F.m_X = 1:size(F.m_Y,2);
			end
			if size(F.m_X,1) ~= ncurves
				if size(F.m_X,1) == 1
					F.m_X = ones(ncurves,1)*F.m_X;
				else
					error('size(m_X,1) must be either 1 or size(m_Y,1)');
				end
			end
			if size(F.m_X,2)~=size(F.m_Y,2)
				error('m_X has %d columns whereas m_Y has %d',size(F.m_X,2),size(F.m_Y,2));
			end
			
			for kr=1:ncurves
				assert(~isempty(F.c_styles{kr}));
				switch(F.ch_plotType2D)
					case 'plot'
						if F.b_logx
							if F.b_logy
								loglog(F.m_X(kr,:),F.m_Y(kr,:),F.c_styles{kr},'LineWidth',2);
							else
								semilogx(F.m_X(kr,:),F.m_Y(kr,:),F.c_styles{kr},'LineWidth',2);
							end
						else
							if F.b_logy
								semilogy(F.m_X(kr,:),F.m_Y(kr,:),F.c_styles{kr},'LineWidth',2);
							else								
								plot(F.m_X(kr,:),F.m_Y(kr,:),F.c_styles{kr},'LineWidth',2);
							end
						end
					case 'stem'
						stem(F.m_X(kr,:),F.m_Y(kr,:),F.c_styles{kr});
					case 'bar'
						bar(F.m_X(kr,:),F.m_Y(kr,:));
					otherwise
						error('unrecognized plot type')
				end
				hold on
			end
			hold off
			if strcmp(F.ch_plotType2D,'bar')  % bar color cannot be set like plot objects
				return
			end
			
			% colors
			chldv=get(gca,'Children');
			for kr=1:size(F.m_Y,1)
				
				get(chldv(length(chldv)-kr+1),'XData');
				try
					set(chldv(length(chldv)-kr+1),'Color',F.m_colorset(mod(kr-1,min(F.colorPeriod,size(F.m_colorset,1)))+1,:));
				catch
					warning('it seems we are trying to set the color of something that is not a curve. Please report this issue.')		
				end
			end
			
			if ~isempty(F.ch_gridStyle)
				grid on
				set(gca,'GridLineStyle',F.ch_gridStyle);
			end
			
			if ~isempty(F.m_errorBars)
				F.plot_error_bars;
			end
			
			set(gcf,'PaperPositionMode','auto');
			
		end
		
		function plotAxis_3D(F)
			if isempty(F.m_Z)
				error('m_Z cannot be empty in 3D figures')
			end
			switch(F.ch_plotType3D)
				% note that some of these figures may not display correctly in
				% pdf. This is a bug and can be found on export_fig documentation.
				% the same happens if we export through matlab standard commands.
				% A solution is to export to eps.
				case 'imagesc'
					pcolor(F.m_X,F.m_Y,F.m_Z);
					colormap(gcf,jet(100));
				case 'surf'
					surf(F.m_X,F.m_Y,F.m_Z);					
				case 'plot3'
					plot3(F.m_X,F.m_Y,F.m_Z,F.c_styles{1},'LineWidth',2);
			end
					
			if F.b_viewFromTop 
				view([0 90]);
			end
			shading interp
			if ~isempty(F.v_caxis)
				caxis(F.v_caxis);
			end
			
		end
			
		function plot_error_bars(F)
						
			hold on
			s_bar_width = (max(F.m_X(:))-min(F.m_X(:)))/100;  % same as in errorbar.m
			
			for m = 1:size(F.m_errorBars,1)
				vx = [];
				vy = [];
				for n = 1:size(F.m_errorBars,2)
					x = F.m_X(m,n);
					y = F.m_Y(m,n);
					slen = F.m_errorBars(m,n);
					vx_now = [x       x       NaN   x-s_bar_width   x+s_bar_width   NaN   x-s_bar_width   x+s_bar_width   NaN];   
					vy_now = [y+slen  y-slen  NaN   y+slen          y+slen          NaN   y-slen          y-slen          NaN];
					
					vx = [vx vx_now];
					vy = [vy vy_now];
				end
				plot(vx,vy,'-k')
			end
			
			hold off

			% colors set as the figures
			chldv=get(gca,'Children');
			for kr=1:size(F.m_Y,1)
				get(chldv(size(F.m_Y,1)-kr+1),'XData');
				set(chldv(size(F.m_Y,1)-kr+1),'Color',F.m_colorset(mod(kr-1,F.colorPeriod)+1,:));
			end
			
		end
		
		function check_parameters(F)
			% this function is for scalar objets of class GFigure		
			assert(numel(F)==1);
			if ~isempty(F.m_multiplot)
				assert(isempty(F.m_X)&&isempty(F.m_Y)&&isempty(F.m_Z),'Whenever the property m_multiplot is not empty, then m_X, m_Y, and m_Z must be empty.');
			end
			
		end
			
		function print_caption( F )
			if ~isempty(F.ch_caption)
				fprintf('%s\n',F.ch_caption);
			end
		end
		
		function save_figure( F , experimentIndex,fig_index)
			
			global displaySettings

			% parsing globals
			if ~isempty(displaySettings)
				if isfield(displaySettings,'b_savePlots')
					b_savePlots= displaySettings.b_savePlots;
				else
					b_savePlots = 0;
				end
				if isfield(displaySettings,'figbasename')
					figbasename= displaySettings.figbasename;
					                   % this should be assigned in
					                   % ExperimentFunctionSet.m
				else
					figbasename = 'figure_';
				end
				if isfield(displaySettings,'b_outputInPdf')
					b_outputInPdf = displaySettings.b_outputInPdf;
				else
					b_outputInPdf = 1;
				end
				
			else
				% we do not support saving files without globals beind defined
				return
			end
			
			if ~b_savePlots
				return
			end
			
			name = sprintf('%s_%d-%d',figbasename,experimentIndex,fig_index);

			set(gcf,'color',[1 1 1]);
			if b_outputInPdf % exporting to pdf
				fprintf('Saving to %s\n', [name '.pdf']);
				export_fig(name,'-pdf');
			else % exporting to eps
				fprintf('Saving to %s\n', [name '.eps']);
				export_fig(name,'-eps');
			end
			set(gcf,'color',0.8*[1.0 1.0 1.0])

			% saving .fig file
			saveas(gcf,[name '.fig']);
			
			% saving caption to txt file
			caption_text = F.ch_caption;			
			if ~isempty(caption_text)
				caption_fn = [name '.txt'];
				fprintf('Saving caption to %s\n',caption_fn);
				fdes=fopen(caption_fn,'w');
				fprintf(fdes,caption_text);
				fclose(fdes);
			end					
			
		end
		
        function F = translateAxis(F)
			
			% local table
			if ~isempty(F.c_localTranslationTable)
				F = F.applyTranslationTableToAxis(F.c_localTranslationTable);
			end
			
			% global table
			global displaySettings
			if (~isempty(displaySettings))&&(isfield(displaySettings,'c_globalTranslationTable'))
				F = F.applyTranslationTableToAxis(displaySettings.c_globalTranslationTable);	
			end
			
		end
		
		function F = applyTranslationTableToAxis(F,translation_tab)
						
			% splitting a possible N x 3 table into two Nx2 tables
			if isempty( translation_tab )
				return;
			elseif size(translation_tab,2) == 2
				title_table = translation_tab;
				common_table = translation_tab;
			elseif size(translation_tab,2) == 3
				rows = size(translation_tab,1);
				title_table = cell(rows,2);
				common_table = cell(rows,2);
				for k=1:rows
					common_table{k,1} = translation_tab{k,1};
					common_table{k,2} = translation_tab{k,2};
					title_table{k,1} = translation_tab{k,1};
					title_table{k,2} = translation_tab{k,3};
				end
			elseif ~isempty(translation_tab)
				error('the translation table must have 2 or 3 columns');
			end
			
			%%%% translating
			
			if (~isempty(F.ch_title))
				F.ch_title = GFigure.translate_string(F.ch_title,title_table);
			end
			if (~isempty(F.ch_xlabel))
				F.ch_xlabel = GFigure.translate_string(F.ch_xlabel,common_table);
			end
			if (~isempty(F.ch_ylabel))
				F.ch_ylabel = GFigure.translate_string(F.ch_ylabel,common_table);
			end
			if (~isempty(F.c_legend))
				F.c_legend = GFigure.translate_string(F.c_legend,common_table);
			end				
			if (~isempty(F.ch_caption))
				F.ch_caption = GFigure.translate_string(F.ch_caption,title_table);
			end	
		end
			
	end	
	
	
	methods(Static)
		
		function gf_out = captureCurrentFigure()
			% F is an object of class GFigure that contains what is
			% displayed in all axes in the current figure.
			
			currentFigure = gcf;
			gf_out=GFigure('axes',currentFigure.Children);
			
		end
		
		function gf_out = captureCurrentAxis()
			% returns a GFigure that contains what is displayed in the
			% current axis. Observe the difference with respect to
			% captureCurrentFigure, which captures all axis in the current
			% figure.
			
			currentAxis = gca;
			gf_out = GFigure('axes', currentAxis);
		end
		
		function Mat = formMatrixByPaddingRows(varargin)
			%
			% function Mat=fill_mat_f(row1, row2, ...)
			%
			%   Arranges all the arguments as rows of a matrix, padding with NaN at the
			%   end if needed. Useful for plotting curves with different lengths in
			%   plot_f
			%
			
			for m=1:length(varargin)
				rows(m)=size( varargin{m} , 1 );
				cols(m)=size( varargin{m} , 2 );
			end
			N=max(cols);
			M=sum(rows);
			Mat=zeros(M,N);
			for m=1:length(varargin)
				row_start = sum( rows(1:m-1) )+1;
				row_end = row_start + rows(m)-1;
				Mat(row_start:row_end,:)=[varargin{m} ones(rows(m),N-cols(m))*NaN];
			end
			
			
		end
		
		function test
			
			% Example of how to plot multiple curves
		    simple_test = 1;			
			if simple_test
				
				Ncurves = 4;
				Npoints = 40;
				
				m_X = linspace(3,10,Npoints);
				m_Y = randn(Ncurves,Npoints);
				
				my_F = GFigure('m_X',m_X,'m_Y',m_Y);
				
				my_F.plot;
			end
			
			% Example of how to create a sequence of figures
			sequence_test = 0;
			if sequence_test
				
				Ncurves = 4;
				Npoints = 40;
				
				m_X = linspace(3,10,Npoints);
				
				
				for k = 1:10
					m_Y = randn(Ncurves,Npoints);
					Fseq(k) = GFigure('m_X',m_X,'m_Y',m_Y);					
				end
				
				my_F = GFigure('m_multiplot',Fseq,'ch_multiplotType','sequence');
				
				my_F.plot;
				
				
			end
			
		end
		
		function str_out = translate_string(str_in,c_translationTable)
			assert( iscell(c_translationTable)&&(size(c_translationTable,2)==2));
			for k = 1:size(c_translationTable,1)
				str_in = regexprep(str_in,c_translationTable{k,1},c_translationTable{k,2});
			end
			str_out = str_in;
			
		end
		
	end
		
end