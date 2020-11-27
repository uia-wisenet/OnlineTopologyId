classdef CopyHandle < handle
	
	methods
		function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this));
 
            % Copy all non-hidden properties.
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
		end
	end
end