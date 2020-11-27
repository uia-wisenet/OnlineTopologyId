function str = print_time(tm)

str = sprintf('%.2f seconds ',tm);
tm = round(tm);
% days hours minutes seconds
days=floor(tm/(3600*24));
tm=tm-days*3600*24;
hours=floor(tm/3600);
tm=tm-hours*3600;
minutes=floor(tm/60);
tm=tm-minutes*60;
seconds=round(tm);
if days
	str = sprintf('%s( %d days, %d h, %d m and %d s )',str,days,hours,minutes,seconds);
elseif hours
	str = sprintf('%s( %d h, %d m and %d s )',str,hours,minutes,seconds);
elseif minutes
	str = sprintf('%s( %d m and %d s )',str,minutes,seconds);
end

if nargout==0
	fprintf('%s\n',str);
end
return
