function  gsimStartup

% moving to the folder that contains gsimStartup
ch_thisFile = which('gsimStartup');
ind = strfind(ch_thisFile,'gsimStartup');
ch_thisFolder = ch_thisFile(1:max(1,ind-1));
cd(ch_thisFolder);

if exist('gsim.m','file')
	movefile('gsim.m','gsim-old.m');
	disp('Your old gsim.m file is now renamed as gsim-old.m')
end
copyfile([ch_thisFolder 'Execution/gsim-template.m'],'./gsim.m')	

disp('GSim correctly initialized.')
disp('You can start working with GSim.');
disp('Example:    >> gsim(0)');

end

