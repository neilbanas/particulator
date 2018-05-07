function P = par_concatSteps(A,fields);

% P = par_concatSteps(steps);
% P = par_concatSteps(filebasename);
%                              ..., whichFields);
%
% concatenates steps variable by variable into a single structure P that's 
% easier to analyze.
% either takes the steps from a vector of structures, or from a sequence of
% files filebasename1.mat, filebasename2.mat,...

if isstruct(A) % series of steps stored in a variable --------------------------
	N = length(A);
	if nargin < 2
		fields = fieldnames(A(1));
	end
	for i=1:length(fields)
		if ~strcmpi(fields{i},'profiles')		
			P.(fields{i}) = nan.*ones([N size(A(1).(fields{i}))]);
		end
	end
	for n=1:N
		for i=1:length(fields)
			if ~strcmpi(fields{i},'profiles')
				P.(fields{i})(n,:) = A(n).(fields{i})(:);
			else
				pro = fieldnames(A(n).profiles);
				for j=1:length(pro)
					P.profiles.(pro{j})(n,:,:) = A(n).profiles.(pro{j});
				end
			end
		end
	end
	
elseif ischar(A) % series of files ---------------------------------------------
	thefiles = dir([A '*.mat']);
	filenames = {thefiles.name};
	if isempty(filenames)
		warning(['no files found at ' A '*.mat']);
	end
	N = length(filenames);
	load(filenames{1},'step');
	if nargin < 2
		fields = fieldnames(step);
	end
	for i=1:length(fields)
		if ~strcmpi(fields{i},'profiles')
			P.(fields{i}) = nan.*ones([N size(step.(fields{i}))]);
		end
	end	
	for n=1:N
		load(filenames{n},'step');
		for i=1:length(fields)
			if ~strcmpi(fields{i},'profiles')
				P.(fields{i})(n,:) = step.(fields{i})(:);
			else
				pro = fieldnames(step.profiles);
				for j=1:length(pro)
					P.profiles.(pro{j})(n,:,:) = step.profiles.(pro{j});
				end
			end
		end
	end	
	
else
	error('neither a structure nor a string.');
end