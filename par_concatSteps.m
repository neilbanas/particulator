function P = par_concatSteps(steps);

% P = par_concatSteps(steps);
% P = par_concatSteps(filebasename);
%
% concatenates steps variable by variable into a single structure P that's 
% easier to analyze.
% either takes the steps from a vector of structures, or from a sequence of
% files filebasename1.mat, filebasename2.mat,...
%
% TODO: the file-based version!

N = length(steps);
fields = fieldnames(steps(1));
for i=1:length(fields)
	P.(fields{i}) = nan.*ones([N size(steps(1).(fields{i}))]);
end
for n=1:N
	for i=1:length(fields)
		P.(fields{i})(n,:) = steps(n).(fields{i})(:);
	end
end