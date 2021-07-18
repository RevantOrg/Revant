% [UTILITY OCTAVE SCRIPT] Plots the coverage histogram of every consensus string, in a
% distinct figure.
%
% Remark: the script assumes variable $fileName$ to be already set to the path of a text
% file containing a list of coverage files to plot.
%
list=fopen(fileName);
str=fgetl(list);
nDirs=fskipl(list,Inf);
fclose(list);

list=fopen(fileName);
for i=1:nDirs
	str=fgetl(list);
	A=load(str);
	X=[1:length(A)];
	figure(i); plot(X,A,'-');
	title(sprintf('Module %s',str));
endfor
fclose(list);
