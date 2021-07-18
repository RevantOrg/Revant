% [UTILITY OCTAVE SCRIPT] Plots several statistics on consensus strings.
%
% Remark: the script assumes the following string variables to be already initialized:
%
% consensiStats: directory containing statistics files for consensi;
% dominatorStats: directory containing statistics files for dominators;
% consensiStrings: directory containing consensus strings.
%

% What to print
printExcessSurface=0;
printCorrelations=0;
printQualityOfExcessSurface=0;
printExcessSurfaceIntersectTrack=0;
printCoverageFromWholeBlock=0;
printUtility=0;
errorRateHistogramsNormalized=0;


minAlignmentLength1=500;
minAlignmentLength2=1500;




% ----------------------------------- ERROR RATE -----------------------------------------
currentFigure=1;
figure(currentFigure); hold on;
if errorRateHistogramsNormalized==1
	A=load(sprintf('%s/errorRateHistogram.txt',dominatorStats));
	plot(A(:,1),A(:,2)./sum(A(:,2)),'.-k');
	nFRAlignments=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram.txt',consensiStats));
	plot(A(:,1),A(:,2)./sum(A(:,2)),'.-b');
	nFCAlignments1=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram-shortPeriod.txt',consensiStats));
	plot(A(:,1),A(:,2)./sum(A(:,2)),'.-r');
	nFCAlignments2=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram-longPeriod.txt',consensiStats));
	plot(A(:,1),A(:,2)./sum(A(:,2)),'.-g');
	nFCAlignments3=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram-unknownPeriod.txt',consensiStats));
	plot(A(:,1),A(:,2)./sum(A(:,2)),'.-m');
	nFCAlignments4=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram-cyclic.txt',consensiStats));
	plot(A(:,1),A(:,2)./sum(A(:,2)),'.-c');
	nFCAlignments5=sum(A(:,2));
	xlabel('Error rate'); ylabel('Fraction of alignments');
else
	A=load(sprintf('%s/errorRateHistogram.txt',dominatorStats));
	plot(A(:,1),A(:,2),'.-k');
	nFragmentReferenceAlignments=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram.txt',consensiStats));
	plot(A(:,1),A(:,2),'.-b');
	nFCAlignments1=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram-shortPeriod.txt',consensiStats));
	plot(A(:,1),A(:,2),'.-r');
	nFCAlignments2=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram-longPeriod.txt',consensiStats));
	plot(A(:,1),A(:,2),'.-g');
	nFCAlignments3=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram-unknownPeriod.txt',consensiStats));
	plot(A(:,1),A(:,2),'.-m');
	nFCAlignments4=sum(A(:,2));
	A=load(sprintf('%s/errorRateHistogram-cyclic.txt',consensiStats));
	plot(A(:,1),A(:,2),'.-c');
	nFCAlignments5=sum(A(:,2));
	xlabel('Error rate'); ylabel('N. of alignments');
endif
title(sprintf('Fragment-consensus alignments (same repeat, >=%d bp alignments)',minAlignmentLength1));
legend('FD','FC, non-p.','FC short p.','FC long p.','FC unknown p.','Cyclic','location','northwest');
grid on;

nFragmentReferenceAlignments
nFragmentConsensusAlignment=nFCAlignments1+nFCAlignments2+nFCAlignments3+nFCAlignments4+nFCAlignments5




% -------------------------------------- LENGTH ------------------------------------------
currentFigure=currentFigure+1;
figure(currentFigure);

subplot(1,2,1); hold on;
A=load(sprintf('%s/lengthHistogram.txt',consensiStats)); plot(A(:,1),A(:,2),'.b');
plot(A(:,1),A(:,1),'-k');
A=load(sprintf('%s/lengthHistogram-longPeriod.txt',consensiStats)); plot(A(:,1),A(:,2),'.g');
A=load(sprintf('%s/lengthHistogram-unknownPeriod.txt',consensiStats)); plot(A(:,1),A(:,2),'.m');
xlabel('Dominator length'); ylabel('Consensus length');
title('Effect of consensus on reference length');
grid on; axis square;
legend('Nonperiodic','identity','Long period','Unknown period','location','northwest');

subplot(1,2,2); hold on;
N_BINS=10;  % Arbitrary
A=load(sprintf('%s/lengthHistogram.txt',consensiStats)); 
[Y,X]=hist(100*(A(:,1)-A(:,2))./A(:,1),N_BINS); plot(X,Y./sum(Y),'ob');
A=load(sprintf('%s/lengthHistogram-longPeriod.txt',consensiStats)); 
[Y,X]=hist(100*(A(:,1)-A(:,2))./A(:,1),N_BINS); plot(X,Y./sum(Y),'og');
A=load(sprintf('%s/lengthHistogram-unknownPeriod.txt',consensiStats)); 
[Y,X]=hist(100*(A(:,1)-A(:,2))./A(:,1),N_BINS); plot(X,Y./sum(Y),'om');
xlabel('(|Reference| - |Consensus|) / |Reference| (%)');
ylabel('Fraction of modules'); title('Delta lengths (10 bins)');
grid on; axis square;
legend('Nonperiodic','Long period','Unknown period','location','northwest');




% --------------------------------- FALSE POSITIVES --------------------------------------
errorRates=[0.3,0.23,0.2,0.17];
longPeriod=1000;

% --------------------------------- EXCESS SURFACE ---------------------------------------
if printExcessSurface==1
	currentFigure=currentFigure+1;
	figure(currentFigure);
	A=load(sprintf('%s/blockReferenceHistogram.txt',consensiStats));
	[nRows,nColumns]=size(A);
	B=load(sprintf('%s/blockReferenceHistogram-periodic.txt',consensiStats));
	% 0=unknown period; -1=non-periodic, cyclic; -2=non-periodic.
	% Same basin
	for j=[0:length(errorRates)-1]
		for i=[1:nRows]
			total=A(i,9*j+2)+A(i,9*j+3);  % Basin surface
			for k=[1:9]
				A(i,9*j+k)=A(i,9*j+k)./total;
			endfor
		endfor
		subplot(3,length(errorRates),j+1); hold on;
		X=find(B==-2); plot(log10( A(X,9*j+1) ),A(X,9*j+3)*100,'.b');
		X=find(B<longPeriod & B>0); plot(log10( A(X,9*j+1) ),A(X,9*j+3)*100,'.r');
		X=find(B>=longPeriod); plot(log10( A(X,9*j+1) ),A(X,9*j+3)*100,'.g');
		X=find(B==0); plot(log10( A(X,9*j+1) ),A(X,9*j+3)*100,'.m');
		X=find(B==-1); plot(log10( A(X,9*j+1) ),A(X,9*j+3)*100,'.c');
		xlabel('|A \ Basin| / |Basin| (log10)'); ylabel('|Basin \ A| / |Basin| * 100');
		title(sprintf('Extra surface WRT basin (>=%d bp alignments, <=%f error)',minAlignmentLength2,errorRates(j+1)));
		legend('Nonperiodic','Short period','Long period');
		axis([-3,3,0,100]);
		grid on; axis square;
	endfor
	% All basins
	for j=[0:length(errorRates)-1]
		subplot(3,length(errorRates),length(errorRates)+j+1); hold on;
		X=find(B==-2); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+4) ),'.b');
		X=find(B<longPeriod & B>0); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+4) ),'.r');
		X=find(B>=longPeriod); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+4) ),'.g');
		X=find(B==0); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+4) ),'.m');
		X=find(B==-1); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+4) ),'.c');
		xlabel('|A \ Basin| / |Basin| (log10)');
		ylabel('|A \ All-basins| / |Basin| (log10)');
		title(sprintf('Extra surface WRT all basins (>=%d bp alignments, <=%f error)',minAlignmentLength2,errorRates(j+1)));
		%legend('Nonperiodic','Short period','Long period');
		axis([-3,3,-4,1]);
		grid on; axis square;
	endfor
	% All intervals
	for j=[0:length(errorRates)-1]
		subplot(3,length(errorRates),2*length(errorRates)+j+1); hold on;
		X=find(B==-2); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+7) ),'.b');
		X=find(B<longPeriod & B>0); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+7) ),'.r');
		X=find(B>=longPeriod); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+7) ),'.g');
		X=find(B==0); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+7) ),'.m');
		X=find(B==-1); plot(log10( A(X,9*j+1) ),log10( A(X,9*j+7) ),'.c');
		xlabel('|A \ Basin| / |Basin| (log10)');
		ylabel('|A \ Intervals| / |Basin| (log10)');
		title(sprintf('Extra surface WRT all intervals (>=%d bp alignments, <=%f error)',minAlignmentLength2,errorRates(j+1)));
		%legend('Nonperiodic','Short period','Long period');
		axis([-3,3,-4,1]);
		grid on; axis square;
	endfor

%find( A(:,7)>=1 )
%fabio=find(A(:,3)>=0.07);
%[fabio,A(fabio,:)]
endif



% ------------------------------------ CORRELATIONS --------------------------------------
if printCorrelations==1
	currentFigure=currentFigure+1;
	figure(currentFigure);

	maxError=0.4; minError=0; nBins=100;
	quantum=(maxError-minError)/nBins;
	errorRateMatrix=load(sprintf('%s/errorRateMatrix.txt',consensiStats));
	errorRateID=0;  % Position in array $errorRates$.
	[nRows,nColumns]=size(errorRateMatrix);
	X=zeros(nColumns,1);
	for i=[1:nColumns]
		X(i)=minError+quantum*(i-1);
	endfor
	descriptorIDs=textread(sprintf('%s/blockReferenceHistogram-ids.txt',consensiStats),'%s');
	for i=[1:length(descriptorIDs)]
		try
			consensiLengths=load(sprintf('%s/%s.daccord.consensus-lengths.txt',consensiStrings,strrep(descriptorIDs{i},'-','_')));
		catch
			sprintf('Error?! %s/%s.daccord.consensus-lengths.txt',consensiStrings,descriptorIDs{i})
			continue
		end_try_catch
		totalLength=sum(consensiLengths);
		avgErrorRate=(errorRateMatrix(i,:)./sum(errorRateMatrix(i,:)))*X;
		fractionOfUncoveredBasin=A(i,9*errorRateID+3)*100;
		extraSurface=log10( A(i,9*errorRateID+1) );  % Not belonging to any interval from factorization
	
		subplot(4,4,1); hold on; plot(log10(totalLength),log10(totalLength),'.b');
		subplot(4,4,2); hold on; plot(avgErrorRate,log10(totalLength),'.b');
		subplot(4,4,3); hold on; plot(fractionOfUncoveredBasin,log10(totalLength),'.b');
		subplot(4,4,4); hold on; plot(extraSurface,log10(totalLength),'.b');
	
		subplot(4,4,5); hold on; plot(log10(totalLength),avgErrorRate,'.b');
		subplot(4,4,6); hold on; plot(avgErrorRate,avgErrorRate,'.b');
		subplot(4,4,7); hold on; plot(fractionOfUncoveredBasin,avgErrorRate,'.b');
		subplot(4,4,8); hold on; plot(extraSurface,avgErrorRate,'.b');
	
		subplot(4,4,9); hold on; plot(log10(totalLength),fractionOfUncoveredBasin,'.b');
		subplot(4,4,10); hold on; plot(avgErrorRate,fractionOfUncoveredBasin,'.b');
		subplot(4,4,11); hold on; plot(fractionOfUncoveredBasin,fractionOfUncoveredBasin,'.b');
		subplot(4,4,12); hold on; plot(extraSurface,fractionOfUncoveredBasin,'.b');
	
		subplot(4,4,13); hold on; plot(log10(totalLength),extraSurface,'.b');
		subplot(4,4,14); hold on; plot(avgErrorRate,extraSurface,'.b');
		subplot(4,4,15); hold on; plot(fractionOfUncoveredBasin,extraSurface,'.b');
		subplot(4,4,16); hold on; plot(extraSurface,extraSurface,'.b');
	endfor
	for i=[1:16]
		subplot(4,4,i); axis square; grid on; 
	endfor
	subplot(4,4,1); ylabel('C total length (log10)');
	subplot(4,4,5); ylabel('Avg. FC error rate');
	subplot(4,4,9); ylabel('Fraction of uncovered basin (%)');
	subplot(4,4,13); ylabel('Excess surface (log10)');
	subplot(4,4,13); xlabel('C total length (log10)');
	subplot(4,4,14); xlabel('Avg. FC error rate');
	subplot(4,4,15); xlabel('Fraction of uncovered basin (%)');
	subplot(4,4,16); xlabel('Excess surface (log10)');
endif




% ------------------------------ QUALITY OF EXCESS SURFACE -------------------------------
if printQualityOfExcessSurface==1
	currentFigure=currentFigure+1;
	figure(currentFigure);
	nQualityBins=52;
	C=load(sprintf('%s/intervalQualities.txt',consensiStats));
	C=C./sum(C);
	A=load(sprintf('%s/blockReferenceQualities.txt',consensiStats));
	[nRows,nColumns]=size(A);
	for i=[0:length(errorRates)-1]
		subplot(2,2,i+1); hold on;
		plot([0:nQualityBins-1],C,'-ko');
		% Non-periodic
		X=find(B==-2); 
		if length(X)>1
			totals=sum(A(X,:));
		else
			totals=A(X,:);
		endif
		if length(A(X,:))>0
			denominator=sum(totals(i*nQualityBins+1:i*nQualityBins+nQualityBins));
			plot([0:nQualityBins-1],totals(i*nQualityBins+1:i*nQualityBins+nQualityBins)./denominator,'-b.');
		endif
		% Short period
		X=find(B<longPeriod & B>0); 
		if length(X)>1
			totals=sum(A(X,:));
		else
			totals=A(X,:);
		endif
		if length(A(X,:))>0
			denominator=sum(totals(i*nQualityBins+1:i*nQualityBins+nQualityBins));
			plot([0:nQualityBins-1],totals(i*nQualityBins+1:i*nQualityBins+nQualityBins)./denominator,'-r.');
		endif
		% Long period
		X=find(B>=longPeriod);
		if length(X)>1
			totals=sum(A(X,:));
		else
			totals=A(X,:);
		endif
		if length(A(X,:))>0
			denominator=sum(totals(i*nQualityBins+1:i*nQualityBins+nQualityBins));
			plot([0:nQualityBins-1],totals(i*nQualityBins+1:i*nQualityBins+nQualityBins)./denominator,'-g.');
		endif
		% Unknown period
		X=find(B==0);
		if length(X)>1
			totals=sum(A(X,:));
		else
			totals=A(X,:);
		endif
		if length(A(X,:))>0
			denominator=sum(totals(i*nQualityBins+1:i*nQualityBins+nQualityBins));
			plot([0:nQualityBins-1],totals(i*nQualityBins+1:i*nQualityBins+nQualityBins)./denominator,'-m.');
		endif
		% Cyclic
		X=find(B==-1);
		if length(X)>1
			totals=sum(A(X,:));
		else
			totals=A(X,:);
		endif
		if length(A(X,:))>0
			denominator=sum(totals(i*nQualityBins+1:i*nQualityBins+nQualityBins));
			plot([0:nQualityBins-1],totals(i*nQualityBins+1:i*nQualityBins+nQualityBins)./denominator,'-c.');
		endif
		title(sprintf('Intrinsic quality of (A \\ Intervals) (>=%d bp alignments, <=%f error)',minAlignmentLength2,errorRates(i+1)));
		xlabel('phred score (lower is better)'); ylabel('Fraction of bases at score');
		grid on; axis square;
	endfor
	subplot(2,2,1); legend('Intervals','Nonperiodic','Short period','Long period','Unknown period','Cyclic','location','northeast');
endif




% ---------------------------- EXCESS SURFACE INTERSECT DUST -----------------------------
if printExcessSurfaceIntersectTrack==1
	currentFigure=currentFigure+1;
	figure(currentFigure);
	Z=load(sprintf('%s/trackIntersection-dust.txt',consensiStats));
	%Z=load(sprintf('%s/trackIntersection-tandem.txt',consensiStats));
	trackID='DUST'; %'TANDEM';

	nBins=100;
	for i=[0:length(errorRates)-1]
		subplot(2,2,i+1); hold on;
		[Y,X]=hist(100*(Z(:,4*i+1)./Z(:,4*i+4)),nBins);
		plot(X,Y./sum(Y),'g.');
		[Y,X]=hist(100*(Z(:,4*i+2)./Z(:,4*i+4)),nBins);
		plot(X,Y./sum(Y),'r.');
		[Y,X]=hist(100*(Z(:,4*i+3)./Z(:,4*i+2)),nBins);
		plot(X,Y./sum(Y),'k.');
		title(sprintf('Excess surface and %s (alignments with <=%f error)',trackID,errorRates(i+1)));
		xlabel('%'); ylabel('Fraction of kernels at %'); grid on; axis square;
		legend(sprintf('percent %s in excess surface',trackID),'percent low-quality in excess surface',sprintf('percent %s in low-quality excess surface',trackID));
	endfor
endif






% --------------------------- COVERAGE FROM THE WHOLE BLOCK ------------------------------
if printCoverageFromWholeBlock==1
	for i=[1:length(descriptorIDs)]
		try
			consensiLengths=load(sprintf('%s/%s.daccord.consensus-lengths.txt',consensiStrings,strrep(descriptorIDs{i},'-','_')));
		catch
			continue
		end_try_catch
		currentFigure=currentFigure+1;
		figure(currentFigure); hold on;
		for j=[1:length(consensiLengths)]
			try
				histogram=load(sprintf('%s/blockCoverage-%s-%d.txt',consensiStats,descriptorIDs{i},j-1));
			catch
				continue
			end_try_catch
			X=[1:consensiLengths(j)];
			plot(X./consensiLengths(j),histogram./sum(histogram),'-');
		endfor
	endfor
endif




% -------------------------------------- UTILITY -----------------------------------------
% Remark: utility is not seen to correlate with length.
%
if printUtility==1
	currentFigure=currentFigure+1;
	figure(currentFigure);

	nBins=20;
	utilitiesConsensus=load(sprintf('%s/utilities-consensus.txt',consensiStats));
	utilitiesRepeat=load(sprintf('%s/utilities-repeat.txt',consensiStats));
	consensiLengths=[];
	repeatLengths=zeros(1,length(descriptorIDs));
	for i=[1:length(descriptorIDs)]
		try
			newLengths=load(sprintf('%s/%s.daccord.consensus-lengths.txt',consensiStrings,strrep(descriptorIDs{i},'-','_')));
		catch
			newLengths=[0];
			continue
		end_try_catch
		consensiLengths=[consensiLengths;newLengths];
		repeatLengths(i)=sum(newLengths);
	endfor
	for r=[1:length(errorRates)]
		subplot(2,2,r); hold on;
		[Y,X]=hist(log10( utilitiesConsensus(r,:) ),nBins);
		plot(X,Y./sum(Y),'-b.');
		[Y,X]=hist(log10( utilitiesRepeat(r,:) ),nBins);
		plot(X,Y./sum(Y),'-r.');
		xlabel('utility bin (log10)'); ylabel('Fraction of elements');
		title(sprintf('Utility (error <=%f)',errorRates(r)));
		axis square; grid on;
	endfor
	subplot(2,2,1); legend('Consensi','Repeats','location','northwest');
endif










currentFigure=currentFigure+1;
figure(currentFigure);
subplot(1,2,1);
A=load(sprintf('%s/fragmentReferenceHistogram.txt',consensiStats));
plot(A(:,1)*100,A(:,2)*100,'ob');
xlabel('Fragments aligning only to consensus (%)'); 
ylabel('Fragments aligning to other consensi (%)');
title(sprintf('Specificity of consensi (same repeat, >=%d bp alignments)',minAlignmentLength1));
grid on; axis square;





% Plotting the error rate distribution of every consensus
% maxError=0.4; minError=0; nBins=100;
% quantum=(maxError-minError)/nBins;
% A=load(sprintf('%s/errorRateMatrix.txt',consensiStats));
% [nRows,nColumns]=size(A);
% X=zeros(nColumns,1);
% for i=[1:nColumns]
% 	X(i)=minError+quantum*(i-1);
% endfor
% for i=[1:nRows]
% 	figure(); plot(X,A(i,:)./sum(A(i,:)),'.-k');
% endfor