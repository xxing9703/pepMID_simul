% a prior version of batch process to simulate unlabeled peptide MID
% batch process peptide listed in the inputfile column 1. Export calculated natural isotopes
clear
inputfile='../data/1036peptidenatural.xlsx';
m=8;  % maximum M
T=readtable(inputfile);
n=size(T,1); 
for i=1:n
    i
    pep=T{i,1}{1};
    tp=peptide_iso(pep,m+1);
    tp=[tp.pct];
    tp=tp(1:m+1);
    output(i,:)=tp;
end
T1 = array2table(output,'VariableNames',cellstr(num2str([0:m]',['M','%01d'])));
writetable([T,T1],'../results/output10531.csv');