%% read aa and peptide measurement table
function Main_fitting(varargin)
folder_out='../results';
fn_aa=varargin{1};
fn_pep=varargin{2};
isotype=find(contains({'13C','15N','2H'},varargin{3})); % 13C, 15N or 2H
tb_aa=tb_fix(readtable(fn_aa)); % input curation
tb_pep=tb_fix(readtable(fn_pep)); % input curation
%%
for i=1:size(tb_pep,1)
  str_pep=tb_pep{i,1};
  str_pep=str_pep{1};
  M=[];
  for k=1:length(str_pep)
    id=find(strcmp(tb_aa{:,1},str_pep(k)));
    M(k,:)=tb_aa{id,2:end};
  end
  D_ab=mergeM(M);
  out0=peptide_mid(str_pep,8);
  mid0=[out0.pct]/100;   %MID for unlabeled natural abundance
  mid1=conv(mid0,D_ab);   %MID for labeled
  mid2=tb_pep{i,2:end};mid2=mid2/sum(mid2);  %MID for measured peptide
  s=min(length(mid0),length(mid1));s=min(s,length(mid2)); %find the shortest length of 3 mids  
  mid0=mid0(1:s); mid1=mid1(1:s); mid2=mid2(1:s); %trancate all 3 to the shortest length
  
  ct=0;para=[];error=[];
  for x=0:0.01:1  % x is turn over rate, parameter fitting
    ct=ct+1;
    fit=mid0*(1-x)+mid1*x;
    fit=fit/sum(fit);
    para(ct)=x;
    error(ct)=sqrt(sum((fit-mid2).^2)/length(fit));
  end  
  %figure,plot(para,error,'.-')
  [err,ind]=min(error);
  x=para(ind);
  fit=mid0*(1-x)+mid1*x;
  fit=fit/sum(fit);
  %figure,plot(fit,mid2,'*')
  
  %output to table
  output(i).peptide=str_pep;
  %output(i).mid_0=fit(1);
  %output(i).mid=fit(2:end);
  output(i).mid_natural_0=mid0(1);
  output(i).mid_natural=mid0(2:end); 
  output(i).mid_labeled_0=mid1(1);
  output(i).mid_labeled=mid1(2:end); 
  output(i).turnover=x;
  output(i).err=err;
   fprintf([num2str(i),'/',num2str(size(tb_pep,1)),' MID(natural) = ',num2str(mid0,'% .4f'),' MID(labeled) = ',num2str(mid1,'% .4f'),' turnover = ',num2str(x,'%.2f'),' ',str_pep]);
   fprintf('\n');
end
% export table to file "ex_*.csv"
T=struct2table(output);
[~,fname,ext]=fileparts(fn_pep);
writetable(T,fullfile(folder_out,['ex_',fname,ext]));


 