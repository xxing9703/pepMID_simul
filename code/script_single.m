% a script version for single peptide query
clear
peptide='ACDDMIPWSHFACDDMIPWSHFACDDMIPWSHF';  %enter peptide here
%-----------
[mass,atoms]=pep2mass(peptide); %ordering: C,N,H,O,S
str={'13C','15N','2H','17O','33S','','','','18O','34S','','','','','36S'};
n=5;
A0=[0.98893,0.996337,0.99985,0.9975904,0.9502]; %abundance for M0
A1=[0.0110694,0.003663,0.00015,3.7387E-04,0.0074962]; %abundance type A for +1
B1=[0,0,0,0.0020358,0.042099]; %abundance type B for +2
C1=[0,0,0,0,0.00020529]; %abundance type C for +4

mzA=[1.00335,0.99703,1.00630,1.00422,0.99940]; %mz diff for +1
mzB=[0,0,0,2.00424,1.99580]; %mz diff for +2
mzC=[0,0,0,0,3.9950]; %mz diff for +4

%------------------------------------------

ab_0=prod(A0.^atoms(1:n)); %abundance for M0
% calculate abundance up to m
m=8;
for i=1:n
   for j=1:m
    ab_A(i,j)=ab_0*(A1(i)/A0(i))^j*nck(atoms(i),j);
    ab_B(i,j)=ab_0*(B1(i)/A0(i))^j*nck(atoms(i),j);
    ab_C(i,j)=ab_0*(C1(i)/A0(i))^j*nck(atoms(i),j);
    dmz_A(i,j)=mzA(i)*j;
    dmz_B(i,j)=mzB(i)*j;  
    dmz_C(i,j)=mzC(i)*j; 
   end
 end
M=[ab_A;ab_B;ab_C];ct=0;
N=[dmz_A;dmz_B;dmz_C];
for i=1:size(M,1)
  for j=1:size(M,2)
      if M(i,j)>0
       ct=ct+1;
       tb(ct).ab=M(i,j);
       tb(ct).str=[str{i},num2str(j)];
       tb(ct).dmz=N(i,j); 
       tb(ct).mz=N(i,j)+mass; 
       tb(ct).type=i;
      end
  end
end

[~,ind]=sort([tb.ab],'descend'); %sorting
tb=tb(ind);
cutoff=0.0000001;%apply cutoff
ind=find([tb.ab]>cutoff);  
tb=tb(ind);
tb_save1=tb;
%% ----tb append
ct=length(tb);
for i=1:length(tb)
    for j=i:length(tb)
        if isempty(intersect(tb(i).type,tb(j).type))
            tp=tb(i).ab*tb(j).ab/ab_0;
          if tp>cutoff
              ct=ct+1;
              tb(ct).ab=tp;
               str_in=[tb(i).str,' ',tb(j).str];
              tb(ct).str=strjoin( sort(strsplit(str_in,' ')));
              tb(ct).dmz=tb(i).dmz+tb(j).dmz; 
              tb(ct).mz=tb(i).dmz+tb(j).dmz+mass; 
              tb(ct).type=union(tb(i).type,tb(j).type);
          end              
        end
    end
end
 [~,ind]=unique({tb.str});  %remove dup
 tb=tb(ind);

 tb=decouple(atoms,tb); %decouple
 
 [~,ind]=find([tb.ab]>cutoff); %remove zero
 tb=tb(ind);
%% ----------
 [~,ind]=sort([tb.ab],'descend'); %sorting by abudance
tb=tb(ind);
tb_save2=tb;
%% ------------------
 [~,ind]=sort([tb.dmz],'ascend'); %sorting by dmz
tb=tb(ind);
tb_save3=tb;
%%
out(1).IsotopeNumber=0;
out(1).mz=mass;
out(1).pct=ab_0*100;
out(1).pctMax=100;
for i=1:round(max([tb.dmz]))
  out(i+1).IsotopeNumber=i;
  [~,ind]=find(round([tb.dmz])==i);
  out(i+1).mz=mean([tb(ind).mz]);
  out(i+1).mz=sum([tb(ind).mz].*[tb(ind).ab])/sum([tb(ind).ab]);
  out(i+1).pct=sum([tb(ind).ab])*100;
  for j=1:length(ind)
    tb(ind(j)).pct=tb(ind(j)).ab/max([tb(ind).ab])*100;
  end  
end
for i=1:round(max([tb.dmz]))+1
   out(i).pctMax=out(i).pct/max([out.pct])*100;
end

writetable(struct2table(out),'/results/output1.csv')
writetable(struct2table(tb),'/results/output2.csv')

%%
resolution=100000;
mz=[tb.mz];
intens=[tb.ab];
sampling=0.0005;
x=(mz(1)-2):sampling:(mz(end)+2); 
c=mz(1)/resolution;
y=ab_0*exp(-(x-mass).^2/c^2);  % M0 peak, gaussian shape
for i=1:length(mz)
    y=y+intens(i)*exp(-(x-mz(i)).^2/c^2);    
end
y=y/max(y)*100; %normalize
f=figure, plot(x,y,'-');
xlabel('m/z');

saveas(f, '../results/figure1.png')



