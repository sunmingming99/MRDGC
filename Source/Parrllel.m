% userpath('\Users\venhe\Documents\MATLAB\top1000')


tic

%%%%%%%%%%%%%%%%���ݵ���%%%%%%%%%%%%%%%%
%gene_difforder = importdata('C:\Users\Sunshine\Desktop\test.txt') % Differentially expressed dataset
gene_difforder =  readtable('GENE_EXP.txt','Delimiter','\t','ReadRowNames',1);
disp('�������ݵ������');
gene_difforder(1:5,1:5);
[numGenes,numSample_Gene]=size(gene_difforder);


%TF_EXP = importdata('TF_EXP.txt')
TF_EXP = readtable('TF_EXP.txt','Delimiter','\t','ReadRowNames',1);
disp('ת¼�������ݵ������');
TF_EXP(1:5,1:5);
[numTFs,numSample_TF]=size(TF_EXP);

%miRNA_EXP <- importdata('miRNA_EXP.txt')
miRNA_EXP = readtable('miRNA_EXP.txt','Delimiter','\t','ReadRowNames',1);
disp('miRNA���ݵ������');
miRNA_EXP(1:5,1:5);
[nummiRNAs,numSample_miRNA]=size(miRNA_EXP);

%%%��Ҫ�����Լ���ʵ������ȷ��case�� control��ֵ
numCase=778 ;
numControl=100;


case_new_gene = gene_difforder{:,1:numCase} ;  % case group   %{}ȡ����Ӧ��ֵ%������ʾȡ�����Ļ���table����
ctl_new_gene = gene_difforder{:,numCase+1:numSample_Gene} ; % control group
case_new_TF = TF_EXP{:,1:numCase};   % case group
ctl_new_TF = TF_EXP{:,numCase+1:numSample_TF};  % control group


genes_names = gene_difforder.Properties.RowNames;
TFs_names = TF_EXP.Properties.RowNames;
miRNAs_names = miRNA_EXP.Properties.RowNames;

%%%%%%%%%%%%%%%%%%%%����TF�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
disp('��ʼ����');
disp('��ʼ�������-TF�����');

s1=corr((case_new_gene)',(case_new_TF)','type', 'Pearson');
s2=corr((ctl_new_gene)',(ctl_new_TF)','type', 'Pearson');
s = abs(s1-s2);

str = array2table(s,'VariableNames',TFs_names,'RowNames',genes_names) ; %���ﻹ����

writetable(str,'gene-TF_cor.csv','WriteVariableNames',true,'WriteRowNames',true,'Delimiter',',');
disp('����TF-GENE���');




%%%%%%%%%%%%%%%%%%%%ȥ��miRNA�ͻ��������в�ͬ������%%%%%%%%%%%%%%%%%%%%%%%%%%


samName_comm = intersect(gene_difforder.Properties.VariableNames(1:numSample_Gene),miRNA_EXP.Properties.VariableNames(1:numSample_miRNA));
samName_comm_case = intersect(samName_comm,gene_difforder.Properties.VariableNames(1:numCase));
samName_comm_ctl = intersect(samName_comm,gene_difforder.Properties.VariableNames(numCase+1:numSample_Gene));


case_new_gene1 = gene_difforder{:,gene_difforder.Properties.VariableNames(samName_comm_case)};      % case group
ctl_new_gene1 = gene_difforder{:,gene_difforder.Properties.VariableNames(samName_comm_ctl)} ;   % control group
case_new_miRNA1= miRNA_EXP{:,miRNA_EXP.Properties.VariableNames(samName_comm_case)}  ;  % case group
ctl_new_miRNA1 =miRNA_EXP{:,miRNA_EXP.Properties.VariableNames(samName_comm_ctl)} ; % control group

disp('���ݴ������');

%%%%%%%%%%%%%%%%%%%%����miRNA�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('��ʼ����miRNA-gene�����');    
s11=corr((case_new_gene1)',(case_new_miRNA1)','type', 'Pearson');
s21=corr((ctl_new_gene1)',(ctl_new_miRNA1)','type', 'Pearson');
    
sw = abs(s11-s21);

t_sw = array2table(sw,'VariableNames',miRNAs_names,'RowNames',genes_names);
writetable(t_sw,'gene-miRNA_cor.csv','WriteVariableNames',true,'WriteRowNames',true,'Delimiter',',');
disp('����miRNA-gene���');


gene_TF_miRNA = [s,sw];


[Y,new_gene_TF_miRNA] = sort(gene_TF_miRNA,2,'descend');% DIM=2��ʾ�Ծ���ĸ����е�Ԫ�����У�MODEĬ������
csvwrite( 'new_use_gene_TF_miRNA.csv',new_gene_TF_miRNA) ; 

[nrow,ncol]=size(new_gene_TF_miRNA);
L=reshape(new_gene_TF_miRNA,nrow*ncol,1);
csvwrite( 'col_sort_L.csv',L);



%%%%%%%%%%%%% ����ÿ���������ӵĸ�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%���м��㷽��

CoreNum=4; %�趨����CPU���������õķ��������������ĺ�Ϊ12������������CoreNum=12

% if parpool('size')<=0 %�жϲ��м��㻷���Ƿ���Ȼ���� 
%     parpool('open','local',CoreNum); %����δ���������������л��� 
% else
%     disp('Already initialized'); %˵�����л����Ѿ�������
% end
 
regulators=zeros(ncol,1);
parfor regulator = 1:ncol
%   disp('calculate i th, i=');
%   disp(regulator);
 
  data1 = find(L==regulator);
  data2 = find(L~=regulator);
  
  [pvalue,H] = ranksum(data1,data2);% �ɼ�alpha����
  regulators(regulator) = pvalue;
end

TF_miRNA_names = [TFs_names;miRNAs_names];

t_regulators = array2table(regulators,'VariableNames',{'pvalue'},'RowNames',TF_miRNA_names);
writetable(t_regulators,'regulators_pvalue.csv','WriteVariableNames',true,'WriteRowNames',true,'Delimiter',',');



%%%%%%%%%%%%%����p-value�������У��õ��ؼ���������%%%%%%%%%%%%%%%%%%%%%%%%
res_candite = sort(regulators);
csvwrite( 'res_candite_totalgenes.csv',res_candite);

toc
disp(['����ʱ��: ',num2str(toc)]);
csvwrite( 'runtime.csv',num2str(toc));

%%%%�رղ��г�
delete(gcp('nocreate'));


save data;