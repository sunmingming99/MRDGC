% userpath('\Users\venhe\Documents\MATLAB\top1000')


tic

%%%%%%%%%%%%%%%%数据导入%%%%%%%%%%%%%%%%
%gene_difforder = importdata('C:\Users\Sunshine\Desktop\test.txt') % Differentially expressed dataset
gene_difforder =  readtable('GENE_EXP.txt','Delimiter','\t','ReadRowNames',1);
disp('基因数据导入完毕');
gene_difforder(1:5,1:5);
[numGenes,numSample_Gene]=size(gene_difforder);


%TF_EXP = importdata('TF_EXP.txt')
TF_EXP = readtable('TF_EXP.txt','Delimiter','\t','ReadRowNames',1);
disp('转录因子数据导入完毕');
TF_EXP(1:5,1:5);
[numTFs,numSample_TF]=size(TF_EXP);

%miRNA_EXP <- importdata('miRNA_EXP.txt')
miRNA_EXP = readtable('miRNA_EXP.txt','Delimiter','\t','ReadRowNames',1);
disp('miRNA数据导入完毕');
miRNA_EXP(1:5,1:5);
[nummiRNAs,numSample_miRNA]=size(miRNA_EXP);

%%%需要根据自己的实验数据确定case和 control的值
numCase=778 ;
numControl=100;


case_new_gene = gene_difforder{:,1:numCase} ;  % case group   %{}取出相应的值%（）表示取出来的还是table类型
ctl_new_gene = gene_difforder{:,numCase+1:numSample_Gene} ; % control group
case_new_TF = TF_EXP{:,1:numCase};   % case group
ctl_new_TF = TF_EXP{:,numCase+1:numSample_TF};  % control group


genes_names = gene_difforder.Properties.RowNames;
TFs_names = TF_EXP.Properties.RowNames;
miRNAs_names = miRNA_EXP.Properties.RowNames;

%%%%%%%%%%%%%%%%%%%%计算TF与基因的相关性%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
disp('开始计算');
disp('开始计算基因-TF相关性');

s1=corr((case_new_gene)',(case_new_TF)','type', 'Pearson');
s2=corr((ctl_new_gene)',(ctl_new_TF)','type', 'Pearson');
s = abs(s1-s2);

str = array2table(s,'VariableNames',TFs_names,'RowNames',genes_names) ; %这里还不对

writetable(str,'gene-TF_cor.csv','WriteVariableNames',true,'WriteRowNames',true,'Delimiter',',');
disp('计算TF-GENE完成');




%%%%%%%%%%%%%%%%%%%%去掉miRNA和基因表达谱中不同的样本%%%%%%%%%%%%%%%%%%%%%%%%%%


samName_comm = intersect(gene_difforder.Properties.VariableNames(1:numSample_Gene),miRNA_EXP.Properties.VariableNames(1:numSample_miRNA));
samName_comm_case = intersect(samName_comm,gene_difforder.Properties.VariableNames(1:numCase));
samName_comm_ctl = intersect(samName_comm,gene_difforder.Properties.VariableNames(numCase+1:numSample_Gene));


case_new_gene1 = gene_difforder{:,gene_difforder.Properties.VariableNames(samName_comm_case)};      % case group
ctl_new_gene1 = gene_difforder{:,gene_difforder.Properties.VariableNames(samName_comm_ctl)} ;   % control group
case_new_miRNA1= miRNA_EXP{:,miRNA_EXP.Properties.VariableNames(samName_comm_case)}  ;  % case group
ctl_new_miRNA1 =miRNA_EXP{:,miRNA_EXP.Properties.VariableNames(samName_comm_ctl)} ; % control group

disp('数据处理完成');

%%%%%%%%%%%%%%%%%%%%计算miRNA与基因的相关性%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('开始计算miRNA-gene相关性');    
s11=corr((case_new_gene1)',(case_new_miRNA1)','type', 'Pearson');
s21=corr((ctl_new_gene1)',(ctl_new_miRNA1)','type', 'Pearson');
    
sw = abs(s11-s21);

t_sw = array2table(sw,'VariableNames',miRNAs_names,'RowNames',genes_names);
writetable(t_sw,'gene-miRNA_cor.csv','WriteVariableNames',true,'WriteRowNames',true,'Delimiter',',');
disp('计算miRNA-gene完成');


gene_TF_miRNA = [s,sw];


[Y,new_gene_TF_miRNA] = sort(gene_TF_miRNA,2,'descend');% DIM=2表示对矩阵的各行中的元素排列，MODE默认升序，
csvwrite( 'new_use_gene_TF_miRNA.csv',new_gene_TF_miRNA) ; 

[nrow,ncol]=size(new_gene_TF_miRNA);
L=reshape(new_gene_TF_miRNA,nrow*ncol,1);
csvwrite( 'col_sort_L.csv',L);



%%%%%%%%%%%%% 计算每个调控因子的富集度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%并行计算方法

CoreNum=4; %设定机器CPU核数，我用的服务器，服务器的核为12，所以在这里CoreNum=12

% if parpool('size')<=0 %判断并行计算环境是否已然启动 
%     parpool('open','local',CoreNum); %若尚未启动，则启动并行环境 
% else
%     disp('Already initialized'); %说明并行环境已经启动。
% end
 
regulators=zeros(ncol,1);
parfor regulator = 1:ncol
%   disp('calculate i th, i=');
%   disp(regulator);
 
  data1 = find(L==regulator);
  data2 = find(L~=regulator);
  
  [pvalue,H] = ranksum(data1,data2);% 可加alpha参数
  regulators(regulator) = pvalue;
end

TF_miRNA_names = [TFs_names;miRNAs_names];

t_regulators = array2table(regulators,'VariableNames',{'pvalue'},'RowNames',TF_miRNA_names);
writetable(t_regulators,'regulators_pvalue.csv','WriteVariableNames',true,'WriteRowNames',true,'Delimiter',',');



%%%%%%%%%%%%%根据p-value升序排列，得到关键调控因子%%%%%%%%%%%%%%%%%%%%%%%%
res_candite = sort(regulators);
csvwrite( 'res_candite_totalgenes.csv',res_candite);

toc
disp(['运行时间: ',num2str(toc)]);
csvwrite( 'runtime.csv',num2str(toc));

%%%%关闭并行池
delete(gcp('nocreate'));


save data;