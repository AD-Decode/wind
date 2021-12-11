%make degree and clustering columns
addpath(genpath('/Users/alex/Documents/MATLAB/BCT/2019_03_03_BCT/'));


runno={'N57496',
'N57500',
'N57694',
'N57702',
'N57709',
'N57498',
'N57502',
'N57513',
'N57552',
'N57692',
'N57700',
'N57582',
'N57584',
'N57590',
'N57546',
'N57548',
'N57550',
'N57580',
'N57587',
'N57442',
'N57447',
'N57451',
'N57518',
'N57522',
'N57437',
'N57446',
'N57449',
'N57515',
'N57520',
'N57554',
'N57559'};

%read connectomes

mypathin='/Users/alex/AlexBadea_MyPapers/DavidDunson/CaseStudyConnectomics/Connectomes_all/' ; %N57496_all_connectomes.xlsx 
for i=1:numel(runno)
myfilein=char(strcat(mypathin,char(runno{i}),'_all_connectomes.xlsx'))
myconnectome=xlsread(myfilein);
 for j=1:size(myconnectome(:,1))
     myconnectome(j,j)=0;
     end
my_deg(i,:)=sum(myconnectome, 1);
my_clus(i,:)=clustering_coef_wu(myconnectome);
end

csvwrite('/Users/alex/AlexBadea_MyPapers/DavidDunson/CaseStudyConnectomics/Connectomes_all/Deg.csv', my_deg);
csvwrite('/Users/alex/AlexBadea_MyPapers/DavidDunson/CaseStudyConnectomics/Connectomes_all/Clus.csv', my_clus);
