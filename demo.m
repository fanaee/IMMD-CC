clear;
load dataset
% data from  Padilha, V.A., Campello, R.J.: A systematic comparative evaluation of biclustering techniques. BMC Bioinformatics 18(1), 55 (2017)
% 5 co-clusters
[ba,noclust]=IMMDCC(x,0.3,0.3,2,Inf,9);
