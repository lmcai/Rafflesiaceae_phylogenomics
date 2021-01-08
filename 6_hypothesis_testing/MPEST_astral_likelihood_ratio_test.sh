#conduct constrained species tree search in astral

java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H1.astral.tre -t 2 -j ../../8_GGI/H1.ref.tre >astral.H1.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H2.astral.tre -t 2 -j ../../8_GGI/H2.ref.tre >astral.H2.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H3.astral.tre -t 2 -j ../../8_GGI/H3.ref.tre >astral.H3.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H4.astral.tre -t 2 -j ../../8_GGI/H4.ref.tre >astral.H4.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H5.astral.tre -t 2 -j ../../8_GGI/H5.ref.tre >astral.H5.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H6.astral.tre -t 2 -j ../../8_GGI/H6.ref.tre >astral.H6.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H7.astral.tre -t 2 -j ../../8_GGI/H7.ref.tre >astral.H7.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H8.astral.tre -t 2 -j ../../8_GGI/H8.ref.tre >astral.H8.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H9.astral.tre -t 2 -j ../../8_GGI/H9.ref.tre >astral.H9.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H10.astral.tre -t 2 -j ../../8_GGI/H10.ref.tre >astral.H10.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H11.astral.tre -t 2 -j ../../8_GGI/H11.ref.tre >astral.H11.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H12.astral.tre -t 2 -j ../../8_GGI/H12.ref.tre >astral.H12.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H13.astral.tre -t 2 -j ../../8_GGI/H13.ref.tre >astral.H13.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H14.astral.tre -t 2 -j ../../8_GGI/H14.ref.tre >astral.H14.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H15.astral.tre -t 2 -j ../../8_GGI/H15.ref.tre >astral.H15.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H16.astral.tre -t 2 -j ../../8_GGI/H16.ref.tre >astral.H16.log
java -jar /n/home08/lmcai/programs/Constrained-search/astral.5.6.9.jar -i G2141.gene.trees -o G2141.H17.astral.tre -t 2 -j ../../8_GGI/H17.ref.tre >astral.H17.log

#calculate tree likelihood in MPEST
mpest mpest.LL.ctl

#in Phybase, conduct pairwise likelihood ratio test
