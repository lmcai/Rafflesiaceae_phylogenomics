# gCF and sCF characterization

We used the species tree, gene trees, and DNA alignments from the coalescent-based Analysis 15 (table S3). 

For gCF:
```
iqtree2 -te G2141.wastral.hybrid.tre --gcf G2141.MFP.trees --prefix G2135.gCF
```
For sCF with the complete DNA alignment:
```
iqtree2 -te ../coalescent/G2141.wastral.hybrid.tre --scfl 200 -s G2135.cds.fas --prefix G2135.full.sCF
```
For sCF with three codon positions seperately
```
iqtree2 -te ../coalescent/G2141.wastral.hybrid.tre --scfl 200 -s G2135.codon1.fas --prefix G2135.codon1st.sCF
iqtree2 -te ../coalescent/G2141.wastral.hybrid.tre --scfl 200 -s G2135.codon2.fas --prefix G2135.codon2nd.sCF
iqtree2 -te ../coalescent/G2141.wastral.hybrid.tre --scfl 200 -s G2135.codon3.fas --prefix G2135.codon3rd.sCF
```

# Gene- and site-wise support for alternative topologies
