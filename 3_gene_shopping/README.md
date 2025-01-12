# Gene tree phylogenetic property summarization

Each gene was quantified for the following phylogenetic properties: number of species, alignment length, average gene tree support, total tree length, root-to-tip branch-length variance, and gene tree–species tree congruence.

1. To get basic metrics like alignment length, number of species, and compositional bias, many of which are available from iqtree log files, apply the python script `compositional_bias_from_iqtree_log.py`. This will generate a `na_aln_chisquare_compositional_bias.tsv` spreadsheet.

2. To get Rafflesiaceae and Apodanthaceae specific metrics (e.g., branch length, node support), follow the steps below:

   a. Root all gene trees using Notung with the following command:
```
java -jar ~/programs/Notung-2.8.1.3-beta.jar -s sp.tre -b gene.list --root --treeoutput newick --speciestag prefix --silent --progressbar --outputdir rooted_geneTr --nolosses
```
**Caveat:** The terminal branch length will be removed. This is only good for getting the correct branch support

   b. Use the python script `geneTr_reroot_branchlen.py` to extract branch length related metrics in `G2141.node_support_branlen.tsv`

3. To obtain root-to-tip branch-length variance, and gene tree–species tree congruence for free-living species using [SortaDate](https://github.com/FePhyFoFum/SortaDate), the following steps were used:

   a. Remove all Rafflesiaceae and Apodanthaceae in the gene trees before calculating the gene tree species tree congruence (such that the choice of the species tree won't biasedly impact the score)

   b. Run SorteDate using `sortdate.sh`

4. To explore the pairwise correlation between tree-wise and lineage-specific phylogenetic properties, consolidate all features into one spreadsheet (Table S2 in the manuscript). Then apply the R script `Gene_feature_pairwise_comparison_plot_v2.R` to test for correlation and plot the result.
