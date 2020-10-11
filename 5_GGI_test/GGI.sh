#module load gcc/7.1.0-fasrc01 RAxML/8.2.11-fasrc01

export ID=$1

mkdir $ID.CONSEL
cd $ID.CONSEL
#prepare constrained tree, remove missing taxa
python ../constrained_tree_prep.py $ID 

#run raxml on partially constrained tree
raxmlHPC-SSE3 -p12345 -g $ID.H1.tre -m GTRGAMMA -s ../../4_na_aln_1to1/$ID.aln -n $ID.H1
raxmlHPC-SSE3 -p12345 -g $ID.H2.tre -m GTRGAMMA -s ../../4_na_aln_1to1/$ID.aln -n $ID.H2
raxmlHPC-SSE3 -p12345 -g $ID.H3.tre -m GTRGAMMA -s ../../4_na_aln_1to1/$ID.aln -n $ID.H3
raxmlHPC-SSE3 -p12345 -g $ID.H4.tre -m GTRGAMMA -s ../../4_na_aln_1to1/$ID.aln -n $ID.H4
raxmlHPC-SSE3 -p12345 -g $ID.H5.tre -m GTRGAMMA -s ../../4_na_aln_1to1/$ID.aln -n $ID.H5
raxmlHPC-SSE3 -p12345 -g $ID.H6.tre -m GTRGAMMA -s ../../4_na_aln_1to1/$ID.aln -n $ID.H6
raxmlHPC-SSE3 -p12345 -g $ID.H7.tre -m GTRGAMMA -s ../../4_na_aln_1to1/$ID.aln -n $ID.H7

rm RAxML_info.$ID.H*
rm RAxML_log.$ID.H*
rm RAxML_result.$ID.H*

#estimate per site log-likelihood in raxml
cat RAxML_bestTree.$ID.* > RAxML_result.$ID.all

if [ -f ../../4_na_aln_1to1/$ID.aln.reduced ]
then
	raxmlHPC -f g -s ../../4_na_aln_1to1/$ID.aln.reduced -z RAxML_result.$ID.all -n $ID.likelihood -m GTRGAMMA	
else
	raxmlHPC -f g -s ../../4_na_aln_1to1/$ID.aln -z RAxML_result.$ID.all -n $ID.likelihood -m GTRGAMMA
fi

rm RAxML_info.$ID.likelihood
#RAxML_perSiteLLs.*

#convert raxml tree-puzzle output for CONSEL
/n/home08/lmcai/programs/consel/bin/makermt --puzzle -b 10 RAxML_perSiteLLs.$ID.likelihood $ID

#get p-value
/n/home08/lmcai/programs/consel/bin/consel $ID
/n/home08/lmcai/programs/consel/bin/catpv $ID >../$ID.CONSEL.result

###############
#The following are hypotheses to be tested
#Euphorbiaceae
#(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora,Drypetes,Clutia);
#Euphorbiaceae+Peraceae
#(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha,Clutia)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora,Drypetes);
#Euphorbiaceae+Peraceae+Putranjivaceae
#(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha,Clutia,Drypetes)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora);
#Euphorbiaceae+Peraceae+Putranjivaceae+Pandaceae
#(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha,Clutia,Drypetes,Galearia)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Erythroxylum,Rhizophora);
#Euphorbiaceae+Peraceae+Putranjivaceae+Pandaceae+ER
#(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha,Clutia,Drypetes,Galearia,Erythroxylum,Rhizophora)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu);
#Passiflora
#(((Sap,Rhi,Rca,Rtu),Passiflora),Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora,Drypetes,Clutia,Ricinus,Hevea,Manihot,Endospermum,Jatropha);
#Ixonanthes
#(((Sap,Rhi,Rca,Rtu),Ixonanthes),Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Passiflora,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora,Drypetes,Clutia,Ricinus,Hevea,Manihot,Endospermum,Jatropha);
