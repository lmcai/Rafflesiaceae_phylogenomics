export ID=$1

#prepare constrained tree, remove missing taxa
python constrained_tree_prep.py [input_tree] EuphoPeraPutraPandaER.ref.tre
python constrained_tree_prep.py [input_tree] EuphoPeraPutraPanda.ref.tre
python constrained_tree_prep.py [input_tree] EuphoPeraPutra.ref.tre
python constrained_tree_prep.py [input_tree] EuphoPera.ref.tre
python constrained_tree_prep.py [input_tree] Eupho.ref.tre
python constrained_tree_prep.py [input_tree] Ixon.ref.tre
python constrained_tree_prep.py [input_tree] Passi.ref.tre

#run raxml on partially constrained tree
raxmlHPC-SSE3 -p12345 -g $ID.Eupho.tre -m GTRGAMMA -s $ID.aln -n $ID.Eupho
raxmlHPC-SSE3 -p12345 -g $ID.EuphoPera.tre -m GTRGAMMA -s 1002.aln -n $ID.EuphoPera
raxmlHPC-SSE3 -p12345 -g $ID.EuphoPeraPutra.tre -m GTRGAMMA -s 1002.aln -n $ID.EuphoPeraPutra
raxmlHPC-SSE3 -p12345 -g $ID.EuphoPeraPutraPanda.tre -m GTRGAMMA -s 1002.aln -n $ID.EuphoPeraPutraPanda
raxmlHPC-SSE3 -p12345 -g $ID.EuphoPeraPutraPandaER.tre -m GTRGAMMA -s 1002.aln -n $ID.EuphoPeraPutraPandaER
raxmlHPC-SSE3 -p12345 -g $ID.Ixon.tre -m GTRGAMMA -s 1002.aln -n $ID.Ixon
raxmlHPC-SSE3 -p12345 -g $ID.Passi.tre -m GTRGAMMA -s 1002.aln -n $ID.Passi

#estimate per site log-likelihood in raxml
cat RAxML_result.$ID.* > RAxML_result.$ID.all
raxmlHPC -f g -s $ID.aln.reduced -z RAxML_result.$ID.all -n $ID.likelihood -m GTRGAMMA

#RAxML_perSiteLLs.*

#convert raxml tree-puzzle output for CONSEL
/n/home08/lmcai/programs/consel/bin/makermt --puzzle -b 10 RAxML_perSiteLLs.raff.likelihood [ouput_prefix]

#get p-value
/n/home08/lmcai/programs/consel/bin/consel RAxML_perSiteLLs.raff [ouput_prefix]
/n/home08/lmcai/programs/consel/bin/catpv [ouput_prefix]

###############
#The following are hypotheses to be tested
#Euphorbiaceae
(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora,Drypetes,Clutia);
#Euphorbiaceae+Peraceae
(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha,Clutia)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora,Drypetes);
#Euphorbiaceae+Peraceae+Putranjivaceae
(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha,Clutia,Drypetes)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora);
#Euphorbiaceae+Peraceae+Putranjivaceae+Pandaceae
(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha,Clutia,Drypetes,Galearia)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Erythroxylum,Rhizophora);
#Euphorbiaceae+Peraceae+Putranjivaceae+Pandaceae+ER
(((Sap,Rhi,Rca,Rtu),(Ricinus,Hevea,Manihot,Endospermum,Jatropha,Clutia,Drypetes,Galearia,Erythroxylum,Rhizophora)),Passiflora,Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu);
#Passiflora
(((Sap,Rhi,Rca,Rtu),Passiflora),Ixonanthes,Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora,Drypetes,Clutia,Ricinus,Hevea,Manihot,Endospermum,Jatropha);
#Ixonanthes
(((Sap,Rhi,Rca,Rtu),Ixonanthes),Sauropus,Bischofia,Chrysobalanus,Linum,Ochna,Garcinia,Clusia,Mammea,Calophyllum,Hypericum,Podostemum,Galphimia,Tristellateia,Bergia,Elatine,Bhesa,Crossopetalum,Oxalis,Elaeocarpus,Viola,Rinorea,Malesherbia,Passiflora,Casearia,Flacourtia,Populus,Salix,SalixSu,Galearia,Erythroxylum,Rhizophora,Drypetes,Clutia,Ricinus,Hevea,Manihot,Endospermum,Jatropha);