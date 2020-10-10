#run raxml on partially constrained tree
raxmlHPC-SSE3 -p12345 -g Eupho_Raff.tre -m GTRGAMMA -s 1002.aln -n Eupho_Raff
raxmlHPC-SSE3 -p12345 -g EuphoPera_Raff.tre -m GTRGAMMA -s 1002.aln -n EuphoPera_Raff
raxmlHPC-SSE3 -p12345 -g EuphoPeraPutra_Raff.tre -m GTRGAMMA -s 1002.aln -n EuphoPeraPutra_Raff
raxmlHPC-SSE3 -p12345 -g EuphoPeraPutraER_Raff.tre -m GTRGAMMA -s 1002.aln -n EuphoPeraPutraER_Raff

#estimate per site log-likelihood in raxml
cat RAxML_result.* > RAxML_result.all
raxmlHPC -f g -s 1002.aln.reduced -z RAxML_result.all -n raff.likelihood -m GTRGAMMA

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