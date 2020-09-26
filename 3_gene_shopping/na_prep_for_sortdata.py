from ete3 import Tree
from Bio import SeqIO
import os

for f in os.listdir('.'):
	t=Tree(f,format=1)
	sp=[leaf.name for leaf in t]
	raff=[i for i in sp if i.startswith(('Sap','Rca','Rhi','Rtu'))]
	if len(raff)>0:
		seqs=SeqIO.index('../na_aln/'+f.split('.')[0]+'.na.aln.trimmed.fas','fasta')
		out=open(f.split('.')[0]+'.aln','a')
		for i in sp:
			d=SeqIO.write(seqs[i],out,'fasta')
		out.close()
	else:
		print('rm '+f)
		

rm 2298.inclade1.ortho1.tre
rm 1981.inclade1.ortho1.tre
rm 2931.inclade1.ortho1.tre
rm 796.inclade1.ortho1.tre
rm 3073.inclade1.ortho1.tre
rm 1351.inclade1.ortho1.tre
rm 1795.inclade1.ortho1.tre
rm 2781.inclade1.ortho1.tre
rm 1186.inclade1.ortho1.tre
rm 2122.inclade1.ortho1.tre
rm 2841.inclade1.ortho1.tre
rm 2562.inclade1.ortho1.tre
rm 2098.inclade1.ortho1.tre
rm 2559.inclade1.ortho1.tre
rm 2670.inclade1.ortho1.tre
rm 1882.inclade1.ortho1.tre
rm 1048.inclade1.ortho1.tre
rm 1060.inclade1.ortho1.tre
rm 143.inclade1.ortho1.tre
rm 2724.inclade1.ortho1.tre
rm 1605.inclade1.ortho1.tre
rm 2108.inclade1.ortho1.tre
rm 1007.inclade1.ortho1.tre
rm 2263.inclade1.ortho1.tre
rm 1545.inclade1.ortho1.tre
rm 2261.inclade1.ortho1.tre
rm 2022.inclade1.ortho1.tre
rm 2517.inclade1.ortho1.tre
rm 573.inclade1.ortho1.tre
rm 358.inclade1.ortho1.tre
rm 1302.inclade1.ortho1.tre
rm 1321.inclade1.ortho1.tre
rm 1829.inclade1.ortho1.tre
rm 585.inclade1.ortho1.tre
rm 557.inclade1.ortho1.tre
rm 809.inclade1.ortho1.tre
rm 2715.inclade1.ortho1.tre
rm 2997.inclade1.ortho1.tre
rm 3044.inclade1.ortho1.tre
rm 3067.inclade1.ortho1.tre
rm 1586.inclade1.ortho1.tre
rm 1328.inclade1.ortho1.tre
rm 1794.inclade1.ortho1.tre
rm 1201.inclade1.ortho1.tre
rm 1313.inclade1.ortho1.tre
rm 2335.inclade1.ortho1.tre
rm 1540.inclade1.ortho1.tre
rm 3049.inclade1.ortho1.tre
rm 3087.inclade1.ortho1.tre
rm 2291.inclade1.ortho1.tre
rm 2306.inclade1.ortho1.tre
rm 2656.inclade1.ortho1.tre
rm 1026.inclade1.ortho1.tre
rm 533.inclade1.ortho1.tre
rm 1778.inclade1.ortho1.tre
rm 2682.inclade1.ortho1.tre
rm 2597.inclade1.ortho1.tre
rm 980.inclade1.ortho1.tre
rm 2976.inclade1.ortho1.tre
rm 572.inclade1.ortho1.tre
rm 2678.inclade1.ortho1.tre
rm 2328.inclade1.ortho1.tre
rm 2869.inclade1.ortho1.tre
rm 2259.inclade1.ortho1.tre
rm 1772.inclade1.ortho1.tre
rm 2317.inclade1.ortho1.tre
rm 2265.inclade1.ortho1.tre
rm 1285.inclade1.ortho1.tre
rm 2520.inclade1.ortho1.tre
rm 2522.inclade1.ortho1.tre
rm 1531.inclade1.ortho1.tre
rm 581.inclade1.ortho1.tre
rm 2426.inclade1.ortho1.tre
rm 2856.inclade1.ortho1.tre
rm 2297.inclade1.ortho1.tre
rm 784.inclade1.ortho1.tre
rm 1039.inclade1.ortho1.tre
rm 2681.inclade1.ortho1.tre