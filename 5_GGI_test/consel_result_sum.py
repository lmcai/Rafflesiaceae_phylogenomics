import os,fnmatch

filename=[]
for fn in os.listdir('./'):
    if fnmatch.fnmatch(fn,'*.CONSEL.result'):filename.append(fn)


output=open('AU_sum.tsv','a')
output.write('\t'.join(['Gene','H1','H2','H4','H5','H7','H9','H12','H15','H16','H17'])+'\n')

output=open('WSH_sum.tsv','a')
output.write('\t'.join(['Gene','H1','H2','H4','H5','H7','H9','H12','H15','H16','H17'])+'\n')

au_results={}
wsh_result={}

for fn in filename:
	x=open(fn).readlines()
	gene=fn.split('.')[0]
	a1=[]
	for i in os.listdir(gene+'.CONSEL'):
		if fnmatch.fnmatch(i,gene+'.H??.tre'):a1.append(i)
	a2=[]
	for i in os.listdir(gene+'.CONSEL'):
		if fnmatch.fnmatch(i,gene+'.H?.tre'):a2.append(i)
	a1.sort()
	a2.sort()
	y=a1+a2
	y=[i.split('.')[1] for i in y]
	res={}
	for l in x[3:-1]:
		#e.g. res['H2']=0.98
		res[y[int(l.split()[2])-1]]=l.split()[4]
	au_results[gene]=res


for key in au_results.keys():
	result=key+'\t'
	for i in ['H1','H2','H4','H5','H7','H9','H12','H15','H16','H17']:
		try:
			result=result+au_results[key][i]+'\t'
		except KeyError:
			result=result+'NA'+'\t'
	output.write(result+'\n')


output.close()
