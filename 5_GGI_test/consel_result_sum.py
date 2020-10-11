import os,fnmatch

filename=[]
for fn in os.listdir('./'):
    if fnmatch.fnmatch(fn,'*.CONSEL.result'):filename.append(fn)


au_output=open('AU_sum.tsv','a')
au_output.write('\t'.join(['Gene','H1','H2','H4','H5','H7','H9','H12','H15','H16','H17'])+'\n')

pp_output=open('PP_sum.tsv','a')
pp_output.write('\t'.join(['Gene','H1','H2','H4','H5','H7','H9','H12','H15','H16','H17'])+'\n')

sh_output=open('SH_sum.tsv','a')
sh_output.write('\t'.join(['Gene','H1','H2','H4','H5','H7','H9','H12','H15','H16','H17'])+'\n')

kh_output=open('KH_sum.tsv','a')
kh_output.write('\t'.join(['Gene','H1','H2','H4','H5','H7','H9','H12','H15','H16','H17'])+'\n')

au_results={}
pp_results={}
sh_results={}
kh_results={}

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
	au_res={}
	pp_res={}
	sh_res={}
	kh_res={}
	for l in x[3:-1]:
		#e.g. res['H2']=0.98
		au_res[y[int(l.split()[2])-1]]=l.split()[4]
		pp_res[y[int(l.split()[2])-1]]=l.split()[8]
		sh_res[y[int(l.split()[2])-1]]=l.split()[10]
		kh_res[y[int(l.split()[2])-1]]=l.split()[9]
	au_results[gene]=au_res
	pp_results[gene]=pp_res
	sh_results[gene]=sh_res
	kh_results[gene]=kh_res


for key in au_results.keys():
	result=key+'\t'
	for i in ['H1','H2','H4','H5','H7','H9','H12','H15','H16','H17']:
		try:
			result=result+au_results[key][i]+'\t'
		except KeyError:
			result=result+'NA'+'\t'
	au_output.write(result+'\n')
	result=key+'\t'
	for i in ['H1','H2','H4','H5','H7','H9','H12','H15','H16','H17']:
		try:
			result=result+pp_results[key][i]+'\t'
		except KeyError:
			result=result+'NA'+'\t'
	pp_output.write(result+'\n')
	result=key+'\t'
	for i in ['H1','H2','H4','H5','H7','H9','H12','H15','H16','H17']:
		try:
			result=result+sh_results[key][i]+'\t'
		except KeyError:
			result=result+'NA'+'\t'
	sh_output.write(result+'\n')
	result=key+'\t'
	for i in ['H1','H2','H4','H5','H7','H9','H12','H15','H16','H17']:
		try:
			result=result+kh_results[key][i]+'\t'
		except KeyError:
			result=result+'NA'+'\t'
	kh_output.write(result+'\n')

au_output.close()
pp_output.close()
sh_output.close()
kh_output.close()