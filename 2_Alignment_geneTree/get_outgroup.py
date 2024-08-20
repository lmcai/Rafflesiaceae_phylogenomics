import sys
x=open(sys.argv[1]).readlines()
x=[l[1:].strip() for l in x if l.startswith('>')]
outgr=[]
for candidate in ['Vitis','Elaeocarpus','Crossopetalum','Oxalis']:
	outgr = outgr + [leaf for leaf in x if leaf.startswith(candidate)]
print(outgr[0])
