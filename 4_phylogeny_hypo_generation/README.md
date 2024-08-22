```
../ASTER-Linux/bin/astral4 -i G2141.aster.genetree.trees -o G2141.aster4.tre
../../ASTER-Linux/bin/wastral -i G2141.aster.genetree.trees -o G2141.wastral.support.tre -r 16 -s 16 -x 100 -n 0 --mode 2
../../ASTER-Linux/bin/wastral -i G2141.aster.genetree.trees -o G2141.wastral.brlen.tre -r 16 -s 16 -x 100 -n 0 --mode 3
../../ASTER-Linux/bin/wastral -i G2141.aster.genetree.trees -o G2141.wastral.hybrid.tre -r 16 -s 16 -x 100 -n 0

#full node annotation and polytomy test
java -jar ../../Astral/astral.5.7.8.jar -q G2141.aster4.tre -i G2141.aster.genetree.trees -t 2 -o G2141.astral_anno.tre

java -jar ../../Astral/astral.5.7.8.jar -q G2141.aster4.tre -i G2141.aster.genetree.trees -t 10 -o G2141.astral_polytomy.tre
```


MPEST rooted gene tree preparation
```
../../../RangerDTL/CorePrograms/OptRoot.linux -i input.tre -o test.tre
```
