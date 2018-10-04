import sys, os

def main():
   i=0
#    for kern in ['linear','rbf']:
   for kern in ['tree100','tree4k','linear','rbf']:
       for rseed in [0]:
#            for icv in range(5):
           for icv in range(1):
               for fset in ['epw','ep','w']:
                   for shuf in [0,1]:
                       ofile=open('q'+str(i),'w')
                       ofile.write('#!/bin/sh\n')
                       ofile.write('python3 reproduce1_k562.py {0} {1} {2} {3} {4}\n'.format(icv,fset,shuf,kern,rseed))
                       ofile.close()
                       os.system('sbatch q'+str(i))
                       i=i+1

main()
