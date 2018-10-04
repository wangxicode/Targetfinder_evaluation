import sys, os

def main():
   i=0
   for kern in ['tree100','tree4k','linear','rbf']:
               for fset in ['epw','ep','w']:
                       ofile=open('qgm'+str(i),'w')
                       ofile.write('#!/bin/sh\n')
                       ofile.write('#SBATCH --time=48:0:0\n\n')
                       ofile.write('python3 reproduce2_gm.py {0} {1}\n'.format(kern,fset))
                       ofile.close()
                       os.system('sbatch qgm'+str(i))
                       i=i+1

main()
