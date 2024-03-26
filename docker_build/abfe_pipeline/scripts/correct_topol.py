import sys

f = open(sys.argv[1]).read()
pattern1 = f'#include "{sys.argv[2]}"\n#endif\n'
pattern2 = f'#include "{sys.argv[2]}"\n#endif\n; Include ligand topology\n#include "{sys.argv[3]}"\n'

f = f.replace(pattern1, pattern2)
out=open(sys.argv[1], 'w')
out.write(f)
out.close()
