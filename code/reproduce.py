from pathlib import Path
import subprocess,sys,os
root=Path(__file__).resolve().parent
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for script in ['study.py','extensions.py','test_checks.py','cache_screen.py','heldout.py','analyze.py','summarize_extra.py']:
    print('Running',script,flush=True)
    subprocess.run([sys.executable,str(root/script)],cwd=root,env=env,check=True)
print('All stages completed.')
