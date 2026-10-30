"""Run all experiments in a new output directory without overwriting recorded data."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def main():
    target = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'results/reproduced'
    if target.exists():
        raise SystemExit('Output directory already exists. Choose a new path to preserve prior results.')
    target.mkdir(parents=True)
    for source in (ROOT/'code').glob('*.py'):
        shutil.copy2(source,target/source.name)
    shutil.copy2(ROOT/'requirements.txt',target/'requirements.txt')
    env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    print('Running all stages; log:',target/'reproduction.log',flush=True)
    with (target/'reproduction.log').open('w') as log:
        subprocess.run([sys.executable,str(target/'reproduce.py')],cwd=target,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    subprocess.run([sys.executable,str(ROOT/'verify_results.py'),str(ROOT/'results/recorded'),str(target)],check=True)
    print('Completed:',target)

if __name__=='__main__':
    main()
