from pathlib import Path
import subprocess,sys
base=Path(__file__).resolve().parent
for name in ['experiments.py','extended_experiments.py','archive_diagnostics.py','make_tables.py','verify_native.py']:
    subprocess.run([sys.executable,str(base/name)],check=True)
print('All experiments, figures, tables and native verification completed')
