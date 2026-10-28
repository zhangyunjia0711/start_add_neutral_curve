"""Isolated additional-frequency context; original matrices and CFD stay unchanged."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
import run_adapter as driver
driver.ROOT=ROOT;driver.CASE=ROOT/'cases'
run,table,sha=driver.run,driver.table,driver.sha
SOLVER=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028/work/solver_source/bin/START_LST')
SOLVER_HASH='e66c15b63c0a2dfaf74450e712e4c62f42a6a20d50b7612efc057855f3659893'
CASES=json.loads((ROOT/'outputs/cases.json').read_text())
def snapshot():
    data=json.loads((ROOT/'outputs/input_manifest.json').read_text())
    for r in data:assert sha(r['path'])==r['sha256'],r['path']
    return data
assert sha(SOLVER)==SOLVER_HASH
