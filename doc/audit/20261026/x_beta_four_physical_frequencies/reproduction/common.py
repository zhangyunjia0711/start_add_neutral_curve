"""Immutable, isolated fixed-omega run context."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
import run_adapter as driver
driver.ROOT=ROOT;driver.CASE=ROOT/'cases'
run,table,sha=driver.run,driver.table,driver.sha
SOLVER=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/outputs/bin/START_LST_wall_local_v10')
SOLVER_HASH='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
CASES=json.loads((ROOT/'outputs/cases.json').read_text())
def snapshot():
    data=json.loads((ROOT/'outputs/input_manifest.json').read_text())
    for r in data:assert sha(r['path'])==r['sha256'],r['path']
    return data
assert sha(SOLVER)==SOLVER_HASH
