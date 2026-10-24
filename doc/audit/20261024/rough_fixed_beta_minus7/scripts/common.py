"""Isolated rough-case run adapter with immutable source guards."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
PREVIOUS=ROOT/'reference'
import run_adapter as driver
driver.ROOT=ROOT;driver.CASE=ROOT
run,table,sha=driver.run,driver.table,driver.sha
SOLVER=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/outputs/bin/START_LST_wall_local_v10')
SOLVER_HASH='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
OPTIONS=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8','-st_pc_factor_mat_solver_type','superlu']
assert sha(SOLVER)==SOLVER_HASH
def snapshot():
    result=json.loads((ROOT/'outputs/input_manifest.json').read_text())
    for e in result:assert sha(e['path'])==e['sha256'],e['path']
    assert sha(SOLVER)==SOLVER_HASH
    return result
