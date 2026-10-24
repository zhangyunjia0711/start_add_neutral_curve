"""Isolated fixed-beta M3 executor, reusing the frozen v10 run adapter."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
PREVIOUS=ROOT.parent/'neutral_beta_full_20261022'
sys.path.insert(0,str(PREVIOUS/'work'))
import full_driver as driver
driver.ROOT=ROOT
run,table,sha,snapshot=driver.run,driver.table,driver.sha,driver.snapshot
SOLVER=PREVIOUS/'outputs/bin/START_LST_wall_local_v10'
SOLVER_HASH='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
OPTIONS=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8','-st_pc_factor_mat_solver_type','superlu']
assert sha(SOLVER)==SOLVER_HASH
