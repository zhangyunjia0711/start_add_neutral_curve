#!MC 1410
$!READDATASET '"../../outputs/variable_beta_full_20261024/tecplot/neutral_variable_beta_all_candidates.plt"'
  READDATAOPTION = NEW
  RESETSTYLE = YES
$!WRITEDATASET "neutral_variable_beta_all_candidates_roundtrip.dat"
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 17
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
$!READDATASET '"../../outputs/variable_beta_full_20261024/tecplot/neutral_variable_beta_beta.plt"'
  READDATAOPTION = NEW
  RESETSTYLE = YES
$!WRITEDATASET "neutral_variable_beta_beta_roundtrip.dat"
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 17
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
$!READDATASET '"../../outputs/variable_beta_full_20261024/tecplot/neutral_variable_beta_nose.plt"'
  READDATAOPTION = NEW
  RESETSTYLE = YES
$!WRITEDATASET "neutral_variable_beta_nose_roundtrip.dat"
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 17
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
$!READDATASET '"../../outputs/variable_beta_full_20261024/tecplot/neutral_variable_beta_stage.plt"'
  READDATAOPTION = NEW
  RESETSTYLE = YES
$!WRITEDATASET "neutral_variable_beta_stage_roundtrip.dat"
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 17
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
$!READDATASET '"../../outputs/variable_beta_full_20261024/tecplot/neutral_variable_beta_tail.plt"'
  READDATAOPTION = NEW
  RESETSTYLE = YES
$!WRITEDATASET "neutral_variable_beta_tail_roundtrip.dat"
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 17
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
$!QUIT
