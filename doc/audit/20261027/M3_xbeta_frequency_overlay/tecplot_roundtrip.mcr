#!MC 1410
$!INCLUDEMACRO "M3_xbeta_five_frequencies_same_branch.mcr"
$!WRITEDATASET "M3_xbeta_five_frequencies_same_branch_roundtrip.dat"
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 17
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
$!INCLUDEMACRO "bluntplate_M3_xbeta_five_frequencies.mcr"
$!WRITEDATASET "bluntplate_M3_xbeta_five_frequencies_roundtrip.dat"
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 17
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
$!INCLUDEMACRO "rough_M3_xbeta_five_frequencies.mcr"
$!WRITEDATASET "rough_M3_xbeta_five_frequencies_roundtrip.dat"
  BINARY = NO
  USEPOINTFORMAT = YES
  PRECISION = 17
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
$!QUIT
