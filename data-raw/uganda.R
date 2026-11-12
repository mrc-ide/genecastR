## The shipped example panel.
##
## inst/extdata/uganda_c469y.csv is the C469Y rows of the Ugandan surveillance
## extract, four columns, no filtering beyond the locus.  Rows with n == 0 are
## KEPT deliberately: the package's own drop rule is then exercised by its
## example data rather than only by a synthetic test.
##
## Regenerate with (from the repository root, outside R):
##
##   python -c "
##   import pandas as pd
##   d = pd.read_csv('<extract>.txt', sep=' ', quotechar='\"')
##   d = d[d.Locus == 'C469Y'][['District','year','x','n']]
##   d.columns = ['district','year','x','n']
##   d.sort_values(['district','year']).to_csv(
##       'inst/extdata/uganda_c469y.csv', index=False)"
##
## Expected after GCData's drop rule: 16 demes, 96 surveys, 2016-2022,
## t0 = 2015, T = 112 at M = 16.  tests/testthat/test-GCData.R asserts this.
