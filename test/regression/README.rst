How to run regression tests
===========================

To run regression tests with default parameters, simply run:

  cd <ffcdir>/tests/regression/
  python test.py

Look at test.py for more options.


How to update references
========================

To update the references for the FFC regression tests, first run the
regression test (to generate the new references) and then run the
script upload:

  cd <ffcdir>/tests/regression/
  python test.py
  ./scripts/upload

Note: You may be asked for your *Bitbucket* username and password when
uploading the reference data, if use of ssh keys fails.

Note: The upload script will push the new references to the
ffc-reference-data repository. This is harmless even if these
references are not needed later.

Note: The upload script will update the file ffc-regression-data-id
and commit this change to the currently active branch, remember to
include this commit when merging or pushing your changes elsewhere.

Note: You can cherry-pick the commit that updated
ffc-regression-data-id into another branch to use the same set of
references there.

Note: If you ever get merge conflicts in the ffc-regression-data-id,
always pick one version of the file. Most likely you'll need to update
the references again.


How to run regression tests against a different set of regression data
======================================================================

To run regression tests and compare to a different set of regression
data, perhaps to see what has changed in generated code since a
certain version, check out the ffc-regression-data-id file you want
and run tests as usual

  cd <ffcdir>/tests/regression/
  git checkout <ffc-commit-id> ffc-regression-data-id
  python test.py

The test.py script will run scripts/download which will check out the
regression data with the commit id from ffc-regression-data-id in
ffc-regression-data/.
