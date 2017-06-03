How to run regression tests
===========================

To run regression tests with default parameters, simply run::

  cd <ffcdir>/tests/regression/
  python test.py

Look at ``test.py`` for more options.


How to update references
========================

To update the references for the FFC regression tests, first commit
your changes, then run the regression test (to generate the new
references) and finally run the script upload::

  <commit your changes>
  cd <ffcdir>/tests/regression/
  python test.py [--use-tsfc]
  ./scripts/upload

Note: Contributors are encouraged to install also TSFC stack and update
references including ``tsfc`` representation using ``--use-tsfc`` flag.
For the installation instructions see ``doc/sphinx/source/installation.rst``.
Note that ``tsfc`` regression test are run in separate plans on the Bamboo
CI system.

Note: You may be asked for your *Bitbucket* username and password when
uploading the reference data, if use of ssh keys fails.

Note: The upload script will push the new references to the
``ffc-reference-data`` repository. This is harmless even if these
references are not needed later.

Note: The upload script will update the file ``ffc-regression-data-id``
and commit this change to the currently active branch, remember to
include this commit when merging or pushing your changes elsewhere.

Note: You can cherry-pick the commit that updated
``ffc-regression-data-id`` into another branch to use the same set of
references there.

Note: If you ever get merge conflicts in the ``ffc-regression-data-id``,
always pick one version of the file. Most likely you'll need to update
the references again.


How to run regression tests against a different set of regression data
======================================================================

To run regression tests and compare to a different set of regression
data, perhaps to see what has changed in generated code since a
certain version, check out the ``ffc-regression-data-id`` file you want
and run tests as usual::

  cd <ffcdir>/tests/regression/
  git checkout <ffc-commit-id> ffc-regression-data-id
  python test.py

The test.py script will run scripts/download which will check out the
regression data with the commit id from ``ffc-regression-data-id`` in
``ffc-regression-data/``.


How to inspect diff in output from executed generated code
==========================================================

Say you have differences in the output of PoissonDG, you can diff the
``.json`` files (NB! within some tolerance on the floating point
accuracy!) like this::

  python recdiff.py ffc-reference-data/r_uflacs/PoissonDG.json output/r_uflacs/PoissonDG.json
  python recdiff.py output/r_uflacs/PoissonDG.json output/r_tensor/PoissonDG.json

Pick any combination of ``ffc-reference-data | output`` and ``r_foo |
r_bar`` you want to compare.


How to manually delete old reference data
=========================================

If you update the tests such that some data should no longer be kept
in the data repository, this approach allows deleting reference data::

  cd ffc-regression-data
  git checkout master
  git pull
  git rm -rf <old-directory-or-files>
  git commit -a -m"Manually updated data because ..."
  git rev-parse HEAD > ../ffc-regression-data-id
  cd ..
  git commit ffc-regression-data-id -m"Manually updated reference id."

This is not automated because it happens rarely.  Probably a good idea
to coordinate with other devs so they don't introduce the deleted
files with another branch.
