===========================
Changes in the next release
===========================


Summary of changes
==================

.. note:: Developers should use this page to track and list changes
          during development. At the time of release, this page should
          be published (and renamed) to list the most important
          changes in the new release.

- Generalize ufc interface to non-affine parameterized coordinates
- Add ``ufc::coordinate_mapping`` class
- Make ufc interface depend on C++11 features requiring gcc version >= 4.8
- Change the mapping ``pullback as metric`` to ``double covariant piola``
- Include comment with effective representation and integral metadata
  to generated ``tabulate_tensor`` code


Detailed changes
================

.. note:: At the time of release, make a verbatim copy of the
          ChangeLog here (and remove this note).
