# Release 20261019 — Changelog and release notes

Summary
- Development release updating build/CI, adding PARI/FLINT/NTL conversion
	helpers, and migrating internal big-integer usage to NTL::ZZ and FLINT::fmpz_t. Includes many refactors, new tests/programs, and small bugfixes.

Highlights
- New conversion layer: eclib/convert.h + libsrc/convert.cc — conversions between PARI (GEN/INT), FLINT (fmpz_t) and NTL::ZZ.
- Bigint migration: wide replacement of the old bigint typedef with NTL::ZZ across public headers and sources.
- Build & CI: require C++17, updated configure.ac, libtool version bump (LT_CURRENT=16, LT_AGE=0), improved GitHub Actions / conda workflows and macOS fixes; FLINT detection/handling reworked.
- Tests & tools: added tests (e.g., tests/flint_test, tests/tconvert, tests/newhecke) and new programs (progs/makensp, progs/listnsp, progs/newspaces).
- Code modernisation: includes cleanup of includes, use of modern C++ idioms, and many small robustness/bug fixes.
- New classes Newspace and Newform for newforms of dimension >=1 (still only trivial character)

Breaking changes
- Public API change: many headers now use NTL::ZZ — downstream code using the old bigint typedef must be updated or use the new conversion helpers.
- Build requirement: a C++17-capable compiler is now required.
- Possible SONAME change: libtool version bump may change the shared library version/SONAME.

Upgrade & packaging notes
- Ensure NTL, FLINT and PARI are available (or pass --with-flint/--with-ntl to configure).
- Update packager/build scripts for C++17 and the new libtool version.
- Verify FLINT detection/pathing if you relied on system FLINT versions.

Testing & release checklist
- Run the full test suite on supported platforms (Linux/macOS).
- Validate CI workflows (especially FLINT build and conda/macOS steps).
- Verify downstream projects compile against updated headers (NTL::ZZ).
- Check the produced library SONAME/version for packaging impact.

One-line summary:

- Release 20261019: modernized the build and CI, added	PARI/FLINT/NTL conversion support, expanded tests and utility programs, extended modular	symbol code to compute newforms of arbitrary dimension, and included a broad set of arithmetic and robustness updates with	important C++17 and API compatibility changes.
