# First-piece endpoint-prior test results

The correction restores the expected first-piece mean and covariance in all
16 deterministic cases tested. Later-piece results are unchanged. Six existing
default tests and the additional extension recipe pass with and without the fix.
The additional gradient-field recipe fails in both versions and gives mixed
results, so these checks do not establish that the correction improves every
fitting configuration.

These are standalone KinKal software tests using constructed inputs and toy
tracks, not detector reconstruction results. Tests and this report are retained
on a separate reference branch; they are not part of the production-only PR.

## Exactly what was compared

The runs completed on October 22, 2026, using ROOT 6.38.00, GCC 13.3.0,
CMake 3.30.5 and Release/C++20 builds.

| Source | Commit |
| --- | --- |
| Upstream base | `eed6815ee4db8f023e3fd9b5c5afca5484b05ed7` |
| A: unchanged production code plus regression | `26d55b1f1c5484289d7fb98c26bbd0007df086c4` |
| B: identical regression plus endpoint-prior fix | `b7c8901ff55c12dda0da13657e40e78d7ac282fc` |

The only production difference between A and B is the five-line change in
`Fit/Track.hh`: save the original forward endpoint weights before filtering,
then add them once to the completed backward state before constructing the
first piece. The regression source is
[FirstPiecePrior_unit.cc](https://github.com/r-hensley/KinKal/blob/b7c8901ff55c12dda0da13657e40e78d7ac282fc/Tests/FirstPiecePrior_unit.cc).
Adding this report does not change the tested code or represent a new run.

## The calculation with a known answer

### What this test asks

Does the first segment of the fitted track include the same information as
the later segments? We check this using a simple fit with a known answer.

A **piece** is one segment of the trajectory stored by the fitter. Each piece
has six fitted parameter values and a **covariance matrix** describing their
uncertainties. The six diagonal entries are variances (squared uncertainties);
the other entries describe how errors in different parameters vary together.
Here, **means** means the six fitted values, not an average over many tracks.

The fitter starts its forward and backward passes with weak constraints around
the starting track estimate: the **endpoint priors**, or regularizers. Its
existing convention includes both in each final piece. Before the fix, the
first piece includes only one; later pieces already include both. These priors
are not two independent detector measurements.

### A one-parameter example: adding weights, not variances

In the simplest case (uncorrelated uncertainties, temperature 0), the parameter
`d0` starts at 0. Two synthetic measurements each give 1 with variance 1. Each
endpoint prior is centered at 0 with variance 4.

For this Gaussian calculation, **weight = 1 / variance**. Each measurement gets
weight 1; each prior gets only 1/4 because its uncertainty is larger. We add
these weights, not the variances themselves.

1. Add the weights: `1 + 1 + 0.25 + 0.25 = 2.5`.
2. Multiply each value by its weight and add: `1*1 + 1*1 + 0*0.25 + 0*0.25 = 2`.
3. Divide to get the fitted value: `2 / 2.5 = 0.8`.
4. Invert the total weight to get the fitted variance: `1 / 2.5 = 0.4`.

The measurements alone would give value 1 and variance 0.5. The weak priors
pull the answer toward 0 and reduce its variance slightly. Omitting one prior
leaves total weight 2.25, giving value `2 / 2.25 = 8/9` (about 0.889) and
variance `1 / 2.25 = 4/9` (about 0.444). Those are the first-piece values
observed before the fix; afterward they are 0.8 and 0.4 as expected.

### Where the 16 cases and 28 pieces come from

Each case is a constructed fit, not a detector event. We use one iteration of
`KinematicLine`, zero magnetic field, and two measurements that directly
constrain all six parameters. Field correction and end processing are disabled.

We test every combination of:

- **Two uncertainty choices:** diagonal covariance (uncorrelated errors) or
  correlated covariance (errors in different parameters are coupled).
- **Two temperatures, 0 and 1:** fitter settings that change the measurement
  and prior weights, not physical temperatures.
- **Four boundary choices:** none, one active, two active, or one inactive.
  These produce 1, 2, 3 and 1 pieces, respectively. Active boundaries split
  the track; an inactive one does not.

That gives `2 * 2 * 4 = 16` fits, each with one first piece. The four boundary
choices contribute `1 + 2 + 3 + 1 = 7` pieces for each uncertainty/temperature
combination: `4 * 7 = 28` pieces total, comprising 16 first and 12 later pieces.

The artificial boundaries cause no energy loss, scattering or parameter change.
All pieces in a given fit should therefore have the same values and covariance.

### How to read the results

For each piece, the expected six parameter values and 36 covariance entries
are calculated separately from the fitter using Gaussian weighted-fit
equations. An entry passes if
`abs(actual - expected) <= 1e-10 * (1 + abs(expected))`.

| Check | A: without fix | B: with fix |
| --- | ---: | ---: |
| First piece: all six fitted values correct | 0/16 | 16/16 |
| First piece: full covariance correct | 0/16 | 16/16 |
| Later pieces: fitted values and covariance correct | 12/12 | 12/12 |
| All pieces: fitted values and covariance correct | 12/28 | 28/28 |
| Test runner (CTest) exit code | 8: fail | 0: pass |

- **First row:** count a fit only if all six values in its first piece match
  the expected answer. Before the fix, none do; afterward, all 16 do.
- **Second row:** check all 36 uncertainty/correlation entries for those same
  16 first pieces. These are not 16 additional fits.
- **Third row:** check both values and covariance for each of the 12 later
  pieces. All pass before and after the fix.
- **Fourth row:** count all 28 pieces, requiring both checks to pass. Before
  the fix, only the 12 later pieces pass; afterward, all 28 pass.
- **Last row:** CTest's overall result. Code 8 means a test failed, not eight
  cases failed. Failure before the fix is expected: the test detects the bug.

The later pieces also agree directly between A and B at the saved 17-digit
precision: `12 pieces * (6 fitted values + 36 covariance entries) = 504`
unchanged numbers. The largest normalized element error in B is `8.29e-15`,
well below the `1e-10` tolerance.

This establishes the missing-prior correction for the tested calculation.
The artificial boundaries do not test real energy loss or scattering, and a
single iteration does not establish uncertainty calibration after a general
iterative fit.

## Checks exercising more of the fitter

Each fit recipe uses the same seed (`123421`), settings and schedules in A and
B. Counts refer to the 1,000 generated toys after the program's preliminary
example. Existing tolerances and the 250-second execution limit were retained.

| Existing test or recipe | A | B |
| --- | --- | --- |
| Default KinematicLineFit | Pass; 1,000 converged | Pass; 1,000 converged |
| Default CentralHelixFit | Pass; 846 converged | Pass; 847 converged |
| Default LoopHelixFit | Pass; 970 converged | Pass; 970 converged |
| KinematicLine, CentralHelix and LoopHelix BField tests | All three pass | All three pass |
| KinematicLine extension/end processing | Pass; 952 converged, 48 gap-diverged | Pass; 952 converged, 48 gap-diverged |
| CentralHelix gradient-field recipe | Fail; program status -3 | Fail; program status -3 |

"Pass" means the existing executable returned success. Both default
CentralHelix logs still print parameter/covariance mismatch diagnostics, so
successful exit is not a claim that every diagnostic was clean. These were six
selected default tests, not the entire upstream test suite.

The extension recipe completes in 27.949 seconds for A and 29.692 seconds for B.
An earlier test on the older v3.8 source timed out with the correction; that
timeout does not recur on this upstream pin. Upstream lookup and inversion
changes prevent attributing that difference to a single cause. No separate
domain-growth repair is included here.

## What the gradient failure means

The gradient recipe uses `CentralHelixFit --Schedule MCAmbig.txt --Bgrad -0.01
--seed 123421`. It exercises a spatially varying magnetic field. Both versions
finish execution but fail the existing front, middle and back momentum-error
checks, returning program status -3 (shell status 253).

| Status among 1,000 toys | A | B |
| --- | ---: | ---: |
| Converged | 194 | 275 |
| Failed | 22 | 6 |
| Chi-square diverged | 18 | 30 |
| Gap-diverged | 96 | 12 |
| Parameter-diverged | 670 | 677 |

The reported Gaussian widths of reconstructed-minus-true momentum increase:

| Sampling position | A width (MeV/c) | B width (MeV/c) |
| --- | ---: | ---: |
| Front | 6.10984 | 7.03142 |
| Middle | 17.8414 | 19.0527 |
| Back | 18.0498 | 19.4493 |

These are outputs of an already failing diagnostic, not calibrated track
uncertainties. The histograms include fits marked usable by the existing test;
the comparison is not restricted to the same successful tracks in both runs.
The aggregate counts also do not identify which individual tracks were gained
or lost. Thus more converged fits does not establish better momentum accuracy,
and the width changes do not quantify a controlled resolution change on a
matched set of tracks.

The failure existed before this correction, so its existence cannot be blamed
on the new code. However, the correction changes the outcomes, including some
adverse diagnostic values. We have not established that those changes are
harmless. This is the specific unresolved issue behind the earlier phrase
"broader compatibility remains unresolved"; it does not invalidate the
separate calculation with a known answer.

## Repeating the checks

Use separate source/build directories at commits A and B with the toolchain
above. Name each source directory `KinKal`, as required by the header includes.
For each version, with `KinKal` in the current directory:

```bash
cmake -S KinKal -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_DISABLE_FIND_PACKAGE_Git=TRUE
cmake --build build --parallel 2 --target Test_FirstPiecePrior
ctest --test-dir build -V -R '^FirstPiecePrior$'
```

A is intentionally expected to fail; B should pass. The test prints `STATE`,
`CHECK` and `SUMMARY` records for inspecting every piece. For the additional
recipes, build `Test_KinematicLineFit` and `Test_CentralHelixFit`, set
`KINKAL_SOURCE_DIR` to the absolute source directory, and use each version's
own built libraries. Run in separate output directories with the same limit:

```text
KinematicLineFit --constrainpar 5 --Bz 0.0 --Schedule driftfit.txt --extend driftextend.txt --seed 123421
CentralHelixFit --Schedule MCAmbig.txt --Bgrad -0.01 --seed 123421
```

The fork-only preview PR is based on older fork main `71cc7d9`; these results
belong to the explicitly listed A/B commits on upstream `eed6815`, not a new
qualification of that preview commit.
