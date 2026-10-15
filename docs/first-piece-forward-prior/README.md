# First-piece forward endpoint prior

**Unofficial personal archive of version-specific development code.** This branch
is unmaintained as a distribution, is not an endorsed Mu2e tutorial or result, and
may have known limitations. This correction is a draft pending native compilation
and testing against its current base.

The branch is based on KinKal `71cc7d928d27b04d3475ac47bea29709481f32f4`
(`v3.8.0`). It contains a focused correction in `Track::iterate`, a deterministic
regression test, and this technical description.

## Defect and correction

The first fitted piece is initialized from the completed backward-filtered state.
That state includes the backward endpoint regularizer. Later material-smoothed
pieces combine both directional states and therefore contain both endpoint
regularizers. A material boundary that changes neither parameters nor covariance
can expose this discrepancy even though the physical state should be unchanged.

The correction saves a value copy of the original forward endpoint information
immediately after `initFitState`. After the backward pass, it adds that saved
information once before `initTraj` constructs the first piece. Both the parameter
estimate and covariance use the completed information. The two filtering passes
and material caches retain their existing processing.

Saving the prior before filtering matters: the completed forward state contains
measurements and transported information that must not be counted twice. Material
caches already contain the forward prior, so the correction does not add it there.

## Regression test and execution status

`Tests/FirstPiecePrior_unit.cc` uses two full-rank parameter constraints with
diagonal covariance. It checks one iteration with and without an active synthetic
material boundary whose parameter change and process noise are both zero, at two
temperatures. Every fitted mean and covariance element is compared with the
independent scalar normal equations, including the first and later pieces.

For example, a unit-variance seed at zero, covariance deweight factor four, and
two unit-variance measurements at one imply a fitted mean of `0.8` and variance
`0.4` at temperature zero. Omitting the forward prior gives `8/9` and `4/9` for
the first piece. These are synthetic regression values, not detector results.

Source review and whitespace checks are complete. Native compilation and test
execution are pending; the preparation environment does not provide ROOT.
In a configured KinKal build, run:

```sh
cmake --build <build-dir> --target Test_FirstPiecePrior
ctest --test-dir <build-dir> -R '^FirstPiecePrior$' --output-on-failure
```

## Scope

This restores consistency with the existing two-endpoint numerical regularization
convention. Endpoint priors are derived from the reference trajectory and are not
independent physical measurements. The correction does not establish a calibrated
prior model or guarantee exact geometric continuity across nonlinear transitions.

The regression uses `KinematicLine`, disables field corrections and end processing,
and checks a single iteration. General nonzero-field domain behavior, other
trajectory types, and full reconstruction performance remain outside this test.
