# Accepted fitted-reference context (experimental)

`Track::iterate` assigns a newly accepted fitted trajectory to every effect through
`updateReference`. The same interface is also used during construction and domain
replacement. A consumer that records the state corresponding to a cached
measurement weight must distinguish these operations.

`AcceptedFitReference<KTRAJ>` identifies the trajectory passed during the accepted
attachment loop. Its `matches(trajectory)` method returns true only for that exact
trajectory while the scope is alive. The prior context is restored on normal or
exceptional exit; nested scopes and different threads keep their own context.

This is an experimental integration hook for an optional straw-hit local-information
diagnostic in Offline. It adds no members or virtual functions to `Hit`, `Effect`,
or `Track`. It does not change filter arithmetic, hit weights, accepted-fit status,
the first-piece prior, or closest-approach policy. It has no numerical effect unless
a consumer explicitly uses the context.

The context identifies an accepted attachment; it does not prove that a particular
hit contributed to that posterior. Consumers must still check the hit's activity,
cached weight identity, weight generation, and supported measurement location.

## Integration and validation status

The companion Offline change uses this hook to retain a copy of a local posterior
at the time of the cached measurement. The Offline option remains disabled by
default. This header is not a generic leave-one-hit-out implementation or a claim
of calibrated uncertainty.

This public-source port is awaiting a build and runtime validation. Source review
alone does not establish binary compatibility, control parity, exception behavior,
or reconstruction performance. A consistent rebuild of consumers is required;
do not mix changed template definitions with precompiled implementations.

This draft is experimental, version-specific, and not an endorsed reconstruction
prescription or tutorial.
