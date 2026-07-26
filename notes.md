# IKDH testing notes

## Objective

Evaluate the correctness and robustness of the inverse kinematics solver across the robot library (1178 DH configurations under robots/), find failure modes, and decide whether the algebraic core (HuPf) is worth keeping in its current form.

## Architecture recap

The solver has two layers.

The HuPf layer (src/hupf) computes closed form candidate solutions from the Husty Pfurner algebraic method: Study parameters, a hyperplane construction split into a shoulder part (three joints) and a wrist part (three joints), and a resultant elimination producing a degree 56 univariate polynomial whose roots map back to joint angles through the tangent half angle substitution v = tan(theta/2).

The refinement layer (src/ikdh.cpp, function refineIK) is a damped least squares Newton iteration. It takes every HuPf candidate plus a fixed set of Halton low discrepancy seeds, refines each with an analytic geometric Jacobian, and keeps the distinct converged results. This layer has no dependency on HuPf's internals; it only needs a forward kinematics function and a Jacobian, both computed directly from the DH table.

## Testing methodology

Two kinds of tests were used throughout.

Round trip test: for a given robot and a chosen joint configuration q, compute the end effector pose T = FK(q), then call solve(T) and check that q appears among the returned solutions. This checks completeness (did we miss a real branch) rather than just correctness (every returned solution already satisfies FK by construction of refineIK).

Wide test: run the round trip test across most of the robot library at a small number of joint configurations per robot (home position q = 0 always included, plus two or three configurations spread across the joint limits), capped so a full run stays in the minutes range rather than hours. Two things were tracked for every (robot, configuration) pair: whether q = 0 (or the tested q) was recovered, and whether the number of returned solutions was even. For a generic, non singular pose a 6R solver should return solutions in mirrored pairs (elbow up/down, wrist flip, shoulder left/right), so an odd count is a signal worth investigating rather than a hard error.

Crash survival: because a bad candidate in HuPf can in principle abort the process, the wide test isolates each (robot, configuration) call in its own subprocess so one crash does not kill the whole run and is instead recorded as a data point.

## Findings

### Two real crashes, fixed

Two robots (denso_cobotta, newker_new2006b) crashed the process. Root cause was in Polynomial::operator* and operator+= (src/hupf/Polynomial.h): when either operand was the zero polynomial, the result's realCoeffType flag (univariate vs bivariate representation) was not propagated from the actual operands, it was hardcoded. This produced a polynomial internally tagged as bivariate while holding univariate data, which later corrupted degree bookkeeping deep in the resultant pipeline and crashed. Fixed by propagating the flag correctly in both operators. Verified against the full 1178 robot library: zero crashes, and solution counts unchanged versus the pre fix baseline everywhere the comparison was possible.

### HuPf's algebraic core does not run for most robots

The Calculate::InverseKin entry point (src/hupf/Calculate.cpp) requires the hyperplane construction to produce 8 hyperplanes (4 for the shoulder triple, 4 for the wrist triple); if fewer are produced it returns no candidates at all and the refinement layer runs on Halton seeds alone. A full library scan showed this happens for 875 of 1178 robots (74 percent). The trigger is the classic spherical wrist DH pattern (a4 = a5 = 0, alpha4 = plus or minus 90 degrees, alpha5 = minus or plus 90 degrees, d5 = 0), which is extremely common in industrial arms and for which the wrist side hyperplane construction has no implemented case. In other words, for three out of four robots in the library, HuPf contributes zero candidates today; every solution actually returned comes from the Halton seeded Newton refinement.

### Condition 1: non convergence near the wrist singularity

The reported home position gap (q = 0 not recovered) traces to the wrist self motion singularity: when theta5 (or its DH offset) is a multiple of 180 degrees, joints 4 and 6 become coaxial and the inverse kinematics problem has a one dimensional family of solutions rather than isolated points. This is a structural property of the pose, not a bug in either solver layer: an isolated root finder (HuPf) cannot enumerate a positive dimensional solution set, and a local Newton method can fail to converge or converge to an arbitrary point on the manifold depending on the seed. Checked against the DH theta offset directly across the library: 983 of 985 candidate cases matched this explanation.

### Condition 2: odd solution counts

Majority explanation (79 percent of the odd count cases inspected) is legitimate joint limit filtering: the mirror solution that would restore parity exists mathematically but falls outside that joint's configured range, so it is correctly excluded. This was confirmed, not assumed, by checking the excluded mirror solutions directly against the limits.

One specific alternative hypothesis was tested and rejected: using closed intervals [lo, hi] instead of half open ones for the limit check. Measured violations of the excluded mirror solutions were 2.75 to 4.5 degrees outside the limit, far beyond a boundary or floating point epsilon, so widening the interval boundary would not recover them. The residual 21 percent of odd count cases is attributed to root classification fragility rather than a single identified cause, and was not pursued further.

### Jacobian conditioning at the singularity

Measured cond(J) and cond(J^T J) directly at the exact wrist singularity and at a generic pose. Generic pose: cond(J) about 10, cond(J^T J) about 100, unremarkable. At the exact singularity: cond(J) about 3e8, cond(J^T J) about 8.35e16, i.e. at the edge of double precision's dynamic range (about 1e16). Solving the damped step through normal equations (J^T J + lambda I) dq = J^T err squares the condition number and loses essentially all correct digits in the degenerate direction right where the singularity handling matters most.

Fix: solveDampedLeastSquares (src/ikdh.cpp) solves the same damped least squares problem by folding the augmented system [J; sqrt(lambda) I] into an upper triangular form with sequential Givens rotations, then back substitutes. This never forms J^T J, so the condition number is never squared. It is a numerical improvement to the refinement layer's stability, independent of HuPf; it does not by itself resolve condition 1, since the underlying issue there is a positive dimensional solution set, not just poor conditioning of an otherwise well posed problem.

### Halton only versus the current hybrid

To check how much HuPf's candidates are actually worth given the 74 percent coverage gap above, the refinement layer was run in two modes: current hybrid (HuPf candidates plus Halton seeds) and Halton only (HuPf disabled, same Halton seed set). Compared across 3534 (robot, configuration) test cases: 98.4 percent exact match in the returned solution sets. Of the 11 genuine mismatches, 10 came from the subset of robots where HuPf's hyperplane construction does succeed (h.size() == 8), meaning HuPf does contribute real, occasionally unique candidates there, it is not simply redundant with Halton even when it does run. It costs nothing to drop for the 74 percent where it already contributes nothing.

## Current refinement layer details

Jacobian: computed analytically per Newton iteration from the accumulated DH transforms (T_accum[0..6]), not numerically and not from any pre derived symbolic expression. For a revolute joint j, the position Jacobian column is z_j cross (p_e minus p_j) and each rotation column is z_j cross R_e's corresponding column, where z_j and p_j come directly from the accumulated transform up to joint j. This is standard geometric Jacobian construction, robot agnostic, and requires no prior derivation per robot.

Note: each robot's YAML file also stores a jacobian: block, a symbolic expression generated offline by tools/derive_jacobian.py using sympy. Confirmed by grep that nothing under include/ or src/ reads that key; it is disconnected from the runtime solver, which always computes its own Jacobian as described above.

Halton seeding: 64 seeds per solve call, one Halton point per seed index using bases 2, 3, 5, 7, 11, 13 for the six joints, mapped from [0,1] into each joint's actual limit range. A single Newton pass refines every seed, followed by the flip expansion pass described below.

Damping: lambda is a fixed constant, 1e-6, not adaptively updated between iterations. The QR based solve above was specifically what made a fixed, small lambda safe to use even near a singularity; with the old normal equations solve a fixed lambda that small would have been numerically unusable in the degenerate direction.

## Open questions

Whether to invest in deriving the missing spherical wrist hyperplane case for HuPf (would close the 74 percent coverage gap, a real algebraic derivation task) or to formally adopt Halton only and remove HuPf from the runtime path (simpler, and the comparison above shows the cost is small but not zero). Not yet decided.

An experimental #ifdef IKDH_TEST_HALTON_ONLY toggle is currently sitting uncommitted in ikdh.cpp, used to build the Halton only comparison binary. Needs to be either removed or formalized depending on the decision above.

Condition 1 (singularity convergence) has a proposed seed selection approach from the user, not yet implemented or discussed in this document.

## Decision: HuPf removed

The open question above is resolved. src/hupf/ has been deleted entirely, along with the lastPolynomial() API, the algebraic seeding path in Solver::solve, and the IKDH_TEST_HALTON_ONLY toggle (moot now that Halton is the only path). The solver is now exactly the Jacobian plus Halton multi-start Newton refinement described above, with no algebraic core and no third-party code, which also simplified the license to a single, full repository PolyForm Noncommercial License. Verified equivalent before/after solution counts and solve time across the full robot library with tools/run_benchmark.cpp; see the repository history around this change for the actual before/after numbers.

## Testing whether every remaining step and constant is necessary

Once HuPf was gone, each remaining piece of Solver::solve (flip expansion, damping, Newton iteration budget, Halton seed count, duplicate merge threshold) was tested individually by disabling or perturbing it and rerunning the full library benchmark, rather than assumed necessary.

Flip expansion (post-convergence, +180 degrees on one joint of each already found solution, one more Newton pass) turned out to matter far more than expected: disabling it dropped 38.8 percent of robots to fewer solutions, an aggregate 2.37 percent solution count loss, roughly eight times the cost of removing HuPf. Kept as is.

Damping (fixed lambda = 1e-6): setting it to zero actually improved aggregate solution count slightly (plus 0.84 percent, 13 percent faster), but the robots that got worse without it were concentrated in the same spherical wrist collaborative arm families (AUBO, DUCO, Siasun) already known to be singularity prone, exactly what damping exists to protect. Kept as is, since the aggregate number masks that it is doing real, targeted work on the hard cases.

Newton iteration budget (maxIter = 50) and Halton seed count (previously 32) both showed a clean monotonic cost: cutting either lost real solutions with no upside. Duplicate merge threshold (about 1 degree combined) also showed no benefit from loosening.

One dead constant found in passing: refineIK's own default parameter maxIter = 40 is never actually used, every call site passes an explicit value (50, or solveFromSeed's own default of 100). Harmless but misleading if read as the effective budget.

## Mirror seeding attempt (not adopted)

Explored a proposed alternative to Halton: converge once from q = 0, then generate further seeds by reflecting each joint's angle across the plane through that joint's pivot, its axis, and the target end effector origin, instead of relying on random coverage. The plane based reflection formula was verified correct against a literal 3D point reflection (exact match). Tried across the full library as up to 64 combinations of six independent per joint reflections: 93.6 percent of ground truth solutions recovered in aggregate, but the mapping from combination to actual solution was not clean, most combinations converged back to the same solution rather than a genuinely different one. Directly testing reflections between already known real solutions (not from a single q = 0 base) confirmed why: mirroring one joint alone while holding the other five exactly fixed is nearly always undone by Newton, 27 of 28 tested single joint reflections of real solutions fell straight back to themselves. This is the same underlying limitation as the existing plus 180 degree flip expansion, not an improvement on it. Capping the combinations at a fixed 16 (fewest joints changed first) lost a further 1185 solutions across the library compared to the full 64. Not adopted; the existing Halton plus flip expansion approach already covers this ground at lower complexity.

## Explicit singularity seeds and final seed count decision

Separately confirmed that most residual FK to IK round trip failures are not a seed coverage problem: of 3534 random (robot, pose) tests, only 0.6 percent failed to recover the original joint configuration, and of those, about six times as many were genuine near singularity conditioning failures (Jacobian condition number in the hundreds of millions or worse) as were well conditioned points simply missed by chance (confirmed directly: two such well conditioned misses, on AUCTECH AN-220-170-2.7 and Staubli TX90L, converged instantly once seeded within a couple degrees of the true answer, and were found by simply extending the same Halton sequence from 32 to 64 points, both times at index 49).

Tested two ways to close this specific gap: doubling Halton to 64 seeds (plus 0.21 percent solutions library wide, plus 51 percent time), and adding 16 fixed seeds at known spherical wrist singular joint values (joint 2 and 5 at 0 or 180 degrees, joint 4 at 0, 90, -90, or 180, other joints at 0; plus 0.16 percent solutions, plus 25 percent time, a better cost per solution ratio). Combining both gave the best raw coverage (plus 0.25 percent) but at the worst combined cost (plus 76 percent), and still left most of the near singularity failures (15 of 18 in the fixed test set) unresolved, since those are genuine rank deficient self motion singularities that no finite seed set fully closes.

Decision: keep Halton at 64 seeds, drop the singularity seed experiment. Simpler than maintaining a second, structurally different seeding mechanism for a comparable gain.

## Adaptive Levenberg-Marquardt damping (adopted)

The seeding experiments above treat the symptom (a fixed damping value handles conditioning badly near a singularity) by adding more or better placed starting points. A more direct fix is to stop using a single fixed lambda at all. Near a singularity the Jacobian's smallest singular value goes to zero, so a step computed with a fixed lambda is simultaneously too weak to be safe there (large, unstable corrections) and too strong everywhere else (slower convergence than necessary for well conditioned poses). This is not a coordinate artifact like the earlier tan half angle issue in the removed HuPf core; a rank deficient Jacobian is an intrinsic, coordinate independent fact about the manipulator's geometry at that pose, so no reparametrization removes it, only better handling of the step itself can.

Implemented standard adaptive Levenberg-Marquardt in refineIK: at each iteration, try the damped step at the current lambda, evaluate the trial FK error, accept and shrink lambda by 3x if it improved, otherwise grow lambda by 3x and retry from the same point (up to 12 trials before giving up on that seed). Verified two ways against the 64 Halton seed baseline: full library benchmark, total solutions up 2.22 percent (5976.4 to 6109.2) for 18.65 percent more time, 360 robots improved against only 6 trivially regressed; and the FK to IK round trip coverage gap test, missed solutions down from 18 to 3 out of 3534, with the near singularity failures specifically down from 15 to 1. This is roughly ten times the benefit of every seeding based experiment tried this session combined, at a lower time cost than simply doubling Halton seeds, confirming the seeding experiments were compensating for a solver limitation rather than a genuine coverage gap.

Decision: adopted as the permanent refineIK implementation, replacing the fixed lambda. The specific constants (12 trial cap, 3x growth and shrink factor, lambda bounds 1e-12 to 1e12) are conventional Levenberg-Marquardt defaults, not individually ablation tested the way the other constants above were; a further tuning pass on these is a plausible next step but not done here.
