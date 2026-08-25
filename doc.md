# IKDH

The IKDH solver was optimized into a numerical method based on Newton-Raphson applied to the Jacobian from Halton points. This solver is able to find all solutions for a given pose, subject to good convergence.

## General principle of the solver

The IKDH solver relies on a numerical approach to inverse kinematics, based on cancelling the kinematic error with a Newton-Raphson method applied to the robot's Jacobian. The goal is to find a joint configuration $\mathbf{q} \in \mathbb{R}^6$ such that the pose obtained by forward kinematics matches the desired pose $\mathbf{x}_d$.

**Problem formulation.** The kinematic error is defined as

$$F(\mathbf{q}) = FK(\mathbf{q}) - \mathbf{x}_d \qquad (2.1)$$

where $FK(\mathbf{q})$ is the robot's forward kinematics. Finding an inverse kinematics solution is then equivalent to cancelling this error:

$$F(\mathbf{q}) = \mathbf{0} \qquad (2.2)$$

This can also be seen as minimizing the norm

$$E(\mathbf{q}) = \| F(\mathbf{q}) \| \qquad (2.3)$$

but IKDH uses a root-finding method rather than gradient descent. **Newton-Raphson method.** The method relies on the local linearization of $F$ around $\mathbf{q}_k$:

$$F(\mathbf{q}_k + \Delta\mathbf{q}) \approx F(\mathbf{q}_k) + J(\mathbf{q}_k)\, \Delta\mathbf{q} \qquad (2.4)$$

where $J(\mathbf{q}) = \partial F / \partial \mathbf{q}$ is the robot's Jacobian. Imposing $F(\mathbf{q}_{k+1}) = \mathbf{0}$ gives the Newton increment:

$$\Delta\mathbf{q}_k = -J(\mathbf{q}_k)^{\dagger} F(\mathbf{q}_k) \qquad (2.5)$$

and the iteration:

$$\mathbf{q}_{k+1} = \mathbf{q}_k + \Delta\mathbf{q}_k \qquad (2.6)$$

The pseudo-inverse $J^{\dagger}$ handles non-square cases and situations close to singularities. IKDH uses a regularized version to improve numerical stability.



**Multi-start strategy by Halton sequence.** Since the Newton method is local, it only converges if the initial point $\mathbf{q}_0$ belongs to a solution's basin of attraction. To explore the full set of solution branches, IKDH generates a set of initial points $\mathbf{q}_0^{(i)}$ using a quasi-random Halton sequence, which ensures uniform coverage of the joint space:

$$\mathbf{q}_0^{(i)} \in H_N \subset [\mathbf{q}_{\min}, \mathbf{q}_{\max}]^n \qquad (2.7)$$

For each initial point, the Newton method is applied:

$$\mathbf{q}^{(i)} = \text{Newton-Raphson}\left(\mathbf{q}_0^{(i)}\right) \qquad (2.8)$$

The converged solutions are then filtered, normalized, and grouped by clustering in joint space to identify the different inverse kinematics branches.

## Limits and difficulties encountered

The solver is able to find the solutions for a given pose, but it happens that with 64 Halton points at the start of the optimization, some solutions are not found. To find the missing solution(s), the symmetry of kinematic chain can be used.
