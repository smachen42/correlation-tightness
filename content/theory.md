+++
title = "Theory"
template = "theory.html"
+++

An $n \times n$ correlation matrix $C$ is symmetric, has unit diagonal, and is positive semidefinite (written $C \succeq 0$), which means $x^\top C x \geq 0$ for every $x \in \mathbb{R}^n$. Practically, positive semidefiniteness ensures there are no contradictions across the variables — that some set of random variables could actually exhibit these correlations.

We can learn a lot about the assumptions embedded in a correlation matrix by fixing all entries but one (the correlation between two of the variables) and asking how far that entry can move while $C$ stays a valid correlation matrix. When the rest already explain the two variables well, their correlation is pinned into a narrow range and the entry is *tight*; when the rest barely explain them, the range is wide. This range has a closed form, and is connected to *partial correlation*.

### Closed-form bounds

Fix a pair of indices $i \neq j$, and let $C(\rho)$ denote the matrix obtained from $C$ by setting $C_{ij} = C_{ji} = \rho$. The *feasible range* of entry $(i, j)$ is
$$
\mathcal{F}_{ij}(C) = \lbrace \rho \in \mathbb{R} \mid C(\rho) \succeq 0 \rbrace.
$$
The feasible range is an interval, so the problem reduces to finding its endpoints.

<details>
<summary>Proof</summary>

The condition $C(\rho) \succeq 0$ is equivalent to $\lambda_{\min}(C(\rho)) \geq 0$, where $\lambda_{\min}(M)$ is the smallest eigenvalue of a symmetric matrix $M$. As a pointwise minimum of affine functions of $\rho$, $\lambda_{\min}(C(\rho))$ is concave, so the feasible range is an interval.

</details>

Let $\mathcal{R} = \lbrace 1, \ldots, n \rbrace \setminus \lbrace i, j \rbrace$, let $A = C_{\mathcal{R}, \mathcal{R}}$ be the $(n-2) \times (n-2)$ principal submatrix indexed by $\mathcal{R}$, and let $r_i = C_{\mathcal{R}, i}$ and $r_j = C_{\mathcal{R}, j}$ be the columns of $C$ at indices $i$ and $j$ respectively restricted to $\mathcal{R}$. Assume $A \succ 0$. Then the feasible range is
$$
\mathcal{F}_{ij}(C) = \left[ c - \sqrt{(1 - \alpha)(1 - \gamma)}, c + \sqrt{(1 - \alpha)(1 - \gamma)} \right],
$$
where
$$
\alpha = r_i^\top A^{-1} r_i, \qquad \gamma = r_j^\top A^{-1} r_j, \qquad c = r_i^\top A^{-1} r_j.
$$

<details>
<summary>Proof</summary>

Without loss of generality, assume $i$ and $j$ are the last two indices, so that
$$
C(\rho) = \begin{bmatrix} A & r_i & r_j \\\\ r_i^\top & 1 & \rho \\\\ r_j^\top & \rho & 1 \end{bmatrix}.
$$
Since $A \succ 0$, $C(\rho) \succeq 0$ if and only if the Schur complement of $A$,
$$
S(\rho) = \begin{bmatrix} 1 & \rho \\\\ \rho & 1 \end{bmatrix} - \begin{bmatrix} r_i^\top \\\\ r_j^\top \end{bmatrix} A^{-1} \begin{bmatrix} r_i & r_j \end{bmatrix} = \begin{bmatrix} 1 - \alpha & \rho - c \\\\ \rho - c & 1 - \gamma \end{bmatrix},
$$
is positive semidefinite.
A symmetric $2 \times 2$ matrix is positive semidefinite if and only if its diagonal entries and its determinant are nonnegative. The diagonal conditions $1 - \alpha, 1 - \gamma \geq 0$ already hold, since $C_{\mathcal{R} \cup \lbrace i \rbrace, \mathcal{R} \cup \lbrace i \rbrace}$ and $C_{\mathcal{R} \cup \lbrace j \rbrace, \mathcal{R} \cup \lbrace j \rbrace}$ are positive semidefinite principal submatrices. The remaining determinant condition,
$$
(1 - \alpha)(1 - \gamma) \geq (\rho - c)^2,
$$
gives the stated interval.

</details>

This interval is exactly the blue region on the main range bar shown for the selected entry in the [explorer](/explore/): its center sits at $c$ and its half-width is $\sqrt{(1 - \alpha)(1 - \gamma)}$. The constants $\alpha$, $\gamma$, and $c$ might at first look like opaque quadratic forms, but each has a clean statistical reading that we explore below.

### Connection to partial correlation

To interpret these bounds, treat $C$ as the covariance of standardized random variables $X_1, \ldots, X_n$. The partial correlation between $X_i$ and $X_j$ given $\lbrace X_k \rbrace_{k \in \mathcal{R}}$ (written $\rho_{ij \mid \mathcal{R}}$) is, by definition, the correlation of the residuals from regressing each on $\lbrace X_k \rbrace_{k \in \mathcal{R}}$. This evaluates to
$$
\rho_{ij \mid \mathcal{R}} = \frac{\rho_{ij} - c}{\sqrt{(1 - \alpha)(1 - \gamma)}}.
$$

<details>
<summary>Proof</summary>

The orthogonal projection of $X_i$ onto $\mathrm{span}(X_{\mathcal{R}})$ is $\widehat{X_i} = (A^{-1} r_i)^\top X_{\mathcal{R}}$, with residual $U_i = X_i - \widehat{X_i}$; define $\widehat{X_j}, U_j$ analogously. By definition, $\rho_{ij \mid \mathcal{R}} = \mathrm{Corr}(U_i, U_j)$, which we now compute.

Since $\mathrm{Cov}(X_{\mathcal{R}}) = A$ and $\mathrm{Cov}(X_{\mathcal{R}}, X_j) = r_j$, the linear combination $\widehat{X_i}$ satisfies
$$
\begin{aligned}
\mathrm{Var}(\widehat{X_i}) &= (A^{-1} r_i)^\top A (A^{-1} r_i) = r_i^\top A^{-1} r_i = \alpha, \\\\
\mathrm{Cov}(\widehat{X_i}, X_j) &= (A^{-1} r_i)^\top r_j = r_i^\top A^{-1} r_j = c.
\end{aligned}
$$

By construction $U_i$ is orthogonal to $\mathrm{span}(X_{\mathcal{R}})$, which contains both $\widehat{X_i}$ and $\widehat{X_j}$. Hence the split $X_i = \widehat{X_i} + U_i$ is orthogonal, giving
$$
\mathrm{Var}(U_i) = \mathrm{Var}(X_i) - \mathrm{Var}(\widehat{X_i}) = 1 - \alpha,
$$
and similarly $\mathrm{Var}(U_j) = 1 - \gamma$. Likewise $\mathrm{Cov}(U_i, \widehat{X_j}) = 0$, so
$$
\mathrm{Cov}(U_i, U_j) = \mathrm{Cov}(U_i, X_j) = \mathrm{Cov}(X_i, X_j) - \mathrm{Cov}(\widehat{X_i}, X_j) = \rho_{ij} - c.
$$
Dividing by the residual standard deviations gives the result.

</details>

Comparing with the feasible range, $\rho_{ij}$ is feasible exactly when this partial correlation lies in $[-1, 1]$. The center $c$ is where it vanishes; the endpoints are where it reaches $\pm 1$.

We can also directly interpret the constants $\alpha$, $\gamma$, and $c$. The scalars $\alpha$ and $\gamma$ are the coefficients of determination ($R^2$) from regressing $X_i$ and $X_j$ on the rest: the fraction of each variable's variance the rest explain. The center $c$ can be interpreted as the correlation implied solely by the rest: splitting each variable into its prediction and residual, $X_i = \widehat{X_i} + U_i$, the correlation decomposes orthogonally,
$$
\rho_{ij} = \mathrm{Cov}(\widehat{X_i}, \widehat{X_j}) + \mathrm{Cov}(U_i, U_j) = c + (\rho_{ij} - c),
$$
into this implied part $c$ plus the residual covariance — so $c$ is exactly the value $\rho_{ij}$ would take if the residuals were uncorrelated.

<details>
<summary>Proof</summary>

$R^2 = \mathrm{Var}(\widehat{X_i})/\mathrm{Var}(X_i) = \alpha$, since $\mathrm{Var}(X_i) = 1$ and $\mathrm{Var}(\widehat{X_i}) = \alpha$; likewise $\gamma$ for $X_j$. The covariances in the decomposition were computed above: $\mathrm{Cov}(\widehat{X_i}, \widehat{X_j}) = r_i^\top A^{-1} r_j = c$ and $\mathrm{Cov}(U_i, U_j) = \rho_{ij} - c$.

</details>

When $\alpha$ or $\gamma$ approaches $1$ ($X_i$ or $X_j$ is nearly a linear combination of the rest), the residual variance ($1 - \alpha$ or $1 - \gamma$) collapses to zero, and the entry $\rho_{ij}$ is pinned to its implied value $c$. This is the meaning of *tightness*.
