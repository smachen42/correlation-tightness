// Small dense linear algebra for n <= 100.
// All routines use plain 2D arrays of arrays for readability;
// at this size that's faster to write and plenty fast to run.

function zeros(n, m) {
  const A = new Array(n);
  for (let i = 0; i < n; i++) A[i] = new Array(m).fill(0);
  return A;
}

// ── RNG: Mulberry32 (32-bit seeded PRNG) ──────────────────────────────────
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6D2B79F5) | 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Box-Muller standard normal (returns one sample; the second is discarded).
function randn(rng) {
  let u;
  do { u = rng(); } while (u === 0);
  const v = rng();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

// ── Wishart-style random correlation matrix ──────────────────────────────
// Mirrors covkit.random_matrices.wishart_covariance with scale = I, df = n,
// then normalises diag to 1.
function wishartCorr(n, rng) {
  const df = n;
  const S = zeros(n, n);
  for (let k = 0; k < df; k++) {
    const z = new Array(n);
    for (let i = 0; i < n; i++) z[i] = randn(rng);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) S[i][j] += z[i] * z[j];
    }
  }
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      S[i][j] /= df;
      S[j][i] = S[i][j];
    }
  }
  const d = new Array(n);
  for (let i = 0; i < n; i++) d[i] = Math.sqrt(S[i][i]);
  const corr = zeros(n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) corr[i][j] = S[i][j] / (d[i] * d[j]);
    corr[i][i] = 1;
  }
  return corr;
}

// ── Cholesky factorization: C = L L^T for symmetric strictly-PD C ─────────
// Returns the lower-triangular factor L, or null on a nonpositive pivot
// (C singular or indefinite) — so success doubles as the cheapest
// positive-definiteness test. Does not mutate input.
function cholesky(C) {
  const n = C.length;
  const L = zeros(n, n);
  for (let j = 0; j < n; j++) {
    let s = C[j][j];
    for (let p = 0; p < j; p++) s -= L[j][p] * L[j][p];
    if (s <= 0) return null;
    L[j][j] = Math.sqrt(s);
    for (let i = j + 1; i < n; i++) {
      let s2 = C[i][j];
      for (let p = 0; p < j; p++) s2 -= L[i][p] * L[j][p];
      L[i][j] = s2 / L[j][j];
    }
  }
  return L;
}

// Solve L L^T X = B given the Cholesky factor L: forward then back
// substitution, column by column of the n×k RHS B.
function cholSolve(L, B) {
  const n = L.length;
  const k = B[0] ? B[0].length : 0;
  const Y = zeros(n, k);
  for (let col = 0; col < k; col++) {
    for (let i = 0; i < n; i++) {
      let s = B[i][col];
      for (let p = 0; p < i; p++) s -= L[i][p] * Y[p][col];
      Y[i][col] = s / L[i][i];
    }
  }
  const X = zeros(n, k);
  for (let col = 0; col < k; col++) {
    for (let i = n - 1; i >= 0; i--) {
      let s = Y[i][col];
      for (let p = i + 1; p < n; p++) s -= L[p][i] * X[p][col];
      X[i][col] = s / L[i][i];
    }
  }
  return X;
}

// ── Solve C X = B for symmetric C, n×k RHS B ─────────────────────────────
// Tries Cholesky (fast, requires C strictly PD); on failure (singular or
// indefinite C), falls back to the minimum-norm least-squares solution via
// the symmetric eigendecomposition — no special-case detection needed.
function solveSymmetric(C, B) {
  const n = C.length;
  if (n === 0) return [];
  const L = cholesky(C);
  if (!L) return solvePseudoinverse(C, B);
  return cholSolve(L, B);
}

// Minimum-norm least-squares solve via symmetric eigendecomposition:
// C = Q diag(λ) Q^T  →  X = Q diag(λ⁺) Q^T B
// where λ⁺_i = 1/λ_i if |λ_i| > tol, else 0. Handles singular and indefinite
// inputs uniformly; reduces to the exact inverse when C is nonsingular.
function solvePseudoinverse(C, B) {
  const n = C.length;
  const k = B[0] ? B[0].length : 0;
  const { values, vectors } = eigh(C);
  // Standard SVD-style tolerance: n × machine-eps-ish × largest singular value.
  const tol = n * 1e-12 * Math.max(1e-300, ...values.map(Math.abs));

  // Form (diag(λ⁺) Q^T B) in one pass, then multiply by Q.
  const tmp = zeros(n, k);
  for (let i = 0; i < n; i++) {
    const inv = Math.abs(values[i]) > tol ? 1 / values[i] : 0;
    for (let col = 0; col < k; col++) {
      let s = 0;
      for (let r = 0; r < n; r++) s += vectors[i][r] * B[r][col];
      tmp[i][col] = s * inv;
    }
  }
  const X = zeros(n, k);
  for (let r = 0; r < n; r++) {
    for (let col = 0; col < k; col++) {
      let s = 0;
      for (let i = 0; i < n; i++) s += vectors[i][r] * tmp[i][col];
      X[r][col] = s;
    }
  }
  return X;
}

// ── Symmetric eigenvalues: Householder tridiagonalization + implicit QL ──
// Eigenvalues only (no vectors), sorted descending. Does not mutate input.
// O(n³) overall — the Householder reduction is O(4n³/3) and the QL sweeps add
// O(n²) per eigenvalue — versus the O(n⁴) cyclic Jacobi this replaces; ~5–10×
// faster, which matters for the per-selection validity check and the λ_min
// curve at large n. The 2×2 closed form skips the iterative path for the
// commonest small case. `eigh` (values + vectors) still uses Jacobi; it is
// off the hot path. Algorithm: EISPACK tred1/tql1 (Numerical Recipes
// tred2/tqli, vectorless variant), validated against the Jacobi path.
function eigvalsh(M) {
  const n = M.length;
  if (n === 0) return [];
  if (n === 1) return [M[0][0]];
  if (n === 2) {
    const a = M[0][0], b = M[0][1], c = M[1][1];
    const mid = (a + c) / 2;
    const disc = Math.sqrt(Math.max(0, ((a - c) / 2) ** 2 + b * b));
    return [mid + disc, mid - disc];
  }

  const a = M.map(r => r.slice());
  const d = new Array(n).fill(0), e = new Array(n).fill(0);

  // Householder reduction to symmetric tridiagonal form. After this, d holds
  // the diagonal and e[1..n-1] the subdiagonal of a matrix orthogonally
  // similar to M (so with the same spectrum). Vectorless: the transforms are
  // applied but never accumulated into a Q.
  for (let i = n - 1; i > 0; i--) {
    const l = i - 1;
    let h = 0, scale = 0;
    if (l > 0) {
      for (let k = 0; k <= l; k++) scale += Math.abs(a[i][k]);
      if (scale === 0) {
        e[i] = a[i][l];
      } else {
        for (let k = 0; k <= l; k++) { a[i][k] /= scale; h += a[i][k] * a[i][k]; }
        let f = a[i][l];
        let g = f >= 0 ? -Math.sqrt(h) : Math.sqrt(h);
        e[i] = scale * g;
        h -= f * g;
        a[i][l] = f - g;
        f = 0;
        for (let j = 0; j <= l; j++) {
          g = 0;
          for (let k = 0; k <= j; k++) g += a[j][k] * a[i][k];
          for (let k = j + 1; k <= l; k++) g += a[k][j] * a[i][k];
          e[j] = g / h;
          f += e[j] * a[i][j];
        }
        const hh = f / (h + h);
        for (let j = 0; j <= l; j++) {
          f = a[i][j];
          e[j] = g = e[j] - hh * f;
          for (let k = 0; k <= j; k++) a[j][k] -= f * e[k] + g * a[i][k];
        }
      }
    } else {
      e[i] = a[i][l];
    }
  }
  for (let i = 0; i < n; i++) d[i] = a[i][i];

  // Implicit-shift QL on the tridiagonal (d, e). e is shifted so e[i] is the
  // subdiagonal below d[i]; e[n-1] is a sentinel zero.
  for (let i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
  for (let l = 0; l < n; l++) {
    let iter = 0, m;
    do {
      for (m = l; m < n - 1; m++) {
        const dd = Math.abs(d[m]) + Math.abs(d[m + 1]);
        if (Math.abs(e[m]) <= Number.EPSILON * dd) break;
      }
      if (m !== l) {
        if (iter++ === 50) break;  // fail-safe; never observed in testing
        let g = (d[l + 1] - d[l]) / (2 * e[l]);
        let r = Math.hypot(g, 1);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? Math.abs(r) : -Math.abs(r)));
        let s = 1, c = 1, p = 0, i;
        for (i = m - 1; i >= l; i--) {
          let f = s * e[i];
          const b = c * e[i];
          r = Math.hypot(f, g);
          e[i + 1] = r;
          if (r === 0) { d[i + 1] -= p; e[m] = 0; break; }
          s = f / r;
          c = g / r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
        }
        if (r === 0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0;
      }
    } while (m !== l);
  }

  d.sort((x, y) => y - x);
  return d;
}

// ── Symmetric eigendecomposition: values AND vectors ─────────────────────
// Same cyclic Jacobi as eigvalsh, but each Givens rotation also accumulates
// into Q (initialised to identity). On exit, A is diagonal (eigenvalues on
// the diagonal) and the columns of Q are the corresponding eigenvectors.
// Returns {values, vectors}, sorted by value descending; vectors[k] is the
// k-th eigenvector as an array. Eigenvector signs are canonicalised so the
// first component with |x| > 1e-12 is positive (makes shareable URLs stable).
function eigh(M) {
  const n = M.length;
  if (n === 1) return { values: [M[0][0]], vectors: [[1]] };

  const A = new Array(n);
  for (let i = 0; i < n; i++) A[i] = M[i].slice();
  const Q = zeros(n, n);
  for (let i = 0; i < n; i++) Q[i][i] = 1;

  const maxSweeps = 100;
  const tol = 1e-14;

  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) {
      for (let q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
    }
    if (off < tol) break;

    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = A[p][q];
        if (apq === 0) continue;
        const app = A[p][p], aqq = A[q][q];
        let c, s, t;
        const diff = aqq - app;
        if (diff === 0) {
          c = s = Math.SQRT1_2; t = 1;
        } else {
          const tau = diff / (2 * apq);
          t = tau >= 0
            ? 1 / (tau + Math.sqrt(1 + tau * tau))
            : -1 / (-tau + Math.sqrt(1 + tau * tau));
          c = 1 / Math.sqrt(1 + t * t);
          s = t * c;
        }
        A[p][p] = app - t * apq;
        A[q][q] = aqq + t * apq;
        A[p][q] = A[q][p] = 0;
        for (let r = 0; r < n; r++) {
          if (r === p || r === q) continue;
          const arp = A[r][p], arq = A[r][q];
          A[r][p] = A[p][r] = c * arp - s * arq;
          A[r][q] = A[q][r] = s * arp + c * arq;
        }
        // Apply the same rotation to Q's columns p and q.
        for (let r = 0; r < n; r++) {
          const qrp = Q[r][p], qrq = Q[r][q];
          Q[r][p] = c * qrp - s * qrq;
          Q[r][q] = s * qrp + c * qrq;
        }
      }
    }
  }

  // Sort by eigenvalue descending; permute Q's columns to match.
  const idx = new Array(n);
  for (let i = 0; i < n; i++) idx[i] = i;
  idx.sort((a, b) => A[b][b] - A[a][a]);

  const values = new Array(n);
  const vectors = new Array(n);
  for (let k = 0; k < n; k++) {
    values[k] = A[idx[k]][idx[k]];
    vectors[k] = new Array(n);
    for (let r = 0; r < n; r++) vectors[k][r] = Q[r][idx[k]];
  }

  // Sign convention: first nonzero component positive.
  for (let k = 0; k < n; k++) {
    for (let r = 0; r < n; r++) {
      if (Math.abs(vectors[k][r]) > 1e-12) {
        if (vectors[k][r] < 0) {
          for (let s = 0; s < n; s++) vectors[k][s] = -vectors[k][s];
        }
        break;
      }
    }
  }

  return { values, vectors };
}

// ── Feasible range for off-diagonal entry (i,j) given the rest ───────────
// Direct port of covkit.bounds_analyser.get_feasible_range.
// Returns {lo, hi, alpha, beta, gamma}. The Schur-complement scalars
// (alpha, beta, gamma) are exposed so callers can derive the partial
// correlation and the half-width without re-solving.
function feasibleRange(corr, i, j) {
  const n = corr.length;
  const rest = [];
  for (let k = 0; k < n; k++) if (k !== i && k !== j) rest.push(k);
  const m = rest.length;
  if (m === 0) return { lo: -1, hi: 1, alpha: 0, beta: 0, gamma: 0 };

  const a = new Array(m), b = new Array(m);
  const C = zeros(m, m);
  for (let r = 0; r < m; r++) {
    a[r] = corr[i][rest[r]];
    b[r] = corr[j][rest[r]];
    for (let s = 0; s < m; s++) C[r][s] = corr[rest[r]][rest[s]];
  }
  const rhs = zeros(m, 2);
  for (let r = 0; r < m; r++) { rhs[r][0] = a[r]; rhs[r][1] = b[r]; }
  const X = solveSymmetric(C, rhs);

  let alpha = 0, beta = 0, gamma = 0;
  for (let r = 0; r < m; r++) {
    alpha += a[r] * X[r][0];
    beta  += a[r] * X[r][1];
    gamma += b[r] * X[r][1];
  }
  const disc = Math.max((1 - alpha) * (1 - gamma), 0);
  const width = Math.sqrt(disc);
  return { lo: beta - width, hi: beta + width, alpha, beta, gamma };
}

// ── Precision matrix Ω = C⁻¹ (robust to indefinite / singular C) ──────────
// Well-conditioned PD C inverts via Cholesky. Everything else — indefinite
// (dragged past the boundary), exactly singular (the partial-mode slider
// endpoints land on the feasibility boundary; rank-deficient uploads), or
// near-singular PD (where the pair identities would amplify roundoff by
// 1/λ_min) — inverts via the eigendecomposition with a selective ridge:
// eigenvalues with |λ| below the cut are treated as the null channel and
// replaced by a sign-preserving λ ± μ; all others invert exactly. For PSD C
// the μ → 0 limit of the resulting pair quantities equals the per-pair
// pseudoinverse values (range(B) ⊆ range(A) in any PSD Schur partition), so
// the ridge computes the boundary answer in Θ(n³); μ = 1e-12·λ_max is the
// measured optimum of the truncation-vs-roundoff tradeoff (~3e-6 error, two
// decades below the 3-decimal display). Always returns an n×n matrix.
function precisionMatrix(corr) {
  const n = corr.length;
  if (n === 0) return [];
  const L = cholesky(corr);
  if (L) {
    // Cholesky can "succeed" on a numerically singular matrix with a garbage
    // last pivot; the pivot ratio is a cheap conditioning signal. Reroute to
    // the eigendecomposition when it indicates λ_min/λ_max ≲ n·1e-12.
    let pmin = Infinity, pmax = 0;
    for (let j = 0; j < n; j++) {
      const p = L[j][j];
      if (p < pmin) pmin = p;
      if (p > pmax) pmax = p;
    }
    if (pmin * pmin > n * 1e-12 * pmax * pmax) {
      const I = zeros(n, n);
      for (let i = 0; i < n; i++) I[i][i] = 1;
      return cholSolve(L, I);
    }
  }
  const { values, vectors } = eigh(corr);
  const amax = Math.max(1e-300, ...values.map(Math.abs));
  const cut = 1.5e-8 * amax;
  const mu = 1e-12 * amax;
  const inv = values.map(v =>
    Math.abs(v) >= cut ? 1 / v : 1 / (v >= 0 ? v + mu : v - mu));
  const Om = zeros(n, n);
  for (let a = 0; a < n; a++) {
    for (let b = a; b < n; b++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += vectors[k][a] * vectors[k][b] * inv[k];
      Om[a][b] = Om[b][a] = s;
    }
  }
  return Om;
}

// ── Feasible ranges for all pairs from a precision matrix ─────────────────
// Every pair's Schur scalars come from three entries of Ω = C⁻¹: with
// d = Ω_ii Ω_jj − Ω_ij² (the signed determinant of Ω's 2×2 block),
//   1 − α = Ω_jj / d,   1 − γ = Ω_ii / d,   c = ρ_ij + Ω_ij / d,
// and half-width √(max(Ω_ii Ω_jj, 0)) / |d| — the clamp mirrors the
// disc-clamp in feasibleRange, collapsing the band when it is empty. Derivation
// on the theory page under "The precision matrix". Θ(n²) given Ω. A non-finite
// or zero d (catastrophic cancellation at an exact singularity) falls the
// affected pair back to the per-pair feasibleRange. Returns info with
// info[i][j] = info[j][i] = {lo, hi, alpha, beta, gamma}, null diagonal.
function feasibleRangeFromPrecision(corr, Om) {
  const n = corr.length;
  const info = new Array(n);
  for (let i = 0; i < n; i++) info[i] = new Array(n).fill(null);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const oii = Om[i][i], ojj = Om[j][j], oij = Om[i][j];
      const d = oii * ojj - oij * oij;
      let r = null;
      if (isFinite(d) && d !== 0) {
        const beta = corr[i][j] + oij / d;
        const prod = oii * ojj;
        const width = prod > 0 ? Math.sqrt(prod) / Math.abs(d) : 0;
        r = { lo: beta - width, hi: beta + width,
              alpha: 1 - ojj / d, beta, gamma: 1 - oii / d };
      }
      info[i][j] = info[j][i] = r || feasibleRange(corr, i, j);
    }
  }
  return info;
}

// ── Feasible ranges for all off-diagonal entries at once ─────────────────
// One inverse of the full matrix (Θ(n³)) then every pair in O(1), versus
// Θ(n⁵) for calling feasibleRange on each pair. Used on full-matrix events
// (generate / upload / typed value); drags use the Woodbury path below.
function feasibleRangeAll(corr) {
  return feasibleRangeFromPrecision(corr, precisionMatrix(corr));
}

// ── Rank-2 precision update for a single moved entry (Woodbury) ───────────
// During a drag only entry (i, j) changes: C(δ) = C₀ + δ(eᵢeⱼᵀ + eⱼeᵢᵀ),
// where C₀ is the drag-start matrix with precision Ω₀. The Sherman-Morrison-
// Woodbury identity gives Ω(δ) = C(δ)⁻¹ directly from Ω₀ in Θ(n²), with no
// refactorization:
//   Ω(δ) = Ω₀ − Ω₀U K⁻¹ Uᵀ Ω₀,   U = [eᵢ eⱼ],
//   K = M⁻¹ + Uᵀ Ω₀ U = [[Ω₀_ii, Ω₀_ij + 1/δ], [Ω₀_ij + 1/δ, Ω₀_jj]].
// Ω₀U is just columns i, j of Ω₀. Valid on both sides of the PSD boundary
// (indefinite C(δ) is still invertible), so this also handles red-zone drags
// without an eigendecomposition. Anchoring at Ω₀ each tick (with the running
// total δ) keeps errors from accumulating. The determinant lemma gives
// det C(δ) = −δ²·det C₀·det K, so K is singular exactly at the dragged pair's
// interval endpoints — and the *partial-correlation slider's ±1 ends land
// there exactly*. There C(δ) is singular, Ω(δ) blows up, and the per-pair
// identities lose all precision (other pairs' bounds go haywire); det K loses
// it gradually, so the guard is relative — |det K| small versus the scale of
// its terms — returning null so the caller falls back to the eigendecomposition
// route in feasibleRangeAll, which stays accurate at the boundary. Likewise
// for δ ≈ 0, where the column structure is rank-deficient.
function woodburyPrecision(Om0, i, j, delta) {
  const n = Om0.length;
  if (!(Math.abs(delta) > 1e-12)) return null;  // δ ≈ 0: no-op; let caller reuse Ω₀
  const oii = Om0[i][i], ojj = Om0[j][j], oij = Om0[i][j];
  const kij = oij + 1 / delta;
  const detK = oii * ojj - kij * kij;
  // Relative conditioning guard. Empirically the bounds stay accurate to
  // ~1e-11 while |det K|/scale ≳ 1e-4 and only break when it collapses toward
  // machine zero at the boundary; 1e-9 sits safely in that gap.
  const scale = Math.abs(oii * ojj) + kij * kij;
  if (!isFinite(detK) || Math.abs(detK) < 1e-9 * scale) return null;  // at/near an endpoint
  // K⁻¹ entries.
  const k11 = ojj / detK, k12 = -kij / detK, k22 = oii / detK;
  // Columns i, j of Ω₀.
  const p = new Array(n), q = new Array(n);
  for (let a = 0; a < n; a++) { p[a] = Om0[a][i]; q[a] = Om0[a][j]; }
  const Om = zeros(n, n);
  for (let a = 0; a < n; a++) {
    const pa = p[a], qa = q[a];
    for (let b = a; b < n; b++) {
      const v = Om0[a][b] - (k11 * pa * p[b] + k12 * (pa * q[b] + qa * p[b]) + k22 * qa * q[b]);
      Om[a][b] = Om[b][a] = v;
    }
  }
  return Om;
}

// ── Full spectrum as a function of one off-diagonal entry ─────────────────
// Returns {cs, eigvals} with cs evenly spaced over [-1, 1] and
// eigvals[k] = sorted-descending eigenvalues at c = cs[k].
// Callers wanting just λ_min read eigvals[k][n - 1].
// Default 201 samples → step size 0.01, ~15 ms at n=20.
function eigvalsCurve(corr, i, j, nSamples = 201, lo = -1, hi = 1) {
  const M = corr.map(row => row.slice());
  const cs = new Array(nSamples);
  const eigvals = new Array(nSamples);
  for (let k = 0; k < nSamples; k++) {
    const c = lo + (hi - lo) * k / (nSamples - 1);
    cs[k] = c;
    M[i][j] = M[j][i] = c;
    eigvals[k] = eigvalsh(M);
  }
  return { cs, eigvals };
}
