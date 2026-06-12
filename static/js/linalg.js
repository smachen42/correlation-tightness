// Small dense linear algebra for n <= 20.
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

// ── Solve C X = B for symmetric C, n×k RHS B ─────────────────────────────
// Tries Cholesky (fast, requires C strictly PD); on failure (singular or
// indefinite C), falls back to the minimum-norm least-squares solution via
// the symmetric eigendecomposition — no special-case detection needed.
function solveSymmetric(C, B) {
  const n = C.length;
  if (n === 0) return [];
  const k = B[0] ? B[0].length : 0;

  const L = zeros(n, n);
  let chol_ok = true;
  for (let j = 0; j < n && chol_ok; j++) {
    let s = C[j][j];
    for (let p = 0; p < j; p++) s -= L[j][p] * L[j][p];
    if (s <= 0) { chol_ok = false; break; }
    L[j][j] = Math.sqrt(s);
    for (let i = j + 1; i < n; i++) {
      let s2 = C[i][j];
      for (let p = 0; p < j; p++) s2 -= L[i][p] * L[j][p];
      L[i][j] = s2 / L[j][j];
    }
  }
  if (!chol_ok) return solvePseudoinverse(C, B);

  // Forward solve L Y = B
  const Y = zeros(n, k);
  for (let col = 0; col < k; col++) {
    for (let i = 0; i < n; i++) {
      let s = B[i][col];
      for (let p = 0; p < i; p++) s -= L[i][p] * Y[p][col];
      Y[i][col] = s / L[i][i];
    }
  }
  // Back solve L^T X = Y
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

// ── Symmetric eigendecomposition: cyclic Jacobi rotations ────────────────
// Returns eigenvalues sorted descending. Does not mutate input.
// O(n^4) but for n <= 20 that's a few thousand ops, well under a millisecond.
function eigvalsh(M) {
  const n = M.length;
  if (n === 1) return [M[0][0]];

  const A = new Array(n);
  for (let i = 0; i < n; i++) A[i] = M[i].slice();

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
          // Numerically stable: pick t with smaller magnitude.
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
      }
    }
  }

  const evals = new Array(n);
  for (let i = 0; i < n; i++) evals[i] = A[i][i];
  evals.sort((a, b) => b - a);
  return evals;
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
