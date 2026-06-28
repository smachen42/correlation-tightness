// Arc · step 3 — three free entries → a 3-D body. Same matrix-with-per-cell-bars as
// the other arc pages (each cell shows its 1-D feasible range); click THREE cells and
// their joint feasible region is an orbitable Three.js mesh, built tomographically
// from the 1-free engine. Edit the current point by dragging in the 2-D slice (two
// entries) and the depth slider (the third). The body is invariant under its own three
// values, so it only re-meshes when the REST of the matrix changes; editing the three
// free entries just moves the marker (and updates every cell's bar). n=3 → the elliptope.
(function () {
  'use strict';
  const N = 6, INDEF = -1e-10, NSAMP = 91, NZ = 61;   // fixed 6×6: purely the 3-free geometry
  // Modifier-click a cell to type its value; a modifier keeps plain clicks instant for axis-picking.
  // Lead with ⌘ on Mac (Ctrl-click is the system right-click there) and Ctrl elsewhere; Alt also works.
  const IS_MAC = /Mac|iPhone|iPad|iPod/.test((navigator.platform || '') + ' ' + (navigator.userAgent || ''));
  const EDIT_KEY_HINT = (IS_MAC ? '⌘' : 'Ctrl') + '-click to type a value';
  const RED = 'rgba(180,40,40,0.55)', BAR_COLOR = 'rgba(50,110,200,0.75)', BAR_H = 3;
  const FILL = 'rgba(28,79,150,0.16)', STROKE = 'rgba(28,79,150,0.7)', DARK = '#2d2d2d', BAD = '#b02828', SLICE_LINE = 0x53699a;
  const PW = 300, PH = 300, PPAD = 30, piW = PW - 2 * PPAD, piH = PH - 2 * PPAD;
  const bX = v => PPAD + (v + 1) / 2 * piW;
  const bY = v => PPAD + (1 - (v + 1) / 2) * piH;
  const $ = id => document.getElementById(id);
  const clamp = (lo, hi, x) => Math.max(lo, Math.min(hi, x));

  let corr = null, info = null;
  let cellRefs = null, barGeom = null, barIndef = false, barRO = null;
  // sel is a fixed 3-slot array [x, y, depth]; a slot is an entry or null. Removing an
  // entry nulls ITS slot (the others keep their axis), so deselecting y leaves x and z
  // where they are and frees the y slot for the next pick. selSeq tracks fill order so a
  // 4th pick (all slots full) evicts the oldest in place — again without reshuffling axes.
  let sel = [null, null, null], selSeq = [0, 0, 0], selTick = 0, body = null, drag = null;
  const ready = () => sel.every(Boolean);
  const fillSlot = (slot, e) => { sel[slot] = e; selSeq[slot] = ++selTick; };
  const setSel = es => { sel = [null, null, null]; selSeq = [0, 0, 0]; selTick = 0; es.forEach((e, k) => fillSlot(k, e)); };
  const entry = (i, j) => i > j ? { i, j } : { i: j, j: i };
  const sameE = (a, b) => a && b && a.i === b.i && a.j === b.j;
  const val = e => corr[e.i][e.j];
  const key = e => `${e.i},${e.j}`;
  const lam = () => eigvalsh(corr)[corr.length - 1];

  // ── Region / body math ───────────────────────────────────────────────────
  function sliceInterval(M, Bi, Bj) {
    const n = M.length, rest = [];
    for (let k = 0; k < n; k++) if (k !== Bi && k !== Bj) rest.push(k);
    const m = rest.length; if (m === 0) return { lo: -1, hi: 1 };
    const A = zeros(m, m), rhs = zeros(m, 2);
    for (let r = 0; r < m; r++) { rhs[r][0] = M[Bi][rest[r]]; rhs[r][1] = M[Bj][rest[r]]; for (let s = 0; s < m; s++) A[r][s] = M[rest[r]][rest[s]]; }
    const L = cholesky(A); if (!L) return null;
    const X = cholSolve(L, rhs); let al = 0, be = 0, ga = 0;
    for (let r = 0; r < m; r++) { al += rhs[r][0] * X[r][0]; be += rhs[r][0] * X[r][1]; ga += rhs[r][1] * X[r][1]; }
    if (al > 1 + 1e-12 || ga > 1 + 1e-12) return null;
    const w = Math.sqrt(Math.max(0, (1 - al) * (1 - ga))); const lo = be - w, hi = be + w;
    if (hi < lo) return null; return { lo: Math.max(-1, lo), hi: Math.min(1, hi) };
  }
  function regionOn(M, A, B, refine) {
    const W = M.map(r => r.slice());
    const evalAt = x => { W[A.i][A.j] = W[A.j][A.i] = x; return sliceInterval(W, B.i, B.j); };
    const xs = [];
    for (let s = 0; s < NSAMP; s++) xs.push(-1 + 2 * s / (NSAMP - 1));
    xs.push(M[A.i][A.j]);
    if (refine) {
      const valid = xs.filter(x => evalAt(x)).sort((a, b) => a - b);
      if (valid.length >= 2) {
        const xIn = valid[valid.length >> 1];
        for (const dir of [-1, 1]) {
          if (evalAt(dir)) continue;
          let a = xIn, b = dir;
          for (let it = 0; it < 30; it++) { const mid = (a + b) / 2; if (evalAt(mid)) a = mid; else b = mid; }
          const xE = dir < 0 ? valid[0] : valid[valid.length - 1];
          for (let t = 1; t <= 4; t++) xs.push(xE + (a - xE) * t / 4);
        }
      }
    }
    const order = [...new Set(xs)].sort((a, b) => a - b);
    const top = [], bot = [];
    for (const x of order) { const iv = evalAt(x); if (iv) { top.push({ x, y: iv.hi }); bot.push({ x, y: iv.lo }); } }
    return { top, bot };
  }
  function computeBody(o0, o1, o2) {
    const base = corr.map(r => r.slice());
    const sliceAt = z => { base[o2.i][o2.j] = base[o2.j][o2.i] = z; return regionOn(base, o0, o1, false); };
    const grid = [];
    for (let s = 0; s < NZ; s++) { const z = -1 + 2 * s / (NZ - 1); const reg = sliceAt(z); if (reg.top.length) grid.push({ z, reg }); }
    if (!grid.length) return { slices: [], zlo: 0, zhi: 0 };
    const ne = z => sliceAt(z).top.length > 0, step = 2 / (NZ - 1);
    const refineEnd = (zIn, dir) => {
      const edge = dir < 0 ? -1 : 1; if (ne(edge)) return edge;
      let a = zIn, b = clamp(-1, 1, zIn + dir * step * 1.5); if (ne(b)) return zIn;
      for (let it = 0; it < 22; it++) { const mid = (a + b) / 2; if (ne(mid)) a = mid; else b = mid; }
      return a;
    };
    const zlo = refineEnd(grid[0].z, -1), zhi = refineEnd(grid[grid.length - 1].z, 1);
    const slices = [];
    if (zlo < grid[0].z - 1e-6) slices.push({ z: zlo, ...sliceAt(zlo) });
    for (const g of grid) slices.push({ z: g.z, top: g.reg.top, bot: g.reg.bot });
    if (zhi > grid[grid.length - 1].z + 1e-6) slices.push({ z: zhi, ...sliceAt(zhi) });
    return { slices, zlo, zhi };
  }

  // Raw correlations only: display maps are the identity. (A partial-correlation
  // toggle once reparametrized the axes; it was removed — this explorer shows the
  // correlations directly. `toDisp`/`toRaw` are kept as the seams the plotting/drag
  // code maps through, now no-ops.)
  const toDisp = (e, raw) => raw;
  const toRaw = (e, d) => d;
  const labelFor = e => `&rho;(${e.i + 1},${e.j + 1})`;

  // ── Matrix + per-cell bars (ported from the original explorer) ───────────
  function buildTable() {
    const n = corr.length; let html = '';
    for (let i = 0; i < n; i++) { html += '<tr>'; for (let j = 0; j < n; j++) html += `<td${i > j ? ' class="interactive" title="' + EDIT_KEY_HINT + '"' : ''} data-i="${i}" data-j="${j}">${corr[i][j].toFixed(3)}</td>`; html += '</tr>'; }
    const tbl = $('corr-table'); tbl.classList.toggle('dense', n > 14); tbl.innerHTML = html;
    cellRefs = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => tbl.querySelector(`td[data-i="${i}"][data-j="${j}"]`)));
    barIndef = false; measureBars();
  }
  function measureBars() {
    const canvas = $('bar-canvas'), tbl = $('corr-table');
    if (!canvas || !tbl || !cellRefs) { barGeom = null; return; }
    // Self-heal: the "studio" flex layout settles the matrix's box asynchronously
    // (the body canvas sizes via its own ResizeObserver after buildTable runs), and
    // later reflows — scrollbars, viewport changes — can shift it again. Re-measure
    // whenever the table's rendered box changes, so the bars never drift out of sync
    // with the cells (matches the body canvas's own observer; window-resize handles dpr).
    if (!barRO && window.ResizeObserver) { barRO = new ResizeObserver(() => measureBars()); barRO.observe(tbl); }
    const n = corr.length, cssW = tbl.offsetWidth, cssH = tbl.offsetHeight, dpr = window.devicePixelRatio || 1;
    canvas.style.left = tbl.offsetLeft + 'px'; canvas.style.top = tbl.offsetTop + 'px';
    canvas.style.width = cssW + 'px'; canvas.style.height = cssH + 'px';
    canvas.width = Math.max(1, Math.round(cssW * dpr)); canvas.height = Math.max(1, Math.round(cssH * dpr));
    canvas.getContext('2d').setTransform(dpr, 0, 0, dpr, 0, 0);
    const tRect = tbl.getBoundingClientRect(), geom = [];
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) { const r = cellRefs[i][j].getBoundingClientRect(); geom.push({ i, j, x: r.left - tRect.left, y: r.top - tRect.top, w: r.width, h: r.height }); }
    barGeom = geom; drawBars();
  }
  function drawBars() {
    const canvas = $('bar-canvas'); if (!canvas || !barGeom || !info) return;
    const ctx = canvas.getContext('2d');
    // Clear the whole backing store regardless of the active dpr transform, so a
    // dpr change between measure and draw can't leave stale pixels behind.
    ctx.save(); ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.clearRect(0, 0, canvas.width, canvas.height); ctx.restore();
    if (barIndef) return;
    ctx.fillStyle = BAR_COLOR;
    for (const g of barGeom) { const f = info[g.i] && info[g.i][g.j]; if (!f) continue; const x0 = g.x + (clamp(-1, 1, f.lo) + 1) / 2 * g.w, x1 = g.x + (clamp(-1, 1, f.hi) + 1) / 2 * g.w; ctx.fillRect(x0, g.y + g.h - BAR_H, Math.max(0, x1 - x0), BAR_H); }
  }
  function refreshCells(indef) {
    if (!cellRefs) return;
    if (indef !== barIndef) { barIndef = indef; const n = corr.length; for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) cellRefs[i][j].style.backgroundColor = indef ? RED : (i > j ? '' : '#e6e1dc'); }
    drawBars();
  }
  function markSelected() {
    document.querySelectorAll('#corr-table td.selected').forEach(td => { td.classList.remove('selected'); td.removeAttribute('data-axis'); });
    sel.forEach((e, k) => { if (!e) return; cellRefs[e.i][e.j].classList.add('selected'); cellRefs[e.j][e.i].classList.add('selected'); cellRefs[e.i][e.j].setAttribute('data-axis', ['x', 'y', 'z'][k]); });
  }
  function setEntry(e, v) {
    corr[e.i][e.j] = corr[e.j][e.i] = v;
    if (cellRefs[e.i][e.j]) cellRefs[e.i][e.j].textContent = v.toFixed(3);
    if (cellRefs[e.j][e.i]) cellRefs[e.j][e.i].textContent = v.toFixed(3);
  }

  // ── 2D slice (draggable marker; depth = sel[2]'s value) ──────────────────
  function axesSVG(labels) {
    const X0 = bX(-1), X1 = bX(1), Y0 = bY(-1), Y1 = bY(1), zx = bX(0), zy = bY(0);
    const tk = (v, h) => h
      ? `<text x="${bX(v)}" y="${Y0 + 13}" text-anchor="middle" font-size="11" fill="#888">${v}</text>`
      : `<text x="${X0 - 5}" y="${bY(v) + 4}" text-anchor="end" font-size="11" fill="#888">${v}</text>`;
    return `<rect x="${X0}" y="${Y1}" width="${piW}" height="${piH}" fill="none" stroke="#2d2d2d"/>
<line x1="${zx}" y1="${Y1}" x2="${zx}" y2="${Y0}" stroke="#ddd"/><line x1="${X0}" y1="${zy}" x2="${X1}" y2="${zy}" stroke="#ddd"/>
${tk(-1, true)}${tk(0, true)}${tk(1, true)}${tk(-1, false)}${tk(0, false)}${tk(1, false)}${labels || ''}`;
  }
  function drawSlice() {
    const host = $('slice-svg');
    if (!ready()) { host.innerHTML = `<svg viewBox="0 0 ${PW} ${PH}" width="100%">${axesSVG('')}</svg>`; $('slice-info').textContent = ''; updateSliceLine(null); return; }
    const z = val(sel[2]);
    const lab = `<text x="${bX(0)}" y="${PH - 2}" text-anchor="middle" font-size="11" fill="#555">${labelFor(sel[0])}</text>
<text x="11" y="${bY(0)}" text-anchor="middle" font-size="11" fill="#555" transform="rotate(-90 11 ${bY(0)})">${labelFor(sel[1])}</text>`;
    const M = corr.map(r => r.slice()); M[sel[2].i][sel[2].j] = M[sel[2].j][sel[2].i] = z;
    const reg = regionOn(M, sel[0], sel[1], true);
    updateSliceLine(reg, z);
    let path = '';
    if (reg.top.length) { const pt = p => `${bX(toDisp(sel[0], p.x)).toFixed(1)} ${bY(toDisp(sel[1], p.y)).toFixed(1)}`; path = `<path d="M ${reg.top.map(pt).join(' L ')} L ${reg.bot.slice().reverse().map(pt).join(' L ')} Z" fill="${FILL}" stroke="${STROKE}" stroke-width="1.5"/>`; }
    const ok = lam() >= INDEF, col = ok ? DARK : BAD;
    const mxv = bX(toDisp(sel[0], val(sel[0]))), myv = bY(toDisp(sel[1], val(sel[1])));
    host.innerHTML = `<svg id="sl-el" viewBox="0 0 ${PW} ${PH}" width="100%" style="touch-action:none">
${axesSVG(lab)}${path}
<line id="sl-mv" x1="${mxv}" x2="${mxv}" y1="${bY(1)}" y2="${bY(-1)}" stroke="${col}" stroke-width="1" stroke-dasharray="2,2"/>
<line id="sl-mh" x1="${bX(-1)}" x2="${bX(1)}" y1="${myv}" y2="${myv}" stroke="${col}" stroke-width="1" stroke-dasharray="2,2"/>
<circle id="sl-dot" cx="${mxv}" cy="${myv}" r="6" fill="${col}" style="cursor:grab"/></svg>`;
    $('sl-el').addEventListener('pointerdown', startDrag);
    $('slice-info').innerHTML = `Slice at <b>${labelFor(sel[2])} = ${toDisp(sel[2], z).toFixed(3)}</b>`;
  }

  // ── 3D scene (Three.js) ──────────────────────────────────────────────────
  let three = null;
  function initGL() {
    const canvas = $('gl-canvas');
    if (typeof THREE === 'undefined' || !THREE.WebGLRenderer) { $('gl-msg').hidden = false; $('gl-msg').textContent = 'WebGL / Three.js unavailable in this browser.'; $('spin-toggle').hidden = true; return; }
    let renderer;
    try { renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true }); }
    catch (e) { $('gl-msg').hidden = false; $('gl-msg').textContent = 'Could not create a WebGL context.'; $('spin-toggle').hidden = true; return; }
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 100); camera.up.set(0, 0, 1); camera.position.set(3.25, 3.25, 2.4);
    const controls = new THREE.OrbitControls(camera, canvas); controls.enableDamping = true; controls.dampingFactor = 0.12; controls.target.set(0, 0, 0);
    controls.autoRotate = true; controls.autoRotateSpeed = 1.0;   // gentle idle spin; user drag/zoom still works (and overrides it while interacting)
    scene.add(new THREE.AmbientLight(0xffffff, 0.55));
    const dl = new THREE.DirectionalLight(0xffffff, 0.75); dl.position.set(2, 3, 4); scene.add(dl);
    const dl2 = new THREE.DirectionalLight(0xffffff, 0.3); dl2.position.set(-3, -2, -1); scene.add(dl2);
    scene.add(new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(2, 2, 2)), new THREE.LineBasicMaterial({ color: 0xbbb3aa })));
    const axis = (a, b, col) => scene.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(...a), new THREE.Vector3(...b)]), new THREE.LineBasicMaterial({ color: col })));
    axis([-1, -1, -1], [1, -1, -1], 0xb23b3b); axis([-1, -1, -1], [-1, 1, -1], 0x2a7a2a); axis([-1, -1, -1], [-1, -1, 1], 0x2f5fb0);
    const group = new THREE.Group(); scene.add(group);
    const marker = new THREE.Mesh(new THREE.SphereGeometry(0.045, 16, 16), new THREE.MeshBasicMaterial({ color: 0x222222 })); scene.add(marker);
    const plane = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), new THREE.MeshBasicMaterial({ color: 0x2f5fb0, transparent: true, opacity: 0.10, side: THREE.DoubleSide, depthWrite: false })); scene.add(plane);
    three = { renderer, scene, camera, controls, group, bodyMesh: null, marker, plane, sliceLine: null };
    resizeGL();
    // the canvas is flex-sized (fills the views column height) → keep the renderer in sync
    if (window.ResizeObserver) { try { new ResizeObserver(() => resizeGL()).observe(canvas); } catch (_) {} }
    (function loop() { requestAnimationFrame(loop); controls.update(); renderer.render(scene, camera); })();
  }
  function resizeGL() {
    if (!three) return;
    const c = $('gl-canvas'), w = c.clientWidth || 320, h = c.clientHeight || 320;
    three.renderer.setPixelRatio(window.devicePixelRatio || 1); three.renderer.setSize(w, h, false);
    three.camera.aspect = w / h; three.camera.updateProjectionMatrix();
  }
  // ── Body mesh = the 3-D CONVEX HULL of the slice boundaries ────────────────
  // The body is a convex spectrahedron, so it equals the convex hull of its slice
  // boundary points (each z-slice is a convex 2-D region whose boundary lies on the
  // body's surface). We clean each slice to its exact 2-D convex hull — removing the
  // swept-grid's collinear points and any sub-pixel non-convex noise — then build the
  // 3-D hull with THREE.ConvexGeometry (a robust QuickHull3D). This renders cleanly for
  // boxes, curved/deformed prisms, spheres and cylinders alike. (The earlier angular
  // loft sampled each slice by fixed angles from its centroid and bridged like-indexed
  // points across slices, which smeared the corners of boxy cross-sections → the jagged
  // surfaces the user reported.)
  function hull2D(pts) {                                // Andrew's monotone chain → CCW hull, collinear removed
    const P = pts.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    if (P.length < 3) return P;
    const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    const lo = []; for (const p of P) { while (lo.length >= 2 && cross(lo[lo.length - 2], lo[lo.length - 1], p) <= 1e-9) lo.pop(); lo.push(p); }
    const hi = []; for (let i = P.length - 1; i >= 0; i--) { const p = P[i]; while (hi.length >= 2 && cross(hi[hi.length - 2], hi[hi.length - 1], p) <= 1e-9) hi.pop(); hi.push(p); }
    lo.pop(); hi.pop(); return lo.concat(hi);
  }
  // Resample a convex polygon UNIFORMLY by arc length (so curved slices give a regular
  // triangulation → smooth normals), while snapping in each sharp corner exactly (so
  // boxes keep crisp corners). `regionOn` samples a slice densely on its sides and
  // sparsely top/bottom; feeding that straight in made the hull irregular → mottled.
  function resampleSlice(v, ds) {
    const n = v.length; if (n < 4) return v;
    const seg = [], cum = [0]; let per = 0;
    for (let i = 0; i < n; i++) { const a = v[i], b = v[(i + 1) % n]; const l = Math.hypot(b[0] - a[0], b[1] - a[1]); seg.push(l); per += l; cum.push(per); }
    if (per < 1e-9) return v;
    const m = Math.max(8, Math.round(per / ds)), out = []; let si = 0;
    for (let k = 0; k < m; k++) {
      const target = k * per / m;
      while (si < n - 1 && cum[si + 1] < target) si++;
      const a = v[si], b = v[(si + 1) % n], t = seg[si] > 1e-12 ? (target - cum[si]) / seg[si] : 0;
      out.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]);
    }
    for (let i = 0; i < n; i++) {                       // keep sharp corners exactly (turn > ~25°)
      const a = v[(i - 1 + n) % n], b = v[i], c = v[(i + 1) % n];
      const ux = b[0] - a[0], uy = b[1] - a[1], vx = c[0] - b[0], vy = c[1] - b[1];
      const lu = Math.hypot(ux, uy), lv = Math.hypot(vx, vy);
      if (lu > 1e-9 && lv > 1e-9 && (ux * vx + uy * vy) / (lu * lv) < 0.906) out.push(b);
    }
    return out;
  }
  function buildBodyGeometry() {
    if (!body || body.slices.length < 2 || typeof THREE.ConvexGeometry === 'undefined') return null;
    const pts = [];
    for (const sl of body.slices) {
      const ring = [];
      for (const p of sl.top) ring.push([p.x, p.y]);
      for (const p of sl.bot) ring.push([p.x, p.y]);
      for (const v of resampleSlice(hull2D(ring), 0.07)) pts.push(new THREE.Vector3(v[0], v[1], sl.z));
    }
    if (pts.length < 4) return null;
    try { const g = new THREE.ConvexGeometry(pts); creaseNormals(g, 50); return g; }
    catch (e) { return null; }
  }
  // ConvexGeometry gives FLAT per-face normals → curved faces look patchy (each facet a
  // different shade). Recompute normals with a crease angle: at each welded vertex, a
  // face-corner averages only the incident faces within `deg` of its own normal. Near-
  // coplanar facets (gentle curvature) blend → smooth surfaces; genuine sharp edges
  // (the cube's 90° edges, cylinder rim, elliptope edges) exceed the threshold and stay
  // split → crisp. (Plain smoothing would round the cube; flat shading looks faceted.)
  function creaseNormals(geom, deg) {
    const pos = geom.getAttribute('position'), nv = pos.count, nt = nv / 3, fn = new Float32Array(nt * 3), fa = new Float32Array(nt);
    for (let t = 0; t < nt; t++) {
      const i = 3 * t;
      const ux = pos.getX(i + 1) - pos.getX(i), uy = pos.getY(i + 1) - pos.getY(i), uz = pos.getZ(i + 1) - pos.getZ(i);
      const vx = pos.getX(i + 2) - pos.getX(i), vy = pos.getY(i + 2) - pos.getY(i), vz = pos.getZ(i + 2) - pos.getZ(i);
      const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
      const L = Math.hypot(nx, ny, nz) || 1; fn[3 * t] = nx / L; fn[3 * t + 1] = ny / L; fn[3 * t + 2] = nz / L; fa[t] = 0.5 * L;  // triangle area
    }
    const buckets = new Map();                              // weld by exact position (hull vertices are shared/identical)
    for (let i = 0; i < nv; i++) {
      const k = pos.getX(i).toFixed(5) + ',' + pos.getY(i).toFixed(5) + ',' + pos.getZ(i).toFixed(5);
      let b = buckets.get(k); if (!b) buckets.set(k, b = []); b.push(i);
    }
    const cosT = Math.cos(deg * Math.PI / 180), out = new Float32Array(nv * 3);
    for (const b of buckets.values()) for (const i of b) {
      const t = (i / 3) | 0, mx = fn[3 * t], my = fn[3 * t + 1], mz = fn[3 * t + 2];
      let sx = 0, sy = 0, sz = 0;                          // area-weighted average of incident faces within the crease angle
      for (const j of b) { const tj = (j / 3) | 0, jx = fn[3 * tj], jy = fn[3 * tj + 1], jz = fn[3 * tj + 2]; if (jx * mx + jy * my + jz * mz >= cosT) { const w = fa[tj]; sx += jx * w; sy += jy * w; sz += jz * w; } }
      const L = Math.hypot(sx, sy, sz) || 1; out[3 * i] = sx / L; out[3 * i + 1] = sy / L; out[3 * i + 2] = sz / L;
    }
    geom.setAttribute('normal', new THREE.Float32BufferAttribute(out, 3));
  }
  function rebuildBodyMesh() {
    if (!three) return;
    if (three.bodyMesh) { three.group.remove(three.bodyMesh); three.bodyMesh.geometry.dispose(); three.bodyMesh = null; }
    const g = buildBodyGeometry(); if (!g) return;
    // FrontSide (not DoubleSide): the body is convex, so each ray hits the near surface once.
    // Drawing the back faces too (through the transparency) double-shaded the surface and made it
    // look mottled — "different shades of blue". Front-only gives a clean frosted-glass shell.
    const mat = new THREE.MeshPhongMaterial({ color: 0x2f6fd0, transparent: true, opacity: 0.5, shininess: 30, side: THREE.FrontSide, depthWrite: false });
    three.bodyMesh = new THREE.Mesh(g, mat); three.group.add(three.bodyMesh);
  }
  function updateMarkerPlane() {
    if (!three || !ready()) return;
    const ok = lam() >= INDEF;
    three.marker.position.set(toDisp(sel[0], val(sel[0])), toDisp(sel[1], val(sel[1])), toDisp(sel[2], val(sel[2])));
    three.marker.material.color.set(ok ? 0x222222 : 0xb02828);
    three.plane.position.set(0, 0, toDisp(sel[2], val(sel[2])));
    $('depth-slider').style.setProperty('--thumb', ok ? DARK : BAD);   // thumb mirrors the 2-D dot's valid/invalid colour
  }

  // Where the depth plane cuts the body: the cross-section boundary at the current depth,
  // drawn as a thin closed loop lying on the plane (the 3-D twin of the 2-D slice outline).
  // reg/z are handed straight from drawSlice, which already computed them — so this only
  // fires when the cross-section actually changes (depth or the rest of the matrix), not on
  // x/y drags (the slice shape is invariant under its own two axis values).
  function updateSliceLine(reg, z) {
    if (!three) return;
    if (three.sliceLine) { three.scene.remove(three.sliceLine); three.sliceLine.geometry.dispose(); three.sliceLine.material.dispose(); three.sliceLine = null; }
    if (!reg || !reg.top.length) return;
    const pts = [];
    for (const p of reg.top) pts.push(new THREE.Vector3(p.x, p.y, z));
    for (let k = reg.bot.length - 1; k >= 0; k--) pts.push(new THREE.Vector3(reg.bot[k].x, reg.bot[k].y, z));
    pts.push(pts[0].clone());   // close the loop
    const mat = new THREE.LineBasicMaterial({ color: SLICE_LINE, transparent: true, opacity: 0.7, depthTest: false });
    three.sliceLine = new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), mat);
    three.sliceLine.renderOrder = 3;   // draw over the translucent body so the cross-section reads clearly
    three.scene.add(three.sliceLine);
  }

  // ── Drag the slice marker (edits sel[0], sel[1]) ─────────────────────────
  function startDrag(e) {
    if (!ready()) return; e.preventDefault();
    const svg = e.currentTarget;
    drag = { svg, rect: svg.getBoundingClientRect(), id: e.pointerId };
    try { svg.setPointerCapture(e.pointerId); } catch (_) {}
    svg.addEventListener('pointermove', onDrag); svg.addEventListener('pointerup', endDrag); svg.addEventListener('pointercancel', endDrag);
    document.body.classList.add('dragging'); onDrag(e);
  }
  function onDrag(e) {
    if (!drag) return; const r = drag.rect;
    const dx = clamp(-1, 1, ((e.clientX - r.left) / r.width * PW - PPAD) / piW * 2 - 1);
    const dy = clamp(-1, 1, (1 - ((e.clientY - r.top) / r.height * PH - PPAD) / piH) * 2 - 1);
    setEntry(sel[0], clamp(-1, 1, toRaw(sel[0], dx))); setEntry(sel[1], clamp(-1, 1, toRaw(sel[1], dy)));  // display→raw; body invariant under its own values → no re-mesh
    info = feasibleRangeAll(corr); refreshCells(lam() < INDEF); updateMarkerPlane();
    const ok = lam() >= INDEF, col = ok ? DARK : BAD;
    const set = (id, a) => { const el = $(id); if (el) for (const k in a) el.setAttribute(k, a[k]); };
    set('sl-mv', { x1: bX(dx), x2: bX(dx), stroke: col }); set('sl-mh', { y1: bY(dy), y2: bY(dy), stroke: col });
    set('sl-dot', { cx: bX(dx), cy: bY(dy), fill: col });
  }
  function endDrag() {
    if (!drag) return;
    try { drag.svg.releasePointerCapture(drag.id); } catch (_) {}
    drag.svg.removeEventListener('pointermove', onDrag); drag.svg.removeEventListener('pointerup', endDrag); drag.svg.removeEventListener('pointercancel', endDrag);
    drag = null; document.body.classList.remove('dragging');
  }

  // ── Refresh / selection ──────────────────────────────────────────────────
  function recomputeBody() { body = ready() ? computeBody(sel[0], sel[1], sel[2]) : null; }
  function syncDepth() {
    const s = $('depth-slider');
    if (!ready()) { s.disabled = true; return; }
    s.disabled = false; s.min = -1; s.max = 1; s.step = 0.001; s.value = toDisp(sel[2], val(sel[2]));
  }
  function drawAxesLabel() {
    if (!ready()) { $('gl-axes').innerHTML = ''; return; }
    $('gl-axes').innerHTML = `<span class="ax-x">X ${labelFor(sel[0])}</span> · <span class="ax-y">Y ${labelFor(sel[1])}</span> · <span class="ax-z">Z ${labelFor(sel[2])}</span>`;
  }
  function rebuildAll() {            // selection / reset — re-mesh + bars + active button
    markSelected(); recomputeBody(); rebuildBodyMesh(); updateMarkerPlane(); syncDepth(); drawSlice();
    refreshCells(eigvalsh(corr).some(x => x < INDEF)); syncGraphButtons();
  }
  function onCellClick(i, j) {
    if (i <= j) return;
    const e = entry(i, j), k = sel.findIndex(x => sameE(x, e));
    if (k >= 0) { sel[k] = null; selSeq[k] = 0; }      // deselect → free THIS slot; the other axes stay put
    else {
      let slot = sel.indexOf(null);                    // fill the first empty axis slot
      if (slot < 0) slot = [0, 1, 2].reduce((a, b) => selSeq[b] < selSeq[a] ? b : a, 0);  // all full → evict the oldest in place
      fillSlot(slot, e);
    }
    drawAxesLabel(); rebuildAll();
  }

  // ── Type a value into a cell (⌘/Ctrl-click) — mirrors the 1-D explorer ─────
  // Modifier-click any interactive cell to type its correlation; Enter/blur commits,
  // Esc cancels. Editing a free (selected) entry just moves the marker; editing a
  // fixed entry deforms the body, so it re-meshes.
  function editCell(i, j) {
    if (i <= j || !cellRefs) return;
    const td = cellRefs[i][j];
    if (!td || td.querySelector('input')) return;             // already editing this cell
    const e = entry(i, j);
    const input = document.createElement('input');
    input.type = 'text'; input.className = 'cell-input'; input.size = 1;
    input.value = corr[i][j].toFixed(3);
    input.setAttribute('inputmode', 'decimal');
    input.setAttribute('aria-label', `Entry (${i + 1},${j + 1})`);
    let done = false;
    const finish = commit => {
      if (done) return; done = true;
      input.removeEventListener('blur', onBlur);
      const v = parseFloat(input.value);
      const ok = commit && input.value !== corr[i][j].toFixed(3) && isFinite(v) && v >= -1 && v <= 1;
      if (input.parentNode) input.remove();                   // back to a plain text cell
      if (ok) applyValue(e, v); else restoreText(e);
    };
    const onBlur = () => finish(true);
    input.addEventListener('blur', onBlur);
    input.addEventListener('focus', () => input.select());
    input.addEventListener('keydown', ev => {
      if (ev.key === 'Enter') { ev.preventDefault(); input.blur(); }
      else if (ev.key === 'Escape') { ev.preventDefault(); finish(false); }
    });
    input.addEventListener('input', () => { const m = cellRefs[e.j][e.i]; if (m) m.textContent = input.value; });  // echo into the symmetric cell
    td.textContent = ''; td.appendChild(input); input.focus();
  }
  function restoreText(e) {
    if (cellRefs[e.i][e.j]) cellRefs[e.i][e.j].textContent = corr[e.i][e.j].toFixed(3);
    if (cellRefs[e.j][e.i]) cellRefs[e.j][e.i].textContent = corr[e.j][e.i].toFixed(3);
  }
  function applyValue(e, v) {
    setEntry(e, v);
    if (!sel.some(s => sameE(s, e))) { recomputeBody(); rebuildBodyMesh(); }   // a fixed entry changed → the body deforms
    info = feasibleRangeAll(corr);
    refreshCells(lam() < INDEF);
    updateMarkerPlane(); syncDepth(); drawSlice();
  }

  $('depth-slider').addEventListener('input', e => {
    if (!ready()) return;
    const v = +e.target.value;
    setEntry(sel[2], clamp(-1, 1, toRaw(sel[2], v))); info = feasibleRangeAll(corr);
    refreshCells(lam() < INDEF); updateMarkerPlane(); drawSlice();  // body unchanged
  });

  // ── Graph-structure presets (a button just SELECTS its three cells) ───────
  // The three free entries are three EDGES of a graph on the six indices. Up to
  // isomorphism there are exactly five simple 3-edge graphs with no isolated
  // vertex; on the identity each is a clean convex body (verified vs a λ_min
  // grid). A button selects its cells in [x, y, depth] order — it does NOT touch
  // the matrix values, so the body reflects whatever the matrix currently holds
  // (identity by default → the clean shapes; edit entries to deform them).
  const PRESETS = {
    triangle: [[1, 0], [2, 0], [2, 1]],   // K₃    → elliptope (Cayley's inflated tetrahedron)
    star:     [[1, 0], [2, 0], [3, 0]],   // K₁,₃  → unit ball  x²+y²+z² ≤ 1
    path:     [[1, 0], [3, 2], [2, 1]],   // P₄    → cushion (square girdle + rounded caps; depth = the middle edge)
    pathedge: [[1, 0], [2, 1], [4, 3]],   // P₃∪K₂ → cylinder (depth = the lone edge)
    matching: [[1, 0], [3, 2], [5, 4]],   // 3·K₂  → cube [−1,1]³
  };
  const identity = () => Array.from({ length: N }, (_, i) => Array.from({ length: N }, (_, j) => i === j ? 1 : 0));
  const cellsKey = arr => arr.map(([i, j]) => i > j ? `${i},${j}` : `${j},${i}`).sort().join('|');
  // Highlight the picker to show the current selection's graph type. Two levels:
  //  · `.active`        — the selection is EXACTLY this button's canonical cells (strong).
  //  · `.active-muted`  — the selection is this button's TYPE but on non-canonical cells
  //                       (soft) — so a manual triangle-anywhere still reads as "Triangle".
  function syncGraphButtons() {
    const cur = sel.filter(Boolean).map(key).sort().join('|'), type = classifyGraph(sel);
    document.querySelectorAll('.graph-btn').forEach(b => {
      const g = b.dataset.graph, canonical = cellsKey(PRESETS[g]) === cur;
      b.classList.toggle('active', canonical);
      b.classList.toggle('active-muted', !canonical && g === type);
    });
  }
  // Classify ANY three-edge selection into one of the five isomorphism types — works
  // for non-canonical placements too — keyed by (#vertices : sorted degree sequence).
  const GRAPH_SIG = { '3:2,2,2': 'triangle', '4:3,1,1,1': 'star', '4:2,2,1,1': 'path', '5:2,1,1,1,1': 'pathedge', '6:1,1,1,1,1,1': 'matching' };
  function classifyGraph(entries) {
    const es = entries.filter(Boolean);
    if (es.length < 3) return null;
    const deg = new Map();
    for (const e of es) { deg.set(e.i, (deg.get(e.i) || 0) + 1); deg.set(e.j, (deg.get(e.j) || 0) + 1); }
    return GRAPH_SIG[deg.size + ':' + [...deg.values()].sort((a, b) => b - a).join(',')] || null;
  }
  function selectGraph(name) {
    const cells = PRESETS[name]; if (!cells) return;
    setSel(cells.map(([i, j]) => entry(i, j)));       // select cells only; matrix untouched
    drawAxesLabel(); rebuildAll();
  }
  function resetIdentity() {                           // values → 6×6 identity; selection kept
    corr = identity(); info = feasibleRangeAll(corr);
    for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) cellRefs[i][j].textContent = corr[i][j].toFixed(3);
    rebuildAll();
  }

  // ── Wiring ───────────────────────────────────────────────────────────────
  document.querySelectorAll('.graph-btn').forEach(b => b.addEventListener('click', () => selectGraph(b.dataset.graph)));
  $('reset-btn').addEventListener('click', resetIdentity);
  $('spin-toggle').addEventListener('click', () => {
    if (!three) return;
    const on = three.controls.autoRotate = !three.controls.autoRotate;
    const b = $('spin-toggle');
    b.innerHTML = on ? '&#9208;' : '&#9654;';                          // ⏸ pause : ▶ play
    b.title = b.ariaLabel = on ? 'Pause rotation' : 'Play rotation';
  });
  // Plain click selects an axis (instant); a modifier-click (⌘/Ctrl, or Alt) types a value
  // instead. A distinct gesture means no click/double-click disambiguation — selection stays snappy.
  $('corr-table').addEventListener('click', e => {
    if (e.target.tagName === 'INPUT') return;                 // typing in an open editor — leave to the browser
    const td = e.target.closest('td.interactive'); if (!td) return;
    const i = +td.dataset.i, j = +td.dataset.j;
    if (e.ctrlKey || e.metaKey || e.altKey) { e.preventDefault(); editCell(i, j); return; }   // modifier-click → type a value
    if (e.detail > 1) return;                                 // ignore the 2nd click of a double-click, so a stray double-click can't toggle the axis twice
    onCellClick(i, j);                                        // plain click → pick the axis
  });
  let rRAF = 0; window.addEventListener('resize', () => { if (rRAF) return; rRAF = requestAnimationFrame(() => { rRAF = 0; if (corr) { measureBars(); resizeGL(); } }); });

  // ── Init: 6×6 identity, triangle selected (→ the elliptope on first paint) ─
  corr = identity(); info = feasibleRangeAll(corr);
  const g0 = new URLSearchParams(location.search).get('graph');
  setSel((PRESETS[g0] || PRESETS.triangle).map(([i, j]) => entry(i, j)));
  $('results').hidden = false; buildTable(); initGL();
  drawAxesLabel(); rebuildAll();
})();
