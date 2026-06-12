// Explorer page logic: matrix generation, click-to-inspect, type-or-drag to
// vary, and the λ_min curve. All computation runs synchronously in the
// browser; no network round-trips.
(function () {
  'use strict';

  const INDEFINITE_THRESHOLD = -1e-10;

  // SVG geometry constants (viewBox units; SVG scales to container width).
  const NL_W = 420, NL_H = 36, NL_PAD = 44;
  const CV_W = 420, CV_H = 120, CV_PAD_L = 44, CV_PAD_R = 44, CV_PAD_T = 14, CV_PAD_B = 22;

  // ── State ─────────────────────────────────────────────────────────────────
  let currentCorr    = null;
  let currentInfo    = null;   // info[i][j] = {lo, hi, alpha, beta, gamma} for i != j
  let currentEigvals = null;   // sorted-descending spectrum of currentCorr; refreshed each setSelectedValue
  let pristineCorr   = null;   // snapshot at last Generate / Upload, for Reset
  let isEdited       = false;  // true once user drags or types a value
  let selectedCell   = null;   // {i, j}
  let dragState      = null;   // {i, j, pointerId, axisLeft, axisRight}
  let currentCurve   = null;   // {cs, eigvals} from eigvalsCurve; eigvals[k] is n-vector sorted desc
  let curveTransform = null;   // {toX, toY} for the λ_min curve
  let cellRefs       = null;   // cellRefs[i][j] = td element (null for diagonal). Cached on buildCorrTable.
  let viewMode       = 'raw';  // 'raw' = slider in c ∈ [-1,1] with feasible band; 'partial' = slider in ρ ∈ [-1,1], whole axis feasible

  // Coordinate helpers for the selected cell. In partial mode the slider/curve
  // display ρ ∈ [-1, 1] (always feasible); in raw mode they display c. The
  // cell value itself is always raw c.
  function partialK() {
    if (!selectedCell) return 1;
    const { i, j } = selectedCell;
    const info = currentInfo[i][j];
    return Math.sqrt(Math.max(0, (1 - info.alpha) * (1 - info.gamma)));
  }
  function cToAxis(c) {
    if (viewMode !== 'partial' || !selectedCell) return c;
    const { i, j } = selectedCell;
    const K = partialK();
    return K > 1e-12 ? (c - currentInfo[i][j].beta) / K : 0;
  }
  function axisToC(a) {
    if (viewMode !== 'partial' || !selectedCell) return a;
    const { i, j } = selectedCell;
    return currentInfo[i][j].beta + a * partialK();
  }

  // ── Persistence ──────────────────────────────────────────────────────────
  // Save the live matrix + pristine snapshot + selection + view mode to
  // localStorage so the user's work survives a page reload or navigating to
  // /theory and back. Saves are debounced so a 60 Hz drag doesn't hammer
  // synchronous storage. Reads are guarded with try/catch — Safari private
  // mode and quota errors throw and we silently fall through to defaults.
  const STORAGE_KEY = 'correlation-explorer:state';
  let saveTimer = null;
  function saveState() {
    if (!currentCorr) return;
    try {
      localStorage.setItem(STORAGE_KEY, JSON.stringify({
        matrix: currentCorr,
        pristine: pristineCorr,
        selectedCell,
        viewMode,
      }));
    } catch (e) { /* storage unavailable or full — drop silently */ }
  }
  function saveStateThrottled() {
    if (saveTimer) clearTimeout(saveTimer);
    saveTimer = setTimeout(() => { saveTimer = null; saveState(); }, 200);
  }
  function loadState() {
    let state;
    try {
      const raw = localStorage.getItem(STORAGE_KEY);
      if (!raw) return null;
      state = JSON.parse(raw);
    } catch (e) { return null; }
    // Validate matrix shape — same guards as CSV upload, since we trust
    // ourselves to have written it but not whatever else might be in storage.
    if (!state || !Array.isArray(state.matrix)) return null;
    const m = state.matrix, n = m.length;
    if (n < 2 || n > 20) return null;
    for (const row of m) {
      if (!Array.isArray(row) || row.length !== n) return null;
      for (const v of row) if (typeof v !== 'number' || !isFinite(v)) return null;
    }
    return state;
  }

  // ── Compute & render ──────────────────────────────────────────────────────
  function computeInfo(corr) {
    const n = corr.length;
    const info = zeros(n, n).map(row => row.map(_ => null));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const r = feasibleRange(corr, i, j);
        info[i][j] = info[j][i] = r;
      }
    }
    return info;
  }

  function renderResults(corr, opts = {}) {
    currentCorr    = corr;
    currentInfo    = computeInfo(corr);
    currentEigvals = eigvalsh(corr);
    if (!opts.preservePristine) {
      pristineCorr = corr.map(row => row.slice());
    }
    isEdited       = false;
    selectedCell   = null;
    dragState      = null;
    currentCurve   = null;
    curveTransform = null;
    document.getElementById('error').textContent = '';
    buildCorrTable();
    refreshCellColors(currentEigvals.some(x => x < INDEFINITE_THRESHOLD));
    document.getElementById('results').hidden = false;
    renderEmptyPanel();
    updateActionButtons();
    saveState();
  }

  // Slider font size in viewBox units, picked so SVG text renders at the same
  // CSS px size as the matrix cells. Uses the slider container's rendered
  // width — that's the box the SVG actually fills (narrower than #results
  // since the mode toggle takes up a left column in the panel grid).
  function computeSliderFontSize() {
    const refWidth = document.getElementById('nl-svg').clientWidth
                  || document.getElementById('results').clientWidth;
    const cellEl = document.querySelector('#corr-table td');
    if (!refWidth || !cellEl) return 10;
    return parseFloat(getComputedStyle(cellEl).fontSize) * NL_W / refWidth;
  }

  // Placeholder slider + curves (axes only, no data). Keeps the panel's
  // layout stable so cell clicks don't shift the page.
  function renderEmptyPanel() {
    const fontSize = computeSliderFontSize();
    const info = document.getElementById('nl-info');
    info.hidden = false;
    info.innerHTML = '<span style="color:#888">Click a cell to inspect its feasible range.</span>';
    document.getElementById('mode-toggle').hidden = true;
    document.getElementById('nl-svg').innerHTML    = makeEmptyNumberLineSVG(fontSize);
    document.getElementById('nl-curve').innerHTML  = makeEmptyCurveSVG(fontSize);
  }

  function makeEmptyNumberLineSVG(fontSize) {
    const W = NL_W, H = NL_H, padL = NL_PAD, padR = NL_PAD;
    const fs = fontSize.toFixed(1);
    const zeroX = padL + (W - padL - padR) / 2;
    const f = n => n.toFixed(2);
    return `<svg viewBox="0 0 ${W} ${H}" width="100%"
       style="display:block;overflow:visible;height:auto;font-family:inherit">
  <line x1="${padL}" y1="19" x2="${W-padR}" y2="19" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${padL}"   y1="15" x2="${padL}"   y2="23" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${W-padR}" y1="15" x2="${W-padR}" y2="23" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${f(zeroX)}" y1="16" x2="${f(zeroX)}" y2="22" stroke="#aaa" stroke-width="1"/>
  <text x="${padL}"     y="${H}" text-anchor="middle" font-size="${fs}" fill="#888">-1</text>
  <text x="${f(zeroX)}" y="${H}" text-anchor="middle" font-size="${fs}" fill="#888">0</text>
  <text x="${W-padR}"   y="${H}" text-anchor="middle" font-size="${fs}" fill="#888">+1</text>
</svg>`;
  }

  function makeEmptyCurveSVG(fontSize) {
    const W = CV_W, H = CV_H;
    const padL = CV_PAD_L, padR = CV_PAD_R, padT = CV_PAD_T, padB = CV_PAD_B;
    const innerW = W - padL - padR;
    const toX = c => padL + (c + 1) / 2 * innerW;
    const fs = fontSize.toFixed(1);
    return `<svg viewBox="0 0 ${W} ${H}" width="100%"
       style="display:block;height:auto;font-family:inherit">
  <line x1="${padL}" y1="${padT}"  x2="${padL}"   y2="${H-padB}" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${padL}" y1="${H-padB}" x2="${W-padR}" y2="${H-padB}" stroke="#2d2d2d" stroke-width="1"/>
  <text x="${toX(-1)}" y="${H-padB+14}" text-anchor="middle" font-size="${fs}" fill="#888">-1</text>
  <text x="${toX(0)}"  y="${H-padB+14}" text-anchor="middle" font-size="${fs}" fill="#888">0</text>
  <text x="${toX(1)}"  y="${H-padB+14}" text-anchor="middle" font-size="${fs}" fill="#888">+1</text>
</svg>`;
  }

  // Whole-matrix red flash when the matrix becomes indefinite.
  const RED_INDEFINITE = 'rgba(180,40,40,0.55)';

  function buildCorrTable() {
    const tableEl = document.getElementById('corr-table');
    const n = currentCorr.length;
    // Tighter cells at large n so the matrix fits the container width.
    tableEl.classList.toggle('dense', n > 14);
    let html = '';
    for (let i = 0; i < n; i++) {
      html += '<tr>';
      for (let j = 0; j < n; j++) {
        const v = currentCorr[i][j];
        const cls = (i > j) ? ' class="interactive"' : '';
        html += `<td${cls} data-i="${i}" data-j="${j}">${v.toFixed(3)}</td>`;
      }
      html += '</tr>';
    }
    tableEl.innerHTML = html;
    // Cache <td> references for ALL cells (diagonal included, so they can turn
    // red on indefinite). Avoids 400× querySelector at n=20 on every drag tick.
    cellRefs = new Array(n);
    for (let i = 0; i < n; i++) {
      cellRefs[i] = new Array(n);
      for (let j = 0; j < n; j++) {
        cellRefs[i][j] = tableEl.querySelector(`td[data-i="${i}"][data-j="${j}"]`);
      }
    }
  }

  // Repaint every cell. In PSD state: off-diagonal cells set --bar-start /
  // --bar-end CSS vars so the bottom gradient shows their feasible interval
  // on the cell's local [-1, 1] axis. In indefinite state: every cell
  // (diagonal included) is painted uniform red, and the bar is hidden.
  function refreshCellColors(indefinite) {
    if (!cellRefs || !currentCorr) return;
    const n = currentCorr.length;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const td = cellRefs[i][j];
        if (!td) continue;
        if (indefinite) {
          td.style.backgroundColor = RED_INDEFINITE;
          td.style.backgroundImage = 'none';  // hide the bar under the red flash
          continue;
        }
        td.style.backgroundColor = '';
        td.style.backgroundImage = '';  // restore CSS-defined gradient
        if (i > j) {  // lower-triangle only — upper and diagonal don't render a bar
          const info = currentInfo[i][j];
          const startPct = ((info.lo + 1) / 2 * 100).toFixed(1);
          const endPct   = ((info.hi + 1) / 2 * 100).toFixed(1);
          td.style.setProperty('--bar-start', startPct + '%');
          td.style.setProperty('--bar-end',   endPct + '%');
        }
      }
    }
  }

  function updateActionButtons() {
    document.getElementById('reset-btn').disabled = !isEdited;
    document.getElementById('download-btn').disabled = !currentCorr;
  }

  function showError(msg) {
    document.getElementById('error').textContent = msg;
  }

  // ── Cell click → show feasible range + λ_min curve ────────────────────────
  // No toggle-off on re-click: clicking the selected cell is now a "focus the
  // input" gesture (handled by the click listener). Deselection happens by
  // selecting a different cell or by regenerating the matrix.
  function onCellClick(i, j) {
    if (selectedCell && selectedCell.i === i && selectedCell.j === j) return;
    clearSelection();
    const td     = cellRefs[i] && cellRefs[i][j];
    const mirror = cellRefs[j] && cellRefs[j][i];
    if (td)     td.classList.add('selected');
    if (mirror) mirror.classList.add('selected');
    selectedCell = { i, j };
    makeCellEditor(i, j);
    showNumberLine(i, j);
    saveState();
  }

  function clearSelection() {
    document.querySelectorAll('#corr-table td.selected').forEach(el => {
      el.classList.remove('selected');
      const i = +el.dataset.i, j = +el.dataset.j;
      el.textContent = currentCorr[i][j].toFixed(3);  // strips any cell-input child
    });
    selectedCell = null;
    renderEmptyPanel();
  }

  // ── Cell input (type-to-edit) ─────────────────────────────────────────────
  // The selected cell holds a transparent <input> so the user can type a
  // value. Drag-to-vary and typing share setSelectedValue; updateCellText
  // routes to input.value when an input is present, textContent otherwise.
  function makeCellEditor(i, j) {
    const td = cellRefs[i] && cellRefs[i][j];
    if (!td) return;
    const input = document.createElement('input');
    input.type = 'text';
    input.value = currentCorr[i][j].toFixed(3);
    input.className = 'cell-input';
    input.size = 1;  // tiny intrinsic width; CSS width:100% lets it fill the cell
    input.setAttribute('inputmode', 'decimal');
    input.setAttribute('aria-label', `Entry (${i},${j})`);
    input.addEventListener('keydown', e => {
      if (e.key === 'Enter') { e.preventDefault(); input.blur(); }
      else if (e.key === 'Escape') { revertCellInput(input); input.blur(); }
    });
    input.addEventListener('blur', () => commitCellInput(input));
    input.addEventListener('focus', () => input.select());
    // Live-mirror keystrokes into the symmetric cell. Full recompute waits
    // until commit (Enter/blur) — this is just text echo.
    input.addEventListener('input', () => {
      const mirror = cellRefs[j] && cellRefs[j][i];
      if (mirror) mirror.textContent = input.value;
    });
    td.textContent = '';
    td.appendChild(input);
  }

  function commitCellInput(input) {
    if (!selectedCell) return;
    const { i, j } = selectedCell;
    // The input shows v.toFixed(3); drag/toggle/etc. all set it programmatically.
    // On blur, if the displayed text matches what the cell would naturally show,
    // the user didn't actually edit — skip the commit. Otherwise blur (e.g.
    // clicking the mode toggle) would re-parse the rounded display string and
    // overwrite the higher-precision c, which at a feasible boundary can push
    // the matrix slightly indefinite.
    if (input.value === currentCorr[i][j].toFixed(3)) return;
    const v = parseFloat(input.value);
    if (!isFinite(v) || v < -1 || v > 1) { revertCellInput(input); return; }
    setSelectedValue(v);
  }

  function revertCellInput(input) {
    if (!selectedCell) return;
    const { i, j } = selectedCell;
    input.value = currentCorr[i][j].toFixed(3);
  }

  function showNumberLine(i, j) {
    const value = currentCorr[i][j];
    const info = currentInfo[i][j];
    const fontSize = computeSliderFontSize();
    renderModeToggle(i, j);
    document.getElementById('nl-svg').innerHTML = makeNumberLineSVG(value, info.lo, info.hi, fontSize);
    // Curve samples differ per mode: raw sweeps c ∈ [-1, 1] (shows the PSD↔
    // indefinite transition); partial sweeps c ∈ [lo, hi] (all feasible),
    // giving a uniform ρ ∈ [-1, 1] when re-mapped.
    currentCurve = (viewMode === 'partial')
      ? eigvalsCurve(currentCorr, i, j, 201, info.lo, info.hi)
      : eigvalsCurve(currentCorr, i, j);
    drawCurve(fontSize);
    attachDragHandlers();
  }

  // Populate the info area with the symbolic readout and show the toggle.
  // Indices are 1-based (math convention); the comma separator appears only
  // when n ≥ 10 (otherwise single digits are unambiguous).
  function renderModeToggle(i, j) {
    const value = currentCorr[i][j];
    const info  = currentInfo[i][j];
    const K = Math.sqrt(Math.max(0, (1 - info.alpha) * (1 - info.gamma)));
    const partial = K > 1e-12 ? (value - info.beta) / K : NaN;
    const partialStr = isFinite(partial) ? partial.toFixed(4) : '—';
    const n = currentCorr.length;
    const sep = n >= 10 ? ',' : '';
    const sub  = `${i + 1}${sep}${j + 1}`;
    const subR = `${sub}|rest`;
    const infoEl = document.getElementById('nl-info');
    infoEl.hidden = false;
    const subI = `${i + 1}|rest`;
    const subJ = `${j + 1}|rest`;
    const tt = {
      rho:     `Correlation between variables ${i + 1} and ${j + 1}.`,
      partial: `Partial correlation: the correlation of variables ${i + 1} and ${j + 1} after removing the linear influence of the other variables. Lies in [-1, 1] when the matrix is positive semidefinite.`,
      alpha:   `Coefficient of determination (R²) from regressing variable ${i + 1} on the other variables (excluding variable ${j + 1}). Lies in [0, 1]: 0 means uncorrelated with them; 1 means exactly a linear combination of them. (Called α on the theory page.)`,
      gamma:   `Coefficient of determination (R²) from regressing variable ${j + 1} on the other variables (excluding variable ${i + 1}). Lies in [0, 1]: 0 means uncorrelated with them; 1 means exactly a linear combination of them. (Called γ on the theory page.)`,
      c:       `Center of the feasible range. Equivalently, the portion of the correlation already explained by the joint dependence of variables ${i + 1} and ${j + 1} on the other variables.`,
    };
    infoEl.innerHTML =
      `<span title="${tt.rho}">ρ<sub>${sub}</sub> = <strong>${value.toFixed(4)}</strong></span>`
      + ` &nbsp;·&nbsp; <span title="${tt.partial}">ρ<sub>${subR}</sub> = <strong>${partialStr}</strong></span>`
      + ` &nbsp;·&nbsp; <span title="${tt.alpha}">R<sup>2</sup><sub>${subI}</sub> = <strong>${info.alpha.toFixed(4)}</strong></span>`
      + ` &nbsp;·&nbsp; <span title="${tt.gamma}">R<sup>2</sup><sub>${subJ}</sub> = <strong>${info.gamma.toFixed(4)}</strong></span>`
      + ` &nbsp;·&nbsp; <span title="${tt.c}">c = <strong>${info.beta.toFixed(4)}</strong></span>`;
    document.getElementById('mode-toggle').hidden = false;
    syncModeToggle();
  }

  function syncModeToggle() {
    document.querySelectorAll('#mode-toggle .mode-btn').forEach(btn => {
      btn.classList.toggle('active', btn.dataset.mode === viewMode);
    });
  }

  // ── Number line SVG ───────────────────────────────────────────────────────
  // In raw mode the axis is c ∈ [-1, 1] with the feasible band [lo, hi]
  // highlighted; in partial mode the axis is ρ ∈ [-1, 1] (everywhere
  // feasible — band hidden). The marker plots whichever value the axis is in.
  function makeNumberLineSVG(value, lo, hi, fontSize) {
    const W = NL_W, H = NL_H, padL = NL_PAD, padR = NL_PAD;
    const fs = fontSize.toFixed(1);
    const innerW = W - padL - padR;
    const f = n => n.toFixed(2);

    const axisToScreen = a => padL + (a + 1) / 2 * innerW;
    const valX  = axisToScreen(cToAxis(value));
    const zeroX = axisToScreen(0);

    // Band always renders. In raw mode it spans the feasible [lo, hi]; in
    // partial mode the whole ρ axis is feasible, so it spans the full slider.
    const bandLoX = viewMode === 'raw' ? axisToScreen(lo) : axisToScreen(-1);
    const bandHiX = viewMode === 'raw' ? axisToScreen(hi) : axisToScreen(1);
    const bandW = Math.max(1, bandHiX - bandLoX);
    const bandSVG = `<rect x="${f(bandLoX)}" y="14" width="${f(bandW)}" height="10"
        fill="rgba(50,110,200,0.15)" stroke="rgba(50,110,200,0.55)" stroke-width="1"/>`;
    // lo/hi numeric labels are only meaningful in raw mode (in partial mode
    // they would always be -1 and +1, redundant with the axis ticks below).
    let labelSVG = '';
    if (viewMode === 'raw') {
      const loX = axisToScreen(lo), hiX = axisToScreen(hi);
      const narrow = (hi - lo) < 0.1;
      const loAnchor = narrow ? 'end'   : 'middle';
      const hiAnchor = narrow ? 'start' : 'middle';
      labelSVG = `<text x="${f(loX)}" y="12" text-anchor="${loAnchor}" font-size="${fs}" fill="rgba(50,110,200,0.85)">${lo.toFixed(3)}</text>
  <text x="${f(hiX)}" y="12" text-anchor="${hiAnchor}" font-size="${fs}" fill="rgba(50,110,200,0.85)">${hi.toFixed(3)}</text>`;
    }

    return `<svg id="nl-svg-el" viewBox="0 0 ${W} ${H}" width="100%"
       style="display:block;overflow:visible;height:auto;font-family:inherit;touch-action:none"
       data-pad-l="${padL}" data-pad-r="${padR}" data-w="${W}">
  <line x1="${padL}" y1="19" x2="${W-padR}" y2="19" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${padL}"   y1="15" x2="${padL}"   y2="23" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${W-padR}" y1="15" x2="${W-padR}" y2="23" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${f(zeroX)}" y1="16" x2="${f(zeroX)}" y2="22" stroke="#aaa" stroke-width="1"/>
  ${bandSVG}
  <line id="nl-marker-line" x1="${f(valX)}" y1="12" x2="${f(valX)}" y2="26" stroke="#2d2d2d" stroke-width="1.5" style="cursor:grab"/>
  <circle id="nl-marker-dot" cx="${f(valX)}" cy="19" r="3" fill="#2d2d2d" style="cursor:grab"/>
  ${labelSVG}
  <text x="${padL}"     y="${H}" text-anchor="middle" font-size="${fs}" fill="#888">-1</text>
  <text x="${f(zeroX)}" y="${H}" text-anchor="middle" font-size="${fs}" fill="#888">0</text>
  <text x="${W-padR}"   y="${H}" text-anchor="middle" font-size="${fs}" fill="#888">+1</text>
</svg>`;
  }

  // ── Drag the marker along the axis ────────────────────────────────────────
  // Attaches to the slider and the λ_min curve. Both SVGs share viewBox width
  // and padding, so the same startDrag handler works on either; it reads
  // dataset.{w,padL,padR} from whichever SVG fired.
  function attachDragHandlers() {
    document.querySelectorAll('#nl-svg-el, #cv-svg-el').forEach(svg => {
      svg.addEventListener('pointerdown', startDrag);
    });
  }

  function startDrag(e) {
    if (!selectedCell) return;
    e.preventDefault();
    const svgEl = e.currentTarget;
    const rect  = svgEl.getBoundingClientRect();
    const W    = +svgEl.dataset.w;
    const padL = +svgEl.dataset.padL;
    const padR = +svgEl.dataset.padR;
    const scale = rect.width / W;
    const { i, j } = selectedCell;
    dragState = {
      i, j,
      pointerId: e.pointerId,
      axisLeft:  rect.left + padL * scale,
      axisRight: rect.left + (W - padR) * scale,
    };
    svgEl.setPointerCapture(e.pointerId);
    svgEl.addEventListener('pointermove',   onDragMove);
    svgEl.addEventListener('pointerup',     endDrag);
    svgEl.addEventListener('pointercancel', endDrag);
    document.body.classList.add('dragging');
    onDragMove(e);
  }

  function onDragMove(e) {
    if (!dragState) return;
    const { axisLeft, axisRight } = dragState;
    const t = (e.clientX - axisLeft) / (axisRight - axisLeft);
    const axisVal = Math.max(-1, Math.min(1, -1 + 2 * t));
    applyDragValue(axisToC(axisVal));
  }

  function endDrag(e) {
    if (!dragState) return;
    const svgEl = e.currentTarget;
    svgEl.releasePointerCapture(dragState.pointerId);
    svgEl.removeEventListener('pointermove',   onDragMove);
    svgEl.removeEventListener('pointerup',     endDrag);
    svgEl.removeEventListener('pointercancel', endDrag);
    dragState = null;
    document.body.classList.remove('dragging');
  }

  // Set the currently-selected cell's value to v and refresh all dependent UI.
  // Bounds for (i,j) itself are invariant under c, but bounds for other pairs
  // shift — recompute the whole map so a follow-up click sees accurate bounds.
  // ~1 ms at n=20, comfortably under one frame.
  function setSelectedValue(v) {
    if (!selectedCell) return;
    const { i, j } = selectedCell;
    const info = currentInfo[i][j];  // bounds for (i,j) are invariant under c
    currentCorr[i][j] = currentCorr[j][i] = v;
    moveMarker(v);
    moveCurveMarker(v);
    updateCellText(i, j, v);
    updateCellText(j, i, v);
    renderModeToggle(i, j);
    currentInfo = computeInfo(currentCorr);
    currentEigvals = eigvalsh(currentCorr);
    refreshCellColors(currentEigvals.some(x => x < INDEFINITE_THRESHOLD));
    if (!isEdited) {
      isEdited = true;
      updateActionButtons();
    }
    saveStateThrottled();
  }

  function applyDragValue(v) { setSelectedValue(v); }

  function moveMarker(v) {
    const innerW = NL_W - 2 * NL_PAD;
    const axisVal = Math.max(-1, Math.min(1, cToAxis(v)));
    const x = NL_PAD + (axisVal + 1) / 2 * innerW;
    const line = document.getElementById('nl-marker-line');
    const dot  = document.getElementById('nl-marker-dot');
    if (line) { line.setAttribute('x1', x); line.setAttribute('x2', x); }
    if (dot)  { dot.setAttribute('cx', x); }
  }

  function updateCellText(i, j, v) {
    const td = cellRefs && cellRefs[i] && cellRefs[i][j];
    if (!td) return;
    const input = td.firstChild && td.firstChild.tagName === 'INPUT' ? td.firstChild : null;
    if (input) input.value = v.toFixed(3);
    else td.textContent = v.toFixed(3);
    // Background colour comes from refreshCellColors — depends on bounds, not on v.
  }

  // ── λ_min curve (Tab 1 bottom panel) ──────────────────────────────────────
  function drawCurve(fontSize) {
    if (!currentCurve || !selectedCell) return;
    const { cs, eigvals } = currentCurve;
    const n = eigvals[0].length;
    // Project the full spectrum down to λ_min (the last entry of each
    // sorted-descending sample) — the only series this curve plots.
    const lambdas = eigvals.map(ev => ev[n - 1]);
    const W = CV_W, H = CV_H;
    const fs = fontSize.toFixed(1);
    const padL = CV_PAD_L, padR = CV_PAD_R, padT = CV_PAD_T, padB = CV_PAD_B;
    const innerW = W - padL - padR;
    const innerH = H - padT - padB;

    const yMinData = Math.min(0, ...lambdas);
    const yMaxData = Math.max(0, ...lambdas);
    const yPad = (yMaxData - yMinData) * 0.08 || 0.1;
    const yMin = yMinData - yPad, yMax = yMaxData + yPad;

    // toX maps raw c (or ρ) through cToAxis so the curve x-axis matches the
    // slider's axis. axisX places labels by their displayed axis value.
    const toX    = c => padL + (cToAxis(c) + 1) / 2 * innerW;
    const axisX  = a => padL + (a + 1) / 2 * innerW;
    const toY    = y => padT + (yMax - y) / (yMax - yMin) * innerH;
    curveTransform = { toX, toY };

    const points = cs.map((c, k) => `${toX(c).toFixed(1)},${toY(lambdas[k]).toFixed(1)}`).join(' ');
    const zeroY = toY(0);
    const negTop = Math.max(zeroY, padT);
    const negH = (H - padB) - negTop;

    const { i, j } = selectedCell;
    const cVal = currentCorr[i][j];
    const lamVal = interpEigval(cVal, n - 1);
    const markerX = toX(cVal);

    document.getElementById('nl-curve').innerHTML = `<svg id="cv-svg-el" viewBox="0 0 ${W} ${H}" width="100%"
       style="display:block;height:auto;font-family:inherit;touch-action:none"
       data-pad-l="${padL}" data-pad-r="${padR}" data-w="${W}">
  ${negH > 0 ? `<rect x="${padL}" y="${negTop}" width="${innerW}" height="${negH}" fill="rgba(180,40,40,0.07)"/>` : ''}
  <line x1="${padL}" y1="${zeroY}" x2="${W-padR}" y2="${zeroY}" stroke="#bbb" stroke-width="1" stroke-dasharray="3,3"/>
  <line x1="${padL}" y1="${padT}"  x2="${padL}"   y2="${H-padB}" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${padL}" y1="${H-padB}" x2="${W-padR}" y2="${H-padB}" stroke="#2d2d2d" stroke-width="1"/>
  <polyline points="${points}" fill="none" stroke="rgba(50,110,200,0.85)" stroke-width="1.5"/>
  <text x="${padL+2}" y="${padT-4}" font-size="${fs}" fill="#666">λ_min</text>
  <line x1="${padL-3}" y1="${toY(yMax).toFixed(1)}" x2="${padL}" y2="${toY(yMax).toFixed(1)}" stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${padL-3}" y1="${zeroY.toFixed(1)}"     x2="${padL}" y2="${zeroY.toFixed(1)}"     stroke="#2d2d2d" stroke-width="1"/>
  <line x1="${padL-3}" y1="${toY(yMin).toFixed(1)}" x2="${padL}" y2="${toY(yMin).toFixed(1)}" stroke="#2d2d2d" stroke-width="1"/>
  <text x="${padL-5}" y="${(toY(yMax)+5).toFixed(1)}" text-anchor="end" font-size="${fs}" fill="#888">${yMax.toFixed(4)}</text>
  <text x="${padL-5}" y="${(zeroY+3).toFixed(1)}"     text-anchor="end" font-size="${fs}" fill="#888">0</text>
  <text x="${padL-5}" y="${(toY(yMin)+3).toFixed(1)}" text-anchor="end" font-size="${fs}" fill="#888">${yMin.toFixed(4)}</text>
  <text x="${axisX(-1)}" y="${H-padB+14}" text-anchor="middle" font-size="${fs}" fill="#888">-1</text>
  <text x="${axisX(0)}"  y="${H-padB+14}" text-anchor="middle" font-size="${fs}" fill="#888">0</text>
  <text x="${axisX(1)}"  y="${H-padB+14}" text-anchor="middle" font-size="${fs}" fill="#888">+1</text>
  <line   id="cv-marker-line" x1="${markerX}" y1="${padT}" x2="${markerX}" y2="${H-padB}" stroke="#2d2d2d" stroke-width="1" stroke-dasharray="2,2" style="cursor:grab"/>
  <circle id="cv-marker-dot"  cx="${markerX}" cy="${toY(lamVal)}" r="3" fill="#2d2d2d" style="cursor:grab"/>
</svg>`;
  }

  // Linear interp into the sampled curve. `k` selects which eigenvalue series
  // (0 = largest, n-1 = smallest = λ_min).
  function interpEigval(c, k) {
    if (!currentCurve) return 0;
    const { cs, eigvals } = currentCurve;
    if (c <= cs[0]) return eigvals[0][k];
    if (c >= cs[cs.length - 1]) return eigvals[cs.length - 1][k];
    let lo = 0, hi = cs.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (cs[mid] <= c) lo = mid; else hi = mid;
    }
    const t = (c - cs[lo]) / (cs[hi] - cs[lo]);
    return eigvals[lo][k] + t * (eigvals[hi][k] - eigvals[lo][k]);
  }

  function moveCurveMarker(c) {
    if (!curveTransform || !currentCurve) return;
    const { toX, toY } = curveTransform;
    const n = currentCurve.eigvals[0].length;
    const x = toX(c);
    const y = toY(interpEigval(c, n - 1));
    const line = document.getElementById('cv-marker-line');
    const dot  = document.getElementById('cv-marker-dot');
    if (line) { line.setAttribute('x1', x); line.setAttribute('x2', x); }
    if (dot)  { dot.setAttribute('cx', x); dot.setAttribute('cy', y); }
  }

  // ── CSV parsing ───────────────────────────────────────────────────────────
  function parseCSVMatrix(text) {
    const lines = text.split(/\r?\n/).map(l => l.trim()).filter(l => l.length > 0);
    if (lines.length === 0) return null;
    const matrix = lines.map(line => line.split(',').map(s => parseFloat(s.trim())));
    const n = matrix.length;
    for (const row of matrix) {
      if (row.length !== n) return null;
      for (const v of row) if (!isFinite(v)) return null;
    }
    return matrix;
  }

  function validateCorr(corr) {
    if (!corr) return "Could not parse CSV as a numeric matrix";
    const n = corr.length;
    if (n < 2) return "Matrix must be at least 2×2";
    if (n > 20) return "Matrix must be at most 20×20";
    return null;
  }

  // ── Event wiring (runs once on DOMContentLoaded, via defer) ──────────────
  function generateRandom(n) {
    const seed = (Date.now() ^ Math.floor(Math.random() * 0xFFFFFFFF)) >>> 0;
    renderResults(wishartCorr(n, mulberry32(seed)));
  }

  document.getElementById('gen-form').addEventListener('submit', e => {
    e.preventDefault();
    const n = parseInt(document.getElementById('n-input').value, 10);
    if (!Number.isInteger(n) || n < 2 || n > 20) {
      showError("n must be an integer between 2 and 20");
      return;
    }
    generateRandom(n);
  });

  document.getElementById('upload-btn').addEventListener('click', () => {
    document.getElementById('csv-file').click();
  });

  document.getElementById('csv-file').addEventListener('change', async e => {
    const file = e.target.files[0];
    if (!file) return;
    const text = await file.text();
    e.target.value = '';
    const corr = parseCSVMatrix(text);
    const err = validateCorr(corr);
    if (err) { showError(err); return; }
    renderResults(corr);
  });

  document.getElementById('reset-btn').addEventListener('click', () => {
    if (!pristineCorr) return;
    renderResults(pristineCorr.map(row => row.slice()), { preservePristine: true });
  });

  document.getElementById('download-btn').addEventListener('click', () => {
    if (!currentCorr) return;
    const n = currentCorr.length;
    const csv = currentCorr.map(row => row.map(v => v.toFixed(6)).join(',')).join('\n');
    const blob = new Blob([csv + '\n'], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `correlation-${n}x${n}.csv`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  });

  // Event-delegated cell click — survives table rebuilds because the listener
  // is on the table element, which itself is never replaced. Clicks on the
  // selected cell's <input> are left alone so the browser handles focus and
  // cursor placement natively; otherwise we select the cell and focus its
  // input so the user can start typing immediately.
  document.getElementById('corr-table').addEventListener('click', e => {
    if (e.target.tagName === 'INPUT') return;
    const td = e.target.closest('td.interactive');
    if (!td) return;
    e.preventDefault();
    const i = +td.dataset.i, j = +td.dataset.j;
    onCellClick(i, j);
    const input = cellRefs[i] && cellRefs[i][j] && cellRefs[i][j].querySelector('input');
    if (input) input.focus();
  });

  // Flush any pending throttled save on unload so a drag that ended < 200 ms
  // before navigating to /theory still persists.
  window.addEventListener('pagehide', () => {
    if (saveTimer) { clearTimeout(saveTimer); saveTimer = null; saveState(); }
  });

  document.getElementById('mode-toggle').addEventListener('click', e => {
    const btn = e.target.closest('.mode-btn');
    if (!btn || !selectedCell) return;
    const mode = btn.dataset.mode;
    if (mode === viewMode) return;
    viewMode = mode;
    syncModeToggle();
    const { i, j } = selectedCell;
    showNumberLine(i, j);
    saveState();
  });

  // Boot: restore the last session from localStorage if present; otherwise
  // generate a fresh random matrix. Default selection is (1, 0) so the
  // slider/curve panel is populated rather than showing the empty hint.
  const saved = loadState();
  if (saved) {
    pristineCorr = saved.pristine && Array.isArray(saved.pristine)
      ? saved.pristine.map(row => row.slice())
      : saved.matrix.map(row => row.slice());
    renderResults(saved.matrix, { preservePristine: true });
    // Restore the n input + edited state so the toolbar mirrors the session.
    document.getElementById('n-input').value = saved.matrix.length;
    isEdited = JSON.stringify(saved.matrix) !== JSON.stringify(pristineCorr);
    updateActionButtons();
    if (saved.viewMode === 'partial') { viewMode = 'partial'; syncModeToggle(); }
    const sel = saved.selectedCell;
    const n = saved.matrix.length;
    const validSel = sel && Number.isInteger(sel.i) && Number.isInteger(sel.j)
                     && sel.i > sel.j && sel.i < n && sel.j >= 0;
    onCellClick(validSel ? sel.i : 1, validSel ? sel.j : 0);
  } else {
    generateRandom(parseInt(document.getElementById('n-input').value, 10));
    onCellClick(1, 0);
  }
})();
