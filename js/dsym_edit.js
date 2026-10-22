// dsym_edit.js -- the editable model of a tiling's fundamental domain.
//
// WHY THIS EXISTS
// dsym_deform.js's "length cone" was the engine's only model of a type's deformation family, and it
// is far too small: it freezes EVERY corner angle at the equal-share value 2*PI/(v*f*s) -- including
// at regular corners, where nothing forces them -- and it never touches corners interior to the
// fundamental domain. Measured on #8826 (2*22, a 16-gon and a triangle): the cone bottoms out at a
// tile-area ratio of 64.6 and the canonical drawing sits at 133.3, while John's hand-drawn target is
// 9.8. Boundary AND interior together reach 10.5; interior alone only reaches 27.4.
//
// THE STATE IS `P`, AND ONLY `P`
// draw_dsym2.js's develop() destructures {P,E,N} but never reads E (verified: no `E[` in that file);
// N is used only inside side() orientation witnesses. E and N are pure functions of P -- and because
// the three side midpoints sum to the same total as the three corners, buildC's N reduces to the
// chamber centroid.
//
// NEVER CALL relax()
// dsym_disk.js:291-335 is barycentric Gauss-Seidel smoothing whose only pin set is walk corners of
// order > 1, so it averages away every other edit. place_from_points.js:156-157 already records that
// skipping it leaves the canonical realization byte-identical.
//
// VALIDITY IS AN ORIENTATION TEST, NOT verify()
// DrawDsym.verify() is a coverage probe. Measured on #8826: of 3000 random interior perturbations,
// 517 passed verify() and ZERO were actually valid -- every one had flipped chambers, so the derived
// side-pairings were no longer the right group. The sound, O(D.size) gate is: the flag graph of the
// cut-open disk (flags joined by sigma_t across NON-boundary edges) is bipartite, and every chamber
// triangle's signed area must agree with that 2-colouring up to one global sign.
(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports)
    module.exports = factory(require('./dsym_geom.js'), require('./dsym_disk.js'));
  else root.DsymEdit = factory(root.DsymGeom, root.DsymDisk);
})(typeof window !== 'undefined' ? window : this, function (DG, DK) {
  'use strict';
  const PAIR = DG.PAIR;
  const TWO_PI = 2 * Math.PI;

  // ---- combinatorial structure (pure functions of D, I, K; cacheable per type) ----

  // Bipartite 2-colouring of the cut-open disk's flag graph. Adjacent chambers (across a
  // non-boundary edge) are mirror images of one another, so their signed areas differ in sign.
  function chamberSigns(D, I, K) {
    const col = new Int8Array(D.size + 1);
    for (let f0 = 1; f0 <= D.size; f0++) {
      if (col[f0]) continue;
      col[f0] = 1;
      const stack = [f0];
      while (stack.length) {
        const f = stack.pop();
        for (let t = 0; t < 3; t++) {
          if (K.isB[I.nodes[f - 1].ep[t]]) continue;   // boundary: the far side is a DIFFERENT copy
          const g = D.sigma[t][f];
          if (!col[g]) { col[g] = -col[f]; stack.push(g); }
          else if (col[g] !== -col[f]) throw new Error('chamberSigns: flag graph is not bipartite');
        }
      }
    }
    return col;
  }

  // The determinant each side-pairing generator must have, read off the colouring alone. A mirror
  // is orientation-reversing; a cut edge [f,g] reverses exactly when f and g have opposite colours.
  function generatorOrientations(D, I, K) {
    const col = chamberSigns(D, I, K);
    const out = new Map();
    I.edges.forEach((e, ei) => {
      if (K.isMirror[ei]) { out.set(ei, -1); return; }
      if (!K.cut.has(ei) || e.nodes.length !== 2) return;
      const [f, g] = e.nodes;
      out.set(ei, -col[g] / col[f]);
    });
    return out;
  }

  // One entry per orbit on the boundary walk. `targetSum` is Poincare's vertex-cycle condition: the
  // sectors at this orbit's instances glue into one cone of total 2*PI/(v*f). `equalShare` is what
  // dsym_deform froze each instance at -- sufficient, but strictly stronger than necessary, and that
  // gap is precisely the angle freedom the length cone was missing.
  function walkOrbitGroups(D, I, K) {
    const by = new Map();
    K.walk.forEach((w, j) => {
      if (!by.has(w.orbit)) by.set(w.orbit, { orbit: w.orbit, positions: [], corners: [] });
      const g = by.get(w.orbit);
      g.positions.push(j);
      g.corners.push(w.corner);
    });
    const out = [];
    for (const g of by.values()) {
      const o = K.ord[g.orbit];
      g.sCount = g.positions.length;
      g.v = o.v; g.f = o.f;
      g.targetSum = TWO_PI / (o.v * o.f);
      g.equalShare = g.targetSum / g.sCount;
      g.forced = (o.v * o.f) > 1;
      out.push(g);
    }
    return out;
  }

  function interiorCorners(K) {
    const onWalk = new Set(K.walk.map(w => w.corner));
    const out = [];
    for (let i = 0; i < K.nCorners; i++) if (!onWalk.has(i)) out.push(i);
    return out;
  }

  // Walk-index pairs that are the two copies of one cut edge. Their lengths must stay equal, or
  // isometryPair() (draw_dsym2.js:68) returns a SIMILARITY and the group stops being discrete.
  function cutPairs(D, I, K) {
    const by = new Map();
    K.walk.forEach((w, j) => {
      if (!K.cut.has(w.edge)) return;
      if (!by.has(w.edge)) by.set(w.edge, []);
      by.get(w.edge).push(j);
    });
    const out = [];
    for (const [edge, js] of by) for (let k = 1; k < js.length; k++) out.push({ edge, a: js[0], b: js[k] });
    return out;
  }

  // Identifies the corner-instance indexing, so a saved P array can be refused if dsym_disk's greedy
  // chooseCuts() (dsym_disk.js:132-215) ever re-indexes.
  function diskFingerprint(D, I, K) {
    const walk = K.walk.map(w => w.corner);
    const cut = [...K.cut].sort((a, b) => a - b);
    const raw = K.nCorners + '|' + K.nSides + '|' + walk.join(',') + '|' + cut.join(',');
    let h = 2166136261;
    for (let i = 0; i < raw.length; i++) { h ^= raw.charCodeAt(i); h = Math.imul(h, 16777619); }
    return { nCorners: K.nCorners, nSides: K.nSides, walk, cut, hash: (h >>> 0).toString(16) };
  }

  // ---- the model ----

  function buildEditModel(D, I, K, C) {
    const M = {
      D, I, K,
      P: C.P.map(p => (p ? p.slice() : p)),
      corners: C.corners,
      n: K.walk.length,
      walk: K.walk,
      interior: interiorCorners(K),
      orbitGroups: walkOrbitGroups(D, I, K),
      cutPairs: cutPairs(D, I, K),
      chamberSigns: chamberSigns(D, I, K),
      genDet: generatorOrientations(D, I, K),
      fingerprint: diskFingerprint(D, I, K),
    };
    // FLAT CHAMBERS. Some domains have a chamber of zero area (its three corners on one line) -- 61 types, all with glides or only
    // translations (1**, 1*x, 22x, 1xx, o1), e.g. #408, #1110. draw_dsym2 decides which side of a cut edge a chamber lies on from
    // its centre N, and the canonical (relaxed) N sits just off that line, on the right side; the plain average toC() used to
    // compute sits exactly on it, the side test fails, and tiles go missing (found 2026-10-22 from John's #1110 edit, which had
    // holes). So remember, for each chamber, which side of its longest side the canonical N is on, and keep it there.
    M.nSide = [];
    for (let f = 1; f <= D.size; f++) {
      const cs = [0, 1, 2].map(t => C.P[K.cornOf(f, t)]);
      let bi = 0, bl = -1;
      for (let t = 0; t < 3; t++) { const a = cs[t], b = cs[(t + 1) % 3], l = Math.hypot(b[0] - a[0], b[1] - a[1]); if (l > bl) { bl = l; bi = t; } }
      const a = cs[bi], b = cs[(bi + 1) % 3], n = C.N && C.N[f - 1];
      const side = n ? Math.sign((b[0] - a[0]) * (n[1] - a[1]) - (b[1] - a[1]) * (n[0] - a[0])) : 0;
      M.nSide.push({ t: bi, side });
    }
    M.walkOfCorner = new Int32Array(K.nCorners).fill(-1);
    K.walk.forEach((w, j) => { M.walkOfCorner[w.corner] = j; });
    return M;
  }

  function cloneModel(M) {
    const N = Object.assign({}, M);
    N.P = M.P.map(p => (p ? p.slice() : p));
    return N;
  }

  // Interior angle at each walk corner, and the length of walk side j (walk[j] -> walk[j+1]).
  function readBoundary(M) {
    const n = M.n, P = M.P, w = M.walk;
    const ang = new Float64Array(n), len = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      const cur = P[w[j].corner], nxt = P[w[(j + 1) % n].corner], prv = P[w[(j - 1 + n) % n].corner];
      len[j] = Math.hypot(nxt[0] - cur[0], nxt[1] - cur[1]);
      const a1 = Math.atan2(prv[1] - cur[1], prv[0] - cur[0]);
      const a2 = Math.atan2(nxt[1] - cur[1], nxt[0] - cur[0]);
      let d = a1 - a2;
      while (d <= 0) d += TWO_PI;
      while (d > TWO_PI) d -= TWO_PI;
      ang[j] = d;
    }
    return { ang, len };
  }

  // Turtle-walk (ang, len) back into P at the walk corners -- the same recurrence as
  // dsym_deform.placeFromLengths (dsym_deform.js:85-118), but with the ANGLES free rather than
  // frozen at alphOf. anchor:'procrustes' rigidly best-fits the new walk onto the old one, so the
  // drawing deforms in place instead of swinging about walk[0].
  function writeBoundary(M, ang, len, opts) {
    const n = M.n, w = M.walk;
    const anchor = (opts && opts.anchor) || 'procrustes';
    const old = w.map(x => M.P[x.corner].slice());
    const pts = new Array(n);
    let x = 0, y = 0, dir = 0;
    for (let j = 0; j < n; j++) {
      pts[j] = [x, y];
      x += len[j] * Math.cos(dir);
      y += len[j] * Math.sin(dir);
      dir += Math.PI - ang[(j + 1) % n];          // exterior turn at the corner we just reached
    }
    const out = anchor === 'procrustes' ? rigidFit(pts, old) : pts;
    for (let j = 0; j < n; j++) M.P[w[j].corner] = out[j];
    return out;
  }

  // Best rigid (rotation + translation, no scale) placement of `src` onto `dst`.
  function rigidFit(src, dst) {
    const n = src.length;
    let sx = 0, sy = 0, dx = 0, dy = 0;
    for (let i = 0; i < n; i++) { sx += src[i][0]; sy += src[i][1]; dx += dst[i][0]; dy += dst[i][1]; }
    sx /= n; sy /= n; dx /= n; dy /= n;
    let num = 0, den = 0;
    for (let i = 0; i < n; i++) {
      const ax = src[i][0] - sx, ay = src[i][1] - sy, bx = dst[i][0] - dx, by = dst[i][1] - dy;
      num += ax * by - ay * bx;
      den += ax * bx + ay * by;
    }
    const th = Math.atan2(num, den), c = Math.cos(th), s = Math.sin(th);
    return src.map(p => {
      const ax = p[0] - sx, ay = p[1] - sy;
      return [dx + ax * c - ay * s, dy + ax * s + ay * c];
    });
  }

  // ---- realization ----

  // Same formulas as place_from_points.js:159-176 (buildC), with NO relax.
  function toC(M) {
    const { D, K, P } = M;
    const E = new Array(K.nSides), N = new Array(D.size);
    for (let f = 1; f <= D.size; f++) for (let t = 0; t < 3; t++) {
      const [a, b] = PAIR[t];
      const p = P[K.cornOf(f, a)], q = P[K.cornOf(f, b)];
      E[K.sideOf(f, t)] = [(p[0] + q[0]) / 2, (p[1] + q[1]) / 2];
    }
    for (let f = 1; f <= D.size; f++) {
      let x = 0, y = 0;
      for (let t = 0; t < 3; t++) {
        const c = P[K.cornOf(f, t)]; x += c[0]; y += c[1];
        const e = E[K.sideOf(f, t)]; x += e[0]; y += e[1];
      }
      N[f - 1] = [x / 6, y / 6];
      // a flat chamber: nudge N off its line, to the side the canonical N was on (see buildEditModel)
      const ns = M.nSide && M.nSide[f - 1];
      if (ns && ns.side) {
        const cs = [0, 1, 2].map(t => P[K.cornOf(f, t)]);
        const a = cs[ns.t], b = cs[(ns.t + 1) % 3], c = cs[(ns.t + 2) % 3];
        const L = Math.hypot(b[0] - a[0], b[1] - a[1]);
        const area2 = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        if (L > 0 && Math.abs(area2) < 1e-9 * L * L) {
          const nx = -(b[1] - a[1]) / L, ny = (b[0] - a[0]) / L, eps = 1e-4 * L;
          N[f - 1] = [N[f - 1][0] + ns.side * eps * nx, N[f - 1][1] + ns.side * eps * ny];
        }
      }
    }
    return { P: P.map(p => p.slice()), E, N, corners: M.corners };
  }

  function chamberArea(M, f) {
    const { K, P } = M;
    const A = P[K.cornOf(f, 0)], B = P[K.cornOf(f, 1)], C = P[K.cornOf(f, 2)];
    return ((B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])) / 2;
  }

  // THE validity gate: every chamber non-degenerate and agreeing with the 2-colouring under one
  // global sign.
  // A chamber can fail in two DIFFERENT ways, and conflating them made whole types uneditable.
  //   inverted  -- non-zero area of the WRONG sign. Fatal: the derived side-pairings are no longer
  //                the right group.
  //   degenerate-- (near) zero area. The chamber is collapsed, not flipped. 69 of the 1270 types
  //                are drawn this way CANONICALLY (#3218 has one chamber of area 5e-18), so
  //                treating it as fatal meant every drag on those types was rejected and the shape
  //                snapped straight back -- exactly what John reported for #3218. Degeneracy is
  //                reported as a warning; opening it up is precisely what the editor is for.
  function orientationCheck(M, eps) {
    const tol = eps === undefined ? 1e-11 : eps;
    const bad = [], degenerate = [];
    let s0 = 0, minAbs = Infinity;
    for (let f = 1; f <= M.D.size; f++) {
      const ar = chamberArea(M, f);
      if (Math.abs(ar) < minAbs) minAbs = Math.abs(ar);
      if (Math.abs(ar) < tol) { degenerate.push(f); continue; }
      const w = Math.sign(ar) * M.chamberSigns[f];
      if (s0 === 0) s0 = w; else if (w !== s0) bad.push(f);
    }
    return { ok: bad.length === 0, bad, degenerate, minAbsArea: minAbs, globalSign: s0 };
  }

  function segInt(p, q, r, s) {
    const d1 = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
    const d2 = (q[0] - p[0]) * (s[1] - p[1]) - (q[1] - p[1]) * (s[0] - p[0]);
    const d3 = (s[0] - r[0]) * (p[1] - r[1]) - (s[1] - r[1]) * (p[0] - r[0]);
    const d4 = (s[0] - r[0]) * (q[1] - r[1]) - (s[1] - r[1]) * (q[0] - r[0]);
    return ((d1 > 0) !== (d2 > 0)) && ((d3 > 0) !== (d4 > 0));
  }

  // The walk polygon must stay simple: a self-crossing boundary cannot bound an embedded domain.
  function simpleCheck(M) {
    const n = M.n, pts = M.walk.map(w => M.P[w.corner]);
    const crossings = [];
    for (let i = 0; i < n; i++) for (let j = i + 2; j < n; j++) {
      if (i === 0 && j === n - 1) continue;            // adjacent through the wrap
      if (segInt(pts[i], pts[(i + 1) % n], pts[j], pts[(j + 1) % n])) crossings.push([i, j]);
    }
    return { ok: crossings.length === 0, crossings };
  }

  function congruenceCheck(M, tol) {
    const { len } = readBoundary(M);
    let worst = 0;
    for (const cp of M.cutPairs) {
      const m = (len[cp.a] + len[cp.b]) / 2;
      if (m > 0) worst = Math.max(worst, Math.abs(len[cp.a] - len[cp.b]) / m);
    }
    return { ok: worst < (tol === undefined ? 1e-9 : tol), worstRel: worst };
  }

  function angleCheck(M, tol) {
    const { ang } = readBoundary(M);
    let worst = 0;
    const perOrbit = M.orbitGroups.map(g => {
      let sum = 0;
      for (const j of g.positions) sum += ang[j];
      const err = Math.abs(sum - g.targetSum);
      if (err > worst) worst = err;
      return { orbit: g.orbit, sum, target: g.targetSum, err };
    });
    return { ok: worst < (tol === undefined ? 1e-9 : tol), worstAbs: worst, perOrbit };
  }

  function polyArea(p) {
    let A = 0;
    for (let i = 0; i < p.length; i++) { const a = p[i], b = p[(i + 1) % p.length]; A += a[0] * b[1] - b[0] * a[1]; }
    return Math.abs(A) / 2;
  }
  function fatness(pts) {                       // same definition as place_from_points.js:63
    let A = 0, Pm = 0;
    for (let i = 0; i < pts.length; i++) {
      const a = pts[i], b = pts[(i + 1) % pts.length];
      A += a[0] * b[1] - b[0] * a[1];
      Pm += Math.hypot(b[0] - a[0], b[1] - a[1]);
    }
    A = Math.abs(A) / 2;
    return Pm > 0 ? (4 * Math.PI * A) / (Pm * Pm) : 0;
  }
  function metrics(dev) {
    const seen = {}, per = [];
    for (const t of dev.complete) if (seen[t.orbit] === undefined) {
      seen[t.orbit] = 1;
      per.push({ orbit: t.orbit, area: polyArea(t.pts), fat: fatness(t.pts), nVerts: t.pts.length });
    }
    const areas = per.map(p => p.area);
    return {
      areaRatio: areas.length < 2 ? 1 : Math.max(...areas) / Math.min(...areas),
      minFat: per.length ? Math.min(...per.map(p => p.fat)) : 0,
      perOrbit: per,
    };
  }

  // Layered validation. The cheap layers are the trustworthy ones; verify() is only confirmatory --
  // it has false positives, which is why orientation is checked first and independently.
  function validate(M, opts) {
    const o = opts || {};
    const out = {
      orientation: orientationCheck(M),
      simple: simpleCheck(M),
      congruence: congruenceCheck(M),
      angles: angleCheck(M),
      develop: null,
      verify: null,
    };
    out.ok = out.orientation.ok && out.simple.ok && out.congruence.ok && out.angles.ok;
    if (o.cheapOnly || !o.develop || !o.sym) return out;
    let dev = null;
    try {
      dev = o.develop(o.sym, { maxCopies: o.maxCopies || 160, realization: { D: M.D, I: M.I, K: M.K, C: toC(M) } });
    } catch (e) {
      out.develop = { ok: false, reason: e.message };
      out.ok = false;
      return out;
    }
    const orbits = new Set(dev.complete.map(t => t.orbit));
    out.develop = { ok: true, copies: dev.copies.length, complete: dev.complete.length, orbits: orbits.size, dev };
    if (o.verify) {
      let v = false;
      try { v = o.verify(dev, o.nProbe || 1200).ok; } catch (e) { v = false; }
      out.verify = { ok: v };
      if (!v) out.ok = false;
    }
    return out;
  }

  // ---- constrained boundary editing ----------------------------------------------------------
  //
  // Parameters are u = (ang[0..n-1], len[0..n-1]). Positions come from a turtle walk in the local
  // gauge (walk corner 0 at the origin, side 0 along +x):
  //     phi_j = sum_{k=1..j} (PI - ang[k]),   p_{j+1} = p_j + len[j] * (cos phi_j, sin phi_j)
  // The constraints are exactly what makes the derived side-pairings a discrete group:
  //   C1 closure          sum len_j cos phi_j = 0, sum len_j sin phi_j = 0
  //   C2 cut congruence   len_a - len_b = 0   (else isometryPair returns a SIMILARITY)
  //   C3 orbit angle sum  sum_{j in orbit} ang_j = 2*PI/(v*f)   (Poincare vertex cycle)
  // C3 is the freedom dsym_deform threw away: it pinned every instance at the equal share
  // 2*PI/(v*f*s), which is sufficient but strictly stronger than necessary.

  function phases(M, ang) {
    const n = M.n, phi = new Float64Array(n);
    for (let j = 1; j < n; j++) phi[j] = phi[j - 1] + (Math.PI - ang[j]);
    return phi;
  }

  function turtle(M, ang, len) {
    const n = M.n, phi = phases(M, ang), pts = new Array(n);
    let x = 0, y = 0;
    for (let j = 0; j < n; j++) {
      pts[j] = [x, y];
      x += len[j] * Math.cos(phi[j]);
      y += len[j] * Math.sin(phi[j]);
    }
    return { pts, phi };
  }

  // Residual vector g(u) and its Jacobian rows, w.r.t. u = (ang, len).
  function constraintSystem(M, ang, len) {
    const n = M.n, phi = phases(M, ang);
    const rows = [], g = [];
    const zero = () => new Float64Array(2 * n);
    let cx = 0, cy = 0;
    for (let j = 0; j < n; j++) { cx += len[j] * Math.cos(phi[j]); cy += len[j] * Math.sin(phi[j]); }
    const rx = zero(), ry = zero();
    for (let k = 1; k < n; k++) {
      let sx = 0, sy = 0;
      for (let j = k; j < n; j++) { sx += len[j] * Math.sin(phi[j]); sy += len[j] * Math.cos(phi[j]); }
      rx[k] = sx; ry[k] = -sy;
    }
    for (let j = 0; j < n; j++) { rx[n + j] = Math.cos(phi[j]); ry[n + j] = Math.sin(phi[j]); }
    rows.push(rx); g.push(cx);
    rows.push(ry); g.push(cy);
    for (const cp of M.cutPairs) {
      const r = zero();
      r[n + cp.a] = 1; r[n + cp.b] = -1;
      rows.push(r); g.push(len[cp.a] - len[cp.b]);
    }
    for (const grp of M.orbitGroups) {
      const r = zero();
      let sum = 0;
      for (const j of grp.positions) { r[j] = 1; sum += ang[j]; }
      rows.push(r); g.push(sum - grp.targetSum);
    }
    /* C4 -- HOLD THE PINNED POINTS. Rotation centres are what the lattice is built from, so if the
       solve is free to slide them, dragging an ordinary corner silently rescales the whole tiling.
       That is what John hit on #442 (22x, 2026-09-04): "Some make it scale just like the yellow
       point. I don't know that we want that." Measured there: one 0.07 nudge of a free midpoint
       moved all four 2-fold centres, by up to 0.09 -- more than the drag itself.
       The rows say the turtle-frame vector from the first pin to each other pin is unchanged, which
       fixes the pinned points' relative configuration; `placeBoundary` then lands them back on their
       old world positions exactly. Note this is an EDITOR contract, not a theorem: sliding the
       centres still yields a legal tiling, just a differently-proportioned one, and that belongs to
       the dimension handles. Rows are added only when `M.pins` is set (see setPins). */
    if (M.pins && M.pins.length > 1 && M.pinRel) {
      const i0 = M.pins[0], K = M.pins.length;
      const V = [], JX = [], JY = [];
      for (let q = 1; q < K; q++) {
        const i1 = M.pins[q];
        const lo = Math.min(i0, i1), hi = Math.max(i0, i1), sgn = i1 > i0 ? 1 : -1;
        let vx = 0, vy = 0;
        for (let j = lo; j < hi; j++) { vx += len[j] * Math.cos(phi[j]); vy += len[j] * Math.sin(phi[j]); }
        const Jx = zero(), Jy = zero();
        for (let k = 1; k < n; k++) {
          let sx = 0, sy = 0;
          for (let j = Math.max(k, lo); j < hi; j++) { sx += len[j] * Math.sin(phi[j]); sy += len[j] * Math.cos(phi[j]); }
          Jx[k] = sgn * sx; Jy[k] = -sgn * sy;
        }
        for (let j = lo; j < hi; j++) { Jx[n + j] = sgn * Math.cos(phi[j]); Jy[n + j] = sgn * Math.sin(phi[j]); }
        V.push([sgn * vx, sgn * vy]); JX.push(Jx); JY.push(Jy);
      }
      /* Hold the pin configuration up to a UNIFORM SCALE, not rigidly. Fixing it outright also
         froze the overall size of the pinned figure, and on #19 (*632) that is the type's ONLY
         freedom -- its best drawing came from scaling the centre triangle by 0.547, not reshaping
         it, and pinning exactly took that away. So: directions locked (cross product zero) and
         proportions locked (length ratios), leaving one scale free. #4810 is the control -- its old
         drag genuinely RESHAPED the configuration (signature change 0.46), and that stays blocked.
         2K-3 rows rather than 2K-2: the one row we drop is precisely the scale. */
      for (let q = 0; q < V.length; q++) {
        const w = M.pinRel[q + 1];
        if (!w) continue;
        const r = zero();
        for (let k = 0; k < 2 * n; k++) r[k] = JX[q][k] * w[1] - JY[q][k] * w[0];
        rows.push(r); g.push(V[q][0] * w[1] - V[q][1] * w[0]);        // v_q parallel to want_q
      }
      const w1 = M.pinRel[1];
      if (w1) {
        const n1 = w1[0] * w1[0] + w1[1] * w1[1];
        for (let q = 1; q < V.length; q++) {
          const w = M.pinRel[q + 1];
          if (!w) continue;
          const nq = w[0] * w[0] + w[1] * w[1];
          const r = zero();
          for (let k = 0; k < 2 * n; k++) {
            r[k] = 2 * (V[q][0] * JX[q][k] + V[q][1] * JY[q][k]) * n1
                 - 2 * (V[0][0] * JX[0][k] + V[0][1] * JY[0][k]) * nq;
          }
          rows.push(r);
          g.push((V[q][0] * V[q][0] + V[q][1] * V[q][1]) * n1
               - (V[0][0] * V[0][0] + V[0][1] * V[0][1]) * nq);       // |v_q| : |v_1| held
        }
      }
    }
    return { rows, g, phi };
  }

  /* Which walk positions are held fixed by a ROTATION in the developed group. Read off the group,
     not guessed from the drawing -- `constraintAt` returns 'pinned' only for a genuine rotation
     stabiliser (or two crossing mirrors). See reference_fixed_points_rigidity.md. */
  function pinnedWalkPositions(M, dev, constraintOf) {
    const out = [];
    for (let j = 0; j < M.n; j++) {
      const p = M.P[M.walk[j].corner];
      if (!p) continue;
      const c = constraintOf(dev, p);
      if (c && c.kind === 'pinned') out.push(j);
    }
    return out;
  }

  /* Record the pins and the configuration they must keep. The vectors are stored in the TURTLE
     frame, so they are a statement about the SHAPE alone -- independent of where the domain happens
     to sit or which way it points. Call this again whenever the pins are allowed to move (i.e.
     after a dimension drag, which rescales the lattice on purpose); otherwise the solve will keep
     trying to restore a configuration the user deliberately changed. */
  function setPins(M, walkIdxs) {
    M.pins = (walkIdxs || []).slice();
    M.pinRel = null;
    if (M.pins.length < 2) return M.pins;
    const b = readBoundary(M), tw = turtle(M, b.ang, b.len), i0 = M.pins[0];
    M.pinRel = M.pins.map((j, q) => (q === 0 ? [0, 0]
      : [tw.pts[j][0] - tw.pts[i0][0], tw.pts[j][1] - tw.pts[i0][1]]));
    return M.pins;
  }

  // Dense solve A x = b, partial pivoting. Sizes here are tiny (2n + rows, n <= 14).
  function solveDense(A, b) {
    const m = A.length;
    const Mx = A.map((r, i) => { const q = new Float64Array(m + 1); q.set(r); q[m] = b[i]; return q; });
    for (let c = 0; c < m; c++) {
      let piv = c;
      for (let r = c + 1; r < m; r++) if (Math.abs(Mx[r][c]) > Math.abs(Mx[piv][c])) piv = r;
      if (Math.abs(Mx[piv][c]) < 1e-14) continue;
      const t = Mx[c]; Mx[c] = Mx[piv]; Mx[piv] = t;
      const d = Mx[c][c];
      for (let r = 0; r < m; r++) {
        if (r === c) continue;
        const f = Mx[r][c] / d;
        if (!f) continue;
        for (let k = c; k <= m; k++) Mx[r][k] -= f * Mx[c][k];
      }
    }
    const x = new Float64Array(m);
    for (let c = 0; c < m; c++) x[c] = Math.abs(Mx[c][c]) < 1e-14 ? 0 : Mx[c][m] / Mx[c][c];
    return x;
  }

  // Minimum-change Gauss-Newton projection back onto the constraint manifold.
  function projectToConstraints(M, ang, len, opts) {
    const o = opts || {}, n = M.n, tol = o.tol || 1e-12, maxIter = o.maxIter || 12;
    const a = Float64Array.from(ang), L = Float64Array.from(len);
    let res = Infinity, iters = 0;
    for (; iters < maxIter; iters++) {
      const { rows, g } = constraintSystem(M, a, L);
      res = Math.max.apply(null, g.map(Math.abs));
      if (res < tol) break;
      const m = rows.length;
      const JJt = [];
      for (let i = 0; i < m; i++) {
        const r = new Float64Array(m);
        for (let j = 0; j < m; j++) { let s = 0; for (let k = 0; k < 2 * n; k++) s += rows[i][k] * rows[j][k]; r[j] = s; }
        JJt.push(r);
      }
      const lam = solveDense(JJt, g.map(v => -v));
      for (let k = 0; k < 2 * n; k++) {
        let d = 0;
        for (let i = 0; i < m; i++) d += rows[i][k] * lam[i];
        if (k < n) a[k] += d; else L[k - n] += d;
      }
      for (let k = 0; k < n; k++) if (L[k] < 1e-9) L[k] = 1e-9;
    }
    return { ang: a, len: L, ok: res < 1e-8, residual: res, iters };
  }

  // Place a solved boundary into M.P holding `anchorIdx` exactly where it was, and choosing the
  // rotation that moves everything else least. John, 2026-08-29: "I would pin the position as I
  // drag. On release, then I might recenter/resize."
  function placeBoundary(M, ang, len, anchorIdx) {
    const n = M.n, w = M.walk;
    const { pts } = turtle(M, ang, len);
    const old = w.map(x => M.P[x.corner].slice());
    /* With pins, anchor on a PIN and take the rotation from the pins alone. Their relative
       configuration is held exactly by C4, so this lands every pinned point back on its old world
       position; a least-squares fit over all corners would spread the residual and let them creep. */
    const pinned = (M.pins && M.pins.length) ? M.pins : null;
    const ai = pinned ? pinned[0] : (anchorIdx == null ? 0 : anchorIdx);
    const fitOver = (pinned && pinned.length > 1) ? pinned : null;
    const ax = pts[ai][0], ay = pts[ai][1], bx = old[ai][0], by = old[ai][1];
    let num = 0, den = 0;
    for (let ii = 0; ii < (fitOver ? fitOver.length : n); ii++) {
      const i = fitOver ? fitOver[ii] : ii;
      const px = pts[i][0] - ax, py = pts[i][1] - ay;
      const qx = old[i][0] - bx, qy = old[i][1] - by;
      num += px * qy - py * qx;
      den += px * qx + py * qy;
    }
    const th = Math.atan2(num, den), c = Math.cos(th), s = Math.sin(th);
    for (let i = 0; i < n; i++) {
      const px = pts[i][0] - ax, py = pts[i][1] - ay;
      M.P[w[i].corner] = [bx + px * c - py * s, by + px * s + py * c];
    }
  }

  // Interior corners feed no generator, so they carry no constraint -- but they are still bounded
  // by the orientation gate, which is why they must be moved continuously rather than jumped.
  function dragInterior(M, corner, target) {
    if (M.walkOfCorner[corner] >= 0) return { ok: false, reason: 'not an interior corner' };
    M.P[corner] = [target[0], target[1]];
    return { ok: true };
  }

  /* Basis of the constraint null space: the directions in u = (ang, len) that change the shape
     without breaking closure, cut-edge congruence, the per-orbit angle sums, or the pins. */
  function nullSpace(rows, N) {
    const A = rows.map(r => Array.from(r)), piv = [];
    let rk = 0;
    for (let c = 0; c < N && rk < A.length; c++) {
      let b = -1, best = 1e-9;
      for (let r = rk; r < A.length; r++) if (Math.abs(A[r][c]) > best) { best = Math.abs(A[r][c]); b = r; }
      if (b < 0) continue;
      const t = A[rk]; A[rk] = A[b]; A[b] = t;
      const d = A[rk][c];
      for (let k = 0; k < N; k++) A[rk][k] /= d;
      for (let r = 0; r < A.length; r++) {
        if (r === rk) continue;
        const f = A[r][c];
        if (!f) continue;
        for (let k = 0; k < N; k++) A[r][k] -= f * A[rk][k];
      }
      piv.push(c); rk++;
    }
    const out = [];
    for (let c = 0; c < N; c++) {
      if (piv.indexOf(c) >= 0) continue;
      const v = new Float64Array(N); v[c] = 1;
      for (let r = 0; r < piv.length; r++) v[piv[r]] = -A[r][c];
      out.push(v);
    }
    return out;
  }

  // d(position of walk corner j, measured from walk corner i0)/du, in the turtle frame.
  function posJacobian(M, ang, len, phi, i0, j) {
    const n = M.n, N = 2 * n;
    const Jx = new Float64Array(N), Jy = new Float64Array(N);
    const lo = Math.min(i0, j), hi = Math.max(i0, j), sgn = j > i0 ? 1 : -1;
    for (let k = 1; k < n; k++) {
      let sx = 0, sy = 0;
      for (let q = Math.max(k, lo); q < hi; q++) { sx += len[q] * Math.sin(phi[q]); sy += len[q] * Math.cos(phi[q]); }
      Jx[k] = sgn * sx; Jy[k] = -sgn * sy;
    }
    for (let q = lo; q < hi; q++) { Jx[n + q] = sgn * Math.cos(phi[q]); Jy[n + q] = sgn * Math.sin(phi[q]); }
    return { Jx, Jy };
  }

  /* WHERE CAN THIS HANDLE ACTUALLY GO? Contract the constraint null space with the handle's
     position Jacobian: rank 2 means it is free in the plane, rank 1 means it can only travel along
     ONE line, rank 0 means it cannot move at all.

     John asked for this on #442 (2026-09-04): "Some of the points (blue squares) do not seem right
     when moved. Some make it jitter... We could restrict its movement to be along the line it is
     on." He was right, and this says exactly which points and which line -- on #442 two of the four
     midpoints are rank 1 and two are rank 2, which is why only some of them misbehaved. Asking a
     rank-1 handle to chase a 2-D pointer is what produced both the jitter and the handle sitting
     0.36 away from the cursor: the solve was hunting for a target that does not exist. */
  function reachable(M, walkIndex) {
    const b = readBoundary(M), n = M.n, N = 2 * n;
    const sys = constraintSystem(M, b.ang, b.len);
    const ns = nullSpace(sys.rows, N);
    const i0 = (M.pins && M.pins.length) ? M.pins[0] : (walkIndex === 0 ? 1 : 0);
    if (i0 === walkIndex) return { rank: 0, dirs: [], ns, cols: [], i0, base: b, phi: sys.phi };
    const pj = posJacobian(M, b.ang, b.len, sys.phi, i0, walkIndex);
    const cols = ns.map(v => {
      let x = 0, y = 0;
      for (let k = 0; k < N; k++) { x += pj.Jx[k] * v[k]; y += pj.Jy[k] * v[k]; }
      return [x, y];
    });
    const dirs = [];
    for (const c of cols) {
      let w = [c[0], c[1]];
      for (const u of dirs) { const d = w[0] * u[0] + w[1] * u[1]; w = [w[0] - d * u[0], w[1] - d * u[1]]; }
      const L = Math.hypot(w[0], w[1]);
      if (L > 1e-7) dirs.push([w[0] / L, w[1] / L]);
    }
    return { rank: dirs.length, dirs, ns, cols, i0, base: b, phi: sys.phi };
  }

  /* Turtle frame -> world rotation, taken over the pins when there are pins (they are held exactly)
     and over the whole walk otherwise. The same fit `placeBoundary` uses, so the two agree. */
  function frameAngle(M, pts) {
    const n = M.n, w = M.walk;
    const set = (M.pins && M.pins.length > 1) ? M.pins : null;
    const i0 = set ? set[0] : 0;
    const ax = pts[i0][0], ay = pts[i0][1];
    const bx = M.P[w[i0].corner][0], by = M.P[w[i0].corner][1];
    let num = 0, den = 0;
    for (let ii = 0; ii < (set ? set.length : n); ii++) {
      const i = set ? set[ii] : ii;
      const px = pts[i][0] - ax, py = pts[i][1] - ay;
      const q = M.P[w[i].corner];
      const qx = q[0] - bx, qy = q[1] - by;
      num += px * qy - py * qx; den += px * qx + py * qy;
    }
    return Math.atan2(num, den);
  }

  /* Move a boundary handle by the least-norm motion that achieves the REACHABLE part of what the
     pointer asked for -- differential inverse kinematics, rather than handing an arbitrary target
     to a least-squares solve and hoping. The unreachable component is simply dropped, so a rank-1
     handle slides along its line and a frozen one does not move, instead of the whole domain
     lurching to fake a solution. Small steps keep it on the manifold; `projectToConstraints` mops
     up the second-order drift after each one. */
  function dragConstrained(M, walkIndex, target, opts) {
    const o = opts || {};
    const steps = o.steps || 4, maxFrac = o.maxFrac === undefined ? 0.34 : o.maxFrac;
    let moved = 0, rank = 0;
    for (let it = 0; it < steps; it++) {
      const rc = reachable(M, walkIndex);
      rank = rc.rank;
      if (!rc.rank) break;
      const b = rc.base, n = M.n, N = 2 * n;
      const tw = turtle(M, b.ang, b.len);
      const th = frameAngle(M, tw.pts), c = Math.cos(th), sn = Math.sin(th);
      const here = M.P[M.walk[walkIndex].corner];
      const ex = target[0] - here[0], ey = target[1] - here[1];
      const e = [ex * c + ey * sn, -ex * sn + ey * c];   // pointer error, in the turtle frame
      let scale = 0;
      for (let i = 0; i < n; i++) scale += b.len[i];
      scale = scale / n || 1;
      const cap = scale * maxFrac, eL = Math.hypot(e[0], e[1]);
      if (eL < 1e-12) break;
      if (eL > cap) { e[0] *= cap / eL; e[1] *= cap / eL; }
      /* Pick the null-space combination that moves the DRAWING least, not the one that moves the
         parameters least. They are not the same thing: on #442, the parameter-least-norm step moved
         the dragged midpoint 0.028 and swung the other two midpoints 0.275 -- a tenfold overshoot
         that reads as the shape lurching. Minimising position change over every walk corner, with
         the handle's own motion held to what was asked, keeps the rest of the domain still. */
      const d = rc.cols.length;
      const W = 1e6;                        // weight on "actually go where the pointer asked"
      const AtA = [], Atb = new Float64Array(d);
      const other = [];                     // position response of every OTHER walk corner
      for (let q = 0; q < n; q++) {
        if (q === walkIndex || q === rc.i0) continue;
        const pj = posJacobian(M, b.ang, b.len, rc.phi, rc.i0, q);
        other.push(rc.ns.map(v => {
          let x2 = 0, y2 = 0;
          for (let k = 0; k < N; k++) { x2 += pj.Jx[k] * v[k]; y2 += pj.Jy[k] * v[k]; }
          return [x2, y2];
        }));
      }
      for (let i = 0; i < d; i++) {
        const r = new Float64Array(d);
        for (let j = 0; j < d; j++) {
          let v = W * (rc.cols[i][0] * rc.cols[j][0] + rc.cols[i][1] * rc.cols[j][1]);
          for (const oc of other) v += oc[i][0] * oc[j][0] + oc[i][1] * oc[j][1];
          r[j] = v;
        }
        r[i] += 1e-9;                       // ridge, so a rank-deficient direction stays finite
        AtA.push(r);
        Atb[i] = W * (rc.cols[i][0] * e[0] + rc.cols[i][1] * e[1]);
      }
      const x = solveDense(AtA, Atb);
      const a2 = Float64Array.from(b.ang), L2 = Float64Array.from(b.len);
      for (let k = 0; k < N; k++) {
        let dv = 0;
        for (let i = 0; i < d; i++) dv += rc.ns[i][k] * x[i];
        if (k < n) a2[k] += dv; else L2[k - n] += dv;
      }
      for (let k = 0; k < n; k++) if (L2[k] < 1e-9) L2[k] = 1e-9;
      const pr = projectToConstraints(M, a2, L2, {});
      if (!pr.ok) break;
      const before = M.P[M.walk[walkIndex].corner].slice();
      placeBoundary(M, pr.ang, pr.len, rc.i0);
      const after = M.P[M.walk[walkIndex].corner];
      moved += Math.hypot(after[0] - before[0], after[1] - before[1]);
    }
    return { ok: rank > 0, rank, moved };
  }

  // Drag walk corner `walkIndex` toward `target`, staying on the constraint manifold.
  function dragBoundary(M, walkIndex, target, opts) {
    const o = opts || {}, n = M.n;
    const base = readBoundary(M);
    let scale = 0;
    for (let i = 0; i < n; i++) scale += base.len[i];
    scale = scale / n || 1;
    const reg = o.reg === undefined ? 1e-3 : o.reg;
    const anchorIdx = o.anchorIdx == null ? (walkIndex === 0 ? 1 : 0) : o.anchorIdx;
    const wOld = M.walk.map(x => M.P[x.corner]);
    const aOld = wOld[anchorIdx];
    const tgt = [target[0] - aOld[0], target[1] - aOld[1]];
    const a = Float64Array.from(base.ang), L = Float64Array.from(base.len);
    const u0 = new Float64Array(2 * n);
    u0.set(base.ang); u0.set(base.len, n);
    for (let it = 0; it < (o.maxIter || 6); it++) {
      const { pts, phi } = turtle(M, a, L);
      const rel = [pts[walkIndex][0] - pts[anchorIdx][0], pts[walkIndex][1] - pts[anchorIdx][1]];
      let num = 0, den = 0;
      for (let i = 0; i < n; i++) {
        const px = pts[i][0] - pts[anchorIdx][0], py = pts[i][1] - pts[anchorIdx][1];
        const qx = wOld[i][0] - aOld[0], qy = wOld[i][1] - aOld[1];
        num += px * qy - py * qx; den += px * qx + py * qy;
      }
      const th = Math.atan2(num, den), c = Math.cos(th), s = Math.sin(th);
      const tLocal = [tgt[0] * c + tgt[1] * s, -tgt[0] * s + tgt[1] * c];
      const lo = Math.min(anchorIdx, walkIndex), hi = Math.max(anchorIdx, walkIndex);
      const sgn = walkIndex > anchorIdx ? 1 : -1;
      const Jx = new Float64Array(2 * n), Jy = new Float64Array(2 * n);
      for (let k = 1; k < n; k++) {
        let sx = 0, sy = 0;
        for (let j = Math.max(k, lo); j < hi; j++) { sx += L[j] * Math.sin(phi[j]); sy += L[j] * Math.cos(phi[j]); }
        Jx[k] = sgn * sx; Jy[k] = -sgn * sy;
      }
      for (let j = lo; j < hi; j++) { Jx[n + j] = sgn * Math.cos(phi[j]); Jy[n + j] = sgn * Math.sin(phi[j]); }
      const ex = rel[0] - tLocal[0], ey = rel[1] - tLocal[1];
      const { rows, g } = constraintSystem(M, a, L);
      const m = rows.length, N = 2 * n, sz = N + m;
      const A = [];
      for (let i = 0; i < sz; i++) A.push(new Float64Array(sz));
      const bvec = new Float64Array(sz);
      const wgt = k => (k < n ? scale * scale : 1) * reg;
      const u = new Float64Array(N); u.set(a); u.set(L, n);
      for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) A[i][j] = Jx[i] * Jx[j] + Jy[i] * Jy[j];
        A[i][i] += wgt(i);
        bvec[i] = -(Jx[i] * ex + Jy[i] * ey) - wgt(i) * (u[i] - u0[i]);
      }
      for (let r = 0; r < m; r++) for (let k = 0; k < N; k++) { A[N + r][k] = rows[r][k]; A[k][N + r] = rows[r][k]; }
      for (let r = 0; r < m; r++) bvec[N + r] = -g[r];
      const sol = solveDense(A, bvec);
      let step = 0;
      for (let k = 0; k < N; k++) step = Math.max(step, Math.abs(sol[k]));
      for (let k = 0; k < n; k++) a[k] += sol[k];
      for (let k = 0; k < n; k++) L[k] = Math.max(1e-9, L[k] + sol[n + k]);
      if (step < 1e-12) break;
    }
    const proj = projectToConstraints(M, a, L, {});
    const oldWalk = M.walk.map(x => M.P[x.corner].slice());
    placeBoundary(M, proj.ang, proj.len, anchorIdx);
    // Carry the interior along. Dragging a slider is a stretch of the WHOLE domain, so the interior
    // corners must ride the same affine map -- otherwise the boundary moves out from under them and
    // you have to drag every inside vertex first. John, 2026-08-29: "if I drag one of the ones that
    // changes the geometry (the two sliders at 90 degrees) ... it's just an affine transform of the
    // whole thing. Otherwise I have to move all the inside vertices away and then move the slider."
    // When the boundary change really is affine (the slider case) the fit is exact; when it is not,
    // it is still the best global approximation, which keeps the interior sensibly placed.
    if (o.transferInterior !== false && M.interior.length) {
      const newWalk = M.walk.map(x => M.P[x.corner]);
      const T = affineFit(oldWalk, newWalk);
      const res = affineResidual(T, oldWalk, newWalk);
      let scl = 0;
      for (const p of oldWalk) scl = Math.max(scl, Math.hypot(p[0], p[1]));
      for (const c of M.interior) M.P[c] = applyAffine(T, M.P[c]);
      return { ok: proj.ok, residual: proj.residual, iters: proj.iters,
               affineExact: res < 1e-9 * Math.max(scl, 1), affineResidual: res };
    }
    return { ok: proj.ok, residual: proj.residual, iters: proj.iters };
  }

  // Least-squares affine map taking src[] onto dst[] (x and y decouple into two 3x3 solves).
  function affineFit(src, dst) {
    let sxx = 0, sxy = 0, sx = 0, syy = 0, sy = 0, s1 = src.length;
    let tx = 0, ty = 0, txx = 0, txy = 0, tyx = 0, tyy = 0;
    for (let i = 0; i < src.length; i++) {
      const x = src[i][0], y = src[i][1], u = dst[i][0], v = dst[i][1];
      sxx += x * x; sxy += x * y; syy += y * y; sx += x; sy += y;
      txx += x * u; tyx += y * u; tx += u;
      txy += x * v; tyy += y * v; ty += v;
    }
    const A = [Float64Array.from([sxx, sxy, sx]), Float64Array.from([sxy, syy, sy]), Float64Array.from([sx, sy, s1])];
    const c1 = solveDense(A.map(r => Float64Array.from(r)), [txx, tyx, tx]);
    const c2 = solveDense(A.map(r => Float64Array.from(r)), [txy, tyy, ty]);
    return [c1[0], c1[1], c1[2], c2[0], c2[1], c2[2]];   // [a b e; c d f]
  }
  const applyAffine = (T, p) => [T[0] * p[0] + T[1] * p[1] + T[2], T[3] * p[0] + T[4] * p[1] + T[5]];
  function affineResidual(T, src, dst) {
    let w = 0;
    for (let i = 0; i < src.length; i++) {
      const q = applyAffine(T, src[i]);
      w = Math.max(w, Math.hypot(q[0] - dst[i][0], q[1] - dst[i][1]));
    }
    return w;
  }

  // The isometry gluing one cut edge's two boundary copies, read off the CURRENT coordinates. The
  // walk traverses the two copies in opposite senses, so side a's start pairs with side b's end.
  function cutPairTransform(M, cp) {
    const n = M.n, w = M.walk;
    const a0 = M.P[w[cp.a].corner], a1 = M.P[w[(cp.a + 1) % n].corner];
    const b0 = M.P[w[cp.b].corner], b1 = M.P[w[(cp.b + 1) % n].corner];
    return similarityFromPairs(a0, a1, b1, b0);
  }
  // The (orientation-preserving or reversing) map sending p1->q1, p2->q2. Both branches are built
  // and the one that best matches the current partner side is kept.
  function similarityFromPairs(p1, p2, q1, q2) {
    const px = p2[0] - p1[0], py = p2[1] - p1[1];
    const qx = q2[0] - q1[0], qy = q2[1] - q1[1];
    const dp = px * px + py * py;
    if (dp < 1e-18) return [1, 0, 0, 0, 1, 0];
    const c = (qx * px + qy * py) / dp, s = (qy * px - qx * py) / dp;   // rotation+scale
    return [c, -s, q1[0] - (c * p1[0] - s * p1[1]),
            s,  c, q1[1] - (s * p1[0] + c * p1[1])];
  }

  // A "free" boundary corner -- a regular point whose orbit's angle total is a full 2*PI -- can be
  // moved anywhere in the plane, PROVIDED its partner instance moves under the cut-edge pairing so
  // the gluing stays an isometry. John, 2026-08-29: "the other points on the boundary (in blue)
  // should be freely movable without changing the slider endpoints. The other blue point would move
  // in tandem." The pairings are captured BEFORE the move and held fixed, which is what keeps the
  // sliders (and every other corner) exactly where they are.
  function dragFree(M, walkIndex, target) {
    const n = M.n, w = M.walk;
    const pre = M.cutPairs.map(cp => ({ cp, g: cutPairTransform(M, cp) }));
    const moved = w[walkIndex].corner;
    M.P[moved] = [target[0], target[1]];
    for (const { cp, g } of pre) {
      const gi = invertAffine(g);
      if (!gi) continue;
      const aS = w[cp.a].corner, aE = w[(cp.a + 1) % n].corner;
      const bS = w[cp.b].corner, bE = w[(cp.b + 1) % n].corner;
      if (moved === aS) M.P[bE] = applyAffine(g, M.P[aS]);
      else if (moved === aE) M.P[bS] = applyAffine(g, M.P[aE]);
      else if (moved === bE) M.P[aS] = applyAffine(gi, M.P[bE]);
      else if (moved === bS) M.P[aE] = applyAffine(gi, M.P[bS]);
    }
    return { ok: true };
  }
  // Is this boundary corner a SLIDER? -- i.e. joined to a forced-angle corner by a MIRROR edge, so
  // it is confined to that mirror line. For #8826 these are the two pink corners either side of the
  // 90-degree one, and sliding them is the "squash in one direction".
  function sliderInfo(M, walkIndex) {
    const n = M.n, w = M.walk;
    const forcedAt = j => {
      const g = M.orbitGroups.find(x => x.positions.indexOf(j) >= 0);
      return g && g.forced && g.sCount === 1;
    };
    for (const [nb, edgeIdx] of [[(walkIndex - 1 + n) % n, (walkIndex - 1 + n) % n], [(walkIndex + 1) % n, walkIndex]]) {
      if (!forcedAt(nb)) continue;
      if (!M.K.isMirror[w[edgeIdx].edge]) continue;
      return { anchorWalk: nb, anchorCorner: w[nb].corner };
    }
    return null;
  }

  // A pure affine stretch of the whole domain along a mirror line. Kept because it is a legitimate
  // deformation in its own right (an invertible affine map preserves chamber orientation, and this
  // one preserves the mirror lines that force the angles), but it must NOT be the default way to
  // drag a slider -- see dragSlider.
  function affineStretch(M, anchorCorner, dirPoint, k) {
    const A = M.P[anchorCorner], Pi = M.P[dirPoint];
    const dx = Pi[0] - A[0], dy = Pi[1] - A[1], L = Math.hypot(dx, dy);
    if (L < 1e-12) return { ok: false, reason: 'degenerate direction' };
    const u = [dx / L, dy / L], v = [-u[1], u[0]];
    for (let i = 0; i < M.P.length; i++) {
      const p = M.P[i]; if (!p) continue;
      const a = (p[0] - A[0]) * u[0] + (p[1] - A[1]) * u[1];
      const b = (p[0] - A[0]) * v[0] + (p[1] - A[1]) * v[1];
      M.P[i] = [A[0] + a * k * u[0] + b * v[0], A[1] + a * k * u[1] + b * v[1]];
    }
    return { ok: true, k };
  }

  // Drag a slider: confine the target to its mirror line (that part IS forced), then let the
  // CONSTRAINED solve place the rest, and repair the interior if a chamber flipped.
  //
  // It is tempting to implement this as a global affine stretch -- for #8826 that is exactly the
  // motion, and it makes the interior ride along instead of being left behind. But a stretch scales
  // every parallel edge by the same factor, which silently LOCKS lengths the symbol leaves free.
  // John caught this on #4810 (2*22): "the mirror edges seem to be locked to always be the same
  // length. I don't think the DD symbol requires this." He was right -- #4810's constraint system
  // has moduli dimension 7 and ties only the two cut pairs, yet the stretch moved all four mirror
  // sides in lockstep (0.5 -> 0.75 together). The solve lets them differ (0.649/0.476/0.524/0.524)
  // and still verifies, while remaining valid on #8826 too once untangle runs.
  function dragSlider(M, walkIndex, target, opts) {
    const o = opts || {};
    const info = o.info || sliderInfo(M, walkIndex);
    if (!info) return { ok: false, reason: 'not a slider' };
    const A = M.P[info.anchorCorner], Pi = M.P[M.walk[walkIndex].corner];
    const dx = Pi[0] - A[0], dy = Pi[1] - A[1], L = Math.hypot(dx, dy);
    if (L < 1e-12) return { ok: false, reason: 'degenerate slider' };
    const u = [dx / L, dy / L];
    const t = Math.max((o.minK || 0.02) * L, (target[0] - A[0]) * u[0] + (target[1] - A[1]) * u[1]);
    const r = dragBoundary(M, walkIndex, [A[0] + t * u[0], A[1] + t * u[1]],
                           { anchorIdx: info.anchorWalk, transferInterior: o.transferInterior });
    if (!orientationCheck(M).ok) untangle(M, {});
    return r;
  }

  // Repair chamber orientation by moving ONLY interior corners (they are unconstrained). Needed
  // because a boundary move can flip the thinnest chambers before the interior catches up -- on
  // #8826 chambers 9 and 12 have areas ~6e-3, the smallest in the domain, and flip immediately.
  // Laplacian smoothing of the INTERIOR corners with every walk corner held fixed. This is the
  // Tutte-style repair: with the boundary pinned, averaging each interior corner toward its chamber
  // neighbours pulls flipped chambers back. Note this is NOT dsym_disk's relax(), which also moves
  // order-1 WALK corners and would therefore undo the user's boundary edit.
  function smoothInterior(M, iters) {
    if (!M.interior.length) return;
    const nbr = new Map();
    for (const c of M.interior) nbr.set(c, new Set());
    for (let f = 1; f <= M.D.size; f++) {
      const cs = [0, 1, 2].map(t => M.K.cornOf(f, t));
      for (const c of cs) {
        if (!nbr.has(c)) continue;
        for (const d of cs) if (d !== c) nbr.get(c).add(d);
      }
    }
    for (let it = 0; it < (iters || 40); it++) {
      for (const c of M.interior) {
        const ns = nbr.get(c);
        if (!ns || !ns.size) continue;
        let x = 0, y = 0;
        for (const d of ns) { x += M.P[d][0]; y += M.P[d][1]; }
        M.P[c] = [x / ns.size, y / ns.size];
      }
    }
  }

  function untangle(M, opts) {
    const o = opts || {}, maxIter = o.maxIter || 200;
    if (!M.interior.length) return { ok: orientationCheck(M).ok, iters: 0 };
    // Smooth objective: push every chamber's signed area to at least `floor` on its correct side,
    // moving ONLY interior corners. A discrete "count of bad chambers" score stalls immediately
    // (it is flat almost everywhere), which is why the first version of this never converged.
    {
      const idxOf = new Map(); M.interior.forEach((c, i) => idxOf.set(c, i));
      let scl = 0;
      for (const p of M.P) if (p) scl = Math.max(scl, Math.hypot(p[0], p[1]));
      scl = scl || 1;
      const areas0 = [];
      for (let f = 1; f <= M.D.size; f++) areas0.push(Math.abs(chamberArea(M, f)));
      const floor = Math.max(1e-9, 0.25 * Math.min.apply(null, areas0.filter(a => a > 1e-12)));
      const s0 = orientationCheck(M).globalSign || 1;
      const energy = () => {
        let E = 0;
        for (let f = 1; f <= M.D.size; f++) {
          const a = s0 * M.chamberSigns[f] * chamberArea(M, f);
          if (a < floor) { const d = floor - a; E += d * d; }
        }
        return E;
      };
      let step = 0.08 * scl, E = energy();
      for (let it = 0; it < (o.maxIter || 400) && E > 0; it++) {
        const grad = M.interior.map(() => [0, 0]);
        for (let f = 1; f <= M.D.size; f++) {
          const want = s0 * M.chamberSigns[f];
          const a = want * chamberArea(M, f);
          if (a >= floor) continue;
          const cs = [0, 1, 2].map(t => M.K.cornOf(f, t));
          const [A, B, C] = cs.map(c => M.P[c]);
          const d = [[(B[1] - C[1]) / 2, (C[0] - B[0]) / 2],
                     [(C[1] - A[1]) / 2, (A[0] - C[0]) / 2],
                     [(A[1] - B[1]) / 2, (B[0] - A[0]) / 2]];
          const coef = -2 * (floor - a) * want;      // d/dP of (floor - want*area)^2
          cs.forEach((c, i) => {
            const gi = idxOf.get(c); if (gi === undefined) return;
            grad[gi][0] += coef * d[i][0];
            grad[gi][1] += coef * d[i][1];
          });
        }
        let gn = 0;
        for (const g of grad) gn = Math.max(gn, Math.hypot(g[0], g[1]));
        if (gn < 1e-18) break;
        const snap = M.P.map(p => (p ? p.slice() : p));
        M.interior.forEach((c, i) => {
          M.P[c] = [M.P[c][0] - step * grad[i][0] / gn, M.P[c][1] - step * grad[i][1] / gn];
        });
        const E2 = energy();
        if (E2 < E) { E = E2; step *= 1.15; }
        else { M.P = snap; step *= 0.5; if (step < 1e-10 * scl) break; }
      }
      if (orientationCheck(M).ok) return { ok: true, iters: 0, via: 'areaPush' };
    }
    let scale = 0;
    for (const p of M.P) if (p) scale = Math.max(scale, Math.hypot(p[0], p[1]));
    let step = (o.step || 0.05) * (scale || 1);
    const idxOf = new Map(); M.interior.forEach((c, i) => idxOf.set(c, i));
    let iters = 0;
    for (; iters < maxIter; iters++) {
      const oc = orientationCheck(M);
      if (oc.ok) return { ok: true, iters };
      const s0 = oc.globalSign || 1;
      const grad = M.interior.map(() => [0, 0]);
      for (const f of oc.bad) {
        const cs = [0, 1, 2].map(t => M.K.cornOf(f, t));
        const [A, B, C] = cs.map(c => M.P[c]);
        // d(signed area)/d(corner), pushed toward the sign this chamber must have
        const d = [[(B[1] - C[1]) / 2, (C[0] - B[0]) / 2],
                   [(C[1] - A[1]) / 2, (A[0] - C[0]) / 2],
                   [(A[1] - B[1]) / 2, (B[0] - A[0]) / 2]];
        const want = s0 * M.chamberSigns[f];
        cs.forEach((c, i) => {
          const gi = idxOf.get(c);
          if (gi === undefined) return;                  // boundary corners are not ours to move
          grad[gi][0] += want * d[i][0];
          grad[gi][1] += want * d[i][1];
        });
      }
      let gmax = 0;
      for (const g of grad) gmax = Math.max(gmax, Math.hypot(g[0], g[1]));
      if (gmax < 1e-15) { step *= 0.5; if (step < 1e-9 * (scale || 1)) break; continue; }
      const before = oc.bad.length;
      const snap = M.P.map(p => (p ? p.slice() : p));
      M.interior.forEach((c, i) => {
        M.P[c] = [M.P[c][0] + step * grad[i][0] / gmax, M.P[c][1] + step * grad[i][1] / gmax];
      });
      const after = orientationCheck(M).bad.length;
      if (after >= before) { M.P = snap; step *= 0.6; if (step < 1e-9 * (scale || 1)) break; }
    }
    return { ok: orientationCheck(M).ok, iters };
  }

  function invertAffine(T) {
    const det = T[0] * T[4] - T[1] * T[3];
    if (Math.abs(det) < 1e-15) return null;
    const a = T[4] / det, b = -T[1] / det, c = -T[3] / det, d = T[0] / det;
    return [a, b, -(a * T[2] + b * T[5]), c, d, -(c * T[2] + d * T[5])];
  }

  // ---- "dimensions of the tiling" vs "shape of the tiles" -------------------------------------
  //
  // John, 2026-08-30: "I'd like a clearer difference between changing the dimensions of the tiling
  // and moving vertices of tiles ... Right now, sometimes (say for #41) when I move any vertex the
  // dimensions change. I can't always predict what the end result will be."
  //
  // The split is real and it is decided by whether a drag changes the SIDE-PAIRING TRANSFORMS:
  //   * interior corners feed no generator at all               -> lattice untouched
  //   * a free boundary corner dragged with dragFree holds the
  //     pairings fixed by construction                          -> lattice untouched
  //   * a slider / any corner moved by the constrained solve
  //     re-derives the pairings from new coordinates            -> lattice CHANGES
  // So `affects` is 'shape' for the first two and 'dimensions' for the last. The 'dimensions'
  // handles are the skeleton John asked for -- and they really are just a few corners.

  // The translation lattice of the developed tiling: the two shortest independent pure
  // translations among dev.copies. This is what "the dimensions" actually means, so it can be shown
  // as numbers instead of being inferred from the picture.
  function latticeBasis(dev, eps) {
    const tol = eps === undefined ? 1e-6 : eps;
    const tr = [];
    for (const T of dev.copies) {
      if (Math.abs(T[0] - 1) > tol || Math.abs(T[1]) > tol ||
          Math.abs(T[3]) > tol || Math.abs(T[4] - 1) > tol) continue;   // linear part must be I
      const v = [T[2], T[5]];
      const L = Math.hypot(v[0], v[1]);
      if (L > tol) tr.push({ v, L });
    }
    tr.sort((a, b) => a.L - b.L);
    if (!tr.length) return null;
    const v1 = tr[0].v;
    let v2 = null;
    for (const t of tr) {
      const cross = v1[0] * t.v[1] - v1[1] * t.v[0];
      if (Math.abs(cross) > tol * Math.max(1, tr[0].L)) { v2 = t.v; break; }
    }
    if (!v2) return { v1, v2: null, a: Math.hypot(v1[0], v1[1]), b: null, angleDeg: null, area: null };
    const a = Math.hypot(v1[0], v1[1]), b = Math.hypot(v2[0], v2[1]);
    const dot = (v1[0] * v2[0] + v1[1] * v2[1]) / (a * b);
    return { v1, v2, a, b,
             angleDeg: Math.acos(Math.max(-1, Math.min(1, dot))) * 180 / Math.PI,
             area: Math.abs(v1[0] * v2[1] - v1[1] * v2[0]) };
  }

  // ---- the two dimension controls -------------------------------------------------------------
  //
  // John, 2026-08-31, looking at #3218: "The only dimensional change allowed is the shape of this
  // rectangle. we don't need 12 square control points... You could have one control point slide for
  // aspect, and another move anywhere to resize (scaling) and rotate. I can't see any way you'd want
  // more than these two control points."
  //
  // He is right, and the reason the editor was offering twelve is that the BOUNDARY moduli (11 for
  // #3218) is not the LATTICE moduli. Most boundary freedom reshapes tiles while leaving the
  // translation lattice alone. The lattice itself is pinned by the wallpaper group up to
  // similarity + at most one aspect parameter, so two controls is exactly right.
  //
  // Both controls are GLOBAL affine maps of every corner. That is safe here (and is NOT the mistake
  // John caught on #4810): an affine map cannot flip a chamber, and applying one to the whole
  // realization is a genuine dimension change. The #4810 bug was using a stretch for a single
  // corner's drag, which silently tied lengths the symbol left free.

  function applyGlobalAffine(M, T) {
    for (let i = 0; i < M.P.length; i++) if (M.P[i]) M.P[i] = applyAffine(T, M.P[i]);
  }

  // Similarity about `centre`: uniform scale s and rotation th. Never changes the tiling's shape,
  // only how it sits on the page -- always available.
  function setSimilarity(M, centre, s, th) {
    const c = Math.cos(th) * s, sn = Math.sin(th) * s;
    applyGlobalAffine(M, [c, -sn, centre[0] - (c * centre[0] - sn * centre[1]),
                          sn,  c, centre[1] - (sn * centre[0] + c * centre[1])]);
  }

  // Aspect: scale by k along `dir`, leaving the perpendicular alone. This is the one genuine
  // dimensional degree of freedom for a rectangular lattice.
  function setAspect(M, centre, dir, k) {
    const L = Math.hypot(dir[0], dir[1]) || 1;
    const u = [dir[0] / L, dir[1] / L], v = [-u[1], u[0]];
    // R^-1 diag(k,1) R, written directly
    const a = k * u[0] * u[0] + v[0] * v[0], b = k * u[0] * u[1] + v[0] * v[1];
    const d = k * u[1] * u[1] + v[1] * v[1];
    applyGlobalAffine(M, [a, b, centre[0] - (a * centre[0] + b * centre[1]),
                          b, d, centre[1] - (b * centre[0] + d * centre[1])]);
  }

  // How many dimensional degrees of freedom the LATTICE has, read straight off the orbifold symbol.
  // John, 2026-08-31: "Of all the orbifolds, any with a rotation 3 or higher is fixed. No changes.
  // (except rotate and zoom) the others have one dimension of freedom - except that 2222 and o each
  // have two dimensions of freedom."
  // This is the standard lattice classification: an order-3/4/6 rotation forces a hexagonal or
  // square lattice (rigid up to similarity); mirrors/glides alone leave a rectangular or rhombic
  // lattice with one ratio; p2 (2222) and p1 (o) are oblique, so both the ratio and the angle move.
  // Authoritative -- use it instead of probing, which cannot tell "rigid" from "I picked a bad axis".
  function latticeFreedom(orbifold) {
    const orb = String(orbifold || '');
    if (/[3-9]/.test(orb)) return 0;          // any rotation of order >= 3 pins the lattice
    // "o" is written "o1" in this census (p1, oblique) -- both are the two-freedom case.
    if (orb === '2222' || /^o/.test(orb)) return 2;
    return 1;
  }

  // For a 1-freedom lattice, WHICH axis does the single stretch run along? Not always a basis
  // vector: for a rhombic lattice (2*22) it is a DIAGONAL of the rhombus, and stretching along v1
  // breaks the cone angles. Try the natural candidates, keeping one that survives a stretch BOTH
  // ways -- a single small test one way lets near-misses through.
  function aspectAxis(M, lat, eps) {
    if (!lat || !lat.v2) return null;
    const tol = eps === undefined ? 1e-7 : eps;
    const cands = [lat.v1, lat.v2,
                   [lat.v1[0] + lat.v2[0], lat.v1[1] + lat.v2[1]],
                   [lat.v1[0] - lat.v2[0], lat.v1[1] - lat.v2[1]]];
    for (const dir of cands) {
      if (Math.hypot(dir[0], dir[1]) < 1e-12) continue;
      let ok = true;
      for (const k of [1.4, 0.7]) {
        const test = cloneModel(M);
        setAspect(test, [0, 0], dir, k);
        // (paired sides must stay equal too -- a glide's axis is the only diagonal that does; #209389, 2026-10-22)
        if (!angleCheck(test, tol).ok || !orientationCheck(test).ok || !congruenceCheck(test, 1e-7).ok) { ok = false; break; }
      }
      if (ok) return dir;
    }
    return null;
  }
  function aspectIsFree(M, lat, eps) { return !!aspectAxis(M, lat, eps); }

  // Replace lattice vector v1 by w, keeping v2 -- the linear map that does it, applied to every
  // corner. Only meaningful for the OBLIQUE cases: John, 2026-08-31, "2222 and o can be based on a
  // parallelogram", so there both the ratio AND the angle of that parallelogram are free, and any
  // invertible linear map preserves p1/p2. (For the 1-freedom groups the FD is a rectangle or
  // rhombus, and only the stretch along one axis is legal -- use setAspect there.)
  function setLatticeVector(M, centre, v1, v2, w) {
    const det = v1[0] * v2[1] - v1[1] * v2[0];
    if (Math.abs(det) < 1e-12) return { ok: false, reason: 'degenerate lattice' };
    // T maps v1 -> w and v2 -> v2, i.e. T = [w v2] [v1 v2]^-1
    const i00 =  v2[1] / det, i01 = -v2[0] / det;
    const i10 = -v1[1] / det, i11 =  v1[0] / det;
    const a = w[0] * i00 + v2[0] * i10, b = w[0] * i01 + v2[0] * i11;
    const c = w[1] * i00 + v2[1] * i10, d = w[1] * i01 + v2[1] * i11;
    applyGlobalAffine(M, [a, b, centre[0] - (a * centre[0] + b * centre[1]),
                          c, d, centre[1] - (c * centre[0] + d * centre[1])]);
    return { ok: true };
  }

  // What is this point allowed to do, according to the GROUP? A point fixed by a reflection must
  // stay on that mirror; a point fixed by a rotation cannot move at all. John, 2026-09-02, on the
  // un-draggable red point in #19/#20: "It is the midpoint of an edge - so I get why I can't move
  // it. But when I move the internal point, that one could move too. (It would be perpendicular to
  // the FD edge since this is a mirror line)". Exactly right -- it is not pinned, it is confined to
  // a line, and that line is read off the stabiliser rather than guessed from the FD's own edges.
  function constraintAt(dev, p, eps) {
    const tol = eps === undefined ? 1e-6 : eps;
    let line = null;
    for (const T of dev.copies) {
      const q = [T[0] * p[0] + T[1] * p[1] + T[2], T[3] * p[0] + T[4] * p[1] + T[5]];
      if (Math.hypot(q[0] - p[0], q[1] - p[1]) > tol) continue;      // not a stabiliser
      const det = T[0] * T[4] - T[1] * T[3];
      if (Math.abs(det - 1) < 1e-9) {
        // a rotation fixing the point (anything but the identity) pins it outright
        if (Math.abs(T[0] - 1) > 1e-9 || Math.abs(T[1]) > 1e-9) return { kind: 'pinned' };
        continue;                                                     // identity
      }
      // reflection: axis direction is the +1 eigenvector of the linear part
      const th = 0.5 * Math.atan2(T[1], T[0]);
      const d = [Math.cos(th), Math.sin(th)];
      if (!line) line = d;
      else if (Math.abs(line[0] * d[1] - line[1] * d[0]) > 1e-9) return { kind: 'pinned' };  // two mirrors
    }
    return line ? { kind: 'line', dir: line } : { kind: 'free' };
  }

  // Project a proposed target onto whatever the point is allowed to do.
  function applyConstraint(con, from, target) {
    if (!con || con.kind === 'free') return target;
    if (con.kind === 'pinned') return from.slice();
    const d = con.dir;
    const t = (target[0] - from[0]) * d[0] + (target[1] - from[1]) * d[1];
    return [from[0] + t * d[0], from[1] + t * d[1]];
  }

  // The dimension controls, placed ON THE FUNDAMENTAL DOMAIN. John, 2026-08-31: "the control points
  // are outside the FD. It would be possible to have the control points and fixed points coincide
  // with vertices of flags of the FD." They were previously put at origin+v1 / origin+v2, i.e. on
  // the LATTICE cell, which is |group| times bigger than the domain -- so they sat well outside it.
  //
  //   anchor  -- a forced-angle walk corner if the type has one, else walk[0]. Held fixed.
  //   shape   -- freedom 1: sliding it along the stretch axis sets the aspect;
  //              freedom 2: it moves anywhere, setting both ratio and angle (oblique cell).
  //   sizer   -- scale + rotate. Always available, never changes the tiling's shape.
  // A drag of either is a GLOBAL map applied to the drag-start snapshot, never compounded frame to
  // frame (compounding is what made the domain judder and collapse to a line).
  function dimensionControls(M, dev, orbifold) {
    const _dev = dev;
    const freedom = latticeFreedom(orbifold);
    const lat = latticeBasis(dev);
    const axis = freedom === 1 ? aspectAxis(M, lat) : null;
    const n = M.n;
    // Prefer real CORNERS of the domain for the controls. John, 2026-09-02 on #32: "could we have
    // the rotate/resize point be one of the vertices of the FD? When it is in the interior of an
    // edge, then it gets pushed when I slide the movable control point." A ctype-0 corner is a
    // tiling vertex -- an actual corner of the domain -- whereas ctype 1 is an edge midpoint, which
    // sits in the middle of a side and gets carried along by neighbouring edits.
    const con = j => constraintAt(dev, M.P[M.walk[j].corner]);
    const ctype = j => M.walk[j].ctype;

    /* Choose the anchor and sizer TOGETHER, preferring a pair whose PERPENDICULAR is a legal
       stretch axis -- then red->blue and red->yellow meet at a right angle and the three controls
       bound the domain's rectangle. John, 2026-09-03 on #41: "You could make the current blue vertex
       be red (so it's the pivot point) make one of the orange dot vertices the blue vertex and make
       the other orange dot vertex the yellow vertex. It seems to me that these are at 90 degrees,
       while the current placement isn't."  Exactly so: on #41 corner 5 -> corner 2 is (1.707, 0)
       and corner 5 -> corner 4 is (0, -1.707). Choosing anchor and sizer independently could not
       see that, because whether the perpendicular is legal depends on the PAIR. */
    /* Best of all: put ALL THREE controls on real corners of the domain, at right angles. John,
       2026-09-04 on #363: "the yellow point is off the FD. we could put the red point directly above
       the blue point and have the yellow point be where the red point now is. It seems we could
       always have the yellow, red and blue points associated with some point."  #363 admits 56 such
       triples. Attaching the ratio handle to a corner also makes it ride every transform for free,
       because the corner is part of the domain -- no re-derivation needed. Falls back to the
       constructed point when no triple exists (that is what keeps #38 from going collinear). */
    if (freedom >= 1) {
      const stretchLegal = (A0, w) => {
        for (const k of [1.4, 0.7]) {
          const test = cloneModel(M);
          setAspect(test, A0, w, k);
          if (!angleCheck(test, 1e-7).ok || !orientationCheck(test).ok || !congruenceCheck(test, 1e-7).ok) return false;
        }
        return true;
      };
      const tyS = j => (ctype(j) === 0 ? 10 : ctype(j) === 2 ? 8 : 0);
      let bestT = null;
      for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) for (let y = 0; y < n; y++) {
        if (a === b || a === y || b === y) continue;
        const A0 = M.P[M.walk[a].corner], S0 = M.P[M.walk[b].corner], Y0 = M.P[M.walk[y].corner];
        const u = [S0[0] - A0[0], S0[1] - A0[1]], w = [Y0[0] - A0[0], Y0[1] - A0[1]];
        const uL = Math.hypot(u[0], u[1]), wL = Math.hypot(w[0], w[1]);
        if (uL < 1e-9 || wL < 1e-9) continue;
        if (Math.abs((u[0] * w[0] + u[1] * w[1]) / (uL * wL)) > (freedom === 1 ? 1e-7 : 0.02)) continue;   // want a right angle
        // (exact when the stretch is along ONE axis: a nearly-perpendicular direction is not a symmetry axis)
        if (!stretchLegal(A0, [w[0] / wL, w[1] / wL])) continue;
        const score = (con(a).kind === 'pinned' ? 60 : 0) + (con(b).kind === 'pinned' ? 40 : 0)
                    + tyS(a) + tyS(b) + tyS(y) + uL + wL;
        if (!bestT || score > bestT.score) bestT = { a, b, y, score };
      }
      if (bestT) {
        const dcT = {
          freedom, lat, axis,
          anchorWalk: bestT.a, anchorCorner: M.walk[bestT.a].corner,
          shapeWalk: bestT.y, shapeCorner: M.walk[bestT.y].corner,   // a REAL corner this time
          shapeAxis: null, shapePoint: null, shapeLen: 0,
          sizerWalk: bestT.b, sizerCorner: M.walk[bestT.b].corner,
          _dev,
        };
        refreshShapeHandle(M, dcT, true);
        return dcT;
      }
      const perpLegal = (Aj, Sj) => {
        const A0 = M.P[M.walk[Aj].corner], S0 = M.P[M.walk[Sj].corner];
        const u = [S0[0] - A0[0], S0[1] - A0[1]], uL = Math.hypot(u[0], u[1]);
        if (uL < 1e-9) return false;
        const perp = [-u[1] / uL, u[0] / uL];
        for (const k of [1.4, 0.7]) {
          const test = cloneModel(M);
          setAspect(test, A0, perp, k);
          if (!angleCheck(test, 1e-7).ok || !orientationCheck(test).ok || !congruenceCheck(test, 1e-7).ok) return false;
        }
        return true;
      };
      const tyScore = j => (ctype(j) === 0 ? 10 : ctype(j) === 2 ? 8 : 0);
      let best = null;
      for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) {
        if (a === b) continue;
        const A0 = M.P[M.walk[a].corner], S0 = M.P[M.walk[b].corner];
        const len = Math.hypot(S0[0] - A0[0], S0[1] - A0[1]);
        if (len < 1e-9) continue;
        if (!perpLegal(a, b)) continue;
        const score = (con(a).kind === 'pinned' ? 60 : 0) + (con(b).kind === 'pinned' ? 40 : 0)
                    + tyScore(a) + tyScore(b) + len;      // longer spokes give better leverage
        if (!best || score > best.score) best = { a, b, score };
      }
      if (best) {
        const A0 = M.P[M.walk[best.a].corner];
        const dc0 = {
          freedom, lat, axis,
          anchorWalk: best.a, anchorCorner: M.walk[best.a].corner,
          shapeWalk: null, shapeCorner: null,
          shapeAxis: null, shapePoint: null, shapeLen: 0,
          sizerWalk: best.b, sizerCorner: M.walk[best.b].corner,
        };
        dc0._dev = _dev;
        refreshShapeHandle(M, dc0, true);
        return dc0;
      }
    }
    // anchor: something genuinely fixed makes the best origin, and a vertex beats a midpoint
    let anchorW = null, anchorScore = -Infinity;
    for (let j = 0; j < n; j++) {
      const c = con(j);
      const score = (c.kind === 'pinned' ? 100 : c.kind === 'line' ? 40 : 0)
                  + (ctype(j) === 0 ? 10 : ctype(j) === 2 ? 5 : 0);
      if (score > anchorScore) { anchorScore = score; anchorW = j; }
    }
    if (anchorW == null) anchorW = 0;
    const A = M.P[M.walk[anchorW].corner];
    const far = [];
    for (let j = 0; j < n; j++) {
      if (j === anchorW) continue;
      const p = M.P[M.walk[j].corner];
      far.push({ j, d: Math.hypot(p[0] - A[0], p[1] - A[1]), vertex: ctype(j) === 0 });
    }
    if (!far.length) return null;
    /* Ranking for the SIZER (scale + rotate). It must not eat a point the user needs for shaping:
       John, 2026-09-02 on #19, "now #19 has the rotate/resize point right on top of the movable
       vertex. We can't move the movable vertex." A PINNED point costs nothing to use as the sizer,
       because it cannot be dragged as a shape point anyway, and the similarity works from any
       reference. So pinned outranks everything; among equals prefer a real corner (vertex, then
       tile centre) over an edge midpoint -- an edge-interior point gets carried along by
       neighbouring edits, which was the #32 complaint. Movable points are the last resort. */
    for (const f of far) {
      const c = con(f.j), ty = ctype(f.j);
      f.pinned = c.kind === 'pinned';
      f.score = (f.pinned ? 100 : 0) + (ty === 0 ? 12 : ty === 2 ? 8 : 0);
    }
    const rank = (a, b) => (b.score - a.score) || (b.d - a.d);
    far.sort(rank);
    // the shape handle should also lie along the stretch axis where there is one
    let shapeW = far[0].j;
    if (axis) {
      const uL = Math.hypot(axis[0], axis[1]) || 1;
      let best = -Infinity;
      for (const f of far) {
        const p = M.P[M.walk[f.j].corner];
        const d = [p[0] - A[0], p[1] - A[1]], dL = Math.hypot(d[0], d[1]) || 1;
        const align = Math.abs((d[0] * axis[0] + d[1] * axis[1]) / (dL * uL));
        // the shape handle is the aspect control, so it must be a point the user can actually
        // move -- the opposite preference from the sizer
        const score = align + (f.vertex ? 0.35 : 0) - (f.pinned ? 2 : 0);
        if (score > best) { best = score; shapeW = f.j; }
      }
    }
    const sizerW = far[0].j;
    /* The ratio handle is a CONSTRUCTED point on the perpendicular through the anchor, not an FD
       corner. John, 2026-09-03 on #38: "when I moved it, the yellow vertex jumped places. It would
       be nice if the yellow vertex was on the perpendicular through the red vertex. (This would put
       blue-red-yellow as vertices around the surrounding rectangle of the FD.)"
       Picking a corner put all three controls COLLINEAR on #38 -- red (0.5,-0.5), blue (0.5,0.5)
       and yellow (0.5,0.0), with the stretch axis along that same line. Dragging yellow past the
       anchor then flipped the sign of k = t1/t0 and the handle leapt. On the perpendicular, with
       t0 fixed at the domain's own extent, the drag is monotone and blue-red-yellow bound the
       domain's rectangle the way John draws it. */
    const dc = {
      freedom, lat, axis,
      anchorWalk: anchorW, anchorCorner: M.walk[anchorW].corner,
      shapeWalk: null, shapeCorner: null,          // the ratio handle is constructed, not a corner
      shapeAxis: null, shapePoint: null, shapeLen: 0,
      sizerWalk: sizerW, sizerCorner: M.walk[sizerW].corner,
    };
    dc._dev = _dev;
    refreshShapeHandle(M, dc, true);
    return dc;
  }

  /* Recompute ONLY the ratio handle's derived geometry, keeping the anchor/sizer choices as they
     are. The handle is a constructed point rather than an FD corner, so unlike every other control
     it does not move by itself when the domain transforms -- John, 2026-09-03: "in #38, If I move
     the blue vertex, the yellow one stays in place. It should rotate and resize with everything
     else. (In fact, isn't it attached to one of the FD vertices?)" It is not attached, deliberately
     (attaching it to a corner is what made all three collinear on #38), so it has to be re-derived
     after every commit. `testPerp` is the expensive validity probe for the perpendicular direction;
     it only depends on the group, so it is done once and the answer reused. */
  function constraintAtCache(M, dc, i) {
    if (!dc._conCache) dc._conCache = {};
    if (dc._conCache[i] === undefined)
      dc._conCache[i] = dc._dev ? constraintAt(dc._dev, M.P[i]).kind : 'free';
    return dc._conCache[i];
  }

  function refreshShapeHandle(M, dc, testPerp, override) {
    if (dc.freedom < 1) return dc;
    const A = M.P[dc.anchorCorner], S = M.P[dc.sizerCorner];
    const u = [S[0] - A[0], S[1] - A[1]];
    const uL = Math.hypot(u[0], u[1]) || 1;
    const perp = [-u[1] / uL, u[0] / uL];
    /* OBLIQUE types (2222 / o) let the ratio handle move anywhere -- that IS the second degree of
       freedom -- so it must not be snapped back onto a constructed axis. John, 2026-09-04 on
       #1728482 (o1, freedom 2): "When I move the yellow point, it warps the FD correctly, but the
       yellow point doesn't follow my mouse." It warped correctly because the drag maths was right;
       the handle then jumped because this function re-derived it from the perpendicular and the
       domain's extent. Instead, hold it in the anchor->sizer FRAME as coefficients, so it rides
       along under a rotate/scale and stays put under everything else. */
    // Attached to a real corner? Then it already moves with the domain -- just read it off.
    if (dc.shapeCorner != null) {
      const Y = M.P[dc.shapeCorner];
      const d = [Y[0] - A[0], Y[1] - A[1]], dL = Math.hypot(d[0], d[1]) || 1;
      dc.shapeAxis = [d[0] / dL, d[1] / dL];
      dc.shapeLen = dL;
      dc.shapePoint = [Y[0], Y[1]];
      return dc;
    }
    /* Hold the ratio handle in the ANCHOR->SIZER frame, for both freedoms. Re-deriving it from a
       world direction makes it lag the pointer, because an aspect stretch also moves the sizer, so
       the perpendicular it was derived from rotates out from under it (measured: 0.25 on #38, 0.43
       on #41). In the frame it rides rotations and scales for free, and after a ratio drag the
       achieved position is written back as the new coefficients so it sits exactly under the mouse.
       John, 2026-09-04 on #1728482: "it warps the FD correctly, but the yellow point doesn't follow
       my mouse." */
    const v = perp;
    const toRel = q => {
      const d = [q[0] - A[0], q[1] - A[1]];
      return { cu: (d[0] * u[0] + d[1] * u[1]) / (uL * uL),
               cv: (d[0] * v[0] + d[1] * v[1]) / uL };
    };
    if (override) dc.shapeRel = toRel(override);
    if (!dc.shapeRel) {
      let ext = 0;
      for (const p of M.P) if (p) ext = Math.max(ext, (p[0] - A[0]) * v[0] + (p[1] - A[1]) * v[1]);
      dc.shapeRel = { cu: 0, cv: (ext > 1e-9 ? ext : uL) / uL };
    }
    const c = dc.shapeRel;
    dc.shapePoint = [A[0] + c.cu * u[0] + c.cv * uL * v[0],
                     A[1] + c.cu * u[1] + c.cv * uL * v[1]];
    if (dc.freedom === 2) {
      dc.shapeAxis = v;
      dc.shapeLen = Math.hypot(dc.shapePoint[0] - A[0], dc.shapePoint[1] - A[1]) || uL;
      return dc;
    }
    if (testPerp) {
      let ok = true;
      for (const k of [1.4, 0.7]) {
        const test = cloneModel(M);
        setAspect(test, A, perp, k);
        if (!angleCheck(test, 1e-7).ok || !orientationCheck(test).ok) { ok = false; break; }
      }
      dc.perpOK = ok;
      /* Store the fallback axis as an angle RELATIVE to anchor->sizer, not as a world direction.
         A world-fixed axis does not turn when the domain is rotated by the sizer, so the ratio
         handle drifts out of position -- on #41 the red-blue/red-yellow angle slid 45 -> 38.3 deg
         during a single rotate. Held relative, it turns with the domain like every other control. */
      if (!ok && dc.axis) {
        const aL = Math.hypot(dc.axis[0], dc.axis[1]) || 1;
        const ax = [dc.axis[0] / aL, dc.axis[1] / aL];
        const uu = [u[0] / uL, u[1] / uL];
        dc.axisRel = Math.atan2(uu[0] * ax[1] - uu[1] * ax[0], uu[0] * ax[0] + uu[1] * ax[1]);
      }
    }
    let chosen;
    if (dc.perpOK) chosen = perp;
    else if (dc.axisRel != null) {
      const c = Math.cos(dc.axisRel), sn = Math.sin(dc.axisRel);
      const uu = [u[0] / uL, u[1] / uL];
      chosen = [uu[0] * c - uu[1] * sn, uu[0] * sn + uu[1] * c];   // turn with anchor->sizer
    } else chosen = dc.axis || perp;
    const aL = Math.hypot(chosen[0], chosen[1]) || 1;
    const w = [chosen[0] / aL, chosen[1] / aL];
    let ext = 0, scale = 0;
    for (const p of M.P) if (p) {
      ext = Math.max(ext, (p[0] - A[0]) * w[0] + (p[1] - A[1]) * w[1]);
      scale = Math.max(scale, Math.hypot(p[0] - A[0], p[1] - A[1]));
    }
    // the handle's own distance along the axis, so the drag maths and the drawn point agree
    let len = (dc.shapePoint ? (dc.shapePoint[0] - A[0]) * w[0] + (dc.shapePoint[1] - A[1]) * w[1] : 0);
    if (!(len > 1e-9)) len = ext > 1e-9 ? ext : uL;
    /* Do not let the constructed handle land ON an existing corner. On #41 it came out exactly on
       corner 1 (an editable edge midpoint), which buried that point completely. Push it outward
       until it is clear of every corner. */
    const movableCorner = i => {
      const j = M.walkOfCorner[i];
      if (j < 0) return true;                       // interior points are always movable
      return constraintAtCache(M, dc, i) !== 'pinned';
    };
    const clearOf = L => {
      const q = [A[0] + w[0] * L, A[1] + w[1] * L];
      const near = Math.max(0.06 * (scale || 1), 1e-6);
      for (let i = 0; i < M.P.length; i++) {
        const p = M.P[i];
        if (!p || !movableCorner(i)) continue;      // a pinned corner cannot be buried: it is inert
        if (Math.hypot(p[0] - q[0], p[1] - q[1]) < near) return false;
      }
      return true;
    };
    for (let i = 0; i < 8 && !clearOf(len); i++) len *= 1.18;
    dc.shapeAxis = w;
    dc.shapeLen = len;
    dc.shapePoint = [A[0] + w[0] * len, A[1] + w[1] * len];
    return dc;
  }

  // Apply a dimension drag to `P0` (the snapshot taken when the drag began) and write the result
  // into M. `which` is 'shape' or 'sizer'; `target` is where the pointer is now.
  function dragDimension(M, dc, P0, which, target) {
    const A = P0[dc.anchorCorner];
    const restore = () => { for (let i = 0; i < P0.length; i++) if (P0[i]) M.P[i] = P0[i].slice(); };
    restore();
    if (which === 'sizer') {
      const q0 = P0[dc.sizerCorner];
      const d0 = [q0[0] - A[0], q0[1] - A[1]], d1 = [target[0] - A[0], target[1] - A[1]];
      const L0 = Math.hypot(d0[0], d0[1]), L1 = Math.hypot(d1[0], d1[1]);
      if (L0 < 1e-12 || L1 < 1e-12) return { ok: false, reason: 'degenerate' };
      setSimilarity(M, A, L1 / L0, Math.atan2(d1[1], d1[0]) - Math.atan2(d0[1], d0[0]));
      return { ok: true, scale: L1 / L0 };
    }
    if (dc.freedom === 0 || !dc.shapeAxis) return { ok: false, reason: 'lattice is rigid' };
    if (dc.freedom === 1) {
      const w = dc.shapeAxis, t0 = dc.shapeLen;
      const t1 = (target[0] - A[0]) * w[0] + (target[1] - A[1]) * w[1];
      if (!(t0 > 1e-12)) return { ok: false, reason: 'degenerate ratio handle' };
      const k = Math.max(0.05, Math.min(20, t1 / t0));   // t0 is the domain's extent: never tiny
      setAspect(M, A, w, k);
      return { ok: true, k };
    }
    // freedom 2: send the ratio handle wherever the pointer is, holding the sizer direction fixed
    const r0 = [dc.shapePoint[0] - A[0], dc.shapePoint[1] - A[1]];
    const s0 = [P0[dc.sizerCorner][0] - A[0], P0[dc.sizerCorner][1] - A[1]];
    return setLatticeVector(M, A, r0, s0, [target[0] - A[0], target[1] - A[1]]);
  }

  // The two handles, in world coordinates, placed on the lattice skeleton.
  function latticeControls(M, dev, orbifold) {
    const lat = latticeBasis(dev);
    if (!lat || !lat.v2) return null;
    let origin = [0, 0];
    for (const p of M.P) if (p) { origin = p.slice(); break; }
    const freedom = latticeFreedom(orbifold);
    // 0 -> similarity only; 1 -> one stretch axis; 2 -> oblique, the second vector moves freely
    const axis = freedom === 1 ? aspectAxis(M, lat) : null;
    const aLen = axis ? Math.hypot(axis[0], axis[1]) : 0;
    return {
      lat, origin, axis, freedom,
      aspectFree: freedom === 1 && !!axis,
      obliqueFree: freedom === 2,
      // the aspect handle sits along the usable stretch axis; the similarity handle on the other
      aspect: axis ? [origin[0] + axis[0], origin[1] + axis[1]] : null,
      similarity: [origin[0] + lat.v2[0], origin[1] + lat.v2[1]],
      aspectLen: aLen,
    };
  }

  // Handle roles follow John's #8826 analysis: a forced-angle corner, corners that slide to keep
  // that angle, a derived centre, and free interior points.
  function classifyHandles(M) {
    const out = [];
    const byCorner = new Map();
    for (const g of M.orbitGroups) for (const j of g.positions) byCorner.set(M.walk[j].corner, g);
    for (let c = 0; c < M.K.nCorners; c++) {
      const j = M.walkOfCorner[c];
      if (j < 0) {
        out.push({ corner: c, kind: 'interior', walkIndex: -1, role: 'free', affects: 'shape' });
        continue;
      }
      const g = byCorner.get(c), w = M.walk[j];
      const role = g.forced ? (g.sCount > 1 ? 'splitForced' : 'forced') : (g.sCount > 1 ? 'splitFree' : 'flat');
      out.push({
        corner: c, kind: 'boundary', walkIndex: j, orbit: w.orbit, cornerType: w.ctype, role,
        // splitFree corners are dragged by dragFree, which holds the pairings fixed, so they reshape
        // the tiles without touching the lattice. Everything else goes through the constrained solve,
        // which re-derives the pairings -- that is a dimension change.
        affects: role === 'splitFree' ? 'shape' : 'dimensions',
        orbitTargetDeg: g.targetSum * 180 / Math.PI, sCount: g.sCount,
        siblings: g.positions.filter(x => x !== j),
      });
    }
    return out;
  }

  return {
    chamberSigns, generatorOrientations, walkOrbitGroups, interiorCorners, cutPairs, diskFingerprint,
    buildEditModel, cloneModel, readBoundary, writeBoundary, rigidFit,
    toC, chamberArea, orientationCheck, simpleCheck, congruenceCheck, angleCheck,
    validate, metrics, fatness, polyArea, latticeBasis,
    applyGlobalAffine, setSimilarity, setAspect, aspectIsFree, aspectAxis, latticeFreedom, latticeControls, setLatticeVector, dimensionControls, refreshShapeHandle, dragDimension, constraintAt, applyConstraint,
    phases, turtle, constraintSystem, solveDense, projectToConstraints, placeBoundary,
    pinnedWalkPositions, setPins, nullSpace, reachable, dragConstrained, frameAngle,
    dragInterior, dragBoundary, dragFree, dragSlider, affineStretch, sliderInfo, untangle, smoothInterior, classifyHandles,
    affineFit, applyAffine, affineResidual, invertAffine, cutPairTransform, similarityFromPairs,
  };
});
