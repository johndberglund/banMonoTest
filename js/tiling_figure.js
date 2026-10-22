// A type's tiling as a stand-alone SVG figure: a window onto the tiling, with the two tile classes shaded
// and every edge drawn as its orbit's curve, so you can see what the tiling is like at a glance.
//
// John, 2026-10-21: our own drawings in place of the paper's figures -- "something like a crop from the
// netshape app complete with curves. I'd want to be able to see enough to tell what the tiling is like."
//
// The curve construction is the builder's (classifyEdgeType / buildEdgeCurves), moved here so the
// builder and the figure script draw the same curves. Works in the browser and in node.
//
//   TilingFigure.svg(typeRow, opts) -> svg string | null
//     typeRow: { sym, ... } from DSYM_TYPES; opts: { id, label, width, height, tilesAcross, colors }
(function (root) {
  'use strict';
  const N = (typeof module !== 'undefined' && module.exports);
  const DD = N ? require('./draw_dsym2.js') : root.DrawDsym;
  const TC = N ? require('./tiler_curves.js') : root.TilerCurves;

  const ptKey = p => Math.round(p[0] * 1e6) + ',' + Math.round(p[1] * 1e6);

  // Which of J/U/S/I an edge orbit gets, from the developed symmetry group (see the builder's comment).
  function classifyEdgeType(dev, P0, P1) {
    const eq = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]) < 1e-6;
    let swapRefl = false, swapRot = false, fixRefl = false;
    for (const T of dev.copies) {
      const det = T[0] * T[4] - T[1] * T[3];
      const Q0 = DD.apply(T, P0), Q1 = DD.apply(T, P1);
      if (eq(Q0, P1) && eq(Q1, P0)) { if (det < 0) swapRefl = true; else swapRot = true; continue; }
      if (det < 0 && eq(Q0, P0) && eq(Q1, P1)) fixRefl = true;
    }
    return fixRefl ? 'I' : swapRefl ? 'U' : swapRot ? 'S' : 'J';
  }

  // One curve per edge orbit, each physical edge once, mirrored per instance by the group (the builder's buildEdgeCurves).
  // grow: optional (orbit length in world units) -> factor >= 1 to enlarge that orbit's motif about the edge's middle
  // choice: optional Map orbit -> { v, k, rev }: which motif variant of the bank (v >= half the count is the flipped side), a
  // height factor k (1, or smaller to keep a corner open), and rev: mirror the curve across the edge's perpendicular bisector
  // (x -> 1-x), which moves a J motif's feature to the other end of its edge (John, 2026-10-22, on #8896: "Reflect that curve
  // over the perpendicular bisector of the edge"); orbits not in it get their default variant at full height
  function buildEdgeCurves(d, grow, choice) {
    if (!TC || !d) return null;
    const orbitRep = new Map();
    for (const t of d.complete) if (t.edgeOrbits) {
      const m = t.pts.length;
      for (let i = 0; i < m; i++) {
        const o = t.edgeOrbits[i];
        if (o != null && !orbitRep.has(o)) orbitRep.set(o, [t.pts[i], t.pts[(i + 1) % m]]);
      }
    }
    if (!orbitRep.size) return null;
    const orbitType = new Map();
    for (const [o, [P0, P1]] of orbitRep) orbitType.set(o, classifyEdgeType(d, P0, P1));
    const orbitIndex = new Map();
    [...orbitRep.keys()].filter(o => orbitType.get(o) !== 'I').sort((a, b) => a - b)
      .forEach((o, i) => orbitIndex.set(o, i));
    const canonCache = new Map(), heights = new Map(), editor = new Map();   // orbit -> motif height (world), editor points
    const canonOf = o => { if (!canonCache.has(o)) {
        const ty = orbitType.get(o), ch = choice && choice.get(o);
        // ch.pts: the curve as EDITOR points (a J's whole curve, a U's or S's half up to the fold) -- a hand-edited curve
        let c = ch && ch.pts ? TC.buildCanon(ty, ch.pts) : TC.testCurve(ty, ch && ch.v != null ? ch.v : (orbitIndex.get(o) || 0));
        if (ch && ch.rev) c = c.map(p => [1 - p[0], p[1]]).reverse();
        // slide a J motif along its edge (John, 2026-10-22: "We can also try sliding the motif left or right along the edge");
        // U and S keep their symmetry about the midpoint, so they never slide
        if (ch && ch.dx) c = c.map((p, i) => (i === 0 || i === c.length - 1) ? p : [Math.min(0.96, Math.max(0.04, p[0] + ch.dx)), p[1]]);
        const [R0, R1] = orbitRep.get(o), g = (grow ? grow(Math.hypot(R1[0] - R0[0], R1[1] - R0[1])) : 1) * (ch && ch.k != null ? ch.k : 1);
        if (g !== 1) c = enlarge(c, g);
        heights.set(o, Math.max(...c.map(p => Math.abs(p[1]))) * Math.hypot(R1[0] - R0[0], R1[1] - R0[1]));
        // the finished curve as EDITOR points (a J's whole curve; a U's or S's half up to the fold -- every transform above keeps
        // their symmetry, so the first half is the editable half), so the builder can open a type looking like its figure
        editor.set(o, ty === 'J' || ty === 'I' ? c.map(p => p.slice()) : c.slice(0, (c.length + 1) / 2).map(p => p.slice()));
        canonCache.set(o, c);
      }
      return canonCache.get(o); };
    const kp = p => Math.round(p[0] / 1e-6) + ',' + Math.round(p[1] / 1e-6);
    const ek = (a, b) => { const A = kp(a), B = kp(b); return A < B ? A + '|' + B : B + '|' + A; };
    const codeOf = new Map();
    for (const [o, [R0, R1]] of orbitRep) {
      if (orbitType.get(o) === 'I') continue;
      const seenAt = new Map();
      for (const T of d.copies) {
        const Q0 = DD.apply(T, R0), Q1 = DD.apply(T, R1), key = ek(Q0, Q1);
        if (seenAt.has(key)) continue;
        const det = T[0] * T[4] - T[1] * T[3];
        seenAt.set(key, { P0: Q0, P1: Q1, code: det < 0 ? 2 : 0 });
      }
      codeOf.set(o, seenAt);
    }
    const seen = new Set(), edges = [];
    for (const t of d.complete) {
      if (!t.edgeOrbits) continue;
      const m = t.pts.length;
      for (let i = 0; i < m; i++) {
        const P0 = t.pts[i], P1 = t.pts[(i + 1) % m], orb = t.edgeOrbits[i];
        if (orb == null || orbitType.get(orb) === 'I') continue;
        const key = ek(P0, P1);
        if (seen.has(key)) continue;
        seen.add(key);
        const info = codeOf.get(orb).get(key) || { P0, P1, code: 0 };
        const mapped = TC.applyMap(canonOf(orb), info.code);
        // P0, P1, code: where this instance sits and which Klein image of the orbit's curve it is -- what an editor needs to
        // turn a drag on this edge back into an edit of the curve itself
        edges.push({ orbit: orb, pts: TC.mapToEdge(mapped, info.P0, info.P1, false), P0: info.P0, P1: info.P1, code: info.code });
      }
    }
    for (const o of orbitIndex.keys()) canonOf(o);      // every orbit's height, even one with no edge in view
    return { edges, orbitIndex, orbitType, heights, editor };
  }

  // How close a tile's edge curves come to one another, in world units: the least distance between points of two DIFFERENT
  // edges of one tile, ignoring points near the corner two edges share (they meet there by design). One tile per class is
  // enough, since every edge belongs to a tile of each class it borders. 0 when two curves cross.
  function clearance(dev, cmap) {
    const reps = new Map();
    for (const t of dev.complete) if (!reps.has(t.orbit)) reps.set(t.orbit, t);
    let worst = Infinity, crossings = 0;
    for (const t of reps.values()) {
      const m = t.pts.length, paths = [];
      for (let i = 0; i < m; i++) {
        const a = t.pts[i], b = t.pts[(i + 1) % m];
        const pth = cmap.get(ptKey(a) + '|' + ptKey(b)) || [a, b], dense = [];
        for (let k = 0; k + 1 < pth.length; k++) for (let q = 0; q < 4; q++) {
          const u = q / 4; dense.push([pth[k][0] + u * (pth[k + 1][0] - pth[k][0]), pth[k][1] + u * (pth[k + 1][1] - pth[k][1])]); }
        dense.push(pth[pth.length - 1]);
        paths.push({ a, b, len: Math.hypot(b[0] - a[0], b[1] - a[1]), pts: dense, raw: pth });
      }
      // a true crossing anywhere, corners included (the distance test below skips the corners, which is where #137's curves
      // cross): segments of two different edges that intersect other than at the corner they share
      const cross = (p, q, r, t) => { const d = (q[0] - p[0]) * (t[1] - r[1]) - (q[1] - p[1]) * (t[0] - r[0]); if (Math.abs(d) < 1e-15) return false;
        const u = ((r[0] - p[0]) * (t[1] - r[1]) - (r[1] - p[1]) * (t[0] - r[0])) / d, v = ((r[0] - p[0]) * (q[1] - p[1]) - (r[1] - p[1]) * (q[0] - p[0])) / d;
        return u > 1e-9 && u < 1 - 1e-9 && v > 1e-9 && v < 1 - 1e-9; };
      // graded, not just 0: each crossing counts, so a search can remove them one at a time (#8996 has several)
      let nx = 0;
      for (let i = 0; i < m; i++) for (let j = i + 1; j < m; j++) {
        const P = paths[i].raw, Q = paths[j].raw;
        for (let a = 0; a + 1 < P.length; a++) for (let b = 0; b + 1 < Q.length; b++)
          if (cross(P[a], P[a + 1], Q[b], Q[b + 1])) nx++;
      }
      if (nx) crossings += nx;
      for (let i = 0; i < m; i++) for (let j = i + 1; j < m; j++) {
        const P = paths[i], Q = paths[j], r = 0.2 * Math.min(P.len, Q.len);
        const corners = [P.a, P.b].filter(c => Math.hypot(c[0] - Q.a[0], c[1] - Q.a[1]) < 1e-9 || Math.hypot(c[0] - Q.b[0], c[1] - Q.b[1]) < 1e-9);
        const nearCorner = p => corners.some(c => Math.hypot(p[0] - c[0], p[1] - c[1]) < r);
        for (const p of P.pts) { if (nearCorner(p)) continue;
          for (const q of Q.pts) { if (nearCorner(q)) continue;
            const d = Math.hypot(p[0] - q[0], p[1] - q[1]); if (d < worst) worst = d; } }
      }
    }
    return crossings ? -crossings : worst;      // negative: that many crossings
  }
  // The most-narrowed VISIBLE corner of any tile, as a fraction of its true angle: at each corner, the least angle between the
  // two boundary paths leaving it, measured from the corner to points along the first half of each edge (curves included),
  // over the angle between the straight edges. A fraction, not degrees: a 45-degree corner that stays 45 is fine, a 90-degree
  // corner drawn as 66 is not. A motif can leave a 90-degree
  // corner looking like a needle (#34's notched squares: John, 2026-10-22, "the angle is so small I can't tell" whether two
  // tiles touch).
  function wedge(dev, cmap) {
    const reps = new Map();
    for (const t of dev.complete) if (!reps.has(t.orbit)) reps.set(t.orbit, t);
    let worst = 1;
    const along = (path, frac) => { const out = []; let L = 0; for (let k = 0; k + 1 < path.length; k++) L += Math.hypot(path[k + 1][0] - path[k][0], path[k + 1][1] - path[k][1]);
      let acc = 0; for (let k = 0; k + 1 < path.length && acc < frac * L; k++) { const a = path[k], b = path[k + 1], d = Math.hypot(b[0] - a[0], b[1] - a[1]);
        for (let q = 1; q <= 4; q++) if (acc + d * q / 4 <= frac * L) out.push([a[0] + (b[0] - a[0]) * q / 4, a[1] + (b[1] - a[1]) * q / 4]); acc += d; }
      return out; };
    for (const t of reps.values()) {
      const m = t.pts.length;
      const pathOf = (a, b) => cmap.get(ptKey(a) + '|' + ptKey(b)) || [a, b];
      // interior is on the left when the outline runs anticlockwise
      let area = 0; for (let i = 0; i < m; i++) { const a = t.pts[i], b = t.pts[(i + 1) % m]; area += a[0] * b[1] - b[0] * a[1]; }
      const ccw = area > 0;
      for (let i = 0; i < m; i++) {
        const c = t.pts[i], nx = t.pts[(i + 1) % m], pv = t.pts[(i - 1 + m) % m];
        const out = along(pathOf(c, nx), 0.5), back = along(pathOf(c, pv), 0.5);
        let geo = (ccw ? Math.atan2(pv[1] - c[1], pv[0] - c[0]) - Math.atan2(nx[1] - c[1], nx[0] - c[0])
                       : Math.atan2(nx[1] - c[1], nx[0] - c[0]) - Math.atan2(pv[1] - c[1], pv[0] - c[0])) * 180 / Math.PI;
        geo = ((geo % 360) + 360) % 360 || 360;
        for (const p of out) for (const q of back) {
          const a1 = Math.atan2(p[1] - c[1], p[0] - c[0]), a2 = Math.atan2(q[1] - c[1], q[0] - c[0]);
          let d = (ccw ? a2 - a1 : a1 - a2) * 180 / Math.PI; d = ((d % 360) + 360) % 360;
          if (d / geo < worst) worst = d / geo;
        }
      }
    }
    return worst;
  }
  // Flip motifs to the other side of their edge, one orbit at a time, keeping each flip that opens up room (John, 2026-10-22,
  // on #1117 and #1347: "the edge curves come too close to each other. Maybe flip them inside out").
  // Each choice is rated by the worse of its curve clearance (in px, against CLEAR_OK) and its narrowest visible corner (against
  // WEDGE_OK degrees). Flipping alone could not open #34's or #77's corners (one class or the other got the needle), so each
  // edge orbit may also take another of the bank's motifs -- one no other orbit of its curve type uses, so matched edges stay
  // unambiguous -- and a lower height (60%). Coordinate search over the orbits, two passes.
  // VISIBILITY counts too: every motif at least MIN_PX tall, or the rating drops -- otherwise, once curves clear and corners
  // are open, shrinking costs nothing and the search shrinks every motif to a speck (#72693, #72712 on 2026-10-22).
  function chooseFlips(dev, grow, s, MIN_PX) {
    const CLEAR_OK = 6, WEDGE_OK = 0.8;
    // rated against 8px, below the 12px a motif is enlarged to: a modest shrink that opens a needle corner (#34) must still win
    const VIS_OK = Math.min(8, MIN_PX);
    // the AVERAGE over orbits, each capped at 1 -- a minimum let one tiny motif block raising any other (#72712 stuck at 1.8px)
    const visible = e => { let t = 0, n = 0; for (const h of e.heights.values()) { t += Math.min(1, h * s / VIS_OK); n++; } return n ? t / n : 1; };
    // crossings rate below any clearance (r = -count), so removing one is always an improvement
    const rate = e => { const cm = curveEdgeMap(e), c = clearance(dev, cm);
      if (c < 0) return { r: c, c: 0 };
      const w = wedge(dev, cm);
      return { r: Math.min(c * s / CLEAR_OK, w / WEDGE_OK, visible(e), 1), c }; };
    let ec = buildEdgeCurves(dev, grow);
    if (!ec) return { ec: null };
    let choice = new Map();
    for (const [o, i] of ec.orbitIndex) choice.set(o, { v: i, k: 1 });
    let best = rate(ec);
    // also try every motif at its lowest height as a start: crossings on several orbits at once can't be undone one by one
    { const low = new Map(); for (const [o, i] of ec.orbitIndex) low.set(o, { v: i, k: 0.15 });
      const el = buildEdgeCurves(dev, grow, low), rl = rate(el);
      if (rl.r > best.r + 0.02) { best = rl; choice = low; ec = el; } }
    const orbs = [...ec.orbitIndex.keys()];
    const half = o => TC.motifCount(ec.orbitType.get(o)) / 2;
    // The simplest motifs to the most crowded edges (John, 2026-10-22, on #72693 and #72712: "We can put the simple motif where
    // it is more crowded and it will fit better"). Crowding is taken as edge length; simplicity as the motif's point count.
    // This is a second starting assignment, kept only if it rates better than the default one.
    {
      const lenOf = new Map(), byType = new Map();
      for (const e of ec.edges) if (!lenOf.has(e.orbit)) { const a = e.pts[0], b = e.pts[e.pts.length - 1]; lenOf.set(e.orbit, Math.hypot(b[0] - a[0], b[1] - a[1])); }
      for (const o of orbs) { const ty = ec.orbitType.get(o); if (!byType.has(ty)) byType.set(ty, []); byType.get(ty).push(o); }
      const assign = new Map(choice);
      for (const [ty, os] of byType) {
        const n = TC.motifCount(ty) / 2;
        const simple = [...Array(n).keys()].sort((a, b) => TC.testCurve(ty, a).length - TC.testCurve(ty, b).length || a - b);
        os.slice().sort((a, b) => (lenOf.get(a) || 0) - (lenOf.get(b) || 0)).forEach((o, i) => assign.set(o, { v: simple[i % n], k: 1 }));
      }
      const ea = buildEdgeCurves(dev, grow, assign), ra = rate(ea);
      if (ra.r > best.r + 1e-9) { best = ra; choice = assign; ec = ea; }
    }
    // swap two orbits' motifs (same curve type), keeping each one's orientation and size: a move the one-orbit search can't make
    const trySwaps = () => {
      for (let i = 0; i < orbs.length; i++) for (let j = i + 1; j < orbs.length; j++) {
        const a = orbs[i], b = orbs[j];
        if (ec.orbitType.get(a) !== ec.orbitType.get(b)) continue;
        const n = half(a), ca = choice.get(a), cb = choice.get(b);
        if (ca.v % n === cb.v % n) continue;
        const c2 = new Map(choice);
        c2.set(a, Object.assign({}, ca, { v: (cb.v % n) + (ca.v >= n ? n : 0) }));
        c2.set(b, Object.assign({}, cb, { v: (ca.v % n) + (cb.v >= n ? n : 0) }));
        const e2 = buildEdgeCurves(dev, grow, c2), r2 = rate(e2);
        if (r2.r > best.r + 0.02) { best = r2; choice = c2; ec = e2; }
      }
    };
    for (let pass = 0; pass < 2; pass++) { trySwaps(); for (const o of orbs) {
      const n = half(o), ty = ec.orbitType.get(o);
      const taken = new Set(orbs.filter(q => q !== o && ec.orbitType.get(q) === ty).map(q => choice.get(q).v % n));
      // Stage 1: which motif, in which of the edge's four Klein-group images (identity, other side, mirrored across the
      // perpendicular bisector, rotated 180 = both), full size, centred. For U and S two of the four are the curve itself.
      // (John, 2026-10-22, on #72677: "rotate the curve ... by 180 degrees about its midpoint. In general, we can try the 4
      // klein group reflection/rotations on the edge.")
      const revs = ty === 'J' ? [false, true] : [false];
      const tryIt = (cand, margin) => {
        const c2 = new Map(choice); c2.set(o, cand);
        const e2 = buildEdgeCurves(dev, grow, c2), r2 = rate(e2);
        if (r2.r > best.r + margin) { best = r2; choice.set(o, cand); ec = e2; return true; }
        return false;
      };
      for (let base = 0; base < n; base++) {
        if (taken.has(base)) continue;
        for (const side of [0, 1]) for (const rev of revs) {
          const cur = choice.get(o), v = base + side * n;
          if (cur.v === v && !!cur.rev === rev && cur.k === 1 && !cur.dx) continue;
          tryIt({ v, k: 1, rev }, 0.02);
        }
      }
      // Stage 2: tune the chosen motif -- lower it (down to 15%: on a tiny tile the enlarged motif can be bigger than the tile,
      // #8996), and slide a J motif along its edge. A lower or off-centre motif has to earn its place clearly.
      const pick = choice.get(o);
      for (const k of [1, 0.6, 0.4, 0.25, 0.15]) for (const dx of (ty === 'J' ? [0, -0.15, 0.15, -0.25, 0.25] : [0])) {
        if (k === pick.k && dx === (pick.dx || 0)) continue;
        tryIt({ v: pick.v, rev: pick.rev, k, dx }, (k < 1 || dx) ? 0.05 : 0.02);
      }
    } }
    if (typeof process !== "undefined" && process.env.FIGDEBUG) console.log("choice", JSON.stringify([...choice]), "rate", JSON.stringify(best), "heights(px)", JSON.stringify([...ec.heights].map(([o, h]) => [o, ec.orbitType.get(o), +(h * s).toFixed(1)])));
    return { ec, choice, clearance: best.c };
  }
  // A short edge's motif is too small to read (John, 2026-10-22: "If the marking is very tiny, we should go to a bigger
  // motif"). Scale it by g about the edge's midpoint (0.5, 0): the end points stay put, inner points move out along the edge
  // (kept inside [0.04, 0.96]) and away from it. Scaling about the midpoint keeps S curves point-symmetric and U curves
  // mirror-symmetric, and every edge of an orbit gets the same g, so matched edges still carry identical marks.
  function enlarge(c, g) {
    // shrinking (g < 1) only lowers the motif; it keeps its spread along the edge
    return c.map((p, i) => (i === 0 || i === c.length - 1) ? p.slice()
      : [Math.min(0.96, Math.max(0.04, 0.5 + (p[0] - 0.5) * Math.max(1, g))), p[1] * g]);
  }
  function curveEdgeMap(ec) {
    const map = new Map();
    if (!ec) return map;
    for (const e of ec.edges) {
      const pts = e.pts, k0 = ptKey(pts[0]), k1 = ptKey(pts[pts.length - 1]);
      if (!map.has(k0 + '|' + k1)) map.set(k0 + '|' + k1, pts);
      if (!map.has(k1 + '|' + k0)) map.set(k1 + '|' + k0, pts.slice().reverse());
    }
    return map;
  }
  function tileBoundaryWorldPts(pts, curveMap) {
    const m = pts.length, out = [];
    for (let i = 0; i < m; i++) {
      const a = pts[i], b = pts[(i + 1) % m];
      const path = curveMap.get(ptKey(a) + '|' + ptKey(b)) || [a, b];
      for (let k = (i > 0 ? 1 : 0); k < path.length; k++) out.push(path[k]);
    }
    return out;
  }
  const polyArea = pts => { let a = 0;
    for (let i = 0; i < pts.length; i++) { const j = (i + 1) % pts.length; a += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]; }
    return Math.abs(a) / 2; };
  const inPoly = (p, pts) => { let c = false;
    for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
      const [xi, yi] = pts[i], [xj, yj] = pts[j];
      if ((yi > p[1]) !== (yj > p[1]) && p[0] < (xj - xi) * (p[1] - yi) / (yj - yi) + xi) c = !c;
    } return c; };

  // Develop enough of the tiling to cover a window of the requested size around the origin.
  function develop(sym, want, realization) {
    let dev = null;
    for (const maxCopies of [160, 400, 1000, 2400, 6000]) {
      try { dev = DD.develop(sym, realization ? { maxCopies, realization } : { maxCopies }); } catch (e) { return null; }
      if (!dev || !dev.complete.length) return null;
      const win = want(dev);
      const covered = (x, y) => dev.complete.some(t => inPoly([x, y], t.pts));
      let ok = true;
      // sample points off any regular grid, so none sits exactly on a tile edge or vertex -- and right out to the
      // window's edges and corners, where a missing tile shows as a white notch
      const fr = [0.0013, 0.0617, 0.2371, 0.4129, 0.5883, 0.7637, 0.9391, 0.9987];
      for (const fx of fr) for (const fy of fr) {
        if (!ok) break;
        ok = covered(win.x0 + (win.x1 - win.x0) * fx, win.y0 + (win.y1 - win.y0) * fy);
      }
      if (ok) return { dev, win };
    }
    return dev ? { dev, win: want(dev), partial: true } : null;
  }

  // Area of the translation lattice's cell: the two shortest independent translations in the developed group.
  function cellArea(dev) {
    const tr = [];
    for (const T of dev.copies) {
      if (Math.abs(T[0] - 1) > 1e-6 || Math.abs(T[1]) > 1e-6 || Math.abs(T[3]) > 1e-6 || Math.abs(T[4] - 1) > 1e-6) continue;
      const L = Math.hypot(T[2], T[5]); if (L > 1e-6) tr.push({ v: [T[2], T[5]], L });
    }
    tr.sort((a, b) => a.L - b.L);
    if (!tr.length) return null;
    const v1 = tr[0].v;
    for (const t of tr) { const c = Math.abs(v1[0] * t.v[1] - v1[1] * t.v[0]); if (c > 1e-6 * Math.max(1, tr[0].L)) return c; }
    return null;
  }

  const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;');

  function figure(row, opts) {
    opts = opts || {};
    const W = opts.width || 300, H = opts.height || 300, LABEL = opts.label ? 26 : 0;
    const scarce = opts.scarce || 7, minCells = opts.minCells || 2;
    const col = Object.assign({ bg: '#ffffff', A: '#c9d6f2', B: '#f2dcb8', line: '#1d2230', text: '#1d2230' }, opts.colors || {});
    // Window: room for about `scarce` tiles of the RARER class (with their share of the other class), and never less than
    // `minCells` repeat cells of the lattice -- so a type with one huge tile and many small ones still shows several of each.
    const want = dev => {
      const r = (H - LABEL) / W, ca = cellArea(dev);
      const count = new Map(), areaOf = new Map();
      for (const t of dev.complete) { count.set(t.orbit, (count.get(t.orbit) || 0) + 1); areaOf.set(t.orbit, (areaOf.get(t.orbit) || 0) + polyArea(t.pts)); }
      const total = [...areaOf.values()].reduce((a, b) => a + b, 0), fewest = Math.min(...count.values());
      const area = Math.max(scarce * total / fewest, ca ? minCells * ca : 0);
      const w = Math.sqrt(area / r), h = w * r;
      // centre on the tile nearest the origin, so the window sits in the developed patch's middle
      const c = dev.complete.reduce((b, t) => Math.hypot(...t.ctr) < Math.hypot(...b.ctr) ? t : b).ctr;
      return { x0: c[0] - w / 2, x1: c[0] + w / 2, y0: c[1] - h / 2, y1: c[1] + h / 2 };
    };
    const got = develop(row.sym, want, opts.realization);
    if (!got) return null;
    const { dev, win, partial } = got;
    const s = W / (win.x1 - win.x0);
    // every motif at least MIN_PX tall (the stock motifs are 0.107 of their edge), up to 3.5x; 12px so it still reads in the
    // 184px-wide review sheet (John, 2026-10-22: "#1106 to 1110 still have very small edge patterns")
    const MIN_PX = opts.motifPx || 12;
    const grow = len => Math.min(3.5, Math.max(1, MIN_PX / (0.107 * len * s)));
    const pick = opts.noFlip ? { ec: buildEdgeCurves(dev, grow) } : chooseFlips(dev, grow, s, MIN_PX);
    const ec = pick.ec, cmap = curveEdgeMap(ec);
    const clearPx = (pick.clearance != null ? pick.clearance : clearance(dev, cmap)) * s;
    const orbits = [...new Set(dev.complete.map(t => t.orbit))].sort((a, b) => a - b);
    const X = x => ((x - win.x0) * s).toFixed(1), Y = y => (LABEL + (win.y1 - y) * s).toFixed(1);
    const pad = 0.5 * (win.x1 - win.x0);
    const near = t => t.pts.some(p => p[0] > win.x0 - pad && p[0] < win.x1 + pad && p[1] > win.y0 - pad && p[1] < win.y1 + pad);
    const paths = { A: [], B: [] };
    for (const t of dev.complete) {
      if (!near(t)) continue;
      const b = tileBoundaryWorldPts(t.pts, cmap);
      paths[orbits.indexOf(t.orbit) === 0 ? 'A' : 'B'].push('M' + b.map(p => X(p[0]) + ' ' + Y(p[1])).join('L') + 'Z');
    }
    const id = 'c' + (opts.id || 'fig');
    let out = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`
      + `<rect width="${W}" height="${H}" fill="${col.bg}"/>`
      + `<clipPath id="${id}"><rect x="0" y="${LABEL}" width="${W}" height="${H - LABEL}"/></clipPath>`
      + `<g clip-path="url(#${id})" stroke="${col.line}" stroke-width="1.3" stroke-linejoin="round">`
      + `<path fill="${col.A}" d="${paths.A.join('')}"/><path fill="${col.B}" d="${paths.B.join('')}"/></g>`;
    if (opts.label) out += `<text x="${W / 2}" y="18" text-anchor="middle" font-family="system-ui,Segoe UI,sans-serif" `
      + `font-size="13" fill="${col.text}">${esc(opts.label)}</text>`;
    const curves = {};
    for (const [o, pts] of ec.editor) if (ec.orbitType.get(o) !== 'I') curves[o] = pts.map(p => [+p[0].toFixed(5), +p[1].toFixed(5)]);
    return { svg: out + `</svg>`, partial: !!partial, clearPx, curves };
  }
  function svg(row, opts) { const f = figure(row, opts); return f ? f.svg : null; }

  const api = { svg, figure, buildEdgeCurves, classifyEdgeType, curveEdgeMap, tileBoundaryWorldPts, clearance, wedge };
  if (N) module.exports = api; else root.TilingFigure = api;
})(typeof window !== 'undefined' ? window : this);
