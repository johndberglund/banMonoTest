// Fatter drawings for the figures: pick, from a type's family of shapes, the one whose THINNEST tile class is as fat as
// possible, so a figure never shows a class as a sliver (John, 2026-10-21: "Is there a way we can make those fatter?").
//
// The family is dsym_deform's: the domain's corner angles are forced, its side lengths are free subject to closure and side
// pairing, a polyhedral cone. Fatness of a tile is the isoperimetric quotient 4 pi A / P^2 of its straight outline (1 for a
// circle, 0.79 for a square, near 0 for a sliver). Fattening a sliver can squeeze its neighbours' edges very short (John,
// 2026-10-21: "be careful not to have any edges get too small"), so edges count too: the shortest tile edge over the longest.
// A drawing's score is the weaker of the two, each scaled so 1 is comfortable -- fatness / FAT_OK (a square), edge ratio /
// EDGE_OK. Scoring alone still let the search trade edge for fatness until the two balanced (#72841: edge 0.44 -> 0.10), so the
// edge is also a hard FLOOR: a reshaped drawing's shortest edge may not fall below the default's, or below EDGE_OK when the
// default is already above it. Search: start at the canonical lengths, then random steps along the cone's null-space
// directions, keeping the best that respects the floor (a simple hill climb).
//
//   Fatten.score(dev)                  -> { score, fat, edge, byOrbit }
//   Fatten.best(sym, opts)             -> { realization | null, score, fat, edge, canonical, canFat, canEdge, evals }
//     opts: { evals (default 120), seed }
(function (root) {
  'use strict';
  const N = (typeof module !== 'undefined' && module.exports);
  const DD = N ? require('./draw_dsym2.js') : root.DrawDsym;
  const DF = N ? require('./dsym_deform.js') : root.DsymDeform;
  const DE = N ? require('./dsym_edit.js') : root.DsymEdit;

  function iq(pts) {
    let a = 0, p = 0;
    for (let i = 0; i < pts.length; i++) {
      const q = pts[i], r = pts[(i + 1) % pts.length];
      a += q[0] * r[1] - r[0] * q[1]; p += Math.hypot(r[0] - q[0], r[1] - q[1]);
    }
    return p > 0 ? 4 * Math.PI * (Math.abs(a) / 2) / (p * p) : 0;
  }
  // SIZE: a fat tile can still be a speck beside the other class (#8996 after fattening), so the smaller class's area over
  // the larger's counts as a third term. SIZE_OK = 1/2: the search keeps growing the small class until it is half the big one
  // (1/4 was too little -- John, 2026-10-22, on #25476, #1134, #1136, #8985: "make the smaller tile bigger").
  const FAT_OK = Math.PI / 4, EDGE_OK = 0.3, SIZE_OK = 0.5;
  function score(dev) {
    if (!dev || !dev.complete.length) return { score: 0, fat: 0, edge: 0, size: 0, byOrbit: {} };
    const byOrbit = {}, rep = {};
    for (const t of dev.complete) if (!(t.orbit in byOrbit)) { byOrbit[t.orbit] = iq(t.pts); rep[t.orbit] = t; }
    const v = Object.values(byOrbit);
    if (v.length < 2) return { score: 0, fat: 0, edge: 0, size: 0, byOrbit };
    let lo = Infinity, hi = 0;
    for (const t of Object.values(rep)) for (let i = 0; i < t.pts.length; i++) {
      const a = t.pts[i], b = t.pts[(i + 1) % t.pts.length], L = Math.hypot(b[0] - a[0], b[1] - a[1]);
      lo = Math.min(lo, L); hi = Math.max(hi, L);
    }
    const areas = Object.values(rep).map(t => Math.abs(polyA(t.pts)));
    const fat = Math.min(...v), edge = hi > 0 ? lo / hi : 0, size = Math.min(...areas) / Math.max(...areas);
    // geometric mean of the three, each capped at 1: taking the MINIMUM let the search push one term and drop the others
    // (#8996: size up, fatness 0.605 -> 0.118); the mean rewards improving all three, and best*() floors keep any from falling
    const c = x => Math.max(1e-6, Math.min(1, x));
    return { score: Math.cbrt(c(fat / FAT_OK) * c(edge / EDGE_OK) * c(size / SIZE_OK)), fat, edge, size, byOrbit };
  }
  // small deterministic PRNG, so a figure comes out the same every run
  function rng(seed) { let s = seed >>> 0 || 1; return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296); }

  function best(sym, opts) {
    opts = opts || {};
    const evals = opts.evals || 120, rand = rng(opts.seed || 12345);
    let fam; try { fam = DF.sampleFamily(sym, { steps: 1 }); } catch (e) { fam = null; }
    let canDev = null; try { canDev = DD.develop(sym, { maxCopies: 60 }); } catch (e) { canDev = null; }
    const can = score(canDev), canonical = can.score;
    if (!fam || !fam.sys.basis.length) return { realization: null, score: canonical, fat: can.fat, edge: can.edge,
                                               canonical, canFat: can.fat, canEdge: can.edge, evals: 0 };
    const { D, I, K, sys } = fam, basis = sys.basis, L0 = fam.samples[0].L;
    const floor = Math.min(can.edge, EDGE_OK) - 1e-9;
    const tryL = L => {
      if (Math.min(...L) <= 1e-6) return null;
      let C; try { C = DF.placeFromLengths(D, I, K, L); } catch (e) { C = null; }
      if (!C) return null;
      const R = { D, I, K, C };
      let dev; try { dev = DD.develop(sym, { maxCopies: 60, realization: R }); } catch (e) { dev = null; }
      const sc = score(dev);
      return sc.score > 0 && sc.edge >= floor ? { R, s: sc.score, sc } : null;
    };
    let cur = { L: L0, s: canonical, R: null, sc: can }, n = 0, step = 0.5;
    const scale = Math.max(...L0);
    while (n < evals) {
      const L = cur.L.slice();
      for (const v of basis) { const vn = Math.max(...v.map(Math.abs)) || 1, t = (rand() * 2 - 1) * step * scale / vn;
        for (let i = 0; i < L.length; i++) L[i] += t * v[i]; }
      n++;
      const r = tryL(L);
      if (r && r.s > cur.s + 1e-9) cur = { L, s: r.s, R: r.R, sc: r.sc };
      else if (n % 20 === 0) step *= 0.6;                       // nothing better nearby: look closer
    }
    return { realization: cur.R, score: cur.s, fat: cur.sc.fat, edge: cur.sc.edge,
             canonical, canFat: can.fat, canEdge: can.edge, evals: n };
  }

  // The same search over the EDITOR's model (dsym_edit), which is much larger: it frees every interior corner of the domain
  // and moves boundary corners under the real constraints (closure, cut congruence, per-orbit angle sums) instead of freezing
  // every angle. The length cone cannot move #8891's 3.8/3.8/172 degree triangle at all (place_from_points.js, 2026-08-22);
  // this can. Each step nudges one corner, keeps it only if the domain stays valid (orientation, simple, congruence, angles)
  // and the score improves without the shortest edge going under the floor.
  function bestEdit(sym, opts) {
    opts = opts || {};
    const evals = opts.evals || 400, rand = rng(opts.seed || 12345);
    let dev0 = null; try { dev0 = DD.develop(sym, { maxCopies: 60 }); } catch (e) { dev0 = null; }
    const can = score(dev0);
    const none = { realization: null, score: can.score, fat: can.fat, edge: can.edge, size: can.size,
                   canonical: can.score, canFat: can.fat, canEdge: can.edge, canSize: can.size, evals: 0 };
    if (!dev0 || !dev0.R) return none;
    let M; try { M = DE.buildEditModel(dev0.R.D, dev0.R.I, dev0.R.K, dev0.R.C); } catch (e) { return none; }
    const real = M2 => ({ D: M2.D, I: M2.I, K: M2.K, C: DE.toC(M2) });
    // INTEGRITY: develop() numbers its copies the same way in every realization of a type, so a sound drawing has exactly the
    // default's count of complete tiles per class. For 61 types (all with glides or only translations: 1**, 1*x, 22x, 1xx, o1)
    // even an UNCHANGED edit model fails this -- tiles go missing, leaving holes (#408, #1110; found 2026-10-22). Those keep the
    // default drawing, and every candidate below must pass.
    const counts = d => { const c = {}; for (const t of d.complete) c[t.orbit] = (c[t.orbit] || 0) + 1; return JSON.stringify(c); };
    const want = counts(dev0);
    const sound = d => !!d && counts(d) === want;
    let d1 = null; try { d1 = DD.develop(sym, { maxCopies: 60, realization: real(M) }); } catch (e) { d1 = null; }
    if (!sound(d1)) return Object.assign(none, { unsound: true });
    try { DE.setPins(M, DE.pinnedWalkPositions(M, dev0, DE.constraintAt)); } catch (e) { M.pins = []; }
    // Floors at 80% of where the search starts (capped at the comfortable levels): growing a small tile usually costs a little
    // fatness or edge, and a strict floor froze #1136, #8933 and #25476; 80% still stops any edge collapsing.
    const SLACK = 0.8;
    let floor = SLACK * Math.min(can.edge, EDGE_OK) - 1e-9, fatFloor = SLACK * Math.min(can.fat, FAT_OK) - 1e-9;
    // A lattice free to stretch (no rotation of order >= 3) but drawn nearly square reads as square: John, 2026-10-22, on #408
    // and #983, "aim for 2:3 or something that is easier to tell it's a rectangle". Stretch it to 3:2 first, and keep every
    // later step at least ASPECT_MIN away from square.
    const freedom = opts.orbifold != null ? DE.latticeFreedom(opts.orbifold) : 0;
    let axis = null, stretched = false;
    if (freedom >= 1) {
      const lat = DE.latticeBasis(dev0);
      try { axis = freedom === 1 ? DE.aspectAxis(M, lat) : (lat && lat.v1); } catch (e) { axis = null; }
      const r = axis && aspectOf(lat, axis), c = r && squareNear(r);
      if (c) {
        const T = DE.cloneModel(M);
        DE.setAspect(T, [0, 0], axis, (r >= c ? c * ASPECT : c / ASPECT) / r);
        let d2 = null; try { if (DE.validate(T, { cheapOnly: true }).ok) d2 = DD.develop(sym, { maxCopies: 60, realization: real(T) }); } catch (e) { d2 = null; }
        const r2 = d2 && aspectOf(DE.latticeBasis(d2), axis);
        if (r2 && !squareNear(r2) && sound(d2)) { M = T; stretched = true; }
      }
    }
    const aspectOk = dev => !stretched || (r => r && !squareNear(r * (1 + 1e-6)) && !squareNear(r * (1 - 1e-6)))(aspectOf(DE.latticeBasis(dev), axis));
    const area = dev0.complete.reduce((a, t) => a + Math.abs(polyA(t.pts)), 0) / dev0.complete.length;
    let step = 0.25 * Math.sqrt(area);
    let cur = { M, s: can.score, sc: can, R: null }, n = 0, since = 0;
    const path = [];                                  // every accepted step, so a caller can fall back to an earlier one
    if (stretched) { const R0 = real(M), d = DD.develop(sym, { maxCopies: 60, realization: R0 }), sc = score(d);
      cur = { M, s: sc.score, sc, R: R0 }; path.push({ R: R0, score: sc.score });
      // the floors protect what the search STARTS from, which is now the stretched drawing
      floor = SLACK * Math.min(sc.edge, EDGE_OK) - 1e-9; fatFloor = SLACK * Math.min(sc.fat, FAT_OK) - 1e-9; }
    const nCorners = M.P.length;
    while (n < evals) {
      n++; since++;
      const T = DE.cloneModel(cur.M);
      const c = Math.floor(rand() * nCorners);
      if (!T.P[c]) continue;
      const tgt = [T.P[c][0] + (rand() * 2 - 1) * step, T.P[c][1] + (rand() * 2 - 1) * step];
      let moved;
      try { moved = T.walkOfCorner[c] >= 0 ? DE.dragBoundary(T, T.walkOfCorner[c], tgt) : DE.dragInterior(T, c, tgt); }
      catch (e) { moved = null; }
      if (!moved || moved.ok === false) continue;
      let v; try { v = DE.validate(T, { cheapOnly: true }); } catch (e) { v = null; }
      if (!v || !v.ok) continue;
      const R = real(T);
      let dev; try { dev = DD.develop(sym, { maxCopies: 60, realization: R }); } catch (e) { dev = null; }
      if (!sound(dev)) continue;
      const sc = score(dev);
      if (sc.score > cur.s + 1e-9 && sc.edge >= floor && sc.fat >= fatFloor && aspectOk(dev)) { cur = { M: T, s: sc.score, sc, R }; since = 0; path.push({ R, score: sc.score }); }
      else if (since >= 40) { step *= 0.6; since = 0; }             // nothing better nearby: look closer
    }
    return { realization: cur.R, path, score: cur.s, fat: cur.sc.fat, edge: cur.sc.edge, size: cur.sc.size, stretched,
             canonical: can.score, canFat: can.fat, canEdge: can.edge, canSize: can.size, evals: n };
  }
  // Lattice aspect along `axis`: the shortest lattice vector parallel to it over the shortest perpendicular to it (among v1, v2,
  // v1+v2, v1-v2 -- a rectangle's sides, or a rhombus's diagonals). null when the lattice has no such pair.
  // A cell of aspect 1, 2 or 1/2 reads as squares: 1:1 directly, 1:2 because the domain is often half the cell (#408, 1**:
  // a 1x2 cell of 1x1 squares). squareNear(r) gives the nearest of those within ASPECT_MIN, and a stretch takes it to c * 1.5
  // (a 1x2 cell to 1x3, whose half is 2:3).
  const ASPECT = 1.5, ASPECT_MIN = 1.3;
  const squareNear = r => [1, 2, 0.5].find(c => Math.abs(Math.log(r / c)) < Math.log(ASPECT_MIN)) || null;
  function aspectOf(lat, axis) {
    if (!lat || !lat.v2 || !axis) return null;
    const L = Math.hypot(axis[0], axis[1]), u = [axis[0] / L, axis[1] / L];
    let along = Infinity, perp = Infinity;
    for (const v of [lat.v1, lat.v2, [lat.v1[0] + lat.v2[0], lat.v1[1] + lat.v2[1]], [lat.v1[0] - lat.v2[0], lat.v1[1] - lat.v2[1]]]) {
      const len = Math.hypot(v[0], v[1]); if (len < 1e-9) continue;
      const cr = Math.abs(u[0] * v[1] - u[1] * v[0]) / len, dt = Math.abs(u[0] * v[0] + u[1] * v[1]) / len;
      if (cr < 1e-6) along = Math.min(along, len); else if (dt < 1e-6) perp = Math.min(perp, len);
    }
    return isFinite(along) && isFinite(perp) ? along / perp : null;
  }
  function polyA(pts) { let a = 0; for (let i = 0; i < pts.length; i++) { const q = pts[i], r = pts[(i + 1) % pts.length]; a += q[0] * r[1] - r[0] * q[1]; } return a / 2; }

  const api = { score, best, bestEdit, iq };
  if (N) module.exports = api; else root.Fatten = api;
})(typeof window !== 'undefined' ? window : this);
