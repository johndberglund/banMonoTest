// Edge-curve model for the tiler.  A curve lives per edge-ORBIT in a canonical
// frame: the edge is the segment (0,0)->(1,0) and the curve is a function-like
// polyline [[x,y],...] with x increasing 0..1, endpoints (0,0) and (1,0); y is the
// perpendicular offset in edge-length units.  Symmetry (J/I/U/S) is baked in at
// store time so rendering and glue-consistency are automatic.
(function (root) {
  'use strict';

  // Build the canonical curve from a drawn polyline, enforcing the edge symmetry.
  //
  // The curve is parametrised by POINT ORDER, not by x.  It may loop back, overhang or
  // self-touch — it is a polyline from (0,0) to (1,0), not a function of x.  So nothing
  // here may sort by x: that would reorder the points of a non-monotone curve and quietly
  // destroy it.  (An earlier `enforceSym` did sort, which is why the overhanging motif of
  // curve_examples.svg could not be represented; the editor's own builder never did.)
  //
  //   I: straight.
  //   J: the drawn polyline, endpoints pinned to (0,0) and (1,0).
  //   U: drawn half + its mirror across x=0.5; the LAST drawn point is the fold and is
  //      snapped to x=0.5, keeping its offset.
  //   S: drawn half + its 180-degree rotation about (0.5,0); the last drawn point is the
  //      centre and is snapped to (0.5,0), since the centre lies ON the curve.
  //
  // This is the single source of truth — tilerTest.js delegates to it, so the editor, the
  // demo motifs and the developed tiling cannot drift apart.
  function buildCanon(type, pts) {
    if (type === 'I' || !pts || pts.length < 2) return [[0, 0], [1, 0]];
    const p = pts.map(q => q.slice()); p[0] = [0, 0];
    if (type === 'J') { p[p.length - 1] = [1, 0]; return p; }
    const out = p.slice();
    if (type === 'U') {
      p[p.length - 1] = [0.5, p[p.length - 1][1]]; out[out.length - 1] = p[p.length - 1];
      for (let i = p.length - 2; i >= 0; i--) out.push([1 - p[i][0], p[i][1]]);
    } else {
      p[p.length - 1] = [0.5, 0]; out[out.length - 1] = [0.5, 0];
      for (let i = p.length - 2; i >= 0; i--) out.push([1 - p[i][0], -p[i][1]]);
    }
    return out;
  }
  const enforceSym = (pts, type) => buildCanon(type, pts);      // legacy argument order

  // U/S store only the drawn HALF (up to the fold) as editable points; buildCanon appends
  // the other half by reflecting (U: x->1-x) or point-rotating (S: x->1-x, y->-y) each raw
  // point in reverse order. Both maps are involutions, so the SAME formula both derives a
  // mirror point from a raw point (used above) and recovers a raw point from a drag on the
  // mirror half (used by the editor below) -- one map, read in either direction.
  function mirrorPt(type, p) { return type === 'S' ? [1 - p[0], -p[1]] : [1 - p[0], p[1]]; }
  // Map an index j into the FULL canonical curve (buildCanon's output, length L=2n-1 for
  // U/S, length n for J/I) back to the underlying editable raw index r (0..n-1) plus
  // whether that full-curve point is the mirror copy rather than the raw point itself.
  // j<n is the raw half verbatim (r=j); j>=n is the reflected/rotated second half, whose
  // raw source is r=(2n-2)-j (buildCanon pushes it in reverse order starting at i=n-2).
  function fullIdxToRaw(type, n, j) {
    if (type !== 'U' && type !== 'S' || j < n) return { r: j, mirror: false };
    return { r: (2 * n - 2) - j, mirror: true };
  }

  // transport a canonical curve to a congruent edge by the engine's Klein-4 map
  // code (tileEdges Mapping / netEdgeData row3).  bit&1 = reverse parameter
  // (x->1-x); bit&2 = negate offset (y->-y).  This is what makes congruent edges
  // render as the engine's symmetry demands: e.g. two S edges related by a tile
  // reflection (code 1) come out as proper mirror images, since reversing the
  // parameter of a point-symmetric S curve flips its offset.
  function applyMap(canon, code) {
    let p = canon.map(q => q.slice());
    // reverse parameter = mirror x AND reverse point order (works for freeform, non-monotone
    // curves; do NOT sort by x — that would scramble a general polyline, e.g. break S curves).
    if (code & 1) p = p.map(q => [1 - q[0], q[1]]).reverse();
    if (code & 2) p = p.map(q => [q[0], -q[1]]);
    return p;
  }

  // ---- demo motif library, v5: traced directly from John's colour-coded `JSU.svg` (2026-08-10)
  // -- 4 hand-drawn examples per curve type, no more needed ("I don't think that we will need
  // more than four of any type of curve" -- John). Earlier versions shared ONE "half-motif" pool
  // across J/U/S and built U/S by doubling it up ourselves; that was a misreading of the source
  // drawings. The SVG's colours make the real intent explicit: black paths ARE complete J curves
  // (already asymmetric by design -- use as drawn), red paths are ALREADY full point-symmetric S
  // curves, green paths are ALREADY full mirror-symmetric U curves ("The Red and green already
  // have this built in. So just take half of the pattern if you need to."). So J motifs are used
  // whole; S/U motifs are halved back down to editor points, which buildCanon re-doubles.
  // Coordinates: each path's own start/end normalised to canonical (0,0)/(1,0), uniform scale
  // (same factor on x and y) so the traced silhouette isn't stretched.
  const J_MOTIFS = [
    [[0, 0], [0.427817, 0], [0.427817, 0.107394], [0.613829, 0], [1, 0]],
    [[0, 0], [0.386171, 0], [0.572183, 0.107394], [0.739437, 0.107394], [0.572183, 0], [1, 0]],
    [[0, 0], [0.386171, 0], [0.572183, 0.107394], [0.572183, 0.053697], [0.684859, 0.053697], [0.684859, 0], [1, 0]],
    [[0, 0], [0.357394, 0], [0.357394, 0.107394], [0.5, 0], [0.5, 0.107394], [0.642606, 0], [1, 0]],
  ];
  // S/U store only the EDITABLE HALF (buildCanon appends the rest by rotating/mirroring it) --
  // each ends at the fold. S's fold is forced to (0.5,0) by buildCanon regardless of what's
  // stored here; U's fold keeps its drawn offset. Recovered from the SVG's full (already-
  // symmetric) curve by keeping its first half and, where the source had no vertex exactly on
  // the fold (an even point count -- the fold falls mid-segment by construction of the
  // symmetry), inserting that segment's midpoint as the fold vertex.
  const S_MOTIFS = [
    [[0, 0], [0.357394, 0], [0.5, 0.107394], [0.5, 0]],
    [[0, 0], [0.357394, 0], [0.357394, 0.107394], [0.5, 0.107394], [0.5, 0]],
    [[0, 0], [0.357394, 0], [0.357394, 0.107394], [0.5, 0]],
    [[0, 0], [0.357394, 0], [0.357394, 0.107394], [0.591549, 0.107394], [0.5, 0]],
  ];
  const U_MOTIFS = [
    [[0, 0], [0.357394, 0], [0.5, 0.107394]],
    [[0, 0], [0.357394, 0], [0.357394, 0.107394], [0.5, 0.107394]],
    [[0, 0], [0.408451, 0], [0.357394, 0.107394], [0.5, 0.107394]],
    [[0, 0], [0.357394, 0], [0.357394, 0.107394], [0.5, 0]],
  ];
  const wrap = (v, n) => ((v | 0) % n + n) % n;
  const MOTIF_BANK = { J: J_MOTIFS, S: S_MOTIFS, U: U_MOTIFS };
  // Variant space per type is DOUBLED: indices [0,N) are the motif as drawn, [N,2N) are the
  // SAME motif with y negated -- flipped to the other side of the baseline (out<->in). Kept as
  // a deliberate escape hatch, not just cosmetic variety: at a small enough tile angle, a curve
  // planted toward the interior can still clip a neighbouring edge even after shifting its bump
  // toward the centre of the edge to shrink the angle at the corner; flipping the whole curve to
  // the other side of the baseline is often enough to clear the intersection without hand-
  // editing every point (John, 2026-08-10: "sometimes the curves ran into each other inside the
  // tile... The mirror variant was just the idea that we could flip a curve from out to in").
  // Negating y on the stored RAW points (J's full curve, or U/S's editable half) flips the
  // WHOLE resulting curve, not just the stored half: buildCanon's U mirror keeps y's sign, so a
  // negated half comes out negated on both sides; its S rotation negates y again for the second
  // half, so a negated half comes out negated on both sides there too (double negative).
  function motifCount(type) { const b = MOTIF_BANK[type]; return b ? b.length * 2 : 0; }
  function motifName(type, variant) {
    const b = MOTIF_BANK[type]; if (!b || !b.length) return '';
    const n = b.length, v = wrap(variant, n * 2);
    return type + '-' + (v % n + 1) + (v >= n ? ' (flipped)' : '');
  }
  // The motif as EDITOR points for this edge type — what gets stored in curveEdits and handed
  // to buildCanon, so a demo curve is exactly a drawn curve and stays draggable. J hands its
  // full drawn curve straight through; U/S hand their drawn HALF straight through (buildCanon
  // mirrors/rotates it into the other half).
  function motifPts(type, variant) {
    const b = MOTIF_BANK[type]; if (!b || !b.length) return [];
    const n = b.length, v = wrap(variant, n * 2);
    const pts = b[v % n];
    return v >= n ? pts.map(([x, y]) => [x, -y]) : pts.map(p => p.slice());
  }
  // a demo curve per type, as a finished canonical curve.  `variant` picks the motif, so
  // different edge orbits get different marks — what makes a matched pair identifiable.
  function testCurve(type, variant) { return buildCanon(type, motifPts(type, variant)); }

  // map a canonical curve onto a world edge P0->P1 (+y to the left; offset scales
  // with edge length).  flip reverses the parametrisation (x->1-x) AND the offset,
  // i.e. renders the *same* physical curve from the other end.
  function mapToEdge(canon, P0, P1, flip) {
    const Ex = P1[0] - P0[0], Ey = P1[1] - P0[1];            // edge vector (length 1 in canon units)
    const Px = -Ey, Py = Ex;                                 // +90 deg (left) perpendicular
    const out = [];
    for (let i = 0; i < canon.length; i++) {
      const c = flip ? canon[canon.length - 1 - i] : canon[i];
      const x = flip ? 1 - c[0] : c[0], y = flip ? -c[1] : c[1];
      out.push([P0[0] + x * Ex + y * Px, P0[1] + x * Ey + y * Py]);
    }
    return out;
  }

  // Inverse of mapToEdge's per-point placement: given a world point on (or near) the edge
  // P0->P1, recover its canonical (x,y). E and its +90 perpendicular Px,Py have equal length
  // (|E|), so projecting (W-P0) onto each and dividing by |E|^2 is exact -- no least-squares
  // needed, since mapToEdge is a similarity transform (rotate + uniform scale + translate).
  function unmapFromEdge(P0, P1, W) {
    const Ex = P1[0] - P0[0], Ey = P1[1] - P0[1];
    const Px = -Ey, Py = Ex;
    const L2 = Ex * Ex + Ey * Ey || 1;
    const dx = W[0] - P0[0], dy = W[1] - P0[1];
    return [(dx * Ex + dy * Ey) / L2, (dx * Px + dy * Py) / L2];
  }

  // edge-curve orbits per config: union-find over row0 (congruent) and row2 (glued).
  function edgeOrbits(row0, row2, eshape) {
    const n = row0.length, par = [...Array(n).keys()];
    const find = x => { while (par[x] !== x) { par[x] = par[par[x]]; x = par[x]; } return x; };
    const uni = (a, b) => { par[find(a)] = find(b); };
    for (let g = 0; g < n; g++) { uni(g, row0[g]); uni(g, row2[g]); }
    const rep = g => find(g), type = {}, members = {};
    for (let g = 0; g < n; g++) { (members[rep(g)] = members[rep(g)] || []).push(g); type[rep(g)] = eshape[g]; }
    return { rep, type, members };
  }

  const ap = (T, p) => [T[0] * p[0] + T[1] * p[1] + T[4], T[2] * p[0] + T[3] * p[1] + T[5]];
  const SYM = { 0: 'J', 1: 'U', 2: 'I', 3: 'S' };           // engine edgeSym codes

  // Build curved boundaries for a developed tiling.  Curves live on the m-gon edges
  // (the tile's actual sides) using the engine's per-edge symmetry, NOT the net
  // edges.  `baseMgon` = the m solved tile vertices (verticesOf); `curveOf[i]` = the
  // canonical curve for m-gon edge i.  Each placed tile is T(baseMgon).  Each physical
  // edge gets ONE world curve (first tile to touch it defines it); neighbours reuse
  // the identical points, so glued curves coincide and reflected tiles mirror
  // automatically.  Returns [{orbit, pts:[...closed boundary...]}].
  function buildCurvedTiles(placed, baseMgon, curveOf) {
    const m = baseMgon.length;
    const kp = p => Math.round(p[0] * 1000) + ',' + Math.round(p[1] * 1000);
    const ek = (a, b) => { const a1 = kp(a), b1 = kp(b); return a1 < b1 ? a1 + '|' + b1 : b1 + '|' + a1; };
    const phys = {}, out = [];
    for (const t of placed) {
      const full = baseMgon.map(v => ap(t.T, v)), bd = [];
      for (let i = 0; i < m; i++) {
        const P0 = full[i], P1 = full[(i + 1) % m], key = ek(P0, P1);
        let wc = phys[key];
        if (!wc) { wc = mapToEdge(curveOf[i], baseMgon[i], baseMgon[(i + 1) % m], false).map(p => ap(t.T, p)); phys[key] = wc; }
        const seg = (kp(wc[0]) === kp(P0)) ? wc : wc.slice().reverse();
        for (let j = (i === 0 ? 0 : 1); j < seg.length; j++) bd.push(seg[j]);
      }
      out.push({ orbit: t.orbit, pts: bd });
    }
    return out;
  }

  // curve on a single (untransformed) tile, for the per-orbit view
  function curvedTile(baseMgon, curveOf) {
    const m = baseMgon.length, bd = [];
    for (let i = 0; i < m; i++) {
      const seg = mapToEdge(curveOf[i], baseMgon[i], baseMgon[(i + 1) % m], false);
      for (let j = (i === 0 ? 0 : 1); j < seg.length; j++) bd.push(seg[j]);
    }
    return bd;
  }

  // Interior angle at each corner measured between the two edge CURVES' tangents, not between the
  // chords. A curve leaves its endpoint at an angle to its chord, so a corner the closing family
  // forces flat is a real, open corner once the curves are drawn -- which is John's rule that a
  // forced-flat angle is not fatal, stated as a number. Measured on his #25543 m7 si5: corner 0 is
  // a 360-degree chord cusp and a 37.25-degree corner in the drawing. On #209516 si7 corner 7, 0.00
  // becomes 313.57. Straight (I) edges return exactly the chord angle, so nothing else shifts.
  //
  // Same convention as tiling_checks.interiorAngle: feed it the neighbouring TANGENT points instead
  // of the neighbouring corners and the existing sign and orientation handling carries over.
  function curveCornerAngles(baseMgon, curveOf) {
    const m = baseMgon.length;
    const seg = [];
    for (let i = 0; i < m; i++) seg.push(mapToEdge(curveOf[i], baseMgon[i], baseMgon[(i + 1) % m], false));
    let sgn = 0;
    for (let i = 0; i < m; i++) { const q = baseMgon[(i + 1) % m];
      sgn += baseMgon[i][0] * q[1] - q[0] * baseMgon[i][1]; }
    sgn = sgn >= 0 ? 1 : -1;
    const out = [];
    for (let i = 0; i < m; i++) {
      const p = seg[(i - 1 + m) % m], n = seg[i], at = baseMgon[i];
      // one sample in from each side; a degenerate segment falls back to the chord
      const prev = p.length >= 2 ? p[p.length - 2] : baseMgon[(i - 1 + m) % m];
      const next = n.length >= 2 ? n[1] : baseMgon[(i + 1) % m];
      const v1 = [prev[0] - at[0], prev[1] - at[1]], v2 = [next[0] - at[0], next[1] - at[1]];
      let t = Math.atan2(sgn * (v2[0]*v1[1] - v2[1]*v1[0]), v2[0]*v1[0] + v2[1]*v1[1]) * 180 / Math.PI;
      out.push(t < 0 ? t + 360 : t);
    }
    return out;
  }

  const api = { buildCanon, enforceSym, applyMap, testCurve, mapToEdge, unmapFromEdge, edgeOrbits,
                buildCurvedTiles, curvedTile, curveCornerAngles, SYM, MOTIFS: MOTIF_BANK,
                motifCount, motifName, motifPts, mirrorPt, fullIdxToRaw };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.TilerCurves = api;
})(typeof window !== 'undefined' ? window : this);
