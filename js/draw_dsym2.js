// Draw a tiling from the cut-open disk (dsym_disk). Generators are the reflections in the mirror
// sides plus, for each cut edge, the side-pairing isometry carrying one copy onto the other.
// UMD so the interactive netshape page can run this in the browser as well as Node.
(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports)
    module.exports = factory(require('./dsym_geom.js'), require('./dsym_disk.js'));
  else root.DrawDsym = factory(root.DsymGeom, root.DsymDisk);
})(typeof window !== 'undefined' ? window : this, function (DG, DK) {
'use strict';
const PAIR = DG.PAIR;

const ID = [1, 0, 0, 0, 1, 0];
const apply = (T, p) => [T[0] * p[0] + T[1] * p[1] + T[2], T[3] * p[0] + T[4] * p[1] + T[5]];
const compose = (A, B) => [
  A[0] * B[0] + A[1] * B[3], A[0] * B[1] + A[1] * B[4], A[0] * B[2] + A[1] * B[5] + A[2],
  A[3] * B[0] + A[4] * B[3], A[3] * B[1] + A[4] * B[4], A[3] * B[2] + A[4] * B[5] + A[5]];
const invert = T => {
  const det = T[0] * T[4] - T[1] * T[3];
  const a = T[4] / det, b = -T[1] / det, d = -T[3] / det, e = T[0] / det;
  return [a, b, -(a * T[2] + b * T[5]), d, e, -(d * T[2] + e * T[5])];
};
function reflectionThrough(p, q) {
  const dx = q[0] - p[0], dy = q[1] - p[1], L = Math.hypot(dx, dy);
  const ux = dx / L, uy = dy / L, a = ux * ux - uy * uy, b = 2 * ux * uy;
  return [a, b, p[0] - (a * p[0] + b * p[1]), b, -a, p[1] - (b * p[0] - a * p[1])];
}
// The two isometries taking p1->q1 and p2->q2 (direct, then composed with a flip).
function isometryPair(p1, p2, q1, q2) {
  const vx = p2[0] - p1[0], vy = p2[1] - p1[1];
  const wx = q2[0] - q1[0], wy = q2[1] - q1[1];
  const lv = Math.hypot(vx, vy);
  const cs = (vx * wx + vy * wy) / (lv * lv), sn = (vx * wy - vy * wx) / (lv * lv);
  const direct = [cs, -sn, q1[0] - (cs * p1[0] - sn * p1[1]),
                  sn,  cs, q1[1] - (sn * p1[0] + cs * p1[1])];
  const ux = vx / lv, uy = vy / lv;
  const flip = [ux * ux - uy * uy, 2 * ux * uy, p1[0] - ((ux * ux - uy * uy) * p1[0] + 2 * ux * uy * p1[1]),
                2 * ux * uy, uy * uy - ux * ux, p1[1] - (2 * ux * uy * p1[0] + (uy * uy - ux * ux) * p1[1])];
  return [direct, compose(direct, flip)];
}
const side = (a, b, p) => (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);

function develop(sym, opts) {
  const maxCopies = (opts && opts.maxCopies) || 300;
  // opts.realization lets a DEFORMED domain (dsym_deform) be unfolded instead of the canonical one
  const R = (opts && opts.realization) || DK.realize(sym);
  if (!R) throw new Error('no disk realization');
  const { D, I, K, C } = R;
  const { P, E, N } = C;
  // CHAMBER COLOURING. Chambers either side of an interior edge are mirror images, so the flags 2-colour (+1/-1), and the
  // isometry gluing the next copy across a cut edge must have the handedness that puts the two chambers in opposite colours:
  // det = -col[f] * col[g] for the generator, -det(T) * col[f] * col[g] when crossing from a chamber placed by T. Choosing by
  // that instead of by which side of the edge a chamber's centre N lies on means the domain may contain INVERTED chambers
  // (John, 2026-10-22: "Maybe we can allow the skeleton to have retroflex triangles that go the opposite direction?") -- a
  // tile that isn't a fan from its centre, #420 with the blue tile rotated. The side test stays only as a fallback.
  const col = (() => {
    const c = new Int8Array(D.size + 1);
    for (let f0 = 1; f0 <= D.size; f0++) {
      if (c[f0]) continue;
      c[f0] = 1; const st = [f0];
      while (st.length) {
        const f = st.pop();
        for (let t = 0; t < 3; t++) {
          if (K.isB[I.nodes[f - 1].ep[t]]) continue;
          const g = D.sigma[t][f];
          if (!c[g]) { c[g] = -c[f]; st.push(g); } else if (c[g] !== -c[f]) return null;
        }
      }
    }
    return c;
  })();
  const detSign = T => Math.sign(T[0] * T[4] - T[1] * T[3]);

  const gens = [];
  // mirror sides -> reflections
  for (const w of K.walk) {
    if (!K.isMirror[w.edge]) continue;
    const e = I.edges[w.edge], f = e.nodes[0];
    const t = [0, 1, 2].find(tt => I.nodes[f - 1].ep[tt] === w.edge);
    const [a, b] = PAIR[t];
    gens.push(reflectionThrough(P[K.cornOf(f, a)], P[K.cornOf(f, b)]));
  }
  // cut edges -> side pairings
  for (const ei of K.cut) {
    const e = I.edges[ei];
    if (e.nodes.length !== 2) continue;
    const [f, g] = e.nodes;
    const t = [0, 1, 2].find(tt => I.nodes[f - 1].ep[tt] === ei);
    const [a, b] = PAIR[t];
    const p1 = P[K.cornOf(f, a)], p2 = P[K.cornOf(f, b)];
    const q1 = P[K.cornOf(g, a)], q2 = P[K.cornOf(g, b)];
    const cands = isometryPair(p1, p2, q1, q2);
    // keep the one putting the image of chamber f opposite chamber g across the shared side: by handedness from the colouring
    // (sound with inverted chambers), else by the side of the shared edge chamber f's centre lands on
    let best = null;
    if (col) best = cands.find(T => detSign(T) === -col[f] * col[g]) || null;
    if (!best) {
      const want = -Math.sign(side(q1, q2, N[g - 1]));
      for (const T of cands) {
        const s = Math.sign(side(q1, q2, apply(T, N[f - 1])));
        if (s === want) { best = T; break; }
      }
    }
    if (!best) best = cands[0];
    gens.push(best, invert(best));
  }
  if (!gens.length) throw new Error('no generators');

  const key = T => T.map(v => (Math.abs(v) < 1e-9 ? 0 : v).toFixed(6)).join(',');
  const copies = [ID], seen = new Set([key(ID)]);
  for (let h = 0; h < copies.length && copies.length < maxCopies; h++)
    for (const g of gens) {
      const T = compose(copies[h], g), k = key(T);
      if (!seen.has(k)) { seen.add(k); copies.push(T); }
    }

  const chambers = [];
  copies.forEach((T, h) => { for (let f = 1; f <= D.size; f++)
    chambers.push({ vert: apply(T, P[K.cornOf(f, 0)]), emid: apply(T, P[K.cornOf(f, 1)]),
                    ctr: apply(T, P[K.cornOf(f, 2)]),
                    tileOrbit: I.nodes[f - 1].op[2], emidOrbit: I.nodes[f - 1].op[1], f, T, h }); });

  // Snap targets for the interactive proof tool: every barycentric point, deduped by position and
  // labelled by kind. John's containment arguments need all three -- "tile centers, edge midpoints
  // as well as vertices".
  //
  // addr:{h,f,ty} is the (group-element index into `copies`, chamber, corner-type) that FIRST
  // reached this position -- a purely combinatorial address, not a coordinate. `gens` is built
  // from K.walk/K.cut alone (topology, never C), so the BFS below explores the same sequence of
  // abstract group words regardless of which realization's coordinates were used to number them;
  // copies[h] always denotes the same abstract word for a given h across any realization sharing
  // this D,I,K. That makes addr a stable identifier netshape's family test can use to re-locate
  // THIS SAME clicked point (not just "some point of the same orbit") after DF.placeFromLengths
  // deforms C -- recomputing apply(copies2[addr.h], P2[K.cornOf(addr.f, addr.ty)]) directly,
  // instead of guessing via orbit + spatial sort (which silently collapsed distinct same-orbit
  // points onto one, since an orbit can contain many geometrically different points).
  const KIND = ['vertex', 'edgeMid', 'tileCentre'];
  const snapMap = new Map();
  copies.forEach((T, h) => { for (let f = 1; f <= D.size; f++) for (let ty = 0; ty < 3; ty++) {
    const p = apply(T, P[K.cornOf(f, ty)]);
    const k = `${Math.round(p[0] / 1e-6)},${Math.round(p[1] / 1e-6)}`;
    if (!snapMap.has(k)) snapMap.set(k, { p, type: ty, kind: KIND[ty], orbit: I.nodes[f - 1].op[ty], addr: { h, f, ty } });
  } });
  const snapPoints = [...snapMap.values()];

  const TOL = 1e-6, buckets = new Map();
  const pkey = p => `${Math.round(p[0] / TOL)},${Math.round(p[1] / TOL)}`;
  for (const ch of chambers) {
    const k = pkey(ch.ctr);
    if (!buckets.has(k)) buckets.set(k, { ctr: ch.ctr, orbit: ch.tileOrbit, chs: [] });
    buckets.get(k).chs.push(ch);
  }
  // Trace each tile's own boundary COMBINATORIALLY (sigma1, sigma0, sigma1, ...), not by sorting
  // chambers by angle around the centroid or matching edges to their nearest midpoint -- both of
  // those are geometric proxies that silently break on a non-convex (reflex-angle) tile: angle-
  // around-the-centre isn't boundary order when the tile isn't star-shaped from its own centre,
  // and "nearest" edge-midpoint can belong to the WRONG edge when a short edge sits close to a
  // longer one's line. #357's notched pentagon (edges 0.5, 2.0, 0.5, 1.118, 1.118) hit exactly
  // that: the geometric match occasionally gave its long edge the short edge's orbit, drawing an
  // S motif sized for the short edge onto the long one where a U motif belonged.
  //
  // sigma1 moves a chamber's edge-midpoint corner and fixes its vertex and tile-centre corners --
  // so f and sigma1(f) are the two chambers meeting at one vertex, on either side of it. sigma0 is
  // the same idea the other way: f and sigma0(f) are the two chambers either side of one edge's
  // midpoint. Both preserve op[2] (the tile orbit), so this walk never leaves the current tile.
  //
  // The subtlety is WHERE the other chamber's corner POSITIONS live. K.cornOf gives each flag's
  // corners as DISK-INSTANCE indices, not orbifold-orbit indices -- when dsym_disk cuts the
  // orbifold open, an orbit that straddles the cut is duplicated into separate disk instances with
  // separate coordinates, even though sigma still (correctly, at the orbifold level) claims to fix
  // it. So "does crossing sigma_ty(f) keep the same corner positions" is NOT decided by whether
  // sigma_ty(f) equals f (a mirror fixed point) -- it's decided by whether the EDGE crossed
  // (ep[ty], of exactly the type being moved) is a disk boundary at all, K.isB[ep[ty]], which is
  // true for BOTH mirror edges (single flag) and cut edges (two flags, artificially split to open
  // the disk). An earlier version only special-cased the fixed-point/mirror case, so any chamber
  // reached by crossing a *cut* edge (two genuinely different flags, still boundary) got the wrong
  // T -- #56 and #1146 needed this; #49 has both a mirror and a cut edge cut by the same move type.
  // For a mirror, the far side is flag f again in a new copy T' = T . reflection through this
  // chamber's own two fixed corners (exactly the reflectionThrough gens above, evaluated at
  // whichever chamber the walk currently sits on). For a cut, the far side is the OTHER flag g of
  // the 2-flag edge, at whatever isometry actually carries this chamber's shared corners onto g's
  // own (disk-instance) corners -- computed the same way as the cut-edge gens above, just directly
  // in the walk's current world frame instead of the base frame develop()'s own gens use.
  function crossBoundary(T, f, ty) {
    const [a, b] = PAIR[ty];
    const p1 = apply(T, P[K.cornOf(f, a)]), p2 = apply(T, P[K.cornOf(f, b)]);
    if (K.isMirror[I.nodes[f - 1].ep[ty]])
      return { T: compose(T, reflectionThrough(P[K.cornOf(f, a)], P[K.cornOf(f, b)])), f };
    const g = D.sigma[ty][f];
    const q1 = P[K.cornOf(g, a)], q2 = P[K.cornOf(g, b)];
    const cands = isometryPair(q1, q2, p1, p2);           // isometries taking g's own corners onto here
    let best = null;
    if (col) best = cands.find(c2 => detSign(c2) === -detSign(T) * col[f] * col[g]) || null;   // by handedness, as above
    if (!best) {
      const want = -Math.sign(side(p1, p2, apply(T, N[f - 1])));
      for (const cand of cands) {
        if (Math.sign(side(p1, p2, apply(cand, N[g - 1]))) === want) { best = cand; break; }
      }
    }
    return { T: best || cands[0], f: g };
  }
  const tiles = [];
  for (const b of buckets.values()) {
    const start = b.chs[0];
    let T = start.T, f = start.f;
    // vaddr: for each vertex, the isometry and domain corner it is the image of (vertex = T(P[corner])), so an editor can turn
    // a drag on any drawn tile corner into a drag of the domain corner behind it. tid: the tile's combinatorial identity (its
    // first chamber's copy and flag) -- the same tile in every realization of the type, so a view can follow it while reshaping.
    const pts = [start.vert], edgeOrbits = [], vaddr = [{ T: start.T, corner: K.cornOf(start.f, 0) }];
    let closed = false, guard = b.chs.length + 2;
    while (guard-- > 0) {
      if (K.isB[I.nodes[f - 1].ep[1]]) { const r = crossBoundary(T, f, 1); T = r.T; f = r.f; }
      else f = D.sigma[1][f];
      edgeOrbits.push(I.nodes[f - 1].op[1]);
      if (K.isB[I.nodes[f - 1].ep[0]]) { const r = crossBoundary(T, f, 0); T = r.T; f = r.f; }
      else f = D.sigma[0][f];
      const v = apply(T, P[K.cornOf(f, 0)]);
      pts.push(v); vaddr.push({ T, corner: K.cornOf(f, 0) });
      if (f === start.f && Math.hypot(v[0] - start.vert[0], v[1] - start.vert[1]) < 1e-6) {
        pts.pop(); vaddr.pop(); closed = true; break;
      }
    }
    if (closed && pts.length >= 3 && pts.length * 2 === b.chs.length)
      tiles.push({ ctr: b.ctr, orbit: b.orbit, pts, edgeOrbits, nch: b.chs.length, vaddr, tid: start.h + ':' + start.f });
  }
  const want = new Map();
  for (let f = 1; f <= D.size; f++) { const o = I.nodes[f - 1].op[2]; want.set(o, 2 * K.ord[o].m); }
  const complete = tiles.filter(t => t.nch === want.get(t.orbit));
  return { R, copies, chambers, tiles, complete, gens, snapPoints };
}

const polyArea = pts => { let a = 0;
  for (let i = 0; i < pts.length; i++) { const j = (i + 1) % pts.length;
    a += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]; } return Math.abs(a) / 2; };
const diameter = pts => { let m = 0;
  for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++)
    m = Math.max(m, Math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1])); return m; };

// Self-check: probe the patch and confirm every point is covered by exactly one tile. This is the
// only verification that measures the thing we care about -- two cheaper "invariants" were tried
// and both were wrong (an area identity that assumed equal chamber areas, and a chamber
// orientation test that ignored that barycentric chambers alternate handedness by construction).
// Known failure: #25537 (22x), where the cut splits two regular points three ways and Tegula's
// equal-share alph = 2*PI/(v*f*s) does not hold; ~2% of its area ends up triple covered.
function verify(dev, nPoints) {
  const n = nPoints || 4000;
  const rmax = Math.max(...dev.complete.map(t => Math.hypot(t.ctr[0], t.ctr[1])));
  const R = Math.min(0.9, rmax * 0.28) || 0.4;
  const inPoly = (p, poly) => { let c = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) { const a = poly[i], b = poly[j];
      if (((a[1] > p[1]) !== (b[1] > p[1])) && (p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0])) c = !c; }
    return c; };
  const cov = {};
  for (let k = 1; k <= n; k++) {
    const th = k * 2.399963, r = R * Math.sqrt((k % (n - 11) + 0.5) / (n - 11));
    const p = [r * Math.cos(th), r * Math.sin(th)];
    let c = 0; for (const t of dev.complete) if (inPoly(p, t.pts)) c++;
    cov[c] = (cov[c] || 0) + 1;
  }
  return { ok: cov[1] === n, coverage: cov, probeRadius: R };
}

return { develop, verify, polyArea, diameter, apply, compose, invert };
});
