// Deform the fundamental domain of a type, so a claim can be tested over the WHOLE family of
// shapes rather than at Tegula's single canonical member.
//
// Tegula fixes the domain: corners sit on the tangential polygon (dist = RAD/sin(alph/2)), which
// is one point of the family. But only the corner ANGLES are forced -- they are alph = 2*pi/i,
// the rotation orders. The side LENGTHS are free, subject to:
//    (a) closure       sum_k L_k * (cos th_k, sin th_k) = 0
//    (b) side pairing  L_a = L_b for the two boundary copies of each cut edge, or the pairing
//                      stops being an isometry and the unfolding breaks
// Both are linear, so the admissible lengths form a polyhedral cone: a linear subspace
// intersected with L_k > 0. This module finds that cone and walks around inside it.
//
// CAUTION -- gauge freedom. Not every direction in the cone changes the tiling. A regular orbit
// (v = 1) lying on a mirror can slide along it: that moves the chamber decomposition but leaves
// every tile vertex, and so every tile, exactly where it was. #4 shows this (dim 2, but the tiles
// are rigid equilateral triangles). Do not read the cone's dimension as "the tiling deforms" --
// measure the quantity you care about and see whether it actually moves.
(function (root) {
  'use strict';
  const DG = (typeof require !== 'undefined') ? require('./dsym_geom.js') : root.DsymGeom;
  const DK = (typeof require !== 'undefined') ? require('./dsym_disk.js') : root.DsymDisk;
  const PAIR = DG.PAIR;

  // Reduced row echelon form; returns {rank, basis} where basis spans the null space.
  function nullBasis(rows, n) {
    const A = rows.map(r => r.slice()), m = A.length;
    const pivotOf = [];
    let rank = 0;
    for (let c = 0; c < n && rank < m; c++) {
      let piv = -1, best = 1e-9;
      for (let r = rank; r < m; r++) if (Math.abs(A[r][c]) > best) { best = Math.abs(A[r][c]); piv = r; }
      if (piv < 0) continue;
      [A[rank], A[piv]] = [A[piv], A[rank]];
      const p = A[rank][c];
      for (let j = 0; j < n; j++) A[rank][j] /= p;
      for (let r = 0; r < m; r++) if (r !== rank && Math.abs(A[r][c]) > 1e-13) {
        const f = A[r][c];
        for (let j = 0; j < n; j++) A[r][j] -= f * A[rank][j];
      }
      pivotOf.push(c); rank++;
    }
    const isPivot = new Set(pivotOf), basis = [];
    for (let c = 0; c < n; c++) {
      if (isPivot.has(c)) continue;
      const v = new Array(n).fill(0); v[c] = 1;
      pivotOf.forEach((pc, ri) => { v[pc] = -A[ri][c]; });
      basis.push(v);
    }
    return { rank, basis };
  }

  // The constraint system on boundary side lengths, one variable per walk position.
  function lengthSystem(K) {
    const walk = K.walk, n = walk.length, th = new Array(n);
    let a = 0;
    for (let j = 0; j < n; j++) { a += Math.PI - K.alphOf.get(walk[j].corner); th[j] = a; }
    const rows = [new Array(n).fill(0), new Array(n).fill(0)];
    for (let j = 0; j < n; j++) { rows[0][j] = Math.cos(th[j]); rows[1][j] = Math.sin(th[j]); }
    const byEdge = new Map();
    walk.forEach((w, j) => {
      if (!K.cut.has(w.edge)) return;
      if (!byEdge.has(w.edge)) byEdge.set(w.edge, []);
      byEdge.get(w.edge).push(j);
    });
    for (const [, js] of byEdge) for (let k = 1; k < js.length; k++) {
      const r = new Array(n).fill(0); r[js[0]] = 1; r[js[k]] = -1; rows.push(r);
    }
    const { basis } = nullBasis(rows, n);
    return { th, rows, basis, n };
  }

  // Side lengths of the canonical (tangential) placement -- one interior point of the cone.
  function canonicalLengths(D, I, K) {
    const C = DK.computeCoords(D, I, K);
    if (!C) return null;
    const walk = K.walk, n = walk.length, L = new Array(n);
    for (let j = 0; j < n; j++) {
      const p = C.P[walk[j].corner], q = C.P[walk[(j + 1) % n].corner];
      L[j] = Math.hypot(q[0] - p[0], q[1] - p[1]);
    }
    return { L, C };
  }

  // Place the domain from an explicit length vector, then relax the interior exactly as before.
  function placeFromLengths(D, I, K, L) {
    const { th } = lengthSystem(K);
    const walk = K.walk, n = walk.length;
    const P = new Array(K.nCorners).fill(null);
    let x = 0, y = 0;
    for (let j = 0; j < n; j++) {
      P[walk[j].corner] = [x, y];
      x += L[j] * Math.cos(th[j]); y += L[j] * Math.sin(th[j]);
    }
    if (Math.hypot(x, y) > 1e-6 * Math.max(1, ...L)) return null;   // did not close
    let cx = 0, cy = 0, c = 0;
    for (const p of P) if (p) { cx += p[0]; cy += p[1]; c++; }
    cx /= c; cy /= c;
    for (let i = 0; i < P.length; i++) if (!P[i]) P[i] = [cx, cy];
    // recentre so shapes are comparable between samples
    for (const p of P) { p[0] -= cx; p[1] -= cy; }

    const E = new Array(K.nSides).fill(null);
    for (let f = 1; f <= D.size; f++) for (let t = 0; t < 3; t++) {
      const si = K.sideOf(f, t), [a, b] = PAIR[t];
      const p = P[K.cornOf(f, a)], q = P[K.cornOf(f, b)];
      E[si] = [0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])];
    }
    const N = [];
    for (let f = 1; f <= D.size; f++) {
      let sx = 0, sy = 0;
      for (let t = 0; t < 3; t++) { sx += E[K.sideOf(f, t)][0] + P[K.cornOf(f, t)][0];
                                    sy += E[K.sideOf(f, t)][1] + P[K.cornOf(f, t)][1]; }
      N.push([sx / 6, sy / 6]);
    }
    const C = { P, E, N, corners: walk.map(w => w.corner) };
    DK.relax(D, I, K, C);
    return C;
  }

  // Samples of the cone: the canonical point, then steps along each null-space direction that
  // keep every length strictly positive.
  function sampleFamily(sym, opts) {
    const steps = (opts && opts.steps) || 4, spread = (opts && opts.spread) || 0.6;
    const D = DG.parseSymbol(sym), I = DG.incidence(D);
    const K = DK.buildDisk(D, I);
    if (!K) return null;
    const sys = lengthSystem(K);
    const can = canonicalLengths(D, I, K);
    if (!can) return null;
    const L0 = can.L, out = [{ L: L0.slice(), tag: 'canonical' }];
    const minL = Math.min(...L0);
    for (let b = 0; b < sys.basis.length; b++) {
      const v = sys.basis[b];
      const vn = Math.max(...v.map(Math.abs)) || 1;
      for (let s = 1; s <= steps; s++) for (const sign of [1, -1]) {
        const t = sign * spread * (s / steps) * minL / vn;
        const L = L0.map((x, i) => x + t * v[i]);
        if (Math.min(...L) <= 1e-6) continue;
        out.push({ L, tag: `dir${b}${sign > 0 ? '+' : '-'}${s}` });
      }
    }
    return { D, I, K, sys, samples: out, dim: sys.basis.length };
  }

  const api = { nullBasis, lengthSystem, canonicalLengths, placeFromLengths, sampleFamily };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.DsymDeform = api;
})(typeof window !== 'undefined' ? window : this);
