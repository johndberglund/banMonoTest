// Generalised Euclidean realization: cut the orbifold complex open to a DISK, place it as a
// tangential polygon, relax, and unfold.
//
// dsym_geom's first version only handled complexes that were already a disk with every cone point
// on the mirror boundary. An interior cone point cannot sit inside a planar fundamental domain --
// the rotation about it would map the domain onto overlapping copies -- so the complex must be
// cut open. That is Tegula's `s` (splits), and it is why `alph = 2*PI / (v*f*s)`: s instances of
// the orbit each take an equal share of the total angle.
//
// THE RULE THAT GOVERNS EVERY CUT: a cut may not have a free end in the middle of the surface.
// A cut terminating at a regular point leaves a slit tip whose orbit has i = 1, hence alph = 2*PI
// and dist = RAD / cos(beta/2) diverges. Empirically that failure shows up as a closure residual
// of exactly PI per slit -- never an arbitrary value -- which is how it was diagnosed. So every
// cut is either an ARC between anchors (cone points or mirror boundary) or a CLOSED LOOP.
//
// Cutting means one orbit can own several positions, which is exactly why Tegula separates Orbit
// from OrbitCoordinates. Coordinates live on INSTANCES, found by union-find over chamber corners
// and chamber sides.
(function (root) {
  'use strict';
  const DG = (typeof require !== 'undefined') ? require('./dsym_geom.js') : root.DsymGeom;
  const PAIR = DG.PAIR;
  const RAD = DG.RAD;

  function UF(n) {
    const p = new Int32Array(n); for (let i = 0; i < n; i++) p[i] = i;
    const find = x => { while (p[x] !== x) { p[x] = p[p[x]]; x = p[x]; } return x; };
    return { find, union: (a, b) => { a = find(a); b = find(b); if (a !== b) p[a] = b; } };
  }

  // ---- instance structure for a given cut set ----
  function instancesFor(D, I, cut, isMirror) {
    const n = D.size, sid = (f, t) => (f - 1) * 3 + t;
    const isB = I.edges.map((e, i) => isMirror[i] || cut.has(i));
    const ufS = UF(n * 3), ufC = UF(n * 3);
    for (let f = 1; f <= n; f++) for (let t = 0; t < 3; t++) {
      if (!cut.has(I.nodes[f - 1].ep[t])) ufS.union(sid(f, t), sid(D.sigma[t][f], t));
      for (const d of PAIR[t]) if (!isB[I.nodes[f - 1].ep[d]]) ufC.union(sid(f, t), sid(D.sigma[d][f], t));
    }
    const sMap = new Map(), cMap = new Map();
    const grab = (uf, key, map) => { const r = uf.find(key);
      if (!map.has(r)) map.set(r, map.size); return map.get(r); };
    const sideOf = (f, t) => grab(ufS, sid(f, t), sMap);
    const cornOf = (f, c) => grab(ufC, sid(f, c), cMap);
    for (let f = 1; f <= n; f++) for (let t = 0; t < 3; t++) { sideOf(f, t); cornOf(f, t); }
    return { isB, sideOf, cornOf, nSides: sMap.size, nCorners: cMap.size };
  }

  // Walk one boundary component starting from chamber f0's side of type t0.
  function walkOne(D, I, inst, f0, t0, seenSides) {
    const c0 = PAIR[t0][0], walk = [];
    let f = f0, t = t0, c = c0, guard = 60 * D.size + 60;
    do {
      let d;
      for (;;) {
        const [a, b] = PAIR[c];
        d = (t === a) ? b : a;
        if (inst.isB[I.nodes[f - 1].ep[d]]) break;
        f = D.sigma[d][f]; t = d;
        if (--guard <= 0) return null;
      }
      const si = inst.sideOf(f, d);
      if (seenSides) seenSides.add(si);
      walk.push({ corner: inst.cornOf(f, c), orbit: I.nodes[f - 1].op[c], side: si,
                  edge: I.nodes[f - 1].ep[d], flag: f, ctype: c, stype: d });
      t = d;
      const [a, b] = PAIR[t];
      c = (c === a) ? b : a;
    } while (!(f === f0 && t === t0 && c === c0) && --guard > 0);
    return guard > 0 ? walk : null;
  }

  // All boundary components, plus the Euler characteristic of the cut-open complex.
  function topology(D, I, inst) {
    const seen = new Set(), walks = [];
    for (let f = 1; f <= D.size; f++) for (let t = 0; t < 3; t++) {
      if (!inst.isB[I.nodes[f - 1].ep[t]]) continue;
      if (seen.has(inst.sideOf(f, t))) continue;
      const w = walkOne(D, I, inst, f, t, seen);
      if (!w) return null;
      walks.push(w);
    }
    const chi = inst.nCorners - inst.nSides + D.size;
    return { walks, chi, isDisk: chi === 1 && walks.length === 1 };
  }

  // ---- orbit-graph helpers, used to pick cuts ----
  function orbitAdj(I, cut) {
    const adj = new Map();
    I.edges.forEach((e, i) => {
      if (cut.has(i)) return;
      const [u, v] = e.oc;
      if (!adj.has(u)) adj.set(u, []);
      if (!adj.has(v)) adj.set(v, []);
      adj.get(u).push({ to: v, ei: i });
      adj.get(v).push({ to: u, ei: i });
    });
    return adj;
  }
  function pathBetween(adj, srcSet, dstSet, banEdge) {
    const prev = new Map(), q = [];
    for (const s of srcSet) { prev.set(s, null); q.push(s); }
    while (q.length) {
      const u = q.shift();
      if (dstSet.has(u) && prev.get(u) !== null) {
        const edges = [];
        for (let x = u; prev.get(x); x = prev.get(x).from) edges.push(prev.get(x).ei);
        return edges;
      }
      for (const { to, ei } of (adj.get(u) || [])) {
        if (ei === banEdge || prev.has(to)) continue;
        prev.set(to, { ei, from: u });
        q.push(to);
      }
    }
    return null;
  }
  // shortest cycle through edge ei: path between its endpoints avoiding ei, plus ei
  function cycleThrough(adj, I, ei) {
    const [u, v] = I.edges[ei].oc;
    if (u === v) return [ei];
    const p = pathBetween(adj, new Set([u]), new Set([v]), ei);
    return p ? p.concat([ei]) : null;
  }

  function boundaryOrbits(I, cut, isMirror) {
    const s = new Set();
    I.edges.forEach((e, i) => { if (isMirror[i] || cut.has(i)) for (const o of e.oc) s.add(o); });
    return s;
  }

  function chooseCuts(D, I, ord) {
    const isMirror = I.edges.map(e => e.nodes.length === 1);
    const cut = new Set();
    const cones = I.orbits.map((o, i) => i).filter(i => ord[i].i > 1);

    // Phase 0 -- no boundary and no cone point to anchor to (torus `o`, Klein bottle `xx`).
    // Cut a closed loop: it has no ends, so it cannot leave a slit tip.
    if (!boundaryOrbits(I, cut, isMirror).size && !cones.length) {
      let best = null;
      for (let ei = 0; ei < I.edges.length; ei++) {
        const cyc = cycleThrough(orbitAdj(I, cut), I, ei);
        if (cyc && (!best || cyc.length < best.length)) best = cyc;
      }
      if (!best) return null;
      for (const e of best) cut.add(e);
    }

    // Phase 1 -- anchor every cone point. With no mirror anywhere (sphere orbifolds 2222, 632,
    // 442, 333) seed on a cone point and join the rest to it: the classical way to open a
    // sphere-with-cone-points into a disk.
    if (cones.length) {
      const onB = boundaryOrbits(I, cut, isMirror);
      if (!onB.size) {
        let seed = cones[0];
        for (const c of cones) if (ord[c].i > ord[seed].i) seed = c;
        onB.add(seed);
      }
      for (const src of cones) {
        if (onB.has(src)) continue;
        const p = pathBetween(orbitAdj(I, cut), new Set([src]), onB, -1);
        if (!p) return null;
        for (const e of p) { cut.add(e); for (const o of I.edges[e].oc) onB.add(o); }
      }
    }

    // Phase 2 -- topological cuts until the complex is genuinely a disk. Every candidate is an
    // ARC running boundary-to-boundary (possibly through interior orbits) or a CLOSED LOOP, so
    // no candidate can leave a slit tip. Accept one only if it measurably improves the topology,
    // which keeps the cut set minimal.
    const score = t => Math.abs(t.chi - 1) + (t.walks.length - 1);
    for (let round = 0; round < 10; round++) {
      const topo = topology(D, I, instancesFor(D, I, cut, isMirror));
      if (!topo) return null;
      if (topo.isDisk) return { cut, isMirror };
      const onB = boundaryOrbits(I, cut, isMirror);
      const adj = orbitAdj(I, cut);

      // multi-source BFS out from the boundary, so any orbit can be joined back to it
      const par = new Map();
      { const q = [];
        for (const o of onB) { par.set(o, null); q.push(o); }
        while (q.length) { const u = q.shift();
          for (const { to, ei } of (adj.get(u) || [])) if (!par.has(to)) { par.set(to, { ei, from: u }); q.push(to); } } }
      const toBoundary = o => { const es = [];
        for (let x = o; par.get(x); x = par.get(x).from) es.push(par.get(x).ei);
        return par.has(o) ? es : null; };

      let best = null, bestScore = score(topo);
      const consider = edges => {
        if (!edges || !edges.length) return;
        const trial = new Set(cut); for (const e of edges) trial.add(e);
        if (trial.size === cut.size) return;
        const tt = topology(D, I, instancesFor(D, I, trial, isMirror));
        if (!tt) return;
        const sc = score(tt);
        if (sc < bestScore) { bestScore = sc; best = edges; }
      };
      for (let ei = 0; ei < I.edges.length; ei++) {
        if (cut.has(ei) || isMirror[ei]) continue;
        const [u, v] = I.edges[ei].oc;
        if (onB.size) {
          const a = toBoundary(u), b = toBoundary(v);            // arc: boundary -> u -e- v -> boundary
          if (a && b) consider([...new Set([...a, ei, ...b])]);
        }
        consider(cycleThrough(adj, I, ei));                      // closed loop through ei
        if (bestScore === 0) break;
      }
      if (!best) return null;
      for (const e of best) cut.add(e);
    }
    const inst = instancesFor(D, I, cut, isMirror);
    const topo = topology(D, I, inst);
    return (topo && topo.isDisk) ? { cut, isMirror } : null;
  }

  function buildDisk(D, I) {
    const ord = DG.orbitOrders(D, I);
    const ci = chooseCuts(D, I, ord);
    if (!ci) return null;
    const inst = instancesFor(D, I, ci.cut, ci.isMirror);
    const topo = topology(D, I, inst);
    if (!topo || !topo.isDisk) return null;
    const walk = topo.walks[0];

    const sCount = {};
    for (const w of walk) sCount[w.orbit] = (sCount[w.orbit] || 0) + 1;
    const alphOf = new Map();
    for (const w of walk) {
      const o = ord[w.orbit];
      alphOf.set(w.corner, (2 * Math.PI) / (o.v * o.f * sCount[w.orbit]));
    }
    // A slit tip (i < 2) cannot be placed: dist = RAD/cos(beta/2) diverges.
    for (const w of walk) if ((2 * Math.PI) / alphOf.get(w.corner) < 2 - 1e-9) return null;

    return { ord, cut: ci.cut, isMirror: ci.isMirror, isB: inst.isB, walk, alphOf,
             nSides: inst.nSides, nCorners: inst.nCorners,
             sideOf: inst.sideOf, cornOf: inst.cornOf, sCount, chi: topo.chi };
  }

  function computeCoords(D, I, K) {
    const P = new Array(K.nCorners).fill(null);
    let sum = Math.PI / 2; const corners = [];
    for (const w of K.walk) {
      const alph = K.alphOf.get(w.corner);
      if ((2 * Math.PI) / alph > 2 + 1e-9) {
        const beta = Math.PI - alph;
        sum += 0.5 * beta;
        const dist = RAD / Math.cos(0.5 * beta);
        P[w.corner] = [dist * Math.cos(sum), dist * Math.sin(sum)];
        sum += 0.5 * beta;
        corners.push(w.corner);
      }
    }
    if (!corners.length) return null;
    const nw = K.walk.length;
    const start = K.walk.findIndex(w => P[w.corner] !== null);
    let run = [], prev = K.walk[start].corner;
    for (let k = 1; k <= nw; k++) {
      const w = K.walk[(start + k) % nw];
      if (P[w.corner] !== null && !run.includes(w.corner)) {
        for (let j = 0; j < run.length; j++) {
          const t = (j + 1) / (run.length + 1);
          P[run[j]] = [P[prev][0] + t * (P[w.corner][0] - P[prev][0]),
                       P[prev][1] + t * (P[w.corner][1] - P[prev][1])];
        }
        run = []; prev = w.corner;
      } else if (P[w.corner] === null) run.push(w.corner);
    }
    let cx = 0, cy = 0, c = 0;
    for (const p of P) if (p) { cx += p[0]; cy += p[1]; c++; }
    cx /= c; cy /= c;
    for (let i = 0; i < P.length; i++) if (!P[i]) P[i] = [cx, cy];

    const E = new Array(K.nSides).fill(null);
    for (let f = 1; f <= D.size; f++) for (let t = 0; t < 3; t++) {
      const si = K.sideOf(f, t), [a, b] = PAIR[t];
      const p = P[K.cornOf(f, a)], q = P[K.cornOf(f, b)];
      E[si] = [0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1])];
    }
    const N = [];
    for (let f = 1; f <= D.size; f++) {
      let x = 0, y = 0;
      for (let t = 0; t < 3; t++) { x += E[K.sideOf(f, t)][0] + P[K.cornOf(f, t)][0];
                                    y += E[K.sideOf(f, t)][1] + P[K.cornOf(f, t)][1]; }
      N.push([x / 6, y / 6]);
    }
    return { P, E, N, corners };
  }

  function relax(D, I, K, C, opts) {
    const EPS = (opts && opts.eps) || 1e-11, MAX = (opts && opts.maxIter) || 30000;
    const { P, E, N } = C;
    const pinned = new Set();
    for (const w of K.walk) if ((2 * Math.PI) / K.alphOf.get(w.corner) > 1 + 1e-9) pinned.add(w.corner);
    const cornInc = new Map(), sideInc = new Map();
    for (let f = 1; f <= D.size; f++) for (let t = 0; t < 3; t++) {
      const ck = K.cornOf(f, t), [a, b] = PAIR[t];
      if (!cornInc.has(ck)) cornInc.set(ck, []);
      cornInc.get(ck).push({ f, sides: [K.sideOf(f, a), K.sideOf(f, b)] });
      const sk = K.sideOf(f, t);
      if (!sideInc.has(sk)) sideInc.set(sk, { flags: [], corners: null });
      sideInc.get(sk).flags.push(f);
      sideInc.get(sk).corners = [K.cornOf(f, a), K.cornOf(f, b)];
    }
    let it = 0, moved = true;
    while (moved && it++ < MAX) {
      moved = false;
      for (let f = 1; f <= D.size; f++) {
        let x = 0, y = 0;
        for (let t = 0; t < 3; t++) { x += E[K.sideOf(f, t)][0] + P[K.cornOf(f, t)][0];
                                      y += E[K.sideOf(f, t)][1] + P[K.cornOf(f, t)][1]; }
        x /= 6; y /= 6;
        if (Math.abs(x - N[f - 1][0]) >= EPS || Math.abs(y - N[f - 1][1]) >= EPS) moved = true;
        N[f - 1] = [x, y];
      }
      for (const [k, s] of sideInc) {
        let x = P[s.corners[0]][0] + P[s.corners[1]][0], y = P[s.corners[0]][1] + P[s.corners[1]][1], d = 2;
        for (const f of s.flags) { x += N[f - 1][0]; y += N[f - 1][1]; d++; }
        x /= d; y /= d;
        if (Math.abs(x - E[k][0]) >= EPS || Math.abs(y - E[k][1]) >= EPS) moved = true;
        E[k] = [x, y];
      }
      for (const [k, list] of cornInc) {
        if (pinned.has(k)) continue;
        let x = 0, y = 0, d = 0;
        for (const c of list) { x += N[c.f - 1][0]; y += N[c.f - 1][1]; d++;
          for (const s of c.sides) { x += E[s][0]; y += E[s][1]; d++; } }
        x /= d; y /= d;
        if (Math.abs(x - P[k][0]) >= EPS || Math.abs(y - P[k][1]) >= EPS) moved = true;
        P[k] = [x, y];
      }
    }
    return { iterations: it, converged: !moved };
  }

  function realize(sym) {
    const D = DG.parseSymbol(sym), I = DG.incidence(D);
    const K = buildDisk(D, I);
    if (!K) return null;
    const C = computeCoords(D, I, K);
    if (!C) return null;
    const rel = relax(D, I, K, C);
    return { D, I, K, C, rel };
  }

  const api = { UF, instancesFor, topology, chooseCuts, buildDisk, computeCoords, relax, realize };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.DsymDisk = api;
})(typeof window !== 'undefined' ? window : this);
