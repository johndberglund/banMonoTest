# program9 — exploring impossible types

A standalone, GitHub-hosted (GitHub Pages) site for the 2-isohedral types that cannot be monohedral:
browse the 362 proofs, and build or check a containment proof yourself.

## What program4 already had

| App | Where | What it does |
|---|---|---|
| **netshape** (proof builder) | `program4/engine/netshape.html` (111 KB) | Draws a type's tiling from its Delaney symbol. You click points to build a *configuration* (points forced into one tile orbit) and a *region* (the boundary the other orbit cannot leave; right-click = forbidden point). "Does it fit?" tests whether any isometric copy of the configuration fits in the region; "Test across the deformation family" repeats it as the shape deforms. Also: derive an n-gon from a 2-point configuration by the rotation order at a centre, a shape editor (drag the fundamental domain), JSON proof record export, PNG/SVG export. |
| **Impossibles viewer** | `program4/netshape pix proofs/impossibles_viewer.html` + `impossibles_manifest.js` | Searchable sidebar of the 362 impossible types in 11 groups; each shows its proof picture, orbifold, wallpaper group, paper label and the reason text. |
| Supporting writeups | `RADIUS_ARGUMENT.html`, `netshape pix proofs/PERIMETER_ARGUMENT.md`, `FORCED_SYMMETRY_ARGUMENT.md`, `py/monohedral_impossible.md` | The arguments behind the perimeter, forced-symmetry and radius classes. |
| John's working notes | `impossible work.txt`, `unknowns that are impossible.txt`, `impossibles.txt` | Hand arguments and proof records for types still being worked on. |

netshape's script dependencies (all in `program4/engine/`, ~340 KB total):
`dsym_geom.js`, `dsym_disk.js`, `dsym_deform.js`, `_dev_isohedral/draw_dsym2.js`,
`dsym_types.js` (all 1270 types: symbol, k, n, orbifold, verdict), `tiler_curves.js`, `dsym_edit.js`.

Current counts in `dsym_types.js`: **362 impossible / 409 possible / 499 unknown.**

Proof groups in the manifest: perimeter 173, rigid-orbifold 114, rect3 25, forced-symmetry 14,
n-gon 13, rotate-to-square 10, rect-no-rotation 9, rotate-to-square-trapezoid 6,
half-size-region 5, diameter 2, rotation-forces-shape 2.

Proof pictures: 362 PNGs, 34 MB already copied into `program7/program8/proofs/`
(55 MB in program4 with backups and SVGs).

## Pain points to fix

- netshape is a developer tool: dark dense panel, 1270-type dropdown, hidden modes
  (middle-drag pan, right-click forbidden points), jargon ("realization", "skeleton").
- The viewer and the builder are separate; you can't open a proof in the builder, and the
  builder can't show the published proof for a type.
- Proof records are copied to the clipboard as JSON; there's no way to load one back in.
- Files load from several folders by relative path, so neither app works when copied anywhere else.

## Proposed site

```
program9/
  index.html        landing: what an impossible type is, counts, links into the two tools
  proofs.html       viewer: browse/search by group, orbifold, id; picture + reason + link
                    "open in builder"
  builder.html      netshape, reorganised as a step-by-step workflow (below)
  learn.html        short illustrated guide to each proof class (from the writeups above)
  js/               the 7 engine files, copied in unchanged first, then trimmed
  data/             types.js (1270 types), proofs.js (manifest + saved proof records)
  proofs/           the 362 PNGs
  README.md         GitHub front page
```

### Builder workflow

1. **Pick a type**: search by id, filter by verdict and orbifold, thumbnails.
2. **Configuration**: click the tile, then the points it must contain. A help strip says what
   to click next.
3. **Region**: click the other tile, then its boundary points. A toggle button for
   "forbidden point" (keeps right-click as a shortcut).
4. **Test**: fit / family test, with a plain-language verdict and the best placement drawn.
5. **Save**: download or load a `.json` proof record, and export PNG/SVG.

Also: pan with left-drag on empty space, zoom buttons, touch support, undo everywhere, and a
light/dark theme.

### GitHub

- Its own repo (e.g. `johndberglund/impossible-tilings`), served by GitHub Pages from `main`.
  It's all static files, so there's no build step.
- 34 MB of PNGs is fine for a repo. We could shrink them to WebP or lower resolution later if
  load time matters.
- `gh` is not installed on this machine; git is, with user `johndberglund`. Either install
  `gh` (`winget install GitHub.cli`) and run `gh auth login`, or create the empty repo on
  github.com and I'll add it as a remote and push.

## Build order

1. Copy the engine files and data into program9 and get netshape running unchanged from here
   (to prove it is self-contained).
2. Port the viewer onto `data/proofs.js` and `proofs/` and cross-link it with the builder.
3. Restyle the builder and add the step workflow, load/save, and the pan/touch fixes.
4. Landing and learn pages, then README.
5. `git init`, first commit, push, enable Pages.

## Decisions (John, 2026-10-21)

- The builder is public.
- Repo: `johndberglund/banMonoTest`, public.

## Keeping program8 and program9 in sync (proposal)

Neither site is the source of truth. Both are generated from the research side in program7, so
each generator only has to read the other site's results.

```
program8 finds a tile ──► program7 store ──┐
                                           ├─► app_tools/gen_verdicts.js ──► program8/data/types.js
program9 proof record ──► program7 proofs ─┘                             └─► program9/data/types.js
```

- **New impossible type (found in program9):** the builder saves a proof record (`.json`). A visitor
  attaches it to a GitHub issue; you drop your own into `interact/`. `import_proof.js` re-runs the
  fit and family test in Node and, if it holds, adds the record to `program7/data/proofs/`. That
  type then turns red in program8, with the proof picture.
- **New possible type (found in program8):** it goes into the store as it does now. On the next
  regeneration, program9 shows the type as possible, removes it from the builder's "unknown" list,
  and links to the tile in program8.
- **Safety check:** the generator refuses to run if a type has both a tile and an impossibility
  proof, as `build_ledger.py` does now.
- **Cross-links:** both sites use `#<id>` URLs, so each type page links to the other site.
- **Publishing:** after an import, run one command that regenerates both, then commit and push.
