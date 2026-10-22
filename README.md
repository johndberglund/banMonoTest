# banMonoTest

Tools for exploring the **2-isohedral tiling types that cannot be monohedral**: types where no
single tile shape can fill both tile orbits.

Of the 1,270 Euclidean 2-isohedral types (Delgado, Huson and Zamorzaeva), 362 are proven
impossible. Each proof shows that what one tile is forced to contain cannot fit inside the region
the other tile is confined to, so the two tiles cannot be congruent.

## The pages

- **[Proof builder](builder.html)**: pick a type and see its tiling, with the two tile classes
  drawn on their own (one of each, or in their true proportion). Click points that one tile
  must contain and a region the other tile cannot leave, then test whether the points fit,
  across the whole deformation family. Save the result as a proof record, PNG or SVG.
- **[Proof gallery](proofs.html)**: all 362 impossibility proofs, grouped by argument, with a
  picture and the reason for each. Link to a type with `proofs.html#<id>`.

Everything runs in the browser from static files. Open `builder.html` directly, or serve the
folder with any static web server.

## Found a new proof?

Use **Copy proof record (JSON)** in the builder, then
[send it in with the proof form](https://github.com/johndberglund/banMonoTest/issues/new?template=new-proof.yml).
Each one is re-checked before it is added.

## License

Code: [GNU GPL version 3](LICENSE) -- the drawing engine follows Tegula (Huson, Zeller, Westphal; GPL-3), see [NOTICE.md](NOTICE.md). Proof pictures, data and write-ups: [CC BY 4.0](LICENSE-CONTENT.md).
