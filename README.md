# Greenstudios Maquette

An interactive 3D viewer for a physical scale model (1:100) of the Greenstudios building, plot 2869, Pafos (adp architecture, planning set). It has two tabs:

- **Model**: the building, floor by floor, with a green-wall pull slider, floor explode, section cut, see-through outer walls, dimensions, a ceramic preview and a wall-thickness table at 1:80 / 1:100 / 1:150 / 1:200.
- **Clay toolpaths**: how each part would be printed in clay, layer by layer, with a red nozzle point you can move, a plan of the current layer, and a before/after comparison of wall thickness.

## Open it

Double-click `index.html`. It works offline in any recent browser (Chrome, Edge, Firefox, Safari). Fonts come from Google Fonts when online; without internet it falls back to system fonts.

To put it online, drag the whole `Greenstudios_maquette` folder onto Netlify Drop (app.netlify.com/drop) or any static host. `index.html` must stay at the top of the folder.

## Files

| File | What it is | Edit it? |
|---|---|---|
| `index.html` | The whole page: layout, styles and all the code for both tabs. | Yes, this is where changes go. |
| `model-data.js` | The building geometry (`GEO`) and the plot outline (`BASE`), in real metres. Exported from the Blender model. | Only to replace the building. |
| `toolpaths.js` | The clay toolpaths (`DATA_TP`), in mm at 1:100. | No. Regenerate it with `tools/gen.py`. |
| `lib/three.min.js`, `lib/OrbitControls.js` | three.js r128 (MIT licence, see `lib/three-LICENSE.txt`). | No. |
| `tools/gen.py` | Python script that computes the clay toolpaths from the geometry. | Yes, to change the print rules. |
| `tools/geo.json` | Same geometry as `model-data.js`, as plain JSON for the script. | Keep it in sync with `model-data.js`. |

## Working on it with an AI tool

- Give the tool `index.html` and this README. Leave out `model-data.js` and `toolpaths.js`: they are large data files and the tool does not need to read them to change the page.
- In `index.html` the first `<script>` block is the Model tab, the second is the Clay toolpaths tab (wrapped in `(() => { ... })()` so the two never clash). Sections are marked with `// ----------` comment lines.
- Colours and fonts are CSS variables at the top of the `<style>` block (light and dark theme).
- After a change, reload the page in the browser to see it.

## Units and coordinates

- Real building: metres. In `GEO`, x = east, y = north, z = up. Ground floor at 0.00 m.
- Clay toolpaths: millimetres at the final 1:100 size (1 m real = 10 mm).
- Floor levels (finished floor): ground 0.00, floor 1 3.25, floor 2 6.55, floor 3 9.85, floor 4 (penthouse) 13.15, roof slab top 16.55, parapet top 17.60 m. Slabs 0.50 m.
- Model size at 1:100: about 17.8 × 13.4 × 17.6 cm including the green wall.

## How the clay toolpaths are made

`tools/gen.py` applies these rules:

- Every wall, interior partitions included, becomes one bead 2.5 mm wide (250 mm real) on the wall centreline.
- Nothing overhangs: a point prints only if it is solid in every layer below. Lintels over doors and windows disappear, so openings run up to the slab of the floor above.
- Each part prints its slab first (layers turned 90° each time), then the walls on top. Layer height 1.0 mm.
- The stair and lift block prints as a hollow 2.5 mm shell.
- Each wall is printed once. A stroke that ends against a wall already printed stops half a bead short.
- Each layer reuses the centreline of the layer below, so beads stack.
- The green wall prints flat on its back as a 3.5 mm plate with its openings left as holes.
- Not included: balconies, AC louvres. No G-code export yet.
- Paths are at final size. For clay, scale up by the measured shrinkage (about 112% for 11%).

To regenerate after changing the rules:

```
cd tools
pip install -r requirements.txt
python gen.py
```

This rewrites `toolpaths.js`. Reload `index.html` to see the result.

## Source

Geometry from the adp architecture planning drawings (job 001966), modelled in Blender. The model is an envelope study: no basement, no fences, the stair is a solid block, some interior dimensions are approximate.
