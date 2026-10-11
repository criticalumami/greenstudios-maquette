"""Generate clay toolpaths for the Greenstudios maquette, part by part.

Run:  pip install -r requirements.txt   then   python gen.py
Reads geo.json (the building, real metres) and rewrites ../toolpaths.js.
Main settings are the constants BEAD, H (layer height) and RES below, in real metres.

Rules:
- every wall becomes one bead, 0.25 m real (2.5 mm at 1:100), on the wall centreline
- nothing overhangs: a point prints at layer k only if it is solid in every layer below
  (lintels over doors/windows disappear; openings run up to the slab above)
- each part = slab printed first, walls printed on top of it
- walls: every wall printed once (Euler route), the nozzle lifts only between strokes;
  stroke ends that land on a wall already printed are trimmed half a bead so beads don't double up
"""
import json, math, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, trimesh, networkx as nx
from shapely.geometry import Polygon, MultiPolygon, LineString, Point, box
from shapely.ops import unary_union
from shapely import affinity
from skimage.morphology import skeletonize
from skimage.draw import polygon as rpoly

D = json.load(open(os.path.join(HERE, 'geo.json')))
BEAD = 0.25        # m real, = 2.5 mm model
H = 0.10           # m real layer height = 1.0 mm model
RES = 0.02         # raster resolution for centrelines
S = 10.0           # model mm per real metre (1:100)

def mesh(k):
    v = np.array(D[k]['v']).reshape(-1, 3); f = np.array(D[k]['f']).reshape(-1, 3)
    return trimesh.Trimesh(v, f, process=True)

def section(k, z):
    m = mesh(k)
    if z <= m.bounds[0][2] + 1e-6 or z >= m.bounds[1][2] - 1e-6:
        return Polygon()
    s = m.section(plane_origin=[0, 0, z], plane_normal=[0, 0, 1])
    if s is None:
        return Polygon()
    p2, T = s.to_2D(to_2D=np.eye(4))
    polys = [Polygon(p.exterior.coords, [h.coords for h in p.interiors]).buffer(0) for p in p2.polygons_full]
    return unary_union(polys) if polys else Polygon()

def ring(poly, t=BEAD):
    return poly.difference(poly.buffer(-t, join_style=2))

# ---------- centreline graph from a region ----------
def skeleton_graph(region):
    if region.is_empty:
        return nx.Graph()
    minx, miny, maxx, maxy = region.bounds
    minx -= 0.2; miny -= 0.2
    W = int((maxx - minx) / RES) + 12; Hh = int((maxy - miny) / RES) + 12
    img = np.zeros((Hh, W), bool)
    geoms = region.geoms if hasattr(region, 'geoms') else [region]
    for g in geoms:
        # each piece in its own mask, so one piece's hole never erases another piece
        one = np.zeros_like(img)
        ext = np.array(g.exterior.coords)
        rr, cc = rpoly((ext[:, 1] - miny) / RES, (ext[:, 0] - minx) / RES, img.shape)
        one[rr, cc] = True
        for h in g.interiors:
            a = np.array(h.coords)
            rr, cc = rpoly((a[:, 1] - miny) / RES, (a[:, 0] - minx) / RES, img.shape)
            one[rr, cc] = False
        img |= one
    sk = skeletonize(img)
    pts = set(zip(*np.nonzero(sk)))
    nb = lambda p: [(p[0] + dy, p[1] + dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dy or dx) and (p[0] + dy, p[1] + dx) in pts]
    deg = {p: len(nb(p)) for p in pts}
    nodepx = {p for p in pts if deg[p] != 2}
    # cluster adjacent node pixels
    cl = {}; cid = 0
    for p in nodepx:
        if p in cl: continue
        st = [p]; cl[p] = cid
        while st:
            q = st.pop()
            for r in nb(q):
                if r in nodepx and r not in cl:
                    cl[r] = cid; st.append(r)
        cid += 1
    cent = {}
    for p, c in cl.items():
        cent.setdefault(c, []).append(p)
    w = lambda p: (minx + p[1] * RES, miny + p[0] * RES)
    G = nx.MultiGraph()
    for c, ps in cent.items():
        a = np.array(ps).mean(0); G.add_node(c, xy=w(a))
    seen = set()
    def trace(start_c, p0, p1):
        path = [p0, p1]
        prev, cur = p0, p1
        while cur not in cl:
            nx_ = [r for r in nb(cur) if r != prev and r not in path[-3:]]
            if not nx_: return None
            prev, cur = cur, nx_[0]; path.append(cur)
        return path
    for p in nodepx:
        for r in nb(p):
            if r in nodepx: continue
            key = (p, r)
            if key in seen: continue
            path = trace(cl[p], p, r)
            if path is None: continue
            seen.add((path[-1], path[-2]))
            seen.add(key)
            a, b = cl[p], cl[path[-1]]
            coords = [G.nodes[a]['xy']] + [w(q) for q in path[1:-1]] + [G.nodes[b]['xy']]
            G.add_edge(a, b, line=coords)
    # pure loops (no node pixels)
    rest = pts - nodepx - {q for e in G.edges(data=True) for q in []}
    covered = set()
    for _, _, dd in G.edges(data=True):
        pass
    # find loop components not touched
    touched = set()
    for p in nodepx:
        touched.add(p)
    comp_seen = set()
    for p in pts:
        if p in nodepx or p in comp_seen: continue
        # BFS component
        st = [p]; comp = {p}
        while st:
            q = st.pop()
            for r in nb(q):
                if r not in comp: comp.add(r); st.append(r)
        comp_seen |= comp
        if comp & nodepx: continue
        # ordered loop
        order = [p]; prev = None; cur = p
        while True:
            n = [r for r in nb(cur) if r != prev and r not in order[-2:]]
            if not n: break
            if n[0] == p and len(order) > 3: break
            prev, cur = cur, n[0]
            if cur in order: break
            order.append(cur)
        c = cid; cid += 1
        G.add_node(c, xy=w(order[0]))
        G.add_edge(c, c, line=[w(q) for q in order] + [w(order[0])])
    # prune spurs
    for _ in range(3):
        rm = []
        for a, b, k, dd in G.edges(keys=True, data=True):
            L = LineString(dd['line']).length
            da, db = G.degree(a), G.degree(b)
            # only trim side-spurs off a junction; a short free-standing piece (pier, stub) is kept
            if a != b and min(da, db) == 1 and max(da, db) >= 3 and L < 0.20:
                rm.append((a, b, k))
        if not rm: break
        for e in rm:
            if G.has_edge(*e): G.remove_edge(*e)
        G.remove_nodes_from([n for n in list(G.nodes) if G.degree(n) == 0])
    # merge deg-2 chains
    changed = True
    while changed:
        changed = False
        for n in list(G.nodes):
            if G.degree(n) == 2 and G.number_of_edges(n, n) == 0:
                es = list(G.edges(n, keys=True, data=True))
                if len(es) != 2: continue
                (_, a, k1, d1), (_, b, k2, d2) = es
                l1 = d1['line'] if tuple(d1['line'][-1]) == tuple(G.nodes[n]['xy']) else d1['line'][::-1]
                l2 = d2['line'] if tuple(d2['line'][0]) == tuple(G.nodes[n]['xy']) else d2['line'][::-1]
                G.remove_node(n)
                G.add_edge(a, b, line=l1 + l2[1:])
                changed = True
                break
    # any piece of wall left without a centreline (tiny piers) gets one along its long axis
    nid = max(list(G.nodes) + [0]) + 1
    for g in geoms:
        gb = g.buffer(0.02)
        if any(gb.intersects(LineString(dd['line'])) for _, _, dd in G.edges(data=True)):
            continue
        r = g.minimum_rotated_rectangle
        c = list(r.exterior.coords)[:4]
        e1 = (c[0], c[1]); e2 = (c[1], c[2])
        long_, short_ = (e1, e2) if math.dist(*e1) >= math.dist(*e2) else (e2, e1)
        mx = lambda p, q: ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
        # centreline joins the midpoints of the two short sides, pulled in half a bead
        sa = mx(*short_); sb = mx(c[(c.index(short_[1]) + 1) % 4], c[(c.index(short_[1]) + 2) % 4])
        L = math.dist(sa, sb); pull = min(BEAD / 2, max(0, L / 2 - 0.03))
        ux, uy = (sb[0] - sa[0]) / (L or 1), (sb[1] - sa[1]) / (L or 1)
        p0 = (sa[0] + ux * pull, sa[1] + uy * pull); p1 = (sb[0] - ux * pull, sb[1] - uy * pull)
        G.add_node(nid, xy=p0); G.add_node(nid + 1, xy=p1); G.add_edge(nid, nid + 1, line=[p0, p1]); nid += 2
    for a, b, k, dd in G.edges(keys=True, data=True):
        ls = LineString(dd['line']).simplify(0.04)
        dd['line'] = [tuple(c) for c in ls.coords]
    return G

def graph_from_lines(lines):
    G = nx.MultiGraph()
    k = lambda p: (round(p[0], 3), round(p[1], 3))
    for l in lines:
        c = list(l.coords)
        a, b = k(c[0]), k(c[-1])
        G.add_node(a, xy=a); G.add_node(b, xy=b)
        G.add_edge(a, b, line=[tuple(p) for p in c])
    return G

# ---------- route: each edge once, few nozzle lifts ----------
def route(G, start_xy):
    strokes = []
    comps = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    cur = np.array(start_xy)
    while comps:
        comps.sort(key=lambda C: min(np.hypot(*(np.array(C.nodes[n]['xy']) - cur)) for n in C.nodes))
        C = comps.pop(0)
        H2 = nx.MultiGraph(C)
        odd = [n for n in H2.nodes if H2.degree(n) % 2]
        virt = []
        if odd:
            K = nx.Graph()
            for i, a in enumerate(odd):
                for b in odd[i + 1:]:
                    K.add_edge(a, b, weight=-math.dist(H2.nodes[a]['xy'], H2.nodes[b]['xy']))
            M = nx.max_weight_matching(K, maxcardinality=True)
            for a, b in M:
                H2.add_edge(a, b, virtual=True); virt.append((a, b))
        st = min(H2.nodes, key=lambda n: math.dist(H2.nodes[n]['xy'], cur))
        if virt:
            # start at an odd node near the nozzle so the dropped travel is the start
            st = min([n for v in virt for n in v], key=lambda n: math.dist(H2.nodes[n]['xy'], cur))
        circ = list(nx.eulerian_circuit(H2, source=st, keys=True))
        # rotate so the circuit begins right after a virtual edge
        if virt:
            idx = [i for i, (a, b, k) in enumerate(circ) if H2.edges[a, b, k].get('virtual')]
            best = min(idx, key=lambda i: math.dist(H2.nodes[circ[(i + 1) % len(circ)][0]]['xy'], cur))
            circ = circ[best + 1:] + circ[:best + 1]
        stroke = []
        for a, b, k in circ:
            dd = H2.edges[a, b, k]
            if dd.get('virtual'):
                if stroke: strokes.append(stroke)
                stroke = []; continue
            line = dd['line']
            if math.dist(line[0], H2.nodes[a]['xy']) > math.dist(line[-1], H2.nodes[a]['xy']):
                line = line[::-1]
            seg = [tuple(H2.nodes[a]['xy'])] + list(line[1:-1]) + [tuple(H2.nodes[b]['xy'])]
            if stroke: stroke += seg[1:]
            else: stroke = seg
        if stroke: strokes.append(stroke)
        cur = np.array(strokes[-1][-1])
    return strokes

def trim_ends(strokes):
    """Trim stroke ends that butt into a wall already printed (half a bead), count overlaps."""
    from collections import defaultdict
    key = lambda p: (round(p[0], 3), round(p[1], 3))
    passthru = defaultdict(int)
    for s in strokes:
        for p in s[1:-1]: passthru[key(p)] += 1
    owned = set(k for k, c in passthru.items() if c)
    crossings = sum(1 for c in passthru.values() if c >= 2)
    out = []
    for s in strokes:
        s = list(s)
        for end in (0, -1):
            k = key(s[end])
            if k in owned:
                ls = LineString(s if end == -1 else s[::-1])
                if ls.length > BEAD:
                    ls = LineString(s if end == -1 else s[::-1])
                    cut = ls.interpolate(ls.length - BEAD / 2)
                    coords = list(ls.coords)
                    # rebuild up to cut
                    acc = 0; new = [coords[0]]
                    for i in range(1, len(coords)):
                        d = math.dist(coords[i - 1], coords[i])
                        if acc + d >= ls.length - BEAD / 2:
                            new.append((cut.x, cut.y)); break
                        new.append(coords[i]); acc += d
                    s = new if end == -1 else new[::-1]
            else:
                owned.add(k)
        out.append(s)
    return out, crossings

# ---------- slab / plate infill ----------
def slab_paths(poly, angle):
    inset = poly.buffer(-BEAD / 2, join_style=2)
    strokes = []
    peri = inset if not hasattr(inset, 'geoms') else max(inset.geoms, key=lambda g: g.area)
    strokes.append(list(peri.exterior.coords))
    for h in peri.interiors: strokes.append(list(h.coords))
    fill = poly.buffer(-BEAD, join_style=2)
    if fill.is_empty: return strokes
    rot = affinity.rotate(fill, -angle, origin=(0, 0))
    minx, miny, maxx, maxy = rot.bounds
    segs = []
    y = miny + BEAD / 2; row = 0
    while y < maxy:
        ln = LineString([(minx - 1, y), (maxx + 1, y)]).intersection(rot)
        parts = [ln] if ln.geom_type == 'LineString' else list(getattr(ln, 'geoms', []))
        parts = [p for p in parts if p.length > 0.05]
        parts.sort(key=lambda p: p.coords[0][0])
        for p in parts:
            c = list(p.coords)
            segs.append((row, c if row % 2 == 0 else c[::-1]))
        y += BEAD; row += 1
    # greedy chaining: next segment = nearest free end; short link inside the plate = keep extruding
    segs = [[tuple(affinity.rotate(Point(p), angle, origin=(0, 0)).coords[0]) for p in c] for _, c in segs]
    free = segs[:]
    cur = free.pop(0)
    while free:
        e = cur[-1]
        best = min(range(len(free)), key=lambda i: min(math.dist(e, free[i][0]), math.dist(e, free[i][-1])))
        c = free.pop(best)
        if math.dist(e, c[-1]) < math.dist(e, c[0]): c = c[::-1]
        link = LineString([e, c[0]])
        if link.length < BEAD * 1.6 and poly.buffer(1e-3).contains(link):
            cur += c
        else:
            strokes.append(cur); cur = c
    strokes.append(cur)
    return strokes

# ---------- parts ----------
PARTS = [
    dict(id='gf', name='Ground floor', slab=None, plate=('F1_slab', 0.0, 0.2),
         walls=['GF_walls', 'GF_lintel_ring'] + [f'GF_column_{i}' for i in range(1, 7)],
         interior=['GF_interior_walls'], stair='GF_stair_block', top=2.75,
         orig=['GF_walls', 'GF_lintel_ring', 'GF_interior_walls', 'GF_stair_block'] + [f'GF_column_{i}' for i in range(1, 7)]),
]
for i, (n, zb, zt) in enumerate([(1, 3.25, 6.05), (2, 6.55, 9.35), (3, 9.85, 12.65)]):
    PARTS.append(dict(id=f'f{n}', name=f'Floor {n}', plate=(f'F{n}_slab', zb - 0.5, 0.5),
                      walls=[f'F{n}_walls', f'F{n}_lintel_ring'], interior=[f'F{n}_interior_walls'],
                      stair=f'F{n}_stair_block', top=zt,
                      orig=[f'F{n}_slab', f'F{n}_walls', f'F{n}_lintel_ring', f'F{n}_interior_walls', f'F{n}_stair_block']))
PARTS.append(dict(id='f4', name='Floor 4 (penthouse)', plate=('F4_slab', 12.65, 0.5),
                  walls=['F4_penthouse_walls', 'F4_lintel_ring', 'F4_terrace_parapet'], interior=['F4_interior_walls'],
                  stair='F4_stair_block', top=16.1,
                  orig=['F4_slab', 'F4_penthouse_walls', 'F4_lintel_ring', 'F4_terrace_parapet', 'F4_interior_walls', 'F4_stair_block']))
PARTS.append(dict(id='roof', name='Roof', plate=('Roof_slab', 16.1, 0.45), walls=['Roof_parapet'], interior=[],
                  stair=None, top=17.6, orig=['Roof_slab', 'Roof_parapet']))

def plan_at(P, z, thick):
    regs = []
    for k in P['walls']:
        s = section(k, z)
        if not s.is_empty: regs.append(s.buffer(0.01))
    for k in P['interior']:
        s = section(k, z)
        if not s.is_empty:
            if thick:
                # widen each partition about its centreline to a full bead
                regs.append(s.buffer(max(0, (BEAD - 0.09) / 2), join_style=2))
            else:
                regs.append(s)
    if P['stair']:
        s = section(P['stair'], z)
        if not s.is_empty: regs.append(ring(s))
    return unary_union(regs) if regs else Polygon()

def fmt(strokes, z):
    return [[round(c * S, 2) for p in s for c in (p[0], p[1])] for s in strokes]

out = {'bead': BEAD * S, 'layer': H * S, 'parts': []}
prev_end = (0, 0)
for P in PARTS:
    pk, pz, pt = P['plate']
    slabpoly = section(pk, (D and mesh(pk).bounds[0][2] + mesh(pk).bounds[1][2]) / 2)
    if P['id'] == 'roof':
        # roof slab is printed upside-down? no: printed flat, parapet on top
        pass
    base_z = pz if P['id'] != 'gf' else 0.0
    layers = []
    geoms = []
    n_sl = max(1, round(pt / H)); hs = pt / n_sl
    for i in range(n_sl):
        st = slab_paths(slabpoly, 0 if i % 2 == 0 else 90)
        if i % 2: st = [s[::-1] for s in st[::-1]]
        geoms.append({'s': fmt(st, 0), 'kind': 'slab'})
        layers.append({'z': round((i + 1) * hs * S, 3), 'g': len(geoms) - 1, 'kind': 'slab'})
    wall_bottom = base_z + pt if P['id'] != 'gf' else 0.0
    if P['id'] == 'gf': wall_bottom = 0.0
    nW = int((P['top'] - (pz + pt if P['id'] != 'gf' else 0.0)) / H + 1e-6)
    zb = pz + pt if P['id'] != 'gf' else 0.0
    cum = None; last_key = None; last_g = None; prev_lines = None
    stats = {'crossings': 0}
    for i in range(nW):
        z = zb + (i + 0.5) * H
        r = plan_at(P, z, True)
        cum = r if cum is None else cum.intersection(r)
        cum = cum.buffer(0)
        if hasattr(cum, 'geoms'):
            cum = unary_union([g for g in cum.geoms if g.area > 0.03])
        key = round(cum.area, 4), len(getattr(cum, 'geoms', [cum]))
        if key != last_key:
            Gs = skeleton_graph(cum)
            cur_lines = [LineString(dd['line']) for _, _, dd in Gs.edges(data=True)]
            if prev_lines:
                # reuse the centreline of the layer below, trimmed to where the wall still stands
                reg = cum.buffer(0.01)
                lines = []
                for l in prev_lines:
                    x = l.intersection(reg)
                    for pc in ([x] if x.geom_type == 'LineString' else list(getattr(x, 'geoms', []))):
                        if pc.geom_type == 'LineString' and pc.length > 0.03: lines.append(pc)
                # anything the new plan has that the old centreline doesn't reach
                cover = unary_union(lines).buffer(BEAD * 0.9) if lines else Polygon()
                for l in cur_lines:
                    x = l.difference(cover)
                    for pc in ([x] if x.geom_type == 'LineString' else list(getattr(x, 'geoms', []))):
                        if pc.geom_type == 'LineString' and pc.length > 0.12: lines.append(pc)
            else:
                lines = cur_lines
            prev_lines = lines
            G = graph_from_lines(lines)
            st = route(G, prev_end)
            st, cr = trim_ends(st)
            stats['crossings'] = max(stats['crossings'], cr)
            geoms.append({'s': fmt(st, 0), 'kind': 'wall'})
            last_g = len(geoms) - 1; last_key = key
        layers.append({'z': round((pt + (i + 1) * H) * S if P['id'] != 'gf' else (pt + (i + 1) * H) * S, 3), 'g': last_g, 'kind': 'wall', 'rev': i % 2})
    # ground plate: walls start on top of the 0.2 m plate
    out['parts'].append({'id': P['id'], 'name': P['name'], 'z0': round((pz if P['id'] != 'gf' else -pt) * S, 2),
                         'slab': pt * S, 'wallTop': P['top'], 'geoms': geoms, 'layers': layers,
                         'orig': P['orig'], 'crossings': stats['crossings']})
    print(P['name'], 'layers', len(layers), 'unique wall geoms', sum(1 for g in geoms if g['kind'] == 'wall'), file=sys.stderr)

# green wall printed flat on its back: section in the facade plane, as a plate
gm = mesh('Greenwall_screen')
s = gm.section(plane_origin=[0, -7.3, 0], plane_normal=[0, 1, 0])
p2, T = s.to_2D(to_2D=np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]))
gpoly = unary_union([Polygon(p.exterior.coords, [h.coords for h in p.interiors]).buffer(0) for p in p2.polygons_full])
# to_2D may flip; normalise to x in [-8,8], y in [0,14.25]
mnx, mny, mxx, mxy = gpoly.bounds
gpoly = affinity.translate(gpoly, 0, -mny)
geoms = []; layers = []
n = round(0.35 / H); hs = 0.35 / n
for i in range(n):
    st = slab_paths(gpoly, 0 if i % 2 == 0 else 90)
    if i % 2: st = [x[::-1] for x in st[::-1]]
    geoms.append({'s': fmt(st, 0), 'kind': 'slab'})
    layers.append({'z': round((i + 1) * hs * S, 3), 'g': i, 'kind': 'slab'})
out['parts'].append({'id': 'gw', 'name': 'Green wall (printed flat)', 'z0': 0, 'slab': 3.5, 'geoms': geoms,
                     'layers': layers, 'orig': ['Greenwall_screen'], 'crossings': 0, 'flat': True})
# original meshes for comparison
used = set(k for P in out['parts'] for k in P['orig'])
out['orig'] = {k: {'v': [round(x * S, 2) for x in D[k]['v']], 'f': D[k]['f'], 'm': D[k]['m']} for k in used}
txt = json.dumps(out, separators=(',', ':'))
with open(os.path.join(HERE, '..', 'toolpaths.js'), 'w') as fh:
    fh.write('// Clay toolpaths for the Clay toolpaths tab. Generated by tools/gen.py — edit the script, not this file.\n'
             '// Units: mm at 1:100. parts[].geoms[].s = strokes, each a flat [x,y,x,y,...] list; '
             'parts[].layers[] = {z (top of layer, mm), g (geom index), kind, rev}.\n')
    fh.write('const DATA_TP = ' + txt + ';\n')
print('bytes', len(json.dumps(out, separators=(',', ':'))), file=sys.stderr)
