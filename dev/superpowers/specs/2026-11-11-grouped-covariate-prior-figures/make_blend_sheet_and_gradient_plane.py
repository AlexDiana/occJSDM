"""Draws blend-sheet-and-gradient-plane.svg. Schematic: species' whole-response vectors, the blend sheet, leftovers and the gradient plane.

Made-up points, two panels. Writes an SVG. No data.
"""
import math, sys

W, H = 1240, 640
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
SURFACE = "#fcfcfb"
BLUE, BLUE_FILL = "#2a78d6", "#b7d3f6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"

def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def mul(a, k): return tuple(x * k for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def unit(a): n = math.sqrt(dot(a, a)); return mul(a, 1 / n)

def projector(cx, cy, s):
    # oblique view: x to the right, y receding up and to the right, z up
    def P(p):
        x, y, z = p
        return (cx + s * (x + 0.55 * y), cy - s * (z + 0.32 * y))
    return P

def poly(P, pts, **kw):
    d = " ".join(f"{P(p)[0]:.1f},{P(p)[1]:.1f}" for p in pts)
    a = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in kw.items())
    return f'<polygon points="{d}" {a}/>'

def line(P, a, b, **kw):
    (x1, y1), (x2, y2) = P(a), P(b)
    at = " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in kw.items())
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" {at}/>'

def text(x, y, s, size=15, fill=INK, anchor="start", weight="400", style=""):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return (f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" '
            f'font-weight="{weight}" {style}>{s}</text>')

def sheet(u, v, lo=-0.15, hi=1.35):
    return [add(mul(u, lo), mul(v, lo)), add(mul(u, hi), mul(v, lo)), add(mul(u, hi), mul(v, hi)), add(mul(u, lo), mul(v, hi))]

# The true plane of the species' whole responses (the habitat gradients), shared by both panels.
U = (1.0, 0.0, 0.18)
V = (0.0, 1.0, 0.42)
N = unit(cross(U, V))
coords = [(0.35, 1.05, 0.10), (0.75, 0.85, -0.12), (1.10, 0.55, 0.14), (0.95, 1.20, 0.08),
          (1.25, 1.00, -0.10), (0.55, 0.45, 0.12), (1.20, 0.25, -0.09)]
tips = [add(add(mul(U, a), mul(V, b)), mul(N, e)) for a, b, e in coords]

out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
       f'font-family="-apple-system, Helvetica, Arial, sans-serif">',
       f'<rect width="{W}" height="{H}" fill="{SURFACE}"/>',
       '<defs>'
       f'<marker id="head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{INK2}"/></marker>'
       '</defs>']

def axes(P):
    o = (0, 0, 0)
    r = [line(P, o, (1.75, 0, 0), stroke=MUTED, stroke_width=1),
         line(P, o, (0, 1.65, 0), stroke=MUTED, stroke_width=1),
         line(P, o, (0, 0, 1.25), stroke=MUTED, stroke_width=1)]
    for p, lab, dx in (((1.8, 0, 0), "column 1", 28), ((0, 1.78, 0), "column 2", 30), ((0, 0, 1.31), "column 3", 0)):
        x, y = P(p); r.append(text(x + dx, y + 4, lab, size=12, fill=MUTED, anchor="middle"))
    return r

def species_arrows(P):
    return [line(P, (0, 0, 0), t, stroke=INK2, stroke_width=2, marker_end="url(#head)") for t in tips]

# ---------- Panel A: small leftovers ----------
PA = projector(105, 430, 215)
out.append(text(40, 52, "A. Shared part carries the response (strong signal)", size=19, weight="600"))
out.append(text(40, 80, "The species' arrows lie close to the blend sheet,", size=14, fill=INK2))
out.append(text(40, 99, "so the gradient plane is the same sheet.", size=14, fill=INK2))
out += axes(PA)
out.append(poly(PA, sheet(U, V), fill=BLUE_FILL, fill_opacity="0.75", stroke=BLUE, stroke_width=2))
for (a, b, e), t in zip(coords, tips):
    s = add(mul(U, a), mul(V, b))
    out.append(line(PA, s, t, stroke=ORANGE, stroke_width=2.5))
out += species_arrows(PA)
out.append(poly(PA, sheet(U, V, -0.2, 1.42), fill="none", stroke=AQUA, stroke_width=2.5, stroke_dasharray="9 6"))
x, y = PA(add(mul(U, 1.38), mul(V, 1.38)))
out.append(text(x - 4, y - 12, "blend sheet = gradient plane", size=14, fill=INK, anchor="end"))

# ---------- Panel B: leftovers dominate ----------
PB = projector(735, 430, 215)
out.append(text(660, 52, "B. Leftovers carry the response (weak signal)", size=19, weight="600"))
out.append(text(660, 80, "The fitted blend sheet is small and points the wrong way.", size=14, fill=INK2))
out.append(text(660, 99, "The gradient plane still fits the whole arrows, so it finds the gradients.", size=14, fill=INK2))
out += axes(PB)
# the fitted (wrong) blend sheet: tilted, and the shared parts are short
U2 = unit((0.25, 0.0, 1.0)); V2 = unit((0.15, 1.0, -0.25))
out.append(poly(PB, sheet(U2, V2, -0.05, 0.6), fill=BLUE_FILL, fill_opacity="0.75", stroke=BLUE, stroke_width=2))
for i, t in enumerate(tips):
    if i not in (0, 2, 6):
        continue                      # take three species apart, leave the rest whole
    s = add(mul(U2, 0.5 * dot(t, U2)), mul(V2, 0.5 * dot(t, V2)))
    out.append(line(PB, (0, 0, 0), s, stroke=BLUE, stroke_width=3))
    out.append(line(PB, s, t, stroke=ORANGE, stroke_width=2.5))
out += species_arrows(PB)
out.append(poly(PB, sheet(U, V, -0.2, 1.42), fill="none", stroke=AQUA, stroke_width=2.5, stroke_dasharray="9 6"))
x, y = PB(add(mul(U, 1.38), mul(V, 1.38)))
out.append(text(x - 4, y - 12, "gradient plane", size=14, fill=INK, anchor="end"))
x, y = PB(add(mul(U2, 0.6), mul(V2, -0.05)))
out.append(text(x - 10, y - 8, "fitted blend sheet", size=14, fill=INK, anchor="end"))

# ---------- legend ----------
pos = [(40, 552), (660, 552), (40, 584), (660, 584)]
items = [("arrow", INK2, "a species' whole response (shared part + leftover)"),
         ("fill", BLUE, "blend sheet (shared part)"),
         ("line", ORANGE, "leftover"),
         ("dash", AQUA, "gradient plane (returnEmbeddingGradients)")]
for (x, ly), (kind, col, lab) in zip(pos, items):
    if kind == "arrow":
        out.append(f'<line x1="{x}" y1="{ly}" x2="{x+30}" y2="{ly}" stroke="{col}" stroke-width="2" marker-end="url(#head)"/>')
    elif kind == "fill":
        out.append(f'<rect x="{x}" y="{ly-8}" width="30" height="16" fill="{BLUE_FILL}" stroke="{col}" stroke-width="2"/>')
    elif kind == "line":
        out.append(f'<line x1="{x}" y1="{ly}" x2="{x+30}" y2="{ly}" stroke="{col}" stroke-width="2.5"/>')
    else:
        out.append(f'<line x1="{x}" y1="{ly}" x2="{x+30}" y2="{ly}" stroke="{col}" stroke-width="2.5" stroke-dasharray="9 6"/>')
    out.append(text(x + 40, ly + 5, lab, size=14, fill=INK2))
out.append(text(40, 628, "Schematic with made-up points. The embedding has 64 columns; 3 are drawn. A species' arrow is its 64 coefficients on the embedding.", size=12, fill=MUTED))
out.append("</svg>")
open(sys.argv[1], "w", encoding="utf-8").write("\n".join(out))
