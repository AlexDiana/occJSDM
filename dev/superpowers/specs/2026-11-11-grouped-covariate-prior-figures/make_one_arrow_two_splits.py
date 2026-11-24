"""Draws one-arrow-two-splits.svg. Schematic: one species' arrow, split two ways by two chains.

The data pin down the sum (the arrow), not the split between the shared part
and the leftover. Made-up points. Writes an SVG. No data.
"""
import math, sys

W, H = 1240, 640
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
SURFACE = "#fcfcfb"
BLUE, BLUE_FILL = "#2a78d6", "#b7d3f6"
ORANGE = "#eb6834"

def add(a, b): return tuple(x + y for x, y in zip(a, b))
def mul(a, k): return tuple(x * k for x in a)

def projector(cx, cy, s):
    def P(p):
        x, y, z = p
        return (cx + s * (x + 0.55 * y), cy - s * (z + 0.32 * y))
    return P

def attrs(kw): return " ".join(f'{k.replace("_", "-")}="{v}"' for k, v in kw.items())
def poly(P, pts, **kw):
    return f'<polygon points="{" ".join(f"{P(p)[0]:.1f},{P(p)[1]:.1f}" for p in pts)}" {attrs(kw)}/>'
def line(P, a, b, **kw):
    (x1, y1), (x2, y2) = P(a), P(b)
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" {attrs(kw)}/>'
def text(x, y, s, size=15, fill=INK, anchor="start", weight="400"):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-weight="{weight}">{s}</text>'

U = (1.0, 0.0, 0.18)
V = (0.0, 1.0, 0.42)
def on_sheet(a, b): return add(mul(U, a), mul(V, b))
SHEET = [on_sheet(-0.1, -0.1), on_sheet(1.4, -0.1), on_sheet(1.4, 1.25), on_sheet(-0.1, 1.25)]
TIP = add(on_sheet(1.05, 0.80), (0.0, -0.05, 0.62))       # the species' whole response
SPLITS = [("Chain 1", on_sheet(1.0, 0.75), "large shared part, small leftover"),
          ("Chain 2", on_sheet(0.30, 0.20), "small shared part, large leftover")]

out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
       f'font-family="-apple-system, Helvetica, Arial, sans-serif">',
       f'<rect width="{W}" height="{H}" fill="{SURFACE}"/>',
       '<defs>'
       f'<marker id="head" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{INK2}"/></marker>'
       f'<marker id="headb" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{BLUE}"/></marker>'
       f'<marker id="heado" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{ORANGE}"/></marker>'
       '</defs>',
       text(40, 50, "One species, two chains: the same arrow, split two ways", size=19, weight="600"),
       text(40, 78, "Both chains reach the same arrow, so they agree on everything the arrow sets: predictions, the gradients,", size=14, fill=INK2),
       text(40, 97, "and the measured coefficients. They disagree on the split, and so on the scales that set it.", size=14, fill=INK2)]

for i, (name, shared, how) in enumerate(SPLITS):
    P = projector(110 + 620 * i, 500, 230)
    o = (0, 0, 0)
    for end, lab, dx in (((1.7, 0, 0), "column 1", 30), ((0, 1.55, 0), "column 2", 30), ((0, 0, 1.2), "column 3", 0)):
        out.append(line(P, o, end, stroke=MUTED, stroke_width=1))
        x, y = P(end); out.append(text(x + dx, y + 4, lab, size=12, fill=MUTED, anchor="middle"))
    out.append(poly(P, SHEET, fill=BLUE_FILL, fill_opacity="0.75", stroke=BLUE, stroke_width=2))
    out.append(line(P, o, shared, stroke=BLUE, stroke_width=4, marker_end="url(#headb)"))
    out.append(line(P, shared, TIP, stroke=ORANGE, stroke_width=4, marker_end="url(#heado)"))
    out.append(line(P, o, TIP, stroke=INK2, stroke_width=2, stroke_dasharray="2 5", marker_end="url(#head)"))
    x, y = P(TIP); out.append(text(x + 12, y - 4, "the species' arrow", size=14))
    out.append(text(40 + 620 * i, 145, name, size=17, weight="600"))
    out.append(text(40 + 620 * i, 166, how, size=14, fill=INK2))

ly = 578
legend = [(40, "blue", "shared part (lies in the blend sheet)"), (420, "orange", "leftover (the species' own)"),
          (720, "arrow", "the whole response, the same in both chains")]
for x, kind, lab in legend:
    if kind == "blue":
        out.append(f'<line x1="{x}" y1="{ly}" x2="{x+30}" y2="{ly}" stroke="{BLUE}" stroke-width="4"/>')
    elif kind == "orange":
        out.append(f'<line x1="{x}" y1="{ly}" x2="{x+30}" y2="{ly}" stroke="{ORANGE}" stroke-width="4"/>')
    else:
        out.append(f'<line x1="{x}" y1="{ly}" x2="{x+30}" y2="{ly}" stroke="{INK2}" stroke-width="2" stroke-dasharray="2 5" marker-end="url(#head)"/>')
    out.append(text(x + 40, ly + 5, lab, size=14, fill=INK2))
out.append(text(40, 622, "Schematic with made-up points. Like 7 = 5 + 2 = 1 + 6: the data pin down the sum, not how it is divided.", size=12, fill=MUTED))
out.append("</svg>")
open(sys.argv[1], "w", encoding="utf-8").write("\n".join(out))
