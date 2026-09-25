#!/usr/bin/env python3
"""Build the TrackCluster-RS methods schematic as editable SVG and vector PDF.

Requires reportlab. PNG exports are rendered from the PDF with Poppler.
All depicted tracks, matrices and call states are schematic, not measured data.
"""
from pathlib import Path
from html import escape
import json
import math

from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/figures"
PDF_OUT = ROOT / "output/pdf"
OUT.mkdir(parents=True, exist_ok=True)
PDF_OUT.mkdir(parents=True, exist_ok=True)

W, H = 1260, 1410
SCALE = (183 / 25.4 * 72) / W
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
for name, filename in [("Arial", "Arial.ttf"), ("Arial-Bold", "Arial Bold.ttf")]:
    pdfmetrics.registerFont(TTFont(name, str(FONT_DIR / filename)))

C = {
    "ink": "#202C34", "muted": "#59656D", "rule": "#D9E0E3",
    "gray": "#94A0A6", "pale": "#E3E9EC", "faint": "#F6F8F9",
    "blue": "#387BA6", "blue_light": "#C7DAE7", "blue_pale": "#F0F5F9",
    "teal": "#168777", "teal_light": "#CAE3DC", "teal_pale": "#EFF7F4",
    "orange": "#C58628", "orange_light": "#F0DFBF", "orange_pale": "#FBF6EB",
    "purple": "#82659C", "purple_light": "#E8DFF0", "purple_pale": "#F6F2F8",
    "white": "#FFFFFF",
}
PDF = PDF_OUT / "trackcluster_rs_fig1.pdf"
cv = canvas.Canvas(str(PDF), pagesize=(W * SCALE, H * SCALE), pageCompression=1)
cv.setTitle("TrackCluster-RS: isoform discovery, quantification and RNA modification analysis")
cv.setAuthor("TrackCluster-RS")
cv.setSubject("Original methods schematic; illustrative tracks and values, no benchmark claims")
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="183mm" height="{H/W*183:.3f}mm" viewBox="0 0 {W} {H}">',
       '<title>TrackCluster-RS workflow and core analysis modules</title>',
       '<desc>Five panels show the workflow, junction and terminal evidence, structural classification, multi-sample quantification, and optional isoform-resolved RNA modification analysis. All tracks and values are schematic.</desc>']
audit_text = []


def color(c):
    return C.get(c, c)


def pcolor(c):
    return HexColor(color(c))


def line(x1, y1, x2, y2, c="ink", sw=1.6, dash=None):
    d = f' stroke-dasharray="{dash[0]} {dash[1]}"' if dash else ""
    svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color(c)}" stroke-width="{sw}"{d}/>')
    cv.setStrokeColor(pcolor(c)); cv.setLineWidth(sw*SCALE)
    cv.setDash([v*SCALE for v in dash] if dash else [])
    cv.line(x1*SCALE,(H-y1)*SCALE,x2*SCALE,(H-y2)*SCALE)
    cv.setDash([])


def rect(x, y, w, h, fill="white", stroke=None, sw=1.0, radius=0):
    svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" fill="{color(fill) if fill else "none"}" stroke="{color(stroke) if stroke else "none"}" stroke-width="{sw}"/>')
    if fill: cv.setFillColor(pcolor(fill))
    if stroke: cv.setStrokeColor(pcolor(stroke))
    cv.setLineWidth(sw*SCALE)
    args=(x*SCALE,(H-y-h)*SCALE,w*SCALE,h*SCALE)
    if radius: cv.roundRect(*args,radius*SCALE,stroke=bool(stroke),fill=bool(fill))
    else: cv.rect(*args,stroke=bool(stroke),fill=bool(fill))


def circle(x, y, r, fill="white", stroke=None, sw=1.3):
    svg.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{color(fill) if fill else "none"}" stroke="{color(stroke) if stroke else "none"}" stroke-width="{sw}"/>')
    if fill: cv.setFillColor(pcolor(fill))
    if stroke: cv.setStrokeColor(pcolor(stroke))
    cv.setLineWidth(sw*SCALE)
    cv.circle(x*SCALE,(H-y)*SCALE,r*SCALE,stroke=bool(stroke),fill=bool(fill))


def poly(points, fill="ink", stroke=None, sw=1.3):
    points_s=" ".join(f"{x},{y}" for x,y in points)
    svg.append(f'<polygon points="{points_s}" fill="{color(fill) if fill else "none"}" stroke="{color(stroke) if stroke else "none"}" stroke-width="{sw}"/>')
    p=cv.beginPath();p.moveTo(points[0][0]*SCALE,(H-points[0][1])*SCALE)
    for x,y in points[1:]: p.lineTo(x*SCALE,(H-y)*SCALE)
    p.close()
    if fill: cv.setFillColor(pcolor(fill))
    if stroke: cv.setStrokeColor(pcolor(stroke))
    cv.setLineWidth(sw*SCALE);cv.drawPath(p,stroke=bool(stroke),fill=bool(fill))


def text(x, y, value, size=19, c="ink", bold=False, anchor="start"):
    font="Arial-Bold" if bold else "Arial"
    width=pdfmetrics.stringWidth(value,font,size)
    left=x if anchor=="start" else x-width/2 if anchor=="middle" else x-width
    audit_text.append({"text":value,"x":round(left,2),"y":y,"width":round(width,2),"font_px":size})
    svg.append(f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{700 if bold else 400}" text-anchor="{anchor}" fill="{color(c)}">{escape(value)}</text>')
    cv.setFont(font,size*SCALE);cv.setFillColor(pcolor(c))
    cv.drawString(left*SCALE,(H-y)*SCALE,value)


def lines(x,y,values,size=19,c="ink",bold=False,leading=25,anchor="start"):
    for i,s in enumerate(values):text(x,y+i*leading,s,size,c,bold,anchor)


def arrow(x1,y1,x2,y2,c="muted",sw=1.9,head=7,dash=None):
    line(x1,y1,x2,y2,c,sw,dash)
    a=math.atan2(y2-y1,x2-x1)
    pts=[(x2,y2),(x2-head*math.cos(a)+head*.46*math.sin(a),y2-head*math.sin(a)-head*.46*math.cos(a)),(x2-head*math.cos(a)-head*.46*math.sin(a),y2-head*math.sin(a)+head*.46*math.cos(a))]
    poly(pts,c)


def track(x,y,w,exons,c="gray",height=9,introns=None,sl=False,end=False):
    """Transcript-order plus-strand cartoon with exon intervals on [0,1]."""
    introns=introns or c
    for (_,b),(a,_) in zip(exons,exons[1:]):
        line(x+b*w,y,x+a*w,y,introns,1.5)
        if (a-b)*w>32:
            mid=x+(a+b)/2*w
            line(mid-3,y-3,mid+1,y,introns,1.1)
            line(mid+1,y,mid-3,y+3,introns,1.1)
    for a,b in exons:rect(x+a*w,y-height/2,(b-a)*w,height,c)
    if sl:rect(x+exons[0][0]*w-5,y-height/2-3,5,height+6,"orange")
    if end:circle(x+exons[-1][1]*w,y,4.5,"orange")


def panel(letter,title,x,y,w):
    text(x,y,letter,29,bold=True)
    text(x+31,y-1,title,23,bold=True)
    line(x,y+16,x+w,y+16,"rule",1.2)


def number(x,y,n):
    circle(x,y-6,12,"blue_pale")
    text(x,y,str(n),17,"blue",True,"middle")


def legend_sw(x,y,c,label):
    rect(x,y-9,15,9,c)
    text(x+22,y,label,17,"muted")


E=[(0,.22),(.40,.56),(.75,1)]
SKIP=[(0,.22),(.75,1)]
SHORT5=[(.13,.22),(.40,.56),(.75,1)]
SHORT3=[(0,.22),(.40,.56),(.75,.88)]

# Header: restrained and legible at a two-column journal width.
rect(0,0,W,H,"white")
text(32,43,"TrackCluster-RS",32,bold=True)
text(32,77,"Isoform discovery, quantification and RNA modification analysis",23,"muted")
line(32,97,1228,97,"rule",1.4)

# a | Workflow.
panel("a","From long reads to an isoform-resolved transcriptome",32,133,1196)
for x,heading in [(52,"Long-read inputs"),(351,"Pool & cluster by gene"),(678,"Shared isoform catalog"),(1000,"Biological outputs")]:
    text(x,181,heading,21,bold=True)
text(52,208,"Aligned reads (BAM / BED)",18,"muted")
for i,(exons,cc) in enumerate([(E,"gray"),(SHORT5,"gray"),(SKIP,"gray"),(E,"gray")]):
    track(83,236+20*i,181,exons,cc,8)
text(52,250,"S1",16,"muted")
text(52,290,"S2",16,"muted")
text(52,324,"One or multiple samples",18,"muted")
track(83,355,181,E,"ink",9)
text(52,386,"Reference: GTF / GFF3 / BED",18,"muted")
arrow(286,263,325,263)

for x,label in [(356,"Gene 1"),(499,"Gene 2")]:
    text(x+54,218,label,18,"muted",False,"middle")
    rect(x,230,110,78,"faint",None,radius=4)
    for i,exons in enumerate([E,SHORT5,SKIP]):track(x+12,246+i*22,86,exons,"gray",7)
text(351,338,"Junction chains (default)",19,"blue",True)
text(351,363,"or exon / intron overlap",18,"muted")
text(351,388,"Parallel processing in Rust",18,"muted")
arrow(622,263,657,263)

for i,(exons,cc) in enumerate([(E,"teal"),(SKIP,"blue"),(SHORT3,"orange")]):
    text(678,236+i*30,f"I{i+1}",17,cc,True)
    track(711,230+i*30,204,exons,cc,11)
text(678,323,"Known + novel structures",18,"muted")
arrow(798,335,798,351,"muted",1.6,6)
rect(675,359,247,37,"teal_pale",None,radius=4)
text(798.5,383,"Unique read-to-isoform map",17,"teal",True,"middle")

arrow(938,230,986,230,"muted",1.6,6)
line(922,377,951,377,"teal",1.6)
line(951,296,951,377,"teal",1.6)
for yy in [296,366]:arrow(951,yy,986,yy,"teal",1.6,6)
track(1001,228,75,E,"blue",7)
text(1090,233,"Structures",19)
for i,h in enumerate([15,24,11]):rect(1003+17*i,306-h,10,h,"teal")
lines(1070,291,["Counts &","isoform usage"],19,leading=23)
for i,cc in enumerate(["purple","white","purple"]):circle(1007+20*i,365,6,cc,"purple",1.4)
lines(1080,360,["RNA","modifications*"],19,leading=23)

# b | Core junction clustering, including the two terminal-evidence rules.
panel("b","Junction correction and end-aware clustering",32,454,742)
for x,n,labels in [(53,1,["Correct splice","junctions"]),(306,2,["Collapse compatible","5′ fragments"]),(559,3,["Retain supported","terminal variants"])]:
    number(x+11,498,n)
    lines(x+32,494,labels,19,bold=True,leading=23)

# Correction: mismatching positions converge on supported reference/read sites.
for xx in [53+.22*193,53+.40*193,53+.56*193,53+.75*193]:
    line(xx,541,xx,648,"rule",1.1,(3,3))
track(53,552,193,E,"ink",8)
track(53,579,193,[(0,.24),(.42,.57),(.75,1)],"gray",8)
track(53,606,193,[(.06,.21),(.38,.55),(.74,1)],"gray",8)
track(53,633,193,E,"blue",9)
line(92,579,107,579,"orange",2.8)
line(123,606,137,606,"orange",2.8)
text(53,668,"Read + reference support",17,"muted")
arrow(260,592,286,592)

track(307,552,193,E,"teal",10)
track(307,581,193,SHORT5,"teal_light",8)
track(307,608,193,[(.43,.56),(.75,1)],"teal_light",8)
arrow(403,620,403,636,"teal",1.6,6)
track(307,648,193,E,"teal",11)
text(307,677,"One representative + members",16.5,"muted")
arrow(514,592,540,592)

track(560,552,192,E,"teal",10)
track(560,593,192,SHORT5,"blue",10,sl=True)
track(560,634,192,SHORT3,"orange",10,end=True)
text(560,677,"Distinct supported endpoints",17,"muted")

# Terminal zoom, using repeated independent molecules as visible evidence.
line(52,701,753,701,"rule",1.0)
text(52,730,"5′ protection with optional SL evidence",18.5,bold=True)
text(418,730,"3′ protection within a splice chain",18.5,bold=True)
for xx,yy in [(85,768),(85,787),(85,806)]:
    track(xx,yy,245,SHORT5,"blue_light",7,sl=True)
track(85,832,245,SHORT5,"blue",10,sl=True)
line(116.9,748,116.9,843,"orange",1.1,(3,3))
for yy in [768,787,806]:track(442,yy,254,SHORT3,"orange_light",7,end=True)
track(442,832,254,SHORT3,"orange",10,end=True)
line(665.5,748,665.5,843,"orange",1.1,(3,3))
lines(52,868,["Optional SL evidence protects supported", "alternative starts; SL is off by default."],17,"muted",leading=23)
lines(418,868,["Repeated support protects alternative ends.","Support is fixed before batching."],17,"muted",leading=23)
text(52,928,"Alternative mode: two-pass exon / intron overlap clustering",18,"muted")

# c | Legacy structural interpretation, redrawn as paired exon cartoons.
panel("c","Interpret transcript structures",812,454,416)
legend_sw(814,496,"ink","Reference")
legend_sw(963,496,"blue","Query isoform")
examples=[
    ("Terminal exon gain / loss",E,[(0,.22),(.40,.56)]),
    ("UTR extension / truncation",E,[(.08,.22),(.40,.56),(.75,.88)]),
    ("Exon inclusion / skipping",E,SKIP),
    ("Intron retention",E,[(0,.56),(.75,1)]),
    ("Alternative splice site",E,[(0,.22),(.46,.56),(.75,1)]),
]
for i,(label,ref,query) in enumerate(examples):
    y=529+i*71
    text(814,y,label,18.5,bold=False)
    track(838,y+20,336,ref,"ink",7)
    track(838,y+41,336,query,"blue",8)
text(814,887,"Gene fusion",18.5)
track(838,907,146,[(0,.35),(.64,1)],"ink",7)
track(1030,907,144,[(0,.4),(.68,1)],"ink",7)
track(838,927,336,[(0,.152),(.278,.435),(.571,.743),(.863,1)],"blue",8)
text(814,958,"11 event labels + annotated reference",17,"muted")

# d | Shared-catalog per-sample quantification (schematic matrix, no real data).
panel("d","Quantify across samples",32,1007,535)
text(52,1047,"One shared catalog, separate sample counts",19)
for i,(exons,cc) in enumerate([(E,"teal"),(SKIP,"blue"),(SHORT3,"orange")]):
    y=1146+i*48
    text(53,y+5,f"I{i+1}",18,cc,True)
    track(86,y,166,exons,cc,11)
arrow(270,1194,296,1194,"muted",1.6,6)
for i,label in enumerate(["A1","A2","B1","B2"]):text(329+i*49,1113,label,18,"muted",False,"middle")
for x1,x2,lab in [(309,399,"Condition A"),(407,497,"Condition B")]:
    line(x1,1091,x2,1091,"gray",2.2)
    text((x1+x2)/2,1081,lab,16,"muted",False,"middle")
shades=[['#5A9D91','#438F83','#CDE4DE','#BBDAD2'],['#BCD4E4','#CCDFEC','#477FA6','#387BA6'],['#F3E8D3','#EEDBB9','#DEC18C','#E8D1AA']]
for i in range(3):
    for j in range(4):rect(307+j*49,1126+i*48,43,39,shades[i][j],"white",.8)
text(308,1283,"Isoform × sample count matrix",17,"muted")
rect(52,1296,493,39,"teal_pale",None,radius=4)
text(298.5,1321,"Counts → within-gene usage → group summaries",18,"teal",False,"middle")
text(52,1361,"Each assigned molecule contributes to one isoform.",17,"muted")

# e | Optional modification analysis: normalize callers, then join exact mapping.
panel("e","Resolve RNA modifications by isoform*",617,1007,611)
text(636,1048,"Dorado modBAM / m6Anet calls",19,"purple",True)
arrow(1030,1042,1060,1042,"purple",1.6,6)
text(1073,1048,"Normalize sites",18,"purple")
text(636,1076,"Join unique assignments; retain per-read call states",17.5,"muted")
for x,label in [(837,"Site 1"),(975,"Site 2"),(1120,"Site 3")]:
    text(x,1110,label,17,"muted",False,"middle")
    line(x,1125,x,1246,"rule",1.2,(3,3))

# Place sites at fixed genomic positions. Site 2 is absent in the skipped exon.
site_x=[837,975,1120]
exons_mod=[E,SKIP,SHORT3]
states=[[['mod','mod'],['mod','unmod'],['mod','unmod']],
        [['unmod','unmod'],['absent'],['mod','unmod']],
        [['unknown','unmod'],['mod','unmod'],['unmod','unmod']]]
for i,(exons,cc) in enumerate(zip(exons_mod,["teal","blue","orange"])):
    y=1140+46*i
    text(638,y+6,f"I{i+1}",18,cc,True)
    track(770,y,405,exons,cc,11)
    for xx,observations in zip(site_x,states[i]):
        if observations==['absent']:
            rect(xx-10,y-11,20,22,"white")
            line(xx-6,y-5,xx+6,y+5,"gray",1.5)
            line(xx-6,y+5,xx+6,y-5,"gray",1.5)
            continue
        for j,state in enumerate(observations):
            yy=y-8+j*16
            if state=='mod':circle(xx,yy,6.3,"purple","white",1.0)
            elif state=='unmod':circle(xx,yy,6.3,"white","purple",1.6)
            elif state=='unknown':
                circle(xx,yy,6.7,"white","gray",1.2)
                text(xx,yy+4.3,"?",12,"muted",True,"middle")

legend_y=1280
circle(645,legend_y-5,5.8,"purple")
text(658,legend_y,"Modified",16,"muted")
circle(757,legend_y-5,5.8,"white","purple",1.5)
text(770,legend_y,"Unmodified",16,"muted")
circle(888,legend_y-5,6.5,"white","gray",1.1)
text(888,legend_y-1,"?",12,"muted",True,"middle")
text(902,legend_y,"Unknown",16,"muted")
line(1016,legend_y-10,1026,legend_y,"gray",1.4)
line(1016,legend_y,1026,legend_y-10,"gray",1.4)
text(1036,legend_y,"Site absent",16,"muted")

rect(636,1296,574,39,"purple_pale",None,radius=4)
text(923,1321,"Modified fraction = modified / callable molecules",18,"purple",False,"middle")
text(636,1361,"Coverage QC · Wilson intervals · effect-size contrasts",17,"muted")

# Figure-wide notes distinguish scientific cartoons from measured results.
line(32,1378,1228,1378,"rule",1.1)
text(32,1398,"Schematic examples; tracks and matrix intensities are illustrative.",16,"muted")
text(1228,1398,"* Optional modification analysis; contrasts are descriptive.",16,"muted",False,"end")

cv.showPage();cv.save()
svg.append('</svg>')
(OUT/"trackcluster_rs_fig1.svg").write_text("\n".join(svg)+"\n")
(OUT/"trackcluster_rs_fig1.text-audit.json").write_text(json.dumps(audit_text,indent=2,ensure_ascii=False)+"\n")
outside=[t for t in audit_text if t['x']<0 or t['x']+t['width']>W or t['y']>H]
if outside: raise RuntimeError(f"Text outside canvas: {outside}")
print(f"SVG: {OUT/'trackcluster_rs_fig1.svg'}")
print(f"PDF: {PDF}")
print(f"Page: 183 × {H/W*183:.1f} mm; {len(audit_text)} editable text labels")
