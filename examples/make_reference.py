"""Generate illustrative sealed reference CAD and a matching one-sheet net.
The flat net is schematic: tabs / actual fold kinematics still need manual checking.
"""
from pathlib import Path
import cadquery as cq
import ezdxf
from shapely.geometry import box
from shapely.ops import unary_union

out = Path(__file__).parent
# all values in millimeters; interior 320x320x50, 3-mm board.
outer = cq.Workplane('XY').box(326,326,56,centered=(False,False,False))
inner = cq.Workplane('XY').box(320,320,50,centered=(False,False,False)).translate((3,3,3))
box_solid = outer.cut(inner)
cq.exporters.export(box_solid, str(out / 'reference_box.step'))
print('STEP mass g:', box_solid.val().Volume()*0.001*0.2)
# Net coordinates: 326x326 base, four side flaps 56 high and lid on one flap.
L = 326; H = 56
regions = [box(0,0,L,L), box(-H,0,0,L),box(L,0,L+H,L),
           box(0,-H,L,0),box(0,L,L,L+H),box(0,L+H,L,2*L+H)]
net=unary_union(regions)
net_dxf=ezdxf.new('R2010')
net_dxf.header['$INSUNITS']=4
net_dxf.layers.add('CUT',color=1)
net_dxf.layers.add('FOLD',color=3)
msp=net_dxf.modelspace()
coords=list(net.exterior.coords)
msp.add_lwpolyline(coords, close=True, dxfattribs={'layer':'CUT'})
for p1,p2 in [((0,0),(L,0)),((0,L),(L,L)),((0,0),(0,L)),((L,0),(L,L)),((0,L+H),(L,L+H))]:
 msp.add_line(p1,p2,dxfattribs={'layer':'FOLD'})
net_dxf.saveas(out / 'reference_net.dxf')
print('DXF mass g',net.area * 0.003 * .2)
