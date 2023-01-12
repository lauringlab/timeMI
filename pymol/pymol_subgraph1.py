###### PREAMBLE

import pymol
from pymol import cmd
import sys

structure = '6rr7'

cmd.reinitialize()
cmd.set('bg_rgb','[1,1,1]') # white
cmd.set('antialias','2')
cmd.set('ray_opaque_background','off')
cmd.set('depth_cue', 'off')

### Modify here
url = 'https://files.rcsb.org/download/6RR7.pdb'
cmd.load(url, 'pcs')
cmd.select('PA', 'chain A')
cmd.select('PB1', 'chain B')
cmd.select('PB2', 'chain C')
cmd.select('RNA', 'chain D:F')

### To get PB2 only
cmd.hide('everything')
cmd.show('surface')
cmd.color('white', 'RNA')
cmd.color('0xb3cde3', 'PB1')
cmd.color('grey70', 'PA')
cmd.color('0xfbb4ae', 'PB2')
# cmd.set('transparency', 0.6, 'PB1')
# cmd.set('transparency', 0.3, 'PA')
# cmd.set('transparency', 0.3, 'PB2')
cmd.orient()

###### END PREAMBLE

metric = 'subgraph_1'
# Nter
# for r in range(1, 247):
#     cmd.color('lightpink', 'resi {0} and PB2'.format(r))
# # Mid-link
# for r in [*range(247,319), *range(481,538)]:
#     cmd.color('wheat', 'resi {0} and PB2'.format(r))    
# # Capbinding
# for r in range(319, 481):
#     cmd.color('wheat', 'resi {0} and PB2'.format(r))
# # 627
# for r in range(538, 680):
#     cmd.color('salmon', 'resi {0} and PB2'.format(r))
# # NLS
# for r in range(680, 742):
#     cmd.color('wheat', 'resi {0} and PB2'.format(r))

# knownAdaptive
for r in [194, 227, 338, 569]:
	cmd.color('red', 'resi {0} and PB2'.format(r))
	cmd.show('spheres', 'resi {{0}} and PB2'.format(r))



###### SETVIEW
# cmd.set_view ('\
#     -0.148019463,   -0.358468264,    0.921732306,\
#     -0.189112291,   -0.904540896,   -0.382152349,\
#      0.970734179,   -0.230877280,    0.066099040,\
#      0.000000000,   -0.000000000, -368.771972656,\
#    134.908599854,  136.783447266,  137.801361084,\
#    290.742645264,  446.801300049,  -20.000000000 ')

###### POSTAMBLE
##No ray for faster development
cmd.select(None)
cmd.disable()
cmd.enable()
cmd.reset()
cmd.orient()
cmd.set('specular', "off")
cmd.rotate('y', angle=180)
# cmd.rotate('x', angle=-20)
# cmd.rotate('z', angle=180)
cmd.draw(width=1200, height=1000)
cmd.png('{0}_{1}.png'.format(structure, metric))
cmd.rotate('y', angle=-60)
cmd.draw(width=1200, height=1000)
cmd.png('{0}_{1}_y.png'.format(structure, metric))
cmd.rotate('y', angle=80)
cmd.rotate('x', angle=-90)
cmd.draw(width=1200, height=1000)
cmd.png('{0}_{1}_x.png'.format(structure, metric))

#Ray for publication
# cmd.select(None)
# cmd.set('specular', "off")
# cmd.rotate('y', angle=180)
# cmd.rotate('x', angle=180)
# cmd.draw(width=1000, height=1000)
# cmd.png('{0}_{1}.png'.format(structure, metric), ray=1)
# cmd.rotate('y', angle=180)
# cmd.draw(width=1000, height=1000)
# cmd.png('{0}_{1}_y.png'.format(structure, metric), ray=1)
# cmd.rotate('y', angle=180)
# cmd.rotate('x', angle=90)
# cmd.draw(width=1000, height=1000)
# cmd.png('{0}_{1}_x.png'.format(structure, metric), ray=1)
