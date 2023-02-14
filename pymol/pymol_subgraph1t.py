###### PREAMBLE

import pymol
from pymol import cmd
import sys

structure = '6xzr'

cmd.reinitialize()
cmd.set('bg_rgb','[1,1,1]') # white
cmd.set('antialias','2')
cmd.set('ray_opaque_background','off')
cmd.set('depth_cue', 'off')

### Modify here
url = 'https://files.rcsb.org/download/6XZR.cif'
cmd.load(url, 'pcs')
# cmd.select('PA', 'chain AP1+DP1')
# cmd.select('PB1', 'chain BP1+EP1')
# cmd.select('PB2', 'chain CP1+FP1')
cmd.select('RNA', 'chain IN1')
cmd.select('ANP32A', 'chain GP1')

### To get PB2 only
cmd.hide('everything')
cmd.show('surface')
cmd.color('white', 'RNA')
cmd.color('0xb3cde3', 'chain BP1+EP1')
cmd.color('grey70', 'chain AP1+DP1')
cmd.color('0xfbb4ae', 'chain CP1+FP1')
cmd.color('grey40', 'ANP32A')
cmd.orient()
###### END PREAMBLE

metric = 'subgraph_1'

# knownAdaptive
for r in [194, 227, 338, 569]:
	cmd.color('red', 'resi {0} and chain CP1'.format(r))
	cmd.show('spheres', 'resi {{0}} and chain CP1'.format(r))

# for r in [227]:
# 	cmd.color('orange', 'resi {0} and PB2'.format(r))
# 	cmd.show('spheres', 'resi {{0}} and PB2'.format(r))
# 
# for r in [338]:
# 	cmd.color('yellow', 'resi {0} and PB2'.format(r))
# 	cmd.show('spheres', 'resi {{0}} and PB2'.format(r))
# 
# for r in [569]:
# 	cmd.color('green', 'resi {0} and PB2'.format(r))
# 	cmd.show('spheres', 'resi {{0}} and PB2'.format(r))


###### POSTAMBLE
##No ray for faster development
cmd.select(None)
cmd.disable()
cmd.enable()
cmd.reset()
cmd.orient()
cmd.set('specular', "off")
cmd.rotate('y', angle=180)
cmd.rotate('x', angle=170)
cmd.draw(width=1200, height=1000)
cmd.png('{0}_{1}.png'.format(structure, metric))
cmd.set('transparency', 0.8, 'chain EP1+DP1+FP1')
cmd.draw(width=1200, height=1000)
cmd.png('{0}_{1}_trans.png'.format(structure, metric))



# cmd.rotate('y', angle=80)
# cmd.rotate('x', angle=-90)
# cmd.draw(width=1200, height=1000)
# cmd.png('{0}_{1}_x.png'.format(structure, metric))

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
