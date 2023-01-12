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
cmd.remove('chain AP1+BP1+CP1+IN1')
cmd.select('PA', 'chain DP1')
cmd.select('PB1', 'chain EP1')
cmd.select('PB2', 'chain FP1')
cmd.select('ANP32A', 'chain GP1')

### To get PB2 only
cmd.hide('everything')
cmd.show('surface')
cmd.color('0xb3cde3', 'PB1')
cmd.color('grey70', 'PA')
cmd.color('0xfbb4ae', 'PB2')
cmd.color('grey40', 'ANP32A')
# cmd.set('transparency', 0.6, 'PB1')
cmd.orient()

###### END PREAMBLE

metric = 'subgraph_2b'

# knownAdaptive
for r in [559, 698]:
	cmd.color('red', 'resi {0} and PB2'.format(r))
	cmd.show('spheres', 'resi {{0}} and PB2'.format(r))
	
for r in [312, 343, 557, 573]:
	cmd.color('grey20', 'resi {0} and PA'.format(r))
	cmd.show('spheres', 'resi {{0}} and PA'.format(r))

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
cmd.rotate('z', angle=-90)
cmd.rotate('y', angle=10)
cmd.draw(width=1200, height=1000)
cmd.png('{0}_{1}.png'.format(structure, metric))
# cmd.rotate('y', angle=-80)
# cmd.draw(width=1200, height=1000)
# cmd.png('{0}_{1}_y.png'.format(structure, metric))
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
