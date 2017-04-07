"""
<keywords>
test, python, collision, rates, lipovka, wrathmal, roueff
</keywords>
<description>
read the collision data by flower and roueff that are used by lipovka.

.. todo:: the raw data can be read and parsed but the reduction into matrix
 form needs to be implemented.
</description>
<seealso>
</seealso>
"""
from __future__ import print_function
from read_cr_lipovka import Reader

reader = Reader('../../data/read/lipovka/flower_roueff_data.dat')

reader.read_data()
print('done')