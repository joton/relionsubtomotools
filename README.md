# reliontomotools

Additional tools for subtomogram analysis in Relion tomo.


## Install
```bash
pip install reliontomotools
```

## Scripts

```python
warptomo2relion

motl2relion

relion2motl
```

Converts refinement of deformation models, particle poses and CTF paramaters from WARP/M to Relion tomo.

### Example
```bash
warptomo2relion -i 'WarpXML/TS_*.xml' -s 'tomograms/TS_*/TS_*_aligned.mrc' -d 1800 -o WarpConverted -p Refine3D/job010/run_data.star

motl2relion allmotl_at_17.em -o allmotl_at_17.star -p 1.75
motl2relion --evenodd allmotl_at_even_17.em allmotl_at_odd_17.em -o allmotl_at_17_2.star -p 1.75

relion2motl -i allmotl_at_17.star -o allmotl_at_17_star.em --angpix 1.75
```
