# simsurvey-paper-scripts
Scripts and files required to reproduce the ZTF lightcurve simulations presented in Feindt et al. (2019)

Requirements
------------

- Install `simsurvey` and its required packages.
- Download the data for `sfdmap` from https://github.com/kbarbary/sfddata
- Set `$SFD_DIR` to the `sfddata` download path

Usage Example
-------------

```
python run_sim.py plan_sim_paper.csv Ia salt2
```

The last two arguments can be replaced by the other supernova types simulated in the paper and further options are available, see `python run_sim_paper.py -h`.

The script will save the output in the directory `lcs` and give it a file name that contains SN type and template. New runs with the same type and template will be save with an incremented number at the end. (Note that each file will be several hundred MB in size.)

To load the output, use the following lines:
```
from simsurvey import LigthcurveCollection
lcs = LightcurveCollection(load='/path/to/file')
```
Note that output files generated using Python 3 cannot be loaded with Python 2 and vice versa.