import os
import sys
from argparse import ArgumentParser

import numpy as np
#from astropy.cosmology import Planck15
from astropy.table import Table
import simsurvey
import sncosmo

def load_fields_ccd(fields_file='ztf_fields.txt', ccd_file='ztf_ccd_corners.txt'):
    fields_raw = np.genfromtxt(fields_file, comments='%')

    fields = {'field_id': np.array(fields_raw[:,0], dtype=int),
              'ra': fields_raw[:,1],
              'dec': fields_raw[:,2]}

    ccd_corners = np.genfromtxt(ccd_file, skip_header=1)
    ccds = [ccd_corners[np.array([0,1,3,2])+4*k, :2] for k in range(16)]

    return fields, ccds
    
def load_ztf_bands(bandpass_dir=''):
    bands = {
        'ztfg' : 'ztfg_eff.txt',
        'ztfr' : 'ztfr_eff.txt',
        'ztfi' : 'ztfi_eff.txt',
    }

    for bandname in bands.keys() :
        fname = bands[bandname]
        b = np.loadtxt(os.path.join(bandpass_dir,fname))
        band = sncosmo.Bandpass(b[:,0], b[:,1], name=bandname)
        sncosmo.registry.register(band, force=True)

def main():
    if 'SFD_DIR' not in os.environ.keys():
        raise ValueError('Please set $SFD_DIR to where you downloaded the SFD dust maps.')
    
    parser = ArgumentParser(description='Run lightcurve simulations for ZTF')
    parser.add_argument('plan', type=str,
                        help='CSV file containing survey plan')
    parser.add_argument('transient', type=str,
                        help='Transient type (e.g. "Ia" for SNe Ia)')
    parser.add_argument('template', type=str,
                        help='Transient template (e.g. "salt2" or "nugent")')
    parser.add_argument('-z', '--redshift', default=None, nargs=2,
                        help='redshift boundaries', type=float)
    parser.add_argument('--no-weather', action='store_true',
                        help='do not apply weather loss')
    
    
    args = parser.parse_args()

    if args.redshift is None:
        args.redshift = (0, 0.2)

    obs = Table.read(args.plan, format='ascii.csv')
    fields, ccds = load_fields_ccd()
    load_ztf_bands()
    
    plan = simsurvey.SurveyPlan(time=obs['time'], band=obs['band'], obs_field=obs['field'],
                            skynoise=obs['skynoise'], comment=obs['comment'],
                            fields={k: v for k, v in fields.items()
                                    if k in ['ra', 'dec', 'field_id',
                                             'width', 'height']},
                            ccds=ccds)

    mjd_range = (plan.pointings['time'].min(), plan.pointings['time'].max())
    ra_range = (0, 360)
    dec_range = (-40, 90)

    tr = simsurvey.get_transient_generator((args.redshift[0], args.redshift[1]),
                                           transient=args.transient, 
                                           template=args.template,
                                           ra_range=ra_range,
                                           dec_range=dec_range,
                                           mjd_range=(mjd_range[0],
                                                      mjd_range[1]))

    survey = simsurvey.SimulSurvey(generator=tr, plan=plan)

    lcs = survey.get_lightcurves(progress_bar=True)

    if not os.path.exists('lcs'):
        os.makedirs('lcs')

    prev_files = [fn for fn in os.listdir('lcs')
                  if fn.startswith('lcs_%s_%s'%(args.transient, args.template))]
    if len(prev_files) > 0:
        k = max([int(fn.split('.')[0].split('_')[-1]) for fn in prev_files]) + 1
    else:
        k = 0
            
    outfile = 'lcs/lcs_%s_%s_%06i.pkl'%(args.transient, args.template, k)

    if not args.no_weather:
        real_nights = np.genfromtxt('hours_per_obsnight.dat')
        idx = np.concatenate((range(31,365), range(31)))
        rn_2016 = np.where(real_nights[idx, -1] > 3.5)[0] + 2458151

        def filterfunc(lc):
            mask = np.array([(int(t_) in rn_2016) and (t_ < lc.meta['t0'] + 100) for t_ in lc['time']])
            
            return lc[mask]

        lcs = lcs.filter(filterfunc)
        
    lcs.save(outfile)
    
if __name__ == '__main__':
    main()
