import os
import sys
import configparser
import argparse

import numpy as np

import pycbc.waveform
from pycbc.detector import Detector
from pycbc.io import FieldArray
from pycbc.inject import InjectionSet


"""
This script can be run using the injection ini file to produce simulated
data on H, L, and V. The syntax to run this code is:

>> python create_datafile.py -c inj_config.ini -o injection_new.hdf

which will expect the existance of a configuration script called inj_config.ini
that it will read the required data from, and then create an injection HDF5
file called injection_new.hdf. It will furthermore read this injection file
to produce the frames files in the three interferometers.
"""


parser = argparse.ArgumentParser()
parser.add_argument("-c", "--configfile", action="store",
                    help="Name of the config file")
parser.add_argument("-o", "--outputfile", action="store",
                    help="Name of the output HDF5 file")
args = parser.parse_args()


configParser = configparser.ConfigParser()
configParser.read(args.configfile)

# Intrinsic parameters
mass1 = configParser.getfloat('intrinsic', 'mass1')
mass2 = configParser.getfloat('intrinsic', 'mass2')
spin1x = configParser.getfloat('intrinsic', 'spin1x')
spin1y = configParser.getfloat('intrinsic', 'spin1y')
spin1z = configParser.getfloat('intrinsic', 'spin1z')
spin2x = configParser.getfloat('intrinsic', 'spin2x')
spin2y = configParser.getfloat('intrinsic', 'spin2y')
spin2z = configParser.getfloat('intrinsic', 'spin2z')
lambda1 = configParser.getfloat('intrinsic', 'lambda1')
lambda2 = configParser.getfloat('intrinsic', 'lambda2')

# Extrinsic parameter
inclination = configParser.getfloat('extrinsic', 'inclination')
distance = configParser.getfloat('extrinsic', 'distance')
ra = configParser.getfloat('extrinsic', 'ra')
dec = configParser.getfloat('extrinsic', 'dec')
polarization = configParser.getfloat('extrinsic', 'polarization')
coa_phase = configParser.getfloat('extrinsic', 'coa_phase')
t_coa = configParser.getfloat('extrinsic', 't_coa')

# Other parameters
f_lower = configParser.getfloat('other', 'f_lower')
f_ref = configParser.getfloat('other', 'f_ref')
approximant = configParser.get('other', 'approximant')
srate = configParser.getint('other', 'srate')
asdfile_ligo = configParser.get('other', 'asdfile_ligo')
asdfile_virgo = configParser.get('other', 'asdfile_virgo')


hp, hc = pycbc.waveform.get_td_waveform(approximant='IMRPhenomPv2_NRTidalv2',
                                        mass1=mass1, mass2=mass2,
                                        lambda1=lambda1, lambda2=lambda2,
                                        spin1x=spin1x, spin1y=spin1y, spin1z=spin1z,
                                        spin2x=spin2x, spin2y=spin2y, spin2z=spin2z,
                                        distance=distance, inclination=inclination,
                                        coa_phase=coa_phase, f_lower=f_lower, f_ref=f_ref,
                                        polarization=polarization, delta_t=1/srate)

sig_len = hp.duration
print("Duration of : signal = {}".format(sig_len))

duration=2**(np.ceil(np.log2(sig_len))+1)
print("Duration of : signal + post-trigger-duration  = {}".format(duration))

frame_start_time = int(t_coa - (1/2)*duration-32*duration)   ###psd calculation
frame_end_time = int(t_coa + (1/2)*duration) + 1
frame_duration = frame_end_time - frame_start_time
print("Frame end time = {}".format(frame_end_time))
print("Frame start time = {}".format(frame_start_time))


dtype = [('mass1', float), ('mass2', float),
         ('spin1x', float), ('spin2x', float),
         ('spin1y', float), ('spin2y', float),
         ('spin1z', float), ('spin2z', float),
         ('tc', float), ('distance', float),
         ('ra', float), ('dec', float),
         ('approximant', 'S32'), ('f_lower', float),
         ('f_ref', float), ('inclination', float),
         ('coa_phase', float), ('polarization', float),
         ('lambda1', float), ('lambda2', float)]

static_params = {'taper': 'start',
                 'eccentricity': 0.
                 }

num_inj = 1
samples = FieldArray(num_inj, dtype=dtype)

samples['mass1'] = mass1
samples['mass2'] = mass2
samples['spin1x'] = spin1x
samples['spin2x'] = spin2x
samples['spin1y'] = spin1y
samples['spin2y'] = spin2y
samples['spin1z'] = spin1z
samples['spin2z'] = spin2z
samples['tc'] = t_coa
samples['ra'] = ra
samples['dec'] = dec
samples['distance'] = distance
samples['inclination'] = inclination
samples['polarization'] = polarization
samples['f_ref'] = f_ref
samples['f_lower'] = f_lower
samples['coa_phase'] = coa_phase
samples['lambda1'] = lambda1
samples['lambda2'] = lambda2
samples['approximant'] = approximant

InjectionSet.write(args.outputfile, samples,
                   static_args=static_params, injtype='cbc',
                   cmd=" ".join(sys.argv))


cmd_L = '''pycbc_condition_strain --injection-file {} \
                                  --channel-name L1:SIM-STRAIN \
                                  --output-strain-file L-L1_STRAIN-{}-{}.gwf \
                                  --fake-strain-from-file {} \
                                  --fake-strain-seed 100 --low-frequency-cutoff 15 \
                                  --fake-strain-flow 15 --frame-duration {} \
                                  --gps-start-time {} --gps-end-time {} \
                                  --sample-rate {}'''.format(args.outputfile,
                                                             "{start}",
                                                             "{duration}",
                                                             asdfile_ligo,
                                                             frame_duration,
                                                             frame_start_time,
                                                             frame_end_time,
                                                             srate)

cmd_H = '''pycbc_condition_strain --injection-file {} \
                                  --channel-name H1:SIM-STRAIN \
                                  --output-strain-file H-H1_STRAIN-{}-{}.gwf \
                                  --fake-strain-from-file {} \
                                  --fake-strain-seed 100 --low-frequency-cutoff 15 \
                                  --fake-strain-flow 15 --frame-duration {} \
                                  --gps-start-time {} --gps-end-time {} \
                                  --sample-rate {}'''.format(args.outputfile,
                                                             "{start}",
                                                             "{duration}",
                                                             asdfile_ligo,
                                                             frame_duration,
                                                             frame_start_time,
                                                             frame_end_time,
                                                             srate)

cmd_V = '''pycbc_condition_strain --injection-file {} \
                                  --channel-name V1:SIM-STRAIN \
                                  --output-strain-file V-V1_STRAIN-{}-{}.gwf \
                                  --fake-strain-from-file {} \
                                  --fake-strain-seed 100 --low-frequency-cutoff 15 \
                                  --fake-strain-flow 15 --frame-duration {} \
                                  --gps-start-time {} --gps-end-time {} \
                                  --sample-rate {}'''.format(args.outputfile,
                                                             "{start}",
                                                             "{duration}",
                                                             asdfile_virgo,
                                                             frame_duration,
                                                             frame_start_time,
                                                             frame_end_time,
                                                             srate)

print("Making frames for L1")
os.system(cmd_L)
print("Making frames for H1")
os.system(cmd_H)
print("Making frames for V1")
os.system(cmd_V)













