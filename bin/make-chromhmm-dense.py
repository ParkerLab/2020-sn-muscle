#!/usr/bin/env python

import os
import sys
import pandas as pd

COLORS = {
        'Active_enhancer': '255,195,77',
        'Active_TSS': '255,0,0',
        'Bivalent_poised_TSS': '205,92,92',
        'Flanking_TSS': '255,69,0',
        'Genic_enhancer': '194,225,5',
        'Quiescent_low_signal': '255,255,255',
        'Repressed_polycomb': '128,128,128',
        'Strong_transcription': '0,128,0',
        'Weak_enhancer': '255,255,0',
        'Weak_transcription': '0,100,0',
        'Weak_TSS': '255,69,0'
        }

per_state_bed_files = sys.argv[1:]

bed = pd.concat([pd.read_csv(f, delimiter='\t', header=None) for f in per_state_bed_files])
bed.columns = ['chrom', 'start', 'end', 'state']
bed['score'] = 0
bed['strand'] = '.'
bed['start2'] = bed.start
bed['end2'] = bed.end
bed['rgb'] = bed.state.map(lambda x: COLORS[x])
bed['x'] = 1
bed['length'] = bed.end - bed.start
bed['x2'] = 0

bed.to_csv(sys.stdout, sep='\t', header=False, index=False)
