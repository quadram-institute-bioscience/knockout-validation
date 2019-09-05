#!/usr/bin/env python
"""
RETRIEVED FROM: https://github.com/lpryszcz/bin/blob/master/embl2gtf.py

Convert embl to gtf.

USAGE:
cat in.embl | embl2gtf.py > out.gtf

NOTE:
It's designed to work with gb files coming from EMBL. gene is used as gene_id and transcript_id (locus_tag if gene not present).
Only entries having types in allowedTypes = ['gene','CDS','tRNA','tmRNA','rRNA','ncRNA'] are stored in GTF. Need to include exon processing.
No frame info is processed. Need to be included in order to process genes having introns!
In addition, some embl files misses gene entries - then gene is considered as entire CDS/*RNA

AUTHOR:
Leszek Pryszcz
lpryszcz@crg.eu

Version 0.1

"""

import os, sys, urllib
try:
    from urllib import quote  # Python 2.X
except ImportError:
    from urllib.parse import quote  # Python 3+
from datetime import datetime
from Bio      import SeqIO

def _get_unique_id( gene_id,products,i=1 ):
  """ """
  new_id = '%s.%s' % ( gene_id,i )
  while new_id in products:
    i += 1
    new_id = '%s.%s' % ( gene_id,i )

  return new_id

def embl2gtf( source='embl2gtf',allowedTypes=set(['CDS', 'tRNA', 'tmRNA', 'rRNA', 'ncRNA']) ): #'mRNA', 'gene',
  """
  """
  #make sure no duplicated products
  products = set()

  handle = sys.stdin
  for r in SeqIO.parse( handle,'embl' ):
    acc     = r.name # r.id r.description
    skipped = 0
