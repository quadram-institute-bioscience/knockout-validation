import pandas
import argparse
import sys
import gffutils

def load_gff(filename, key="gene"):
   """
   Load a GFF file and returns a locus_tag -> gene_id dict
   """
   output = {}
   try:
        db = gffutils.create_db(filename, dbfn=filename+'.db')
   except:
        pass

   db = gffutils.FeatureDB(filename+'.db', keep_order=True)
   for line in db.all_features():
      if key in line.attributes and 'locus_tag' in line.attributes:
         id=line.attributes['locus_tag' ][0]
         output[ id ] = line.attributes[key][0]
      else:
         debprint('No {} in {}'.format(key, line.attributes))
   return output


def eprint(*args, **kwargs):
        """print to STDERR"""
        print(*args, file=sys.stderr, **kwargs)

def debprint(msg):
        if opt.debug:
            eprint('# {}'.format(msg))

def verprint(msg):
        if opt.verbose or opt.debug:
            eprint('* {}'.format(msg))

opt_parser = argparse.ArgumentParser(description='Manipulate feature counts table')

opt_parser.add_argument('-i', '--input',
                        help='Path to featureCounts output',
                        required=True)

opt_parser.add_argument('-f', '--filter',
                        type=float,
                        help='Discard genes having median lower than this',
                        default=1e-4)


opt_parser.add_argument('-a', '--annotation',
                        help='GFF annotation to print gene names (optional)')


opt_parser.add_argument('-n', '--normalizefactor',
                        type=float,
                        help='Normalize sum to this total [1]',
                        default=1)

opt_parser.add_argument('-m', '--mincounts',
                        type=float,
                        help='Discard samples with less counts [1000]',
                        default=1000)

opt_parser.add_argument('-t', '--threshold', type=float, help="Fraction of median used as threshold [0.3]", default=0.3)
                    
opt_parser.add_argument('-v', '--verbose',
                        help='Enable verbose output',
                        action='store_true')

opt_parser.add_argument('-d', '--debug',
                        help='Enable very verbose output',
                        action='store_true')

opt = opt_parser.parse_args()

annotation = {}
if opt.annotation:
    annotation = load_gff(opt.annotation)
    verprint("Annotation file ({}) loaded: {} entries".format(opt.annotation, len(annotation)))
else:
    verprint("Annotation file not provided. Skipping")


dataframe = pandas.read_csv(opt.input, skiprows=1, sep='\t')


(features, cols) = dataframe.shape


norm = dataframe.drop(dataframe.columns[[0, 1, 2, 3, 4, 5]], axis=1)
(features, cols) = norm.shape
verprint('Loaded matrix: {} samples, {} genes'.format(cols, features))

norm.drop([col for col, val in norm.sum().iteritems() if val < opt.mincounts], axis=1, inplace=True)
(features, cols) = norm.shape
verprint('Polished matrix: {} samples, {} genes (after removing low coverage samples)'.format(cols, features))

normdiv = norm.div(norm.sum(axis=0), axis=1)
normdiv = normdiv.mul(opt.normalizefactor, axis=1)
median = normdiv.median(axis=1)
sum  = normdiv.sum(axis=1)
max  = normdiv.max(axis=1)

debprint( 'Sum of counts:\n{}'.format(norm.sum(axis = 0, skipna = True)) )

counters = {
  'filtered_gene':   0,
  'discarded_value': 0,
}

for index, row in normdiv.iterrows():
    m = median[index]
    i = 0
    if m <= (opt.normalizefactor * opt.filter):
        counters['filtered_gene'] += 1
        debprint('gene {} has median {} < {}'.format( index, m, opt.normalizefactor * opt.filter))
        continue

    for value in row:
        i += 1
        if value >= (m * opt.threshold):
            counters['discarded_value'] += 1
        else:
            geneID =  dataframe.loc[ index , 'Geneid' ]
            sampleID = normdiv.columns[i-1]
            geneName = ''
            if geneID in annotation:
                geneName = annotation[geneID]
            # Print: GeneID (GeneName) SampleID, NormalizedCounts, MedianCounts, MaxCounts
            print('{}\t{}\t{}\t{}\t{}\t{}'.format( geneID, geneName, sampleID, value, median[index], max[index]))


print('#Discarded genes:\t{} (< {})\n#Discarded values:\t{}'.format( 
	counters['filtered_gene'],
	(opt.normalizefactor * opt.filter),
	counters['discarded_value'] )
)
