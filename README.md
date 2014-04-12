condenser
=========

simple cli RDF to condensed graph of reaction(CGR) SDF converter for chemoinformatics.

also include CSV generator from fragmentor output for data mining.

based on
DOI: 10.1021/ci300149n
DOI: 10.1007/s10822-005-9008-0

usage
=====
condenser [-h] [--version] [--input INPUT] [--output OUTPUT]
                 [--coords COORDS]

optional arguments:

  -h, --help            show this help message and exit

  --version, -v         show program's version number and exit

  --input INPUT, -i INPUT

                        RDF inputfile

  --output OUTPUT, -o OUTPUT

                        SDF outputfile

  --coords COORDS, -c COORDS

                        write to SDF coordinates of:

                        1 - reagents or

                        2 - products or

                        3 - both (changed spec of MOL)

