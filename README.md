# InserSim
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)
## Description 
InserSim is a set of scripts to simulate insertion variants in a reference genome and to compare insertion calls obtained by SV callers to the simulated ones.

The simulation folder contains three scripts : `Novo_dup_generator.py` to simulate de novo or dispersed duplication insertions, `New_Simulation_type.py` to simulate mobile element, tandem repeat and  tandem duplication insertions and `Microhomology.py` to simulate insertions containing microhomology.

## Requirements
- Python3
- biopython (pip3 install biopython)

## Installation
    git clone https://github.com/WesDe/InserSim.git

## Simulating novel and dispersed duplications
### Parameters :
| Option | Description |
| ------------ | ----------------------------------------- |
| -g | reference genome in FASTA format | 
| -b | bed file, to generate insertions in specific regions |
| -n | number of insertions (int) |
| -s | size of insertion (int) |
| -p | induced snp near insertion (True/false, false by default) |
| -o | output file name |
| -c | induced insertion in specific chromosome |
| -a | up or down possibility, position for snp |
| -e | insertion in genes, based on encode format |
| -r | random insertion location |
| -x | for insertion de novo, provide an other fasta. The absence of -x  option will lead to dispersed duplication simulation |
| -v | for insertion de novo, provide bed file of the other fasta, insertion will come from the region provided by the bed file |

### Usage :
    python3 Novo_dup_generator.py -g reference.fasta -b exon_list.bed -n 200 -s 250 -o 200_insertion_250bp -x other_genome.fasta -v exon_other_genome.bed

## Simulating other types of insertions
### Parameters : 
| Option | Description |
| ------------ | ----------------------------------------- |
| -v | vcf file, position from the vcf will be used to generate new insertion type |
| -g | reference genome, use the same version that the vcf file |
| -m | mobile element database in fasta format |
| -t | generate tandem duplications |
| -r | generate tandem repeats |
| -u | seed size for tandem repeats (6 by default) |
| -me | generate mobile elements |
| -s | size of the insertion to generate. In case of mobile element simulation, a size +/- 50 will be used to find suitable mobile elements in the database |

### Usage :
    python3 New_Simulation_type.py -v position_use_to_simulate.vcf -g reference_genome.fa -m mobile_element.fa -t True -s 500

##  Simulating microhomology
### Parameters :
| Option | Description |
| ------------ | ----------------------------------------- |
| -v | vcf reference file, precise location and resolved sequence in ALT column are required |
| -g | reference genome used to obtain the vcf file |
| -m | microhomology size |
| -t | name for the output file |

### Usage :
    python3 Microhomology.py -v vcf_file.vcf -g reference_genome -m 20 -t genome_20mh

Each simulation provides a vcf file and an altered genome. 
Insertion are simulated with an homozyguous genotype (1/1).
One vcf file and one altered genome are generated per insertion type. Thus multiple insertion types are not simulated in a same file.


## Benchmarking
### Parameters :
| Option | Description |
| ------------ | ----------------------------------------- |
| -t | truth ground vcf file (not compressed format) |
| -m | microhomology, windows size to validate an insertion (+/- microhomology size provided) |
| -v | vcf file to compare |
| -f | filter, (PASS and LowAss by default), only use it to not filter by quality |
| -c | compressing format, indicate True if the vcf to compare is gzipped |

### Usage :
    python3 Benchmark_SVcaller.py -t truth.vcf -m 10 -v tested.vcf.gz -c True

The output file provides three informations : the recall based on insertion site locations, the recall based on insertion sites + inserted sequences and the amount of false positive calls.
