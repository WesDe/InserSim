# InserSim
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)
## Description 
InserSim is an ensemble of script to simulate insertion in a reference genome and to benchmark insertion from SV caller.

Simulation folder contains three scripts : Novo_dup_generator.py to simulate de novo or dispersed duplication insertion, New_Simulation_type.py to simulate mobile element, tandem repeat, tandem duplication insertion and Microhomology.py to simulate insertion containing microhomology.

## Requirements
- Python3
- biopython (pip3 install biopython)

## Installation
git clone https://github.com/WesDe/InserSim.git

## Novel and dispersed duplication simulation usage
### Parameters :
| Option | Description |
| -g | reference genome in FASTA format | 
| -b | bed file, to generate insertion in specific regions |
| -n : number of insertion (int) |
| -s | size of insertion (int) |
| -p | induced snp near insertion (True/false, false by default) |
| -o | output file name |
| -c | induced insertion in specific chromosome |
| -a | up or down possibility, position for snp |
| -e | insertion in genes, based on encode format |
| -r | random insertion location |
| -x | for insertion de novo, provide an other fasta. The absence of -x  option will lead to dispersed duplication simulation |
| -v | for insertion de novo, provide bed file of the other fasta, insertion will come from the region provided by the bed file |

### How to use :
python3 Novo_dup_generator.py -g reference.fasta -b exon_list.bed -n 200 -s 250 -o 200_insertion_250bp -x other_genome.fasta -v exon_other_genome.bed

## Other type simulation usage
### Parameters : 
| Option | Description |
| -v | vcf file, position from the vcf will be used to generate new insertion type |
| -g | reference genome, use the same version that the vcf file |
| -m | mobile element database in fasta format |
| -t | generate tandem duplication |
| -r | generate tandem repeat |
| -u | seed size for tandem repeat (6 by default) |
| -me | generate mobile element |
| -s | size of the insertio to generate. In case of mobile element simulation, a size +/- 50 will be used to found potential mobile element in the database |

### How to use :
python3 New_Simulation_type.py -v position_use_to_simulate.vcf -g reference_genome.fa -m mobile_element.fa -t True -s 500

##  Microhomology simulation usage
### Parameters :
| Option | Description |
| -v | vcf reference file, location and sequence resolved in ALT column required |
| -g | reference genome used to obtain the vcf file |
| -m | micrhomology size |
| -t | name for the output file |

### How to use :
python3 Microhomology.py -v vcf_file.vcf -g reference_genome -m 20 -t genome_20mh

Each simulation provide a vcf file and an altered genome. 
Insertion simulated are homozyguous (1/1).
One vcf file and one altered genome is generated per insertion type. Thus multiple insertion type are not simulated in a same file.


## Benchmark
### Parameters :
| Option | Description |
| -t | truth ground vcf file (not compressed format) |
| -m | microhomology, windows size to validate an insertion (+/- microhomology size provided) |
| -v | vcf file to compare |
| -f | filter, (PASS and LowAss by default), only use it to not filter by quality |
| -c | compressing format, indicate True if the vcf to compare is gzipped |

### How to use :
python3 Benchmark_SVcaller.py -t truth.vcf -m 10 -v tested.vcf.gz -c True

The benchmark provide Three informations : recall based on site location, recall based on site + sequence resolved and the amount of false positive discover.