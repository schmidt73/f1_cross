# F1 Cross Analysis

The F1 cross analysis attempts to find sets of guides that are
specific to each genome. The difficult part in the standard analysis
that I performed was not generating these sets, but finding the
corresponding witnesses. So I had to take another approach.

## Step 1: Generating synthetic genome and mapping using MMARGE

I took the set of indels/SNPS (as a VCF file) and generated
a synthetic 129S1 genome along with an associated mapping file
using the tool called MMARGE. The variant calls are drawn
from Ensemble's mm10 genome, so the resultant genome is 

The script `prepare_genome.sh` in the scripts directory will 
perform this entire procedure from scratch given that the MMARGE
tool is installed.

To put the genomes in a standard FASTA format, simply run:
```bash
$ cat 129S1_SVIMJ/*.fa > 129s1.fna
$ cat genomes/*.fa > mm10.fna
```

## Step 2: Generating complete DBs

We generate the sets of gRNAs for both 129S1 and B6 independently
using our `build_dbs` script that parallelizes execution of Guidescan2. 
The results of this generation are located under the experiment14 
directory.

## Step 3: Generating cross DBs

Now that we have both databases and their corresponding
off-targets, we run the 129S1 set of guides against the B6 genome
and vice versa. This will run more quickly, since most guides are
filtered out because of their perfect match in the other genome. To
perform this procedure, we run the `pull_kmers.py` script on the
SAM databases for mm10 and 129S1 to generate a kmer set for both
organisms.

```bash
$ python pull_kmers.py 129s1.sam > 129s1.kmers.csv
$ python pull_kmers.py mm10.sam > mm10.kmers.csv
```

With the kmer sets in hand, we run the 129S1 kmer set against mm10 and
the mm10 kmer set against 129S1 using Guidescan2 with the no
perfect matches option enabled.

## Post-processing: Annotating disruption events and off-targets

At this point, we have a total of 4 databases:
    - A mm10 database with guides targeting mm10
    - A mm10 specific database with guides targeting mm10 and not 129S1
    - A 129S1 database with guides targeting 129S1
    - A 129S1 specific database with guides targeting 129S1 and not mm10

The first two databases are relatively easy to work with as the guideRNA coordinates
are with respect to mm10, it is only the off-targets that are with respect to another
genome. For the latter two we need to convert sgRNA coordinates to mm10. We do this
using MMARGE as follows where `129s1_dir` is the directory in which the `prepare_genomes.sh`
script was executed in the first step.

```bash
$ $MMARGE shift -data_dir 129s1_dir/ -ind 129S1_SVIMJ -sam -files 129s1_x_mm10.sam 
```

In the final step, we parse out off-target information and classify the type of
disruption to create two CSVs for mm10 and 129S1 allele specific gRNAs.

```bash
$ python scripts/classify_disruption.py 129s1_dir/inputs/129S1_SvImJ.mgp.v5.{snps,indels}.*.sorted.gz 129s1_x_mm10.coords_shifted_from_129S1_SVIMJ.sam --annotation-db /data/leslie/schmidth/genomes/grcm38/GRCm38_68.db --chromosome-mapping inputs/mm10_chr_mapping.txt --non-specific-db 129s1.bam.sorted -o 129s1_x_mm10.csv
$ python scripts/classify_disruption.py ../experiment15/inputs/129S1_SvImJ.mgp.v5.{snps,indels}.*.sorted.gz mm10_x_129s1.sam --annotation-db /data/leslie/schmidth/genomes/grcm38/GRCm38_68.db --chromosome-mapping inputs/mm10_chr_mapping.txt --non-specific-db mm10.bam.sorted -o mm10_x_129s1.csv
```

We then merge these CSVs using `awk`.
