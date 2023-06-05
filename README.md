# Universal Pipeline for Viruses
## Configuration:
* Clone this repository.
* In order to install all dependencies using conda, run: \
    `conda env create -f <your/upv/path>/env/environment.yml`
## Usage:
 `conda activate upv` \
 `python3 upv.py [-h] | [-r REFSEQ_PATH] [-i FASTQ_PATH] [options] `
#### -i /input/path : provide input path  (required)
`/input/path` - the path to fastq files.

#### -r|--refseq <refseq/path/> : provide refseq path (required)
Don't worry about indexing the fasta file, it happens automatically.

 `python3 upv.py -i /input/path -r path/to/ref.fa`

#### -p|--process [INT] (optional)
Number of processes to run in parallel. \
Note that only some of the pipeline will run in multiprocessing. \
Be careful not to run too many processes on your CPU. \
 `python3 upv.py -i /input/path -r path/to/ref.fa -p 3`

#### --gb_file (optional)
Parse a gene bank file. \
Gene bank file can be downloaded from ncbi. \
The result is a regions file in the same path of the gb file, containing the genes region. \
 `python3 upv.py -gb path/to/file.gb`

#### -m |--mutations_table (optional)
mutations table analysis. \
Gene bank flag is required. \
The result is a mutation report file containing the nucleotide mutations and their translation according to the gene regions in the gb file.
 `python3 upv.py -i /input/path -r path/to/ref.fa -gb path/to/file.gb --mutation_table`
 
#### --mini (optional)
This flag enables to run a mini version of the pipeline by inserting alignment file as an input . \
The results are the parsed gb file and the mutations report. \
 `python3 upv.py -i /input/path -r path/to/ref.fa -gb path/to/file.gb --mutation_table --mini`

#### -v|--vcf (optional)
Generate a variant calling file using gatk4. \
This step fill the identity column in QC/report.csv file by counting the mutations in the vcf.

#### --flu (optional)
Run the pipeline for influenza viruses. \
The reference should contain the virus segments seperated by fasta headers. \
 `python3 upv.py  -r path/to/influenza_ref.fa -i path/to/fastq/location --flu `

#### --HIV (optional)
Run the pipeline for HIV analysis. \
NO NEED TO PROVIDE REFERENCE SEQUENCE. \
This analysis considers only the genes: PR, RT and Integrase. their reference is in the pipeline. \
 `python3 upv.py -i path/to/fastq/location --HIV `

#### --polio (optional)
Run the pipeline for polio virus. \
The reference should contain the 3 Sabins seperated by fasta headers. \
 `python3 upv.py  -r path/to/polio_refs.fa -i path/to/fastq/location --polio `
 
#### --de_novo (optional)
Run de-novo analysis in addition to the regular pipeline. \
The reference should be a Blast Database. The pipeline will choose a reference by the best match after using Blast.
The output will be seperated to contig based and fastq based analysis. \
 `python3 upv.py  -r path/to/blast_database.fa -i path/to/fastq/location --de_novo ` \
In order to generate blast database from a fasta/multi-fasta file you can use the following command: \
`makeblastdb -in fasta/file -parse_seqids -title "Viral" -dbtype nucl`

#### --CMV (optional)
Run CMV resistance mutations analysis in addition to the regular pipeline . \
This part was developed for human Herpes virus 5 reference - NC_006273.2. \
The output "cmv_resistance.txt" is a report containing the resistance mutations of each sample if exist. \
 `python3 upv.py  -r path/to/NC_006273.2.fa -i path/to/fastq/location --CMV `

 
## Main Steps:

### General pipeline:

* Index the reference sequence. \
`bwa index /refseq/path`

* For each fastq sample - map to reference and keep as bam file. \
For each fastq file run: \
    `bwa mem /refseq/path sample_R1 sample_R2 | samtools view -b - > BAM/sample_name.bam` \
R1 and R2: forward and reverse paired ends.
This command execute in parallel for 5 fastq files.

* Keep only mapped reads of each sample. \
For each bam file run: \
    `samtools view -b -F 260 file_name.bam > BAM/file_name.mapped.bam` 

* Sort and index each sample in BAM. \
For each mapped bam file run:
    `samtools sort BAM/file_name.mapped.bam BAM/file_name.mapped.sorted.bam`
    
* Create consensus files for each sample. \
For each mapped and sorted bam file run:
    `samtools mpileup -A  BAM/file_name.mapped.sorted.bam | ivar consensus -m 5 -p CNS5/file_name`

* Align all consensus sequences and references sequence: 
Gather all fasta sequences (consensus + reference): 
    `cat CNS5/*.fasta /refseq/path > alignment/all_not_aligned.fasta` 
Align with augur (mafft based) and save output in alignment/ directory: 
  ` augur align 
  --sequences alignment/all_not_aligned.fasta 
  --reference-sequence /refseq/path 
  --output alignment/all_aligned.fasta` 

* report.csv: 
    some coverage statistics and number of mapped reads: `samtools coverage -H BAM/sample_name.mapped.sorted.bam`  
    total number of reads in sample: `samtools view -c BAM/sample_name.bam` 
    depth of each position in sample: `samtools depth BAM/sample_name.mapped.sorted.bam`  
Check out pipeline.sh code for specifics.

- parallel commands are executing using python multiprocessing package. 

### Influenza addional steps:

* Splitting the BAM files by the references they were mapped to to multiple files. further analysis is made using the splitted files. 
for each bam file run: 
`bamtools split -in bam_file -reference` 

* concatenate the consensus files by sample. 

* concatenate the consensus files by reference segment. 

* multiple sequnce alignment: 
concatenate the reference's segments into one sequence with no headers. 
concatenate the sample's segmented consensus into one multi-fasta. 
run augur align. 

### de-novo additional steps:

* de-novo assembly using Spades. 
`spades -1 R1.fastq -2 R2.fastq -o %(output_path)s --rna` 

* Blast: 
Run Blast for each sample usind the contigs from Spades output ("transcripts.fasta"). 
`blastn -db blast_database -query transcripts.fasta -out output_file -outfmt \"6 qseqid sseqid stitle qlen qcovs score bitscore \ `

* Choose reference: 
The referece selected reference of each sample is the reference of the longest highest score contig based on blast output.

* Farther analysis is the analysis of the general pipeline executed seperatly contig based and fastq based. 
contig based analysis using Spades output: "transcripts.fasta". 


### polio additional steps: 

* filter non polio reads and reads that are mapped to more than one Sabin.
`bwa mem -v1 -t 32 references R1.fastq R2.fastq | samtools view -@ 32 -b -f 2 - | samtools fastq > sample.fastq`

* bwa mem mapping to each Sabin seperatly.

`bwa mem -v1 -t 32 references R1.fastq R2.fastq | samtools view -@ 32 -b -f 2 - | samtools fastq > sample.fastq`

### CMV additional steps: 

* Based on litrature cytomegalovirus (Human Herpesvirus 5) has known resistance mutations in the genes UL97 and UL54 (cmv_hotspots.csv). 
* Align (MAFFT) the samples to the reference sequence.

`augur align --sequences all_not_align.fasta --reference-sequence reference --output alignment/all_aligned.fasta`

* Translate the resistance mutations positions.
* Compare to the reference.
* Output "CMV_resistance.txt" report.
