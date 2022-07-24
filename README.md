# PloidyFrost
reference-free estimation of ploidy level from whole genome sequencing data based on de Bruijn graph 



## Install

```sh

$ git clone https://github.com/BlueBerrySun/PloidyFrost.git
$ cd PloidyFrost && mkdir build && cd build
$ cmake ..
$ make

```

## Usage


### Compute Lower coverage threshold
```
Usage: PloidyFrost cutoffL
PloidyFrost cutoffL kmer_histogram_file
```
### Compute Upper coverage threshold
```
Usage: PloidyFrost cutoffU
PloidyFrost cutoffU kmer_histogram_file [quantile[<1 ,default:0.998]]
```
### Run algorithms (Superbubble detection-Variant calling-Variant ranking)
```
Usage: PloidyFrost -g <BifrostGraph> -d <KMCDatabase> -o <outfile_prefix> ...
parameters with required argument:

  -g,             Input Bifrost Graph file (GFA format)
  -f,             Input Bifrost color file (BFG_COLORS format,colored CDBG)
  -o,             Prefix for Output files (default : 'output')
  -t,             Number of Threads (default is 1)
  -d,             Load KMC Database
  -h,             kmer Histogram file
  -l,             Lower coverage threshold (default : 10 )
  -u,             Upper coverage threshold (default : 1000 )
  -C,             Input Coverage thresholds file (colored CDBG)
  -z,             Maximum number of unitigs in superbubble (default : 8 )
  -M,             Match score (default : 2 )
  -D,             Mismatch penalty (default : -1 )
  -G,             Gap penalty (default : -3 )
  
parameters with no argument:

  -v,             Print information messages during construction
  -i,             Output Information about Bifrost graph
```
### Compute GMM log-likelihood to estimate ploidy level
```
Usage: PloidyFrost model
GMM model
  -f,             Prefix of coverage files
  -g,             Allele frequency file
  -l,             Minimum ploidy level (default : 2 )
  -u,             Maximum ploidy level (default : 10 )
  -q,             Minimum allele frequency
  -m,             Weigth minimum threshold ( 1/(p-1)/m , default value of m is 5)
  -n,             Weigth minimum threshold ( first_weight > maximum_weight/(p-1)/n , default value of n is 2)
  -k,             Maximum iterations
  -a,             Maximum delta
  -o,             Output prefix
```

### Generate allelic frequency distribution histogram
```
Usage: Rscript Drawfreq.R

        -f FILE, --file=FILE
                allelic frequency file

        -t TITLE, --title=TITLE
                histogram title

        -o OUTPREFIX, --outprefix=OUTPREFIX
                histogram output prefix

        -p PLOIDY, --ploidy=PLOIDY
                ploidy level

        -h, --help
                Show this help message and exit
```

### Filter the result using variant information (single sample)

```
Usage:  Rscript Filter.R

        -S, --simple
                only simple bubble

        -o OUTPREFIX, --outprefix=OUTPREFIX
                output prefix

        -i INPREFIX, --inprefix=INPREFIX
                input prefix

        -l LOW, --low=LOW
                lower coverage cutoff value

        -u UP, --up=UP
                upper coverage cutoff value

        -I, --indel
                filter indel

        -P, --snp
                filter snp

        -n NUM, --num=NUM
                VarNum cutoff value

        -d DISTANCE, --distance=DISTANCE
                VarDistance cutoff value

        -s SIZE, --size=SIZE
                VarSize cutoff value

        -q FREQUENCY, --frequency=FREQUENCY
                frequency range(frequency,1-frequency)

        -h, --help
                Show this help message and exit
```
### Filter the result using variant information (multiple samples)

```
Usage:  Rscript Filter-multi.R 

        -S, --simple
                only simple bubble

        -o OUTPREFIX, --outprefix=OUTPREFIX
                output prefix

        -c COLOR, --color=COLOR
                sample color

        -i INPREFIX, --inprefix=INPREFIX
                input prefix

        -l LOW, --low=LOW
                lower coverage cutoff value

        -u UP, --up=UP
                upper coverage cutoff value

        -I, --indel
                filter indel

        -P, --snp
                filter snp

        -n NUM, --num=NUM
                VarNum cutoff value

        -d DISTANCE, --distance=DISTANCE
                VarDistance cutoff value

        -s SIZE, --size=SIZE
                VarSize cutoff value

        -q FREQUENCY, --frequency=FREQUENCY
                frequency range(frequency,1-frequency)

        -v CRAMER, --cramer=CRAMER
                Cramer'V

        -h, --help
                Show this help message and exit
```
        
## Ploidy level estimation pipeline
### 1. Remove adapters and filter low-quality bases ([Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic))
### 2. Construct k-mer database(s) and mask low-coverage k-mers ([KMC](https://github.com/refresh-bio/KMC))
  
  ```sh
  #perform the following steps for each sample
  $ mkdir kmc_tmp
  #construct k-mer database for each sample
  $ kmc -ci1 -cs10000 -k25 -t${thread} sample.fq kmc_db kmc_tmp
  #generate k-mer coverage histogram
  $ kmc_tools transform kmc_db histogram hist
  #compute Lower coverage threshold (or set manually)
  $ lower_threshold=$(PloidyFrost cutoffL hist)
  #mask low-coverage k-mers in reads
  $ kmc_tools -t${thread} filter -hm kmc_db sample.fq -ci${lower_threshold} sample_filtered.fq 
  ```
### 3. Construct compacted (C)DBG ([Bifrost](https://github.com/pmelsted/bifrost))

  - **single sample**
      ```sh
      Bifrost build -i -d -k 25 -v -r sample_filtered.fq -o dbg -t ${thread}
      ```
  - **multiple samples**
      ```sh
      Bifrost build -c -i -d -k 25 -v -r sample1_filtered.fq -r sample2_filtered.fq -r ... -o cdbg -t ${thread}
      ```
### 4. Execute ploidy estimation algorithms

  - **single sample**
      ```sh
      #compute Lower coverage threshold (or set manually)
      $ lower_threshold=$(PloidyFrost cutoffL hist) 
      #compute Upper coverage threshold (or set manually)
      $ upper_threshold=$(PloidyFrost cutoffU hist 0.998) 
      #A: set threshold values manually
      $ PloidyFrost -g dbg.gfa -d kmc_db -t ${thread} -v -o single -l ${lower_threshold} -u ${upper_threshold}
      #B: give the k-mer histogram file name
      $ PloidyFrost -g dbg.gfa -d kmc_db -t ${thread} -v -o single -h hist
      ```
      **Output files**
      ```
      PloidyFrost_output/
      |-- single_alignseq.txt         | BubbleId | BubbleType | EntranceId | ExitId | AlignedSeq |
      |-- single_allele_frequency.txt | AlleleFrequency |
      |-- single_bicov.txt            | Cov1 | Cov2 | BubbleType | VarType/IndelSize | BubbleId | VarNum | VarDis |
      |-- single_bifre.txt            | BiAlleleFrequency |      
      |-- single_tricov.txt           | Cov1 | Cov2 | Cov3 | BubbleType | VarType/IndelSize | BubbleId | VarNum | VarDis |
      |-- single_trifre.txt           | TriAlleleFrequency |      
      |-- single_tetracov.txt         | Cov1 | Cov2 | Cov3 | Cov4 | BubbleType | VarType/IndelSize | BubbleId | VarNum | VarDis |
      |-- single_tetrafre.txt         | TetraAlleleFrequency |    
      |-- single_pentacov.txt         | Cov1 | Cov2 | Cov3 | Cov4 | Cov5 | BubbleType | VarType/IndelSize | BubbleId | VarNum | VarDis |
      |-- single_pentafre.txt         | PentaAlleleFrequency |    
      |-- single_super_bubble.txt     | BubbleId | EntranceId | EntranceStrand | ExitId | BubbleType | BubbleIsComplex |
      `-- single_Unitig_Id.txt        | UnitigId | UnitigSeq |
      ```
  - **multiple samples**
  
      **1. write the names of all k-mer databases in a text file named kmc_db_file**

      kmc_db_file
      ```
      kmc_db1
      kmc_db2
      kmc_db3
      ...
      ```

      **2(A). set lower and upper thresholds for each sample manually**

      coverage_file
      ```
      lower_threshold1<tab>upper_threshold1
      lower_threshold2<tab>upper_threshold2
      lower_threshold3<tab>upper_threshold3
      ...
      ``` 

      ```sh
      #A: set threshold values manually
      PloidyFrost -g cdbg.gfa -f cdbg.bfg_colors -d kmc_db_file -t ${thread} -v -o multi -C coverage_file 
      ```  
      **2(B). write the names of all k-mer histogram files in a text file named hist_file**

      hist_file
      ```
      hist1
      hist2
      hist3
      ...
      ```  
      ```sh
      #B: give the k-mer histogram file name
      PloidyFrost -g cdbg.gfa -f cdbg.bfg_colors -d kmc_db_file -t ${thread} -v -o multi -h hist_file 

      ```
      **Output files**
      ```
      PloidyFrost_output/
      |-- multi_alignseq.txt         | BubbleId | BubbleType | EntranceId | ExitId | AlignedSeq |
      |-- multi_allele_frequency.txt | AlleleFrequency |
      |-- multi_bicov.txt            | Cov1 | Cov2 | Color | BubbleType | VarType/IndelSize | BubbleId | VarNum | Cramer's V | VarDis |
      |-- multi_bifre.txt            | BiAlleleFrequency |      
      |-- multi_tricov.txt           | Cov1 | Cov2 | Cov3 | Color | BubbleType | VarType/IndelSize | BubbleId | VarNum | Cramer's V | VarDis |
      |-- multi_trifre.txt           | TriAlleleFrequency |      
      |-- multi_tetracov.txt         | Cov1 | Cov2 | Cov3 | Cov4 | Color | BubbleType | VarType/IndelSize | BubbleId | VarNum | Cramer's V | VarDis |
      |-- multi_tetrafre.txt         | TetraAlleleFrequency |    
      |-- multi_pentacov.txt         | Cov1 | Cov2 | Cov3 | Cov4 | Cov5 | Color | BubbleType | VarType/IndelSize | BubbleId | VarNum | Cramer's V | VarDis |
      |-- multi_pentafre.txt         | PentaAlleleFrequency |    
      |-- multi_super_bubble.txt     | BubbleId | EntranceId | EntranceStrand | ExitId | BubbleType | BubbleIsComplex |
      `-- multi_Unitig_Id.txt        | UnitigId | UnitigSeq |
      ```

### 5. Filter the result using variant information
  - **single sample**
  ```sh
    Rscript Filter.R -i single -o single-filtered -n 6 -s 11 -q 0.05 ...
  ```
  - **multiple samples**
   ```sh
    Rscript Filter-multi.R -i multi -o multi-filtered -n 6 -s 11 -q 0.05 -c 0 -v 0.25 ...
   ```
### 6. Generate visual and quantitative results
  - **log-likelihood of Gaussian Mixture Model**
   ```sh
   PloidyFrost model -g single-filtered_allele_frequency.txt -l 2 -u 10 -q 0.05 -o gmm
   ```
  - **allelic frequency distribution histogram**
   ```sh
   Rscript Drawfreq.R -f single-filtered_allele_frequency.txt -t title -p 2 -o histogram
   ```
   


