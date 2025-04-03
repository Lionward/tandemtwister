<!-- Anchor for Top of the README -->
<a name="readme-top"></a>
<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
1. - [Getting Started](#getting-started)
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]
<a href="link_to_demo_gif.gif"><strong>View Demo</strong></a>
  <h1>TandemTwister Tool</h1>
-->

<!-- PROJECT LOGO -->
<br />
<div align="center">
 <img width="443" alt="grafik" src="https://github.com/user-attachments/assets/77776774-3ff6-4b22-8d78-bcb09bf1ceb5">

  </p>
</div>
<br />
<div align="center">

  <p>TandemTwister is a fast tool for tandem repeat genotyping!</p>
  <p>
    ·
    <a href="https://github.com/Lionward/TandemTwist/issues">Report Bug</a>
    ·
    <a href="https://github.com/Lionward/TandemTwist/pulls">Request Feature</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
## Table of Contents
- [TandemTwister](#TandemTwister)
- [Requirements](#Requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Contributions](#contributions)
- [License](#license)
- [Upcoming Features](#UpcomingFeatures)



## Introducing TandemTwister 

<p style="font-family: Arial, sans-serif; font-size: 1.1em; line-height: 1.6; color: #333; margin: 1em 0;">
  TandemTwister is a user-friendly tool for genotyping tandem repeats that can handle
  long-read data from various technologies — like <strong>CLR</strong>, <strong>CCS</strong>,
  and ONT — as well as <strong>Somatic</strong> and <strong>aligned genomes</strong> as input.
</p>


## Key features
1. Versatile Compatibility: TandemTwister supports long-read sequencing data from CLR, CCS, and ONT technologies, ensuring adaptability to diverse genomic datasets.

2. Phasing Capabilities: The tool incorporates phasing algorithm, by leveraging distinctive features within TR regions.

3. Noise Correction for Short Motifs: TandemTwister includes specialized correction mechanisms for short motifs (≤3) in CLR and ONT reads, ensuring robust and accurate genotyping results in the presence of noisy data.

4. Speed and Scalability: Optimized for efficiency, TandemTwister supports multi-processing and can complete genotyping analyses for approximately 1.2 Mio regions in under 20 minutes.

    
<!-- Requirenments -->
## Requirements

TandemTwister requires the following tools:

- `g++` (GCC) version 14.2.0 or later
- `GNU Make` version 4.4 or later
- `conda` version 22.11.1 or later
  
Please ensure that you have these tools installed and available in your PATH before proceeding with the build process.

<!-- Installation -->
## Installation

1. start by cloning the Github repository
```bash
  git clone https://github.com/Lionward/TandemTwist.git
  cd TandemTwister
  ```

2. Create a new Conda environment:
```bash
   mamba create -n TandemTwist
   ```
3. Activate the newly created environment:
```bash
  mamba activate TandemTwist
  ```
3. Install libdeflate:
```bash
  mamba install -c conda-forge libdeflate=1.21
  ```
4. Install HTSLIB
```bash
   mamba install htslib=1.21
  ```
5. Install mlpack
```bash
   mamba  install mlpack==4.2.1
  ```
6. run make  to install TandemTwister in your /usr/local/bin
```bash
   make install
```
  Alternatively you can just run `make`, and the Image file will be created in the same directory.
 
<!-- Usage -->
# Usage

To use the TandemTwister run the desired commands within the activated environment.
<h1>TandemTwister</h1>

<p><strong>TandemTwister</strong> is a tool for genotyping tandem repeats from long reads and aligned genome input.</p>

<h2>Usage</h2>


<hr />

<h2>Required Input</h2>

<ol>
  <li>
    <strong>Analysis Type</strong>
    <ul>
      <li><code>--germline</code> / <code>--somatic</code> / <code>--assembly</code>: Type of analysis to perform.</li>
    </ul>
  </li>
  <li>
    <strong>Arguments</strong>
    <ul>
      <li><code>-b, --bam</code> &nbsp; Path to the BAM file of the aligned reads to the reference genome.</li>
      <li><code>-r, --ref</code> &nbsp; Path to the input reference file (e.g., .fa/.fna).</li>
      <li><code>-m, --motif_file</code> &nbsp; Path to the file containing reference coordinates and motif sequence (BED/TSV/CSV).</li>
      <li><code>-o, --output_file</code> &nbsp; Output file containing region, motif, hap1 and hap2 copy numbers.</li>
      <li><code>-s, --sex</code> &nbsp; Sample sex (0 = female, 1 = male).</li>
      <li><code>-sn, --sample</code> &nbsp; Name of the sample.</li>
      <li><code>-rt, --reads_type</code> &nbsp; Type of reads (Default: CCS).</li>
      <li><code>-bt, --bam_type</code> &nbsp; Type of BAM file (e.g., reads or assembly).</li>
    </ul>
  </li>
</ol>

<hr />

<h2>Dynamic Programming Alignment Parameters</h2>
<ul>
  <li><code>-mml, --min_match_ratio_l</code>: Minimum match ratio for long motifs (Default: 0.5)</li>
</ul>

<hr />

<h2>Germline &amp; Somatic Analysis Options</h2>

<h3>Read Extraction Parameters</h3>
<ul>
  <li><code>-h, --help</code>: Show help message</li>
  <li><code>-s, --output_file_statistics</code>: Output file containing phasing info &amp; consensus CN call for each region</li>
  <li><code>-pad, --padding</code>: Padding around the STR region to extract reads (Default: 0)</li>
  <li><code>-t, --threads</code>: Number of threads to use (Default: 1)</li>
  <li><code>-kpr, --keepPhasingResults</code>: Keep phasing results (Default: false)</li>
  <li><code>-kcr, --keepCutReads</code>: Keep cut reads (Default: false)</li>
  <li><code>-minR, --minReadsInRegion</code>: Minimum number of reads that should span the region (Default: 2)</li>
  <li><code>-btg, --bamIsTagged</code>: Reads in BAM are phased (Default: false)</li>
  <li><code>-qs, --quality_score</code>: Minimum quality score for a read to be considered (Default: 10, Max: 60)</li>
</ul>

<h3>Correction Parameters</h3>
<ul>
  <li><code>-cor, --correct</code>: Perform genotype calling correction (CCS Default: false, CLR/ONT Default: true)</li>
  <li><code>-crs, --consensus_ratio_str</code>: Min fraction of reads in a cluster that must call an interval (STR) (Default: 0.3)</li>
  <li><code>-crv, --consensus_ratio_vntr</code>: Min fraction of reads in a cluster that must call an interval (VNTR) (Default: 0.3)</li>
  <li><code>-roz, --removeOutliersZscore</code>: Remove outliers for phasing (Default: false)</li>
  <!-- <li><code>-rtr, --refineTrRegions</code>: Refine tandem repeat regions (Default: true)</li> -->
</ul>

<h3>Clustering Parameters</h3>
<ul>
  <li><code>-seps, --start_eps_str</code>: Start radian for clustering in STR regions (Default: 0.2)</li>
  <li><code>-sepv, --start_eps_vntr</code>: Start radian for clustering in VNTR regions (Default: 0.2)</li>
  <li><code>-minPF, --minPts_frac</code>: Min fraction of reads that should be in one cluster (Default: 0.3)</li>
  <li><code>-nls, --noise_limit_str</code>: Noise limit for clustering in STR regions (Default: 0.2)</li>
  <li><code>-nlv, --noise_limit_vntr</code>: Noise limit for clustering in VNTR regions (Default: 0.35)</li>
  <li><code>-ci, --cluster_iter</code>: Number of iterations for clustering (Default: 20)</li>
</ul>

<hr />

<h2>Help</h2>
<ul>
  <li><code>-h, --help</code>: Print help message</li>
</ul>

<hr />

<p><strong>Example:</strong></p>
<pre><code>
TandemTwister --germline \
  -b sample.bam \
  -m motifs.bed \
  -r reference.fna \
  -o output.txt \
  -s 1 \
  -sn SampleName \
  -rt CCS \
  -t 4
</code></pre>



## Example Output

Below is an example of the output in VCF format:

```vcf
chr1    60637   chr1:60636-60665        ATTGTAAAGTCAAACAATTATAAGTCAAAC  ATTGTAAAGTCAAACAATTATAAGTCAAAC  1       PASS    TR_type=VNTR;MOTIFS=AATTATAAGTCAAAC,ATTGTAAAGTCAAAC,AATTGTAAGTCAAAC,TTGTAAAGTCAAAC,AATTATAAGTCAAA;MOTIF_IDs_REF=1,0;MOTIF_IDs_H1=1,0;MOTIF_IDs_H2;UNIT_LENGTH=14;CN_H1=2;CN_H2=2;CN_ref=2   GT:DP:SP        0/1:31,0:(1-15)_(16-30),(1-15)_(16-30),(1-15)_(16-30)

```
<!-- Contributions -->

## test data
The test data is available in the test_data folder. The test data is in the form of a bam file, a reference file, and a motif file. The motif file is in the form of a bed file. The test data can be used to test the TandemTwister tool.
run the following command to test the tool. NOTICE: For running please add the reference file to the test folder and name it hg38.fa

```bash
  make test
```

# Contributions

We welcome contributions from the community. If you find any issues or have suggestions for improvement, please open an issue or create a pull request.


# UpcomingFeatures
1.  Implementation of a Lookup Table for ONT Input Acceleration:
      Integrate a lookup table for ONT input to enhance processing speed, optimizing the tool's performance.

2.  Inclusion of Methylation Information:
      Integrate methylation information into the analysis, providing users with additional insights into the epigenetic characteristics of the tandem repeats.


<!-- Acknowledgements -->

# Acknowledgements

We would like to express our appreciation to the IT team at Max Planck Institute for Molecular Genetics for their support with technical aspects related to this project.

<!-- Anchor for Bottom of the README -->

<a name="readme-bottom"></a>

