#!/bin/zsh
conda activate bioinf
INPUTDIR="Fastq-files"
OUTPUTDIR="analysis"
HUMANMIRBASEINDEX="databases/miRBase"
HUMANREFINDEX="databases/human-Grch38-Ref"

# cutadapt parameters
MINLENGTH=18
MAXLENGTH=30

# mkdir $INPUTDIR/QCReports
# mkdir $OUTPUTDIR/maturemiRNAcounts
# mkdir $OUTPUTDIR/genome-miRNAcounts
alias fastqc='/Users/federico/Desktop/Università/COURSE\ -\ Bioinformatics/Progetto/fastqc/fastqc'
# alias bowtie='/Users/federico/Desktop/Università/COURSE\ -\ Bioinformatics/Progetto/bowtie-1.2.2-macos-x86_64/bowtie'
alias bowtie='/Users/federico/Desktop/Università/COURSE\ -\ Bioinformatics/Progetto/bowtie-1.2.2-macos-x86_64/bowtie-align-s --wrapper basic-0'
alias samtools='/usr/local/install/samtools/bin/samtools'
alias sed='gsed'
# for i in "BC31_001"
for i in 70459 70460 70461 70462 74665 74666 74667 74668 74669 74670 74671 74672 74673 74674 74675 74676 74677 74678 74679 74680 74681 74682 74683 74684
do
    echo "Running FastQC on raw sample $i now"
    fastqc -o $INPUTDIR/QCReports $INPUTDIR/$i".fastq.gz" -t 16
    mkdir $OUTPUTDIR/$i
    echo "Sample $i: Running of UMI extraction with regex on raw reads"
    umi_tools extract --stdin=$INPUTDIR/$i".fastq.gz" --log=$OUTPUTDIR/$i/$i"-UMIextraction-fromrawreads.log" --stdout=$OUTPUTDIR/$i/$i"-directUMIextracted.fastq" --extract-method=regex --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)'
    echo "Sample $i: Applying a read length filter to remove too short and long reads"
    cutadapt --minimum-length=$MINLENGTH --maximum-length=$MAXLENGTH -o $OUTPUTDIR/$i/$i"-directUMIextracted-min"$MINLENGTH"max"$MAXLENGTH"L.fastq" $OUTPUTDIR/$i/$i"-directUMIextracted.fastq" > $OUTPUTDIR/$i/$i"-directUMIextracted-readlengthfilter-cutadapt[min"$MINLENGTH"max"$MAXLENGTH"L].log"
    echo "Running of bowtie to mirbase mature seq"
    bowtie -n 0 -l 30 --norc --best --strata -m 1 --threads 16 mature $OUTPUTDIR/$i/$i"-directUMIextracted-min"$MINLENGTH"max"$MAXLENGTH"L.fastq" --un $OUTPUTDIR/$i/$i"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq" -S $OUTPUTDIR/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" 2> $OUTPUTDIR/$i/$i"-bowtie-maturemiRNA-beststratam1.log"
    echo "Sorting and indexing the BAM file for counting the mirna occurences"
    samtools sort $OUTPUTDIR/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.sam" > $OUTPUTDIR/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam"
    samtools index $OUTPUTDIR/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam"
    echo "Counting step"
    umi_tools count --method=unique --per-contig -I $OUTPUTDIR/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam" -L $OUTPUTDIR/$i/$i"_counts-uniquemethod-maturemiRNA.log" -S $OUTPUTDIR/$i/$i"_counts-finaloutput-uniquemethod-maturemiRNA.txt"
    echo "Deduplicating the aligned BAM"
    umi_tools dedup --method=unique -I $OUTPUTDIR/$i/$i"-maturemiRNA-aligned-bowtie1-beststratam1.bam" -S $OUTPUTDIR/$i/$i"_deduplicated-matureMirna-uniquemethod-beststratam1.bam" -L $OUTPUTDIR/$i/$i"-deduplicate-matureMirna-uniquemethod-beststratam1.log"
    echo "Indexing the BAM output for finding counts"
    samtools index $OUTPUTDIR/$i/$i"_deduplicated-matureMirna-uniquemethod-beststratam1.bam"
    echo "Generation of miRNA counts for file $i"
    samtools idxstats $OUTPUTDIR/$i/$i"_deduplicated-matureMirna-uniquemethod-beststratam1.bam" | cut -f1,3 - | sed "1s/^/miRNA\t${i}-miRNAcount\n/" - > $OUTPUTDIR/maturemiRNAcounts/$i"-maturemiRNA-counts.txt"
    echo "Running of bowtie from unaligned mature miRNA reads to genome as reference (Exiqon style)"
    bowtie -n 1 -l 32 --norc --best --strata -m 1 --threads 16 GCA_000001405.15_GRCh38_no_alt_analysis_set $OUTPUTDIR/$i/$i"-maturemiRNA-unalignedReads-bowtie1-beststratam1.fastq" --al $OUTPUTDIR/$i/$i"-genome-alignedReads-bowtie1-beststratam1.fastq" -S $OUTPUTDIR/$i/$i"-genomeaftermiRNA-aligned-bowtie1-beststratam1.sam" 2> $OUTPUTDIR/$i/$i"-bowtie-genomeaftermiRNA-beststratam1.log"
    echo "Sorting and indexing again for downstream analyses"
    samtools sort $OUTPUTDIR/$i/$i"-genomeaftermiRNA-aligned-bowtie1-beststratam1.sam" > $OUTPUTDIR/$i/$i"-genomeaftermiRNA-aligned-bowtie1-beststratam1.bam"
    samtools index $OUTPUTDIR/$i/$i"-genomeaftermiRNA-aligned-bowtie1-beststratam1.bam"
    # rm -rf $OUTPUTDIR/$i/$i"-genomeaftermiRNA-aligned-bowtie1-beststratam1.sam"
    echo "Deduplicating the BAM file before counting through custom scripts"
    umi_tools dedup --method=unique -I $OUTPUTDIR/$i/$i"-genomeaftermiRNA-aligned-bowtie1-beststratam1.bam" -S $OUTPUTDIR/$i/$i"_deduplicated-genomeaftermiRNA-uniquemethod.bam" -L $OUTPUTDIR/$i/$i"-deduplicate-genomeaftermiRNA-uniquemethod.log"
    echo "Index the deduplicated BAM file"
    samtools index $OUTPUTDIR/$i/$i"_deduplicated-genomeaftermiRNA-uniquemethod.bam"
    echo "Actual filtering sam step only for overlapping microRNA locations"
    tagBam -i $OUTPUTDIR/$i/$i"_deduplicated-genomeaftermiRNA-uniquemethod.bam" -files hsa-genome-miRBase22v-onlymiRNAs-convforTagBAM.bed -names -tag XQ > $OUTPUTDIR/$i/$i"-tagged.bam"
done
echo "Pipeline completed"