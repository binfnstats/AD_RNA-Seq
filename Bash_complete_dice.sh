#!/bin/bash

# for i in `cat ./cond1`; \
# do qsub -N des_${i} -l h_vmem=8G,mem_requested=8G,tmpfree=1024G,tmp_requested=1024G \
# -P ParkinsonsDiseaseandNeurodegeneration \
# -cwd -pe smp 6 -V -b y \
# bash ./scripts/complete.sh ${i} 6 ; done

echo "$tmp_requested $TMPDIR"


module load borgue/star/2.5.1a
module load gi/picard-tools/1.95
module load gi/samtools/1.2
module load gi/gcc/4.8.2
module load gi/novosort/1.02.02
module load fabbus/python/2.7.3
module load borgue/RSeQC/2.6.1
module load borgue/R/3.2.0
module load leemar/rsem/1.3.0
module load gi/ucsc_utils/283
module load gi/trimmomatic/0.32
module load marsmi/kallisto/0.42.4
module load gi/fastqc/0.11.3
module load borgue/stringtie/1.2.2 	
module load gi/bedtools/2.22.0
module load fabbus/perl/5.14.2



if [ -z $1 ]; then
	echo "Need to submit name for output files/folders" && exit
fi

if [ -z $2 ]; then
	echo "Need to provide number of processors to be used e.g. 6" && exit
fi

## Define project folder
outDir="/share/ScratchGeneral/borgue/pdrna"
echo $outDir

inDir="/share/ScratchGeneral/borgue/pdrna/fastq/${1}"
echo $inDir

## Cores best 6 so far
numcores=${2}
echo $numcores

## Annotations
indexDir="/share/ClusterShare/biodata/contrib/borgue/annotation/gencode_24/denovo2/"
echo "indexDir is $indexDir"

genomeDir=$indexDir"GRCh38.p5.genome_ercc.fa"
echo $genomeDir

gtfDir=$indexDir"gencode_24_ercc_denovo2.gtf"
echo $gtfDir

rseqcIndexDir=$indexDir"gencode_24_ercc_denovo2.gtf.bed"
echo $rseqcIndexDir

genSize=$indexDir"chrNameLength.txt"
echo $genSize

kalIndexDir=$indexDir"kallisto_index"
echo $kalIndexDir

rsemIndexDir=$indexDir"rsem"
echo $rsemIndexDir

##Global Internal file structure

fastQCDir=${outDir}"/fastQC/"${1}
echo $fastQCDir && mkdir -p $fastQCDir

fastQC_trim_Dir=${outDir}"/fastQC_trim/"${1}
echo $fastQC_trim_Dir && mkdir -p $fastQC_trim_Dir

trim_Dir=$outDir"/fastq/${1}_trim"
echo $trim_Dir && mkdir -p ${trim_Dir}

starDir=${outDir}"/star/"${1}
echo $starDir && mkdir -p ${starDir}

kalDir=${outDir}"/kallisto/"${1}
echo $kalDir && mkdir -p ${kalDir}

rsemDir=${outDir}"/rsem/"${1}
echo $rsemDir && mkdir -p ${rsemDir}

logDir="${outDir}/logs"
echo $logDir

resDir=${outDir}"/splicing/"${1}
echo $resDir && mkdir -p ${resDir}

juncDir=${resDir}"/rseqc"
echo $juncDir && mkdir -p ${juncDir}

read_len="115"
echo $read_len

gownDir=${outDir}"/stringtie/"${1} 		# add to global file structure
echo $gownDir && mkdir -p $gownDir

###Check if files are there

minFileSize="1M"

inFile1=${inDir}"/*R1.fq.gz"
echo $inFile1

inFile2=${inDir}"/*R2.fq.gz"
echo $inFile2

find $inFile1 -type f -size $minFileSize -delete
find $inFile2 -type f -size $minFileSize -delete

if [ ! -f $inFile1 ];
then
	echo "Can't input inFile1"
else
	echo "Found" $inFile1
fi

if [ ! -f $inFile2 ];
then
	echo "Can't input inFile2"
else
	echo "Found" $inFile2
fi

##########################################################################################################################
##########################################################################################################################
###File combination and potential subsetting
##########################################################################################################################
##########################################################################################################################

###Download from DNANexus
#source /share/ClusterShare/software/contrib/borgue/dx-toolkit/environment
#dx login
#dx select --level VIEW
#dx ls


### CHeck md5 checksums
# find -name "*md5" | parallel -j 8 --gnu "cd {//} && md5sum -c {/}"
# seq 40 123 > conditions
# for i in `cat conditions`; do mkdir ${i} ; done
# for i in `cat conditions`; do mv C8KKRANXX_*_151105_AZX--0${i}* ${i} ; done
# for i in `cat ../conditions1`; do mv C8KKRANXX_*_151105_AZX--${i}* ${i} ; done
# for i in `cat conditions1`; do mv ./${i}/*R1.fastq.gz ./${i}/${i}.R1.fq.gz ; done
# for i in `cat conditions`; do mv ./${i}/*R2.fastq.gz ./${i}/${i}.R2.fq.gz ; done
# for j in `cat ../conditions`;
# do
# 	echo ${j}
#     for i in `cat ../${j}`;
#     do
#     	echo ${i}
#         zcat ${i}/${i}.R2.fq.gz >> ${j}/${j}.R2.fq
#         gzip ${j}/${j}.R2.fq;
#     done;
# done &
# for i in `cat ../conditions`; do mv ${i}/${i}.R2.fq.gz ${i}/${i}.R2.fq ; done
# for i in `cat ../conditions`; do mv ${i}/${i}.R1.fq.gz ${i}/${i}.R1.fq ; done
# for i in `cat ../conditions`; do gzip ${i}/${i}.R2.fq ; done &
# for i in `cat ../conditions`; do gzip ${i}/${i}.R1.fq ; done &
# zcat ./fastq/${1}/${1}_L1_R1.fq.gz ./fastq/${1}/${1}_L2_R1.fq.gz | gzip -c > ./fastq/${1}/${1}_comb_R1.fq.gz 
# zcat ./fastq/${1}/${1}_L1_R2.fq.gz ./fastq/${1}/${1}_L2_R2.fq.gz | gzip -c > ./fastq/${1}/${1}_comb_R2.fq.gz

# module load gi/seqtk/1.0

# seqtk sample -s100 /share/Temp/borgue/harvard/fastq/HC_BN00-14_4_trimgalore/HC_BN00-14_4.R1.fastq 10000 > ./test_trimgalore/test.R1.fq &
# seqtk sample -s100 /share/Temp/borgue/harvard/fastq/HC_BN00-14_4_trimgalore/HC_BN00-14_4.R2.fastq 10000 > ./test_trimgalore/test.R2.fq &
 
# for gzip
# seqtk sample -s100 <(zcat ./sam/Sample_C02/C02_ATCACG_L001_R1_001.fastq.gz) 10000 > ./test/test_R1.fq
 



##########################################################################################################################
##########################################################################################################################
###FastQC preTrim
##########################################################################################################################
##########################################################################################################################


if [ ! -f $fastQCDir/${1}_R1.fq_fastqc.zip ];
then
	echo "Running fastQC preTrim"
	
	scp -r $inDir/*fq.gz $TMPDIR
 
 	inFile1=${TMPDIR}"/*R1.fq.gz"
	echo "Infile1"$inFile1

	inFile2=${TMPDIR}"/*R2.fq.gz"
	echo "Infile2"$inFile2
	
	fastqc_pretrim_tmp_Dir=$TMPDIR"_trim"
	echo $fastqc_pretrim_tmp_Dir && mkdir -p $fastqc_pretrim_tmp_Dir

	echo "fastqc -t $numcores --outdir $fastQCDir $inFile1  $inFile2"
	time fastqc \
	-t $numcores \
	--outdir $fastQCDir \
	$inFile1  $inFile2 
else
	echo "Found" $fastQCDir/${1}_R1.fq_fastqc.zip
fi


##########################################################################################################################
##########################################################################################################################
###Trimming
##########################################################################################################################
##########################################################################################################################

# ###Trimming using Trimgalore
# trim_galore \
# --length 20 --fastqc --phred33 --quality 20 \
# --paired --output_dir ${OUTDir} \
# ${INDir}/*comb_R1.fq.gz ${INDir}/*comb_R2.fq.gz



if [ ! -f $trim_Dir/${1}_R1_trim.fq.gz ];
then
	echo "Running trimmomatic"
	
	scp -r $inDir/*fq.gz $tmpDir
 
 	inFile1=${tmpDir}"/*R1.fq.gz"
	echo "Infile1"$inFile1

	inFile2=${tmpDir}"/*R2.fq.gz"
	echo "Infile2"$inFile2
	
	trim_tmp_Dir=$tmpDir"_trim"
	echo $trim_tmp_Dir && mkdir -p $trim_tmp_Dir

	echo "Running java -jar /share/ClusterShare/software/contrib/gi/trimmomatic/0.32/trimmomatic-0.32.jar PE -phred33 -threads ${numcores} \
	$inFile1 $inFile2 \
	$trim_tmp_Dir/${1}_R1_trim.fq.gz $trim_tmp_Dir/${1}_R1_trim_unpaired.fq.gz \
	$trim_tmp_Dir/${1}_R2_trim.fq.gz $trim_tmp_Dir/${1}_R2_trim_unpaired.fq.gz \
	ILLUMINACLIP:/share/ClusterShare/software/contrib/gi/trimmomatic/0.32/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:40"
	time java -jar /share/ClusterShare/software/contrib/gi/trimmomatic/0.32/trimmomatic-0.32.jar PE -phred33 -threads ${numcores} \
	$inFile1 $inFile2 \
	$trim_tmp_Dir/${1}_R1_trim.fq.gz $trim_tmp_Dir/${1}_R1_trim_unpaired.fq.gz \
	$trim_tmp_Dir/${1}_R2_trim.fq.gz $trim_tmp_Dir/${1}_R2_trim_unpaired.fq.gz \
	ILLUMINACLIP:/share/ClusterShare/software/contrib/gi/trimmomatic/0.32/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:40
	#HEADCROP:10
	mv $trim_tmp_Dir/* $trim_Dir
else
	echo "Found" $trim_Dir/${1}_R1_trim.fq.gz
fi

###Reset inFile to the trimmed fastq's

inFile1=$trim_Dir/${1}_R1_trim.fq.gz
echo "New infile1" $inFile1

inFile2=$trim_Dir/${1}_R2_trim.fq.gz
echo "New infile2" $inFile2



##########################################################################################################################
##########################################################################################################################
###Kallisto counting 
##########################################################################################################################
##########################################################################################################################


if [ ! -f $kalDir/abundance.tsv ];
then
	echo "Running Kallisto"

	scp -r $trim_Dir/* $tmpDir 

	inFile1=${tmpDir}"/*R1_trim.fq.gz"
	echo $inFile1

	inFile2=${tmpDir}"/*R2_trim.fq.gz"
	echo $inFile2
	
	kal_tmp_Dir=$tmpDir"kallisto/${1}"
	echo $kal_tmp_Dir && mkdir -p $kal_tmp_Dir

	echo "Running kallisto quant -i $kalIndexDir -b 100 <(zcat $inFile1) <(zcat $inFile2)"
	time kallisto quant \
	-i $kalIndexDir \
	-o ${kal_tmp_Dir} \
	-b 100 -t $numcores --bias \
	<(zcat $inFile1) \
	<(zcat $inFile2)

	mv $kal_tmp_Dir/* $kalDir
	sed '1,1d' $kalDir/abundance.tsv | (echo ${1}; awk '{print $4}') > $kalDir/${1}_kallisto_count_iso
else
	echo "Found" $kalDir/abundance.tsv
fi


##########################################################################################################################
##########################################################################################################################
###FastQC postTrim
##########################################################################################################################
##########################################################################################################################



	if [ ! -f $fastQC_trim_Dir/${1}_R1_trim.fq_fastqc.zip ];
	then
		echo "fastqc -t $numcores --outdir $fastQCDir $inFile1  $inFile2"
		time fastqc \
		-t $numcores \
		--outdir $fastQC_trim_Dir \
		$inFile1  $inFile2 
	else
		echo "Found" $fastQC_trim_Dir/${1}_R1_trim.fq_fastqc.zip
	fi



##########################################################################################################################
##########################################################################################################################
###STAR 
##########################################################################################################################
##########################################################################################################################


star_tmp_Dir=$tmpDir"_star"
echo $star_tmp_Dir && mkdir -p $star_tmp_Dir


find $starDir/${1}Aligned.sortedByCoord.out.bam -type f -size $minFileSize -delete
find ${starDir}/${1}Aligned.toTranscriptome.out.bam -type f -size $minFileSize -delete

if [ ! -f $starDir/${1}Aligned.sortedByCoord.out.bam ];
then
	echo "Running STAR"

	scp -r $trim_Dir/* $tmpDir 

	inFile1=${tmpDir}"/*R1_trim.fq.gz"
	echo "New infile1" $inFile1

	inFile2=${tmpDir}"/*R2_trim.fq.gz"
	echo "New infile2" $inFile2

	echo "STAR --runMode alignReads \
		--readFilesIn $inFile1 $inFile2 \
		--readFilesCommand zcat \
		--outFileNamePrefix ${star_tmp_Dir}/${1} \
		--genomeDir ${indexDir} \
		--sjdbGTFfile ${gtfDir} \
		--outSJfilterReads Unique \
		--sjdbOverhang 124 \
		--twopassMode Basic \
		--runThreadN ${numcores:=6} \
		--genomeLoad NoSharedMemory \
		--outFilterType BySJout \
		--outFilterMultimapNmax 100 \
		--outFilterMismatchNmax 33 \
		--outFilterMatchNminOverLread 0 \
		--outFilterMismatchNoverLmax 0.3 \
		--outFilterScoreMinOverLread 0.3 \
		--limitOutSJcollapsed 1000000 \
		--limitSjdbInsertNsj 1000000 \
		--alignSJoverhangMin 8 \
		--alignEndsType EndToEnd \
		--alignSJDBoverhangMin 3  \
		--alignSJoverhangMin 15 \
		--alignIntronMin 20 \
		--winAnchorMultimapNmax 50 \
		--seedSearchStartLmax 12 \
		--chimSegmentMin 20 \
		--outSAMattributes All \
		--outSAMstrandField intronMotif \
		--quantMode TranscriptomeSAM \
		--outSAMattrIHstart 0 \
		--outSAMunmapped Within \
		--outSAMtype BAM SortedByCoordinate"
	time STAR --runMode alignReads \
		--readFilesIn $inFile1 $inFile2 \
		--readFilesCommand zcat \
		--outFileNamePrefix ${star_tmp_Dir}/${1} \
		--genomeDir ${indexDir} \
		--sjdbGTFfile ${gtfDir} \
		--outSJfilterReads Unique \
		--sjdbOverhang 124 \
		--twopassMode Basic \
		--runThreadN ${numcores:=6} \
		--genomeLoad NoSharedMemory \
		--outFilterType BySJout \
		--outFilterMultimapNmax 100 \
		--outFilterMismatchNmax 33 \
		--outFilterMatchNminOverLread 0 \
		--outFilterMismatchNoverLmax 0.3 \
		--outFilterScoreMinOverLread 0.3 \
		--limitOutSJcollapsed 1000000 \
		--limitSjdbInsertNsj 1000000 \
		--alignSJoverhangMin 8 \
		--alignEndsType EndToEnd \
		--alignSJDBoverhangMin 3  \
		--alignSJoverhangMin 15 \
		--alignIntronMin 20 \
		--winAnchorMultimapNmax 50 \
		--seedSearchStartLmax 12 \
		--chimSegmentMin 20 \
		--outSAMattributes All \
		--outSAMstrandField intronMotif \
		--quantMode TranscriptomeSAM \
		--outSAMattrIHstart 0 \
		--outSAMunmapped Within \
		--outSAMtype BAM SortedByCoordinate

		mv -f $star_tmp_Dir/* $starDir
else
	echo "Found" $starDir/${1}Aligned.sortedByCoord.out.bam
fi

find $starDir/${1}Aligned.sortedByCoord.out.bam -type f -size $minFileSize -delete
find ${starDir}/${1}Aligned.toTranscriptome.out.bam -type f -size $minFileSize -delete

in_Bam=$starDir"/"${1}"Aligned.sortedByCoord.out.bam"
echo $in_Bam


##########################################################################################################################
##########################################################################################################################
### Homebrew
##########################################################################################################################
##########################################################################################################################
##Intial Filtering and creation of Introns bed

# if 	[ ! -f ${resDir}/${1}.spliced.bam ]; then
# 		echo "samtools view -h -@ ${numcores} $in_Bam | awk -F"\t" 'OFS="\t" ((/NH:i:1[\t$]/) && ($6 ~ /N/ || $1 ~ /^@/)) || NF<7' | samtools view -bS -@ ${numcores} - > $resDir/${1}.spliced.bam"
# 		samtools view -h -@ ${numcores} $in_Bam | awk -F"\t" 'OFS="\t" ((/NH:i:1[\t$]/) && ($6 ~ /N/ || $1 ~ /^@/)) || NF<7' | samtools view -bS -@ ${numcores} - > $resDir/${1}.spliced.bam
# 		echo "samtools view $resDir/${1}.spliced.bam | fgrep -vw NH:i:1 > $resDir/${1}testfgrep.txt"
# 		samtools view -@ ${numcores} $resDir/${1}.spliced.bam | fgrep -vw NH:i:1 > $resDir/${1}testfgrep.txt
# else
# 	echo "Found $resDir/${1}.spliced.bam"
# fi

# ###___INSERT FILTER___###
# if [ ! -s $resDir/${1}testfgrep.txt ]; then 
# 	echo "All Good"
# else
# 	echo "Contaminants Detected" && qdel $JOB_ID
# fi

if 	[ ! -f $resDir/${1}.spliced.bed ]; then
	bamToBed -bed12 -i $in_Bam |  awk ' $10 >1 {print $0 }'> $resDir/${1}.spliced.bed
else
	echo "Found $resDir/${1}.spliced.bed"
fi

if [ ! -f $resDir/${1}.8names.txt ]; then
	awk 'OFS="\t" {print $4, $11}' $resDir/${1}.spliced.bed  | sed 's/,/\t/g' | awk '$3 >=8' | awk '$2 >= 8 {print $1}' > $resDir/${1}.8names.txt
else
	echo "Found $resDir/${1}.8names.txt"
fi

if [ ! -f $resDir/${1}.filtered_8.bed ]; then
	LC_ALL=C grep -wF -f $resDir/${1}.8names.txt $resDir/${1}.spliced.bed > $resDir/${1}.filtered_8.bed
else
	echo "Found $resDir/${1}.filtered_8.bed"
fi

if [ ! -f $resDir/${1}.introns.bed ]; then	
	perl /home/borgue/scripts/bed2introns.pl $resDir/${1}.filtered_8.bed $resDir/${1}.introns.bed
else
	echo "Found $resDir/${1}.introns.bed"
fi

##RSeQC Analysis
if [ ! -f $juncDir/${1}filtered_8.bam ]; then
	bedToBam -i $resDir/${1}.filtered_8.bed -g ${genSize} -bed12 > $juncDir/${1}filtered_8.bam
else
	echo "Found $resDir/${1}filtered_8.bam"
fi

if [ ! -f $juncDir/${1}bam_stat.txt ]; then
        echo "Running bam_stat.py -i $juncDir/${1}filtered_8.bam > $juncDir/${1}bam_stat.txt 2>&1"
        bam_stat.py -i $juncDir/${1}filtered_8.bam > $juncDir/${1}bam_stat.txt 2>&1
else
        echo "Found" $juncDir/${1}bam_stat.txt
fi

if [ ! -f $juncDir/${1}read_dis.txt ]; then
        echo "Running read_distribution.py"
        read_distribution.py -r $refDir -i $juncDir/${1}filtered_8.bam > $juncDir/${1}read_dis.txt 2>&1
else
        echo "Found" $juncDir/${1}read_dis.txt
fi

if [ ! -f $juncDir/${1}.splice_junction.pdf ]; then
        echo "Running junction_annotation.py -r $refDir -i $juncDir/${1}filtered_8.bam --out-prefix $juncDir/${1} 2>&1"
        junction_annotation.py -r $refDir -i $juncDir/${1}filtered_8.bam --out-prefix $juncDir/${1} 2>&1
else
        echo "Found" $juncDir/${1}.splice_junction.pdf
fi

if [ ! -f $juncDir/${1}.junctionSaturation_plot.pdf ]; then
        echo "Running junction_saturation.py -r $refDir -i $juncDir/${1}filtered_8.bam --out-prefix $juncDir/${1}"
        junction_saturation.py -r $refDir -i $juncDir/${1}filtered_8.bam --out-prefix $juncDir/${1}
else
        echo "Found" $juncDir/${1}.junctionSaturation_plot.pdf
fi

if [ ! -f $juncDir/${1}.bw ]; then
        echo "Running bam2wig.py --skip-multi-hits -i $juncDir/${1}filtered_8.bam --out-prefix $juncDir/${1} --chromSize /share/Temp/borgue/harvard/rseqc/hg19.chrom.sizes"
        bam2wig.py --skip-multi-hits -i $juncDir/${1}filtered_8.bam --out-prefix $juncDir/${1} --chromSize $genSize #-d "1+-,1-+,2++,2--" --wigsum=TOTAL_WIGSUM
        rm $juncDir/${1}.wig
else
        echo "Found" $juncDir/${1}.bw
fi


##########################################################################################################################
##########################################################################################################################
### Stringtie
##########################################################################################################################
##########################################################################################################################

if [ ! -f $gownDir/${1}.gtf ];
then
        echo "Running stringtie"
        if [ ! -f $tmpDIR/${1}Aligned.sortedByCoord.out.bam ]; 	# re-work around $tmp_dir for new SGE tmp management
        then
           	scp -r ${in_Bam} $tmpDIR 						# where in_bam is aligned to genome $starDir"/"${1}"Aligned.sortedByCoord.out.bam"
           	echo "copying bam for stringtie run 1"
        else
        	echo "Found" 
        fi

        inFile=$tmpDIR"/${1}Aligned.sortedByCoord.out.bam"
        echo "New infile" $inFile

        echo "stringtie ${inFile} -G ${gtfDir} -o ${gownDir}.gtf -a 8 -p ${numcores} -v -e -b ${gownDir}"
        stringtie ${inFile} -G ${gtfDir} -o ${gownDir}.gtf -a 8 -p ${numcores} -v -b ${gownDir} 				# semi-denovo over the gtf annotation
else
        echo "Found"  $gownDir/${1}.gtf
fi


### For Stringite 2-pass, using the merged gtf of all samples form teh previous run as an annotation we must break for a merge of all gtfs

### Merge gtfs with stringtie --merge
ls ${gownDir}.gtf > gtfs.txt

stringtie --merge -G ${gtfDir} -o ${gownDir}_merge gtfs.txt
####write checkpoint for ${gownDir}_merge

# resubmit


##########################################################################################################################
##########################################################################################################################
### Stringtie-2-pass
##########################################################################################################################
##########################################################################################################################

# gtfDir=${gownDir}_merge 			# new gtfDir to direct towards merged assemblies
# gownDir_2=${outDIR}"/stringtie_2/${1}" 

# if [ ! -f $gownDir_2/${1}.gtf ];
# then
#         echo "Running stringtie"
#         if [! -f  $tmpDIR/${1}Aligned.sortedByCoord.out.bam] 	# re-work around $tmp_dir for new SGE tmp management
#         then
#            	scp -r ${starDir}/${1}Aligned.sortedByCoord.out.bam $tmpDIR 							# where in_bam is aligned to genome $starDir"/"${1}"Aligned.sortedByCoord.out.bam"
#            	echo "copying bam for stringtie run 1"
#         else
#         	echo "Found" 
#         fi

#         inFile=$tmpDIR"/${1}Aligned.sortedByCoord.out.bam"
#         echo "New infile" $inFile

#         echo "stringtie ${inFile} -G ${gtfDir} -o ${gownDir_2}.gtf -a 8 -p ${numcores} -v -e -b ${gownDir_2}" 		# -e option forces assembly to only match those transcripts in the merged annotation
#         stringtie ${inFile} -G ${gtfDir} -o ${gownDir_2}.gtf -a 8 -p ${numcores} -v -b ${gownDir_2} 				
# else
#         echo "Found"  $gownDir_2/${1}.gtf
# fi

##########################################################################################################################
##########################################################################################################################
###Indexing using SamTools 
##########################################################################################################################
##########################################################################################################################

if [ ! -f $in_Bam".bai" ]; 
then
	echo "Running samtools index $in_Bam"
	time samtools index $in_Bam
else
	echo "Found" $in_Bam".bai"
fi





##########################################################################################################################
##########################################################################################################################
###Sorting the transcritpome bam for RSEM
##########################################################################################################################
##########################################################################################################################

find ${starDir}/${1}.sorted.bam -type f -size $minFileSize -delete

if [ ! -f ${starDir}/${1}.sorted.bam ]; then

	echo "Sorting the transcriptome bam staroutput for RSEM"
	echo "scp -r ${starDir}/${1}Aligned.toTranscriptome.out.bam $star_tmp_Dir"
	scp -r ${starDir}/${1}Aligned.toTranscriptome.out.bam $star_tmp_Dir
	echo "samtools view -@ ${numcores:=6} $star_tmp_Dir/${1}Aligned.toTranscriptome.out.bam -f 3 -b > ${star_tmp_Dir}/${1}.out.bam"
	time samtools view -@ ${numcores:=6} $star_tmp_Dir/${1}Aligned.toTranscriptome.out.bam -f 3 -b > ${star_tmp_Dir}/${1}.out.bam
	echo "novosort -n -m 16G -c ${numcores:=6} ${star_tmp_Dir}/${1}.out.bam > ${star_tmp_Dir}/${1}.sorted.bam"
	time novosort -n -m 16G -c ${numcores:=6} ${star_tmp_Dir}/${1}.out.bam > ${star_tmp_Dir}/${1}.sorted.bam
	echo "scp ${star_tmp_Dir}/${1}.sorted.bam ${starDir}"
	scp ${star_tmp_Dir}/${1}.sorted.bam ${starDir}
	# rm ${star_tmp_Dir}/${1}.out.bam
else
	echo "Found" ${starDir}/${1}.sorted.bam
fi

find ${starDir}/${1}.sorted.bam -type f -size $minFileSize -delete

##########################################################################################################################
##########################################################################################################################
###Counting the transcritpome bam in RSEM
##########################################################################################################################
##########################################################################################################################

find $rsemDir/${1}.genes.results -type f -size $minFileSize -delete

if [ ! -f ${rsemDir}/${1}.genes.results ]; then
	echo "${starDir}/${1}.sorted.bam ${star_tmp_Dir}"
	scp ${starDir}/${1}.sorted.bam ${star_tmp_Dir}

	echo "time rsem-calculate-expression --paired-end --bam --forward-prob 0 --no-bam-output -p $numcores ${star_tmp_Dir}/${1}.sorted.bam ${rsemIndexDir} ${rsemDir}/${1}"
	time rsem-calculate-expression \
		--paired-end \
		--bam \
		--forward-prob 0 \
		--no-bam-output \
		-p $numcores \
		${star_tmp_Dir}/${1}.sorted.bam \
		${rsemIndexDir} \
		${rsemDir}/${1}
		(echo ${1}; awk '{print $5}' $rsemDir/${1}.isoforms.results) > $rsemDir/${1}_rsem_iso_count
		(echo ${1}; awk '{print $7}' $rsemDir/${1}.isoforms.results) > $rsemDir/${1}_rsem_iso_fpkm
		(echo ${1}; awk '{print $7}' $rsemDir/${1}.isoforms.results) > $rsemDir/${1}_rsem_iso_fpkm
		(echo ${1}; awk '{print $7}' $rsemDir/${1}.isoforms.results) > $rsemDir/${1}_rsem_iso_fpkm
		(echo ${1}; awk '{print $5}' $rsemDir/${1}.genes.results) > $rsemDir/${1}_rsem_gene_count
		(echo ${1}; awk '{print $7}' $rsemDir/${1}.genes.results) > $rsemDir/${1}_rsem_gene_fpkm
else
	echo "${rsemDir}/${1}.genes.results is already there."
fi

find $rsemDir/${1}.genes.results -type f -size $minFileSize -delete


##########################################################################################################################
##########################################################################################################################
## Running RSeQC
##########################################################################################################################
##########################################################################################################################


echo "Running RSeQC"

mkdir -p ${starDir}/rseqc

scp $in_Bam $star_tmp_Dir
scp $in_Bam".bai" $star_tmp_Dir
in_Bam=$star_tmp_Dir"/"${1}"Aligned.sortedByCoord.out.bam"



if [ ! -f ${starDir}/rseqc/${1}bam_stat.txt ]; then
	echo "Running bam_stat.py -i ${in_Bam} > ${starDir}/rseqc/${1}bam_stat.txt 2>&1"
	time bam_stat.py -i ${in_Bam} > ${starDir}/rseqc/${1}bam_stat.txt 2>&1
else
	echo "Found" ${starDir}/rseqc/${1}bam_stat.txt
fi

if [ ! -f ${starDir}/rseqc/${1}.clipping_profile.pdf ]; then
	echo "Running clipping_profile.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time clipping_profile.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.clipping_profile.pdf
fi

if [ ! -f ${starDir}/rseqc/${1}.deletion_profile.pdf ]; then
	echo "Running deletion_profile.py -i ${in_Bam} -l $read_len --out-prefix ${starDir}/rseqc/${1}"
	time deletion_profile.py -i ${in_Bam} -l $read_len --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.deletion_profile.pdf
fi

if [ ! -f ${starDir}/rseqc/${1}.geneBodyCoverage.curves.pdf ]; then
	echo "Running geneBody_coverage.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time geneBody_coverage.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/${1}.geneBodyCoverage.curves.pdf
fi

if [ ! -f ${starDir}/rseqc/${1}infer_exp.txt ]; then
	echo "Running infer_experiment.py"
	time infer_experiment.py -i ${in_Bam} -r ${rseqcIndexDir} > ${starDir}/rseqc/${1}infer_exp.txt
else
	echo "Found" ${starDir}/rseqc/${1}infer_exp.txt
fi

if [ ! -f ${starDir}/rseqc/${1}.inner_distance_plot.pdf ]; then
	echo "Running inner_distance.py"
	time inner_distance.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.inner_distance_plot.pdf
fi

if [ ! -f ${starDir}/${1}read_dis.txt ]; then
	echo "Running read_distribution.py"
	time read_distribution.py -r ${rseqcIndexDir} -i ${in_Bam} > ${starDir}/${1}read_dis.txt 2>&1
else
	echo "Found" ${starDir}/${1}read_dis.txt
fi

if [ ! -f ${starDir}/rseqc/${1}.insertion_profile.txt ]; then
	echo "Running insertion_profile.py"
	time insertion_profile.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -l $read_len
else
	echo "Found" ${starDir}/rseqc/${1}.insertion_profile.txt
fi

if [ ! -f ${starDir}/rseqc/${1}.insertion_profile.pdf ]; then
	echo "Running mismatch_profile.py"
	time mismatch_profile.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -l $read_len
else
	echo "Found" ${starDir}/rseqc/${1}.insertion_profile.pdf
fi

if [ ! -f ${starDir}/rseqc/${1}.DupRate_plot.pdf ]; then
	echo "Running read_duplication.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time read_duplication.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.DupRate_plot.pdf
fi

if [ ! -f ${starDir}/rseqc/${1}.GC_plot.pdf ]; then
	echo "Running read_GC.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time read_GC.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.GC_plot.pdf
fi

if [ ! -f ${starDir}/rseqc/${1}.NVC_plot.pdf ]; then
	echo "Running read_NVC.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time read_NVC.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.NVC_plot.pdf]
fi

if [ ! -f ${starDir}/rseqc/${1}.qual.heatmap.pdf ]; then
	echo "Running read_quality.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time read_quality.py -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.qual.heatmap.pdf
fi

# -d STRAND_RULE, --strand=STRAND_RULE How read(s) were stranded during sequencing. For example: --strand='1++,1--,2+-,2-+' \
# means that this is a pair-end, strand-specific RNA-seq, and the strand rule is: read1 mapped to '+' => parental gene on '+'; \
# read1 mapped to '-' => parental gene on '-'; read2 mapped to '+' => parental gene on '-'; read2 mapped to '-' => parental gene on '+'.
# If you are not sure about the strand rule, run 'infer_experiment.py'

## note build grep from 'infer_experiment.py'
##strand_rule=$(grep 'stranded' ${starDir}/${1}infer_exp.txt )

if [ ! -f ${starDir}/rseqc/${1}_read_count.xls ]; then
	echo "Running RPKM_count.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -d "1+-,1-+,2++,2--""
	time RPKM_count.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -d "1+-,1-+,2++,2--"
else
	echo "Found" ${starDir}/${1}_read_count.xls
fi

if [ ! -f ${starDir}/rseqc/${1}.eRPKM.xls ]; then
	echo "Running RPKM_saturation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -d "1+-,1-+,2++,2--""
	time RPKM_saturation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} -d "1+-,1-+,2++,2--"
else
	echo "Found" ${starDir}/rseqc/${1}.eRPKM.xls
fi

if [ ! -f ${starDir}/rseqc/${1}.splice_junction.pdf ]; then
	echo "Running junction_annotation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} 2>&1"
	time junction_annotation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} 2>&1
else
	echo "Found" ${starDir}/rseqc/${1}.splice_junction.pdf
fi

if [ ! -f ${starDir}/rseqc/${1}.junctionSaturation_plot.pdf ]; then
	echo "Running junction_saturation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}"
	time junction_saturation.py -r ${rseqcIndexDir} -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1}
else
	echo "Found" ${starDir}/rseqc/${1}.junctionSaturation_plot.pdf
fi

cd ${starDir}




if [ ! -f ${starDir}/rseqc/${1}Aligned.sortedByCoord.out.tin.xls ]; then
	echo "Running tin.py -r ${rseqcIndexDir} -i ${in_Bam}"
	time tin.py -r ${rseqcIndexDir} -i ${in_Bam}
	mv ${starDir}/${1}*.tin.xls ${starDir}/rseqc/
else
	echo "Found" ${starDir}/rseqc/${1}Aligned.sortedByCoord.out.tin.xls
fi


if [ ! -f ${starDir}/rseqc/${1}.Forward.bw ]; then
	echo "Running bam2wig.py --skip-multi-hits -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} --chromSize ${indexDir}/genome.chrom.sizes"
	time bam2wig.py --skip-multi-hits -i ${in_Bam} --out-prefix ${starDir}/rseqc/${1} --chromSize ${indexDir}chrNameLength.txt --wigsum=1000000000 --strand='1++,1--,2+-,2-+'
	rm ${starDir}/rseqc/${1}*.wig
else
	echo "Found" ${starDir}/rseqc/${1}.Forward.bw
fi


##########################################################################################################################
##########################################################################################################################
###GATK pipeline
##########################################################################################################################
##########################################################################################################################



# if [ ! -f ${starDir}/${1}final.vcf ]; then
# 	echo " java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.121/AddOrReplaceReadGroups.jar I=$starDir/${1}Aligned.sortedByCoord.out.bam O=${star_tmp_Dir}/${1}rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample"
# 	time java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.121/AddOrReplaceReadGroups.jar I=$starDir/${1}Aligned.sortedByCoord.out.bam O=${star_tmp_Dir}/${1}rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
# 	echo " java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.121/MarkDuplicates.jar I=${star_tmp_Dir}/${1}rg_added_sorted.bam O=${starDir}/${1}dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${starDir}/output.metrics"
# 	time java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.121/MarkDuplicates.jar I=${star_tmp_Dir}/${1}rg_added_sorted.bam O=${star_tmp_Dir}/${1}dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${starDir}/output.metrics
# 	rm ${starDir}/${1}rg_added_sorted.bam
# 	echo " Running java -jar /share/ClusterShare/software/contrib/pethum/gatk/prebuilt/3.4-46/jar/GenomeAnalysisTK.jar -T SplitNCigarReads -R $genomeDir -I ${star_tmp_Dir}/${1}dedupped.bam  -o ${star_tmp_Dir}/${1}split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS"
# 	time java -jar /share/ClusterShare/software/contrib/pethum/gatk/prebuilt/3.4-46/jar/GenomeAnalysisTK.jar -T SplitNCigarReads -R $genomeDir -I ${star_tmp_Dir}/${1}dedupped.bam  -o ${star_tmp_Dir}/${1}split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
# 	echo " Running java -jar /share/ClusterShare/software/contrib/pethum/gatk/prebuilt/3.4-46/jar/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genomeDir -I ${star_tmp_Dir}/${1}split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o ${star_tmp_Dir}/${1}output.vcf"
# 	time java -jar /share/ClusterShare/software/contrib/pethum/gatk/prebuilt/3.4-46/jar/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genomeDir -I ${star_tmp_Dir}/${1}split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o ${star_tmp_Dir}/${1}output.vcf
# 	echo " Running java -jar /share/ClusterShare/software/contrib/pethum/gatk/prebuilt/3.4-46/jar/GenomeAnalysisTK.jar -T VariantFiltration -R $genomeDir -V ${star_tmp_Dir}/${1}output.vcf -window 35 -cluster 3 -filterName FS -filter FS > 30.0 -filterName QD -filter QD < 2.0 -o ${starDir}/${1}final.vcf"
# 	time java -jar /share/ClusterShare/software/contrib/pethum/gatk/prebuilt/3.4-46/jar/GenomeAnalysisTK.jar -T VariantFiltration -R $genomeDir -V ${star_tmp_Dir}/${1}output.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${starDir}/${1}final.vcf
# 	echo "java -jar /share/ClusterShare/software/contrib/pethum/gatk/prebuilt/3.4-46/jar/GenomeAnalysisTK.jar -T ASEReadCounter -R $genomeDir -dt NONE --sitesVCFFiles ${starDir}/${1}final.vcf -I ${star_tmp_Dir}/${1}split.bam -drf DuplicateRead -o ${starDir}/${1}_ase.csv"
# 	time java -jar /share/ClusterShare/software/contrib/pethum/gatk/prebuilt/3.4-46/jar/GenomeAnalysisTK.jar -T ASEReadCounter -R $genomeDir -dt NONE --includeDeletions --sitesVCFFile ${starDir}/${1}final.vcf -I ${star_tmp_Dir}/${1}split.bam -drf DuplicateRead -o ${starDir}/${1}_ase.csv
# 	rm ${starDir}/${1}split.bam
# 	rm ${starDir}/${1}output.vcf
# else
# 	echo "Found" ${starDir}/${1}final.vcf "GATK is done."
# fi

##########################################################################################################################
##########################################################################################################################
###Cleanup
##########################################################################################################################
##########################################################################################################################


if [ ! -z $tmpDIR ]; then
	>&2 echo “Removing contents of directory “$tmpDIR
	rm -rf $tmpDIR/* 
else
	>&2 echo “[ERROR] Attempting to remove an empty variable”
fi

#mv ${outDir}/*.e* ${logDir}
#mv ${outDir}/*.o* ${logDir}


find ./*.po* -type f -size 1M -delete
find ./*.pe* -type f -size 1M -delete


asd() {
cat <<"EOT"



                      _           _        _____                      _      _           _
    /\               | |         (_)      / ____|                    | |    | |         | |
   /  \   _ __   __ _| |_   _ ___ _ ___  | |     ___  _ __ ___  _ __ | | ___| |_ ___  __| |
  / /\ \ | '_ \ / _` | | | | / __| / __| | |    / _ \| '_ ` _ \| '_ \| |/ _ \ __/ _ \/ _` |
 / ____ \| | | | (_| | | |_| \__ \ \__ \ | |___| (_) | | | | | | |_) | |  __/ ||  __/ (_| |
/_/    \_\_| |_|\__,_|_|\__, |___/_|___/  \_____\___/|_| |_| |_| .__/|_|\___|\__\___|\__,_|
                         __/ |                                 | |
                        |___/                                  |_|
EOT
}

asd
echo "Analysis finished. Have a good day"

qstat -j $JOB_ID | grep usage

##########################################################################################################################
##########################################################################################################################
###Script done!
##########################################################################################################################
##########################################################################################################################


