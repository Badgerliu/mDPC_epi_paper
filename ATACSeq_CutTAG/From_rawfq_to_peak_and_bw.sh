#!/usr/bin/env bash

mkdir ./trimmed

for i in $(ls *.fq.gz | rev |cut -c 10- |rev |uniq)
	do
		java -jar trimmomatic-0.39.jar PE -phred33 ${i}_R1.fq.gz ${i}_R2.fq.gz ${i}_R1_paired.fq.gz ${i}_R1_unpaired.fq.gz ${i}_R2_paired.fq.gz ${i}_R2_unpaired.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:7
		mv ${i}_R* ./trimmed
		echo done with ${i}
	done


cd ./trimmed

mkdir ./mapped
for i in $(ls *_paired.fq.gz|rev | cut -c 17- |rev| uniq)
	do
		bowtie2 -x mm10 -1 ${i}_R1_paired.fq.gz -2 ${i}_R2_paired.fq.gz -X 1500 -p 16 -S ${i}.sam

		sed '/chrM/d;/random/d;/chrUn/d' ${i}.sam > ${i}_m.sam

		samtools view -bT mm10.fa ${i}_m.sam > ${i}_m.bam

		rm *.sam

		samtools sort -o ${i}.bam ${i}_m.bam

		java -jar /home/liuhuan/project/software/picardtools/picard.jar MarkDuplicates I=${i}.bam O=${i}_m_d.bam M=dups.txt REMOVE_DUPLICATES=true

		samtools index ${i}_m_d.bam && rm ${i}.bam && rm ${i}_m.bam
		mv ${i}_m_d.bam* ./mapped
		echo OK done with ${i}
	done



mkdir bw

echo "Make sure you start virtualEnv for deeptools and go to the right folder before starting this shell, or you will have trouble!"

counter=0

for sample in *.bam
	do
		echo $sample
		describer=$(echo ${sample} | sed 's/_m_d.bam//') # define the describer with sed function
		echo $describer  
		#Normalization using deeptools
		bamCoverage --bam ${describer}_m_d.bam -o ${describer}.bw --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --ignoreForNormalization chrX chrM --extendReads
		echo Normalization with ${describer} is done
		let counter=counter+1
		echo ${counter} bamfiles have been finished
	done

mv *.bw ./bw
	
echo "Done with normalization! You will find the results in folder 'bw', and you can visualize the data in IGV now! --Huan"

cd ./mapped

echo "Now convert bam to bed"


for zyy in *_m_d.bam
	do
		echo $zyy
		describer=$(echo ${zyy} | sed 's/.bam//')
		echo $describer conversion begins
		bedtools bamtobed -i $describer.bam >$describer.bed
		echo ${describer} converstion finished
	done
	
echo "bam to bed conversion! --Huan"


mkdir shift
for sample in *_m_d.bed
	do
		echo $sample
		describer=$(echo ${sample} | sed 's/.bed//')
		echo $describer  
		#shift the read position
		awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $3 + 4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}' ${describer}.bed >${describer}_shifted.bed
		echo ${describer} shift finished
	done
	
echo "Done with position shift, and now moving shifted files to shift folder! --Huan"

mv *_shifted.bed ./shift

cd ./shift

mkdir peak
echo "Do peak calling"


for foo in *_m_d_shifted.bed
	do
		echo $foo
		bar=$(echo ${foo} | sed 's/_m_d_shifted.bed//') # define the describer with sed function
		echo $bar  
		#peak calling
		macs2 callpeak --nomodel -t ${bar}_m_d_shifted.bed -n ${bar} --nolambda --gsize 2.7e9 --keep-dup all --slocal 10000
		echo peak calling for ${bar} finished
	done

mv *_summits.bed ./peak
mv *.narrowPeak ./peak
mv *.xls ./peak


rm *_m.bed  # (optional) remove the chrM deleted bedfiles just to save disk space.
	
echo "Done with peak calling. You will find results in 'peak' folder.  --Huan"


