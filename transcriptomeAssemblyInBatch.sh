#! /bin/bash


###############################################################################
stringtie=/mnt/bay3/sswang/software/NGS/stringtie-1.3.2b.Linux_x86_64/stringtie
cufflinks=cufflinks

is_stringtie=false
is_cufflinks=false
cpu=1
is_force="false"


###############################################################################
while [ $# -gt 0 ]; do
	case $1 in
		-i)
			infile=$2
			shift
			;;
		--stringtie|--Stringtie)
			is_stringtie=true
			;;
		--cufflinks|--Cufflinks)
			is_cufflinks=true	
			;;
		--gtf)
			gtf_file=$2
			shift
			;;
		--cpu)
			cpu=$2
			shift
			;;
		--outdir)
			outdir=$2
			shift
			;;
		--force)
			is_force="true"
			;;
	esac
	shift
done


if [ $is_stringtie == false -a $is_cufflinks == false ]; then
	echo "stringtie or cufflinks has to be specified!. Exiting ......"
	exit 1
fi


if [ -z $infile ]; then
	echo "infile has not been given! Exiting ......"
	exit 1
elif [ -z $gtf_file ]; then
	echo "gtf_file has not been given! Exiting ......"
	exit 1
elif [ -z $outdir ]; then
	echo "outdir has not been given! Exiting ......"
	exit
fi


if [ $is_force == "true" ]; then
	[ -d $outdir ] && rm -rf $outdir
fi


###############################################################################
function getCorenameOfBamFile(){
	local bam_file
	bam_file=$1
	b=`basename $1`
	a=`grep -oP '^(SRR|DRR|ERR)\d+' <<< $b`
	if [ ! -z $a ]; then
		corename=$a
	else
		corename=`grep -oP '^[^.]+' <<< $b`
	fi
	echo $corename
}


while read -r line; do
	tissue=`awk '{print $1}' <<< $line`
	list_str=`awk '{print $2}' <<< $line`
	echo $tissue
	for i in `awk -F , '{for(i=1;i<=NF;i++){print $i}}' <<< $list_str`; do
		for bam_file in $i/*bam; do
			echo $bam_file
			corename=$(getCorenameOfBamFile $bam_file)
			[ ! -d $outdir/$tissue/$corename ] && mkdir -p $outdir/$tissue/$corename
			if [ $is_stringtie == true ]; then
				$stringtie -G $gtf_file -p $cpu -c 1 -o $outdir/$tissue/$corename/transcripts.gtf $bam_file
			elif [ $is_cufflinks == true ]; then
				$cufflinks -q -p $cpu -o $outdir/$tissue/$corename/ $bam_file
			fi
		done
	done
done < $infile


