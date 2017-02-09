#! /bin/bash


##############################################################
while [ $# -gt 0 ]; do
	case $1 in
		--indir)
			indir=$2
			shift
			;;
	esac
	shift
done


##############################################################
for i in $indir/*; do
	[ ! -d $i ] && continue
	b=`basename $i`
	corename=`echo $b|cut -d . -f 1`
	cd $i >/dev/null
	phastCons --rho 0.4 --target-coverage 0.25 --expected-length 12 --msa-format MAF --seqname $corename --most-conserved most-cons.bed *.maf phyloFit.mod > $corename.wig
	cd - >/dev/null
done


