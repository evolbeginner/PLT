#! /bin/bash


set_EV=~/tools/self_bao_cun/tba_wrapper/scripts/set_EV.sh


###########################################################################
source $set_EV


while [ $# -gt 0 ]; do
	case $1 in
		--indir)
			indir=$2
			shift
			;;
		--tree)
			tree=$2
			shift
			;;
	esac
	shift
done


###########################################################################
for maf in $indir/*maf; do
	b=`basename $maf`
	corename=`echo $b | cut -f 1 -d .`
	echo $corename
	mkdir $corename
	old_PWD=$PWD
	cd $corename >/dev/null
	time phyloFit --msa-format MAF --tree $tree $old_PWD/$b 
	cd - >/dev/null
done


