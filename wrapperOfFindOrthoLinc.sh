#! /bin/bash


##############################################################################
home_dir=~/project/alex/SNP/plants/
scripts_dir=$home_dir/scripts
ortho_dir=$home_dir/ortho
linc_dupli_noBlock=$scripts_dir/linc_dupli_noBlock.rb
generate_pseudo_blastoutput_gff=$scripts_dir/generate_pseudo_blastoutput_gff.sh
calculate_overlaps=$scripts_dir/calculate_overlaps.rb
gff0_2_3=$scripts_dir/gff0_2_3.rb
outdir=.
is_force=false
is_best=true

blast_in=blast
syntenic_linc_in=syntenic_linc
home_in=.


##############################################################################
while [ $# -gt 0 ]; do
	case $1 in
		--taxa)
			taxa=$2
			shift
			;;
		--ref_taxa|--ref)
			ref_taxa=$2
			shift
			;;
		--outdir)
			outdir=$2
			shift
			;;
		--force)
			is_force=true
			;;
		--no_best)
			is_best=false
			;;
		--blast_in)
			blast_in=$2
			shift
			;;
		--syntenic_linc_in)
			syntenic_linc_in=$2
			shift
			;;
		--home_in)
			home_in=$2
			shift
			;;
	esac
	shift
done


if [ -z $taxa -o -z $ref_taxa ]; then
	echo 'taxa and ref_taxa have to be specified!'
	echo 'Exiting ......'
	exit 1
fi


if [ ! -d $outdir ]; then
	mkdir -p $outdir
else
	if [ $is_force == true ]; then
		rm -rf $outdir
		mkdir -p $outdir
	fi
fi


linc_out_saf=$outdir/$taxa.linc.ortho.saf
linc_out_gff=$outdir/$taxa.linc.ortho.gff


if [ $is_best == true ]; then
	ref_arg='--ref query'
else
	ref_arg=''
fi


##############################################################################
bash $generate_pseudo_blastoutput_gff -i $blast_in/$taxa/blast_result --bias suffix --outdir $syntenic_linc_in/$taxa/MM

ruby2.1 $linc_dupli_noBlock --blast $syntenic_linc_in/$taxa/MM/*blast --saf $home_in/linc_data/${ref_taxa}_linc.saf --saf $syntenic_linc_in/$taxa/MM/*gff --gff $home_in/data/$ref_taxa.gene.gff --gff $home_in/gff/$taxa.gene.gff --gff_suffix suffix --pair $ortho_dir/${ref_taxa}-$taxa.prot.ortho --shared_gene 0.5 $ref_arg --evalue 0.001 | \
	ruby -nae 'arr = $_.split("\t"); arr[0]=~/[-]suffix/ ? (puts arr[0]) : (puts arr[1])' | \
	ruby $calculate_overlaps --i1 - --i2 $syntenic_linc_in/$taxa/MM/*gff --f2 2 --show 12 --content 2 --no_report \
	> $linc_out_saf

sed -i 's/[-]suffix//g' $linc_out_saf
ruby $gff0_2_3 --rev_coor -i $linc_out_saf > $linc_out_gff


