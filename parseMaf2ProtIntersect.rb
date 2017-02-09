#! /bin/env ruby


require 'getoptlong'

require 'Dir'


#######################################################
maf_parse = File.expand_path("~/software/genome/phast-1.3/bin/maf_parse")
tmp_dir = File.expand_path("~/tmp")


maf_file = nil
gff_file = nil
outdir = nil
is_force = false
is_tolerate = false

aligned_regions = Hash.new{|h,k|h[k]=[]}


#######################################################
class Maf_item
  attr_accessor :chr, :start, :stop, :strand, :lincRNA
end


#######################################################
def find_aligned_regions(maf_parse, tmp_dir, gff_file, maf_file, lincRNA_maf_outdir, aligned_regions)
  tmp_gff_file = File.join([tmp_dir, 'maf_intersect.gff'])

  in_fh = File.open(gff_file, 'r')
  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    attr_str = line_arr[-1]
    attr_str =~ /ID=(.+)/
    lincRNA = $1
    out_fh = File.open(tmp_gff_file, 'w')
    out_fh.puts line
    out_fh.close

    lincRNA_maf_outfile = File.join([lincRNA_maf_outdir, lincRNA+'.maf'])
    if File.exist?(lincRNA_maf_outfile)
      ;
    else
      `#{maf_parse} -g #{tmp_gff_file} #{maf_file} > #{lincRNA_maf_outfile}`
    end

    count = 0
    File.open(lincRNA_maf_outfile, 'r').each_line do |line|
      line.chomp!
      next if line =~ /^#/
      if line =~ /^$/
        count = 0
      else
        count+=1
        if line =~ /^a score=/ #count==1
          ;
        else
          #s ATH.1                      133 67 + 30427671 TCA--TTGTGTAT-ATAATGATAATTTTATCGTTTTTA------TGTAATTGCTTATTGTTGTGTGTAGATTTTT----T
          #s ALY.scaffold_8        13643265 74 + 22951293 TAGGGTTGTGTCTTCTAGTGACATTCTTAACGTATTTA------TGCGACCTCCTACTTCAGTATGCAGATCTATATATT
          line_arr = line.split(/\s+/)
          chr_full_name, start, length, strand, chr_length, seq = line_arr.values_at(1,2,3,4,5,6)
          chr_full_name =~ /(^[^.]+)\.(.+)/
          taxon = $1
          chr = $2
          length = length.to_i
          chr_length = chr_length.to_i
          if strand == '+'
            start = start.to_i + 1
            stop = start + length - 1
          else
            stop = chr_length - start.to_i
            start = stop - length + 1
          end
          #puts [chr, '.', 'intersect', start, stop, '.', strand, '.', "ID=#{lincRNA}"].map{|i|i.to_s}.join("\t")
          maf_item_obj = Maf_item.new
          maf_item_obj.chr = chr
          maf_item_obj.start = start
          maf_item_obj.stop = stop
          maf_item_obj.strand = strand
          maf_item_obj.lincRNA = lincRNA
          aligned_regions[taxon] << maf_item_obj
        end
      end
    end
  end
  in_fh.close
  return(aligned_regions)
end


def get_new_maf_item_objs(maf_item_objs)
  maf_item_objs.each do |maf_item_obj1|
    maf_item_objs.each do |maf_item_obj2|
      if maf_item_obj2.start == maf_item_obj1.stop + 1
        maf_item_obj1.stop = maf_item_obj2.stop
        maf_item_objs.delete(maf_item_obj2)
      end
    end
  end
  return(maf_item_objs)
end


#######################################################
opts = GetoptLong.new(
  ['--maf', GetoptLong::REQUIRED_ARGUMENT],
  ['--gff', GetoptLong::REQUIRED_ARGUMENT],
  ['--outdir', GetoptLong::REQUIRED_ARGUMENT],
  ['--force', GetoptLong::NO_ARGUMENT],
  ['--tolerate', GetoptLong::NO_ARGUMENT],
  ['--tmp_dir', GetoptLong::REQUIRED_ARGUMENT],
  ['--gene_gff_dir', '--gene_gff_indir', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--maf'
      maf_file = value
    when '--gff'
      gff_file = value
    when '--outdir'
      outdir = value
    when '--force'
      is_force = true
    when '--tolerate'
      is_tolerate = true
    when '--tmp_dir'
      tmp_dir = value
  end
end


mkdir_with_force(outdir, is_force, is_tolerate)

blast_outdir = File.join([outdir, 'blast'])
mkdir_with_force(blast_outdir, is_force, is_tolerate)

lincRNA_maf_outdir = File.join([outdir,'lincRNA_maf'])

mkdir_with_force(lincRNA_maf_outdir, is_force, is_tolerate)


#######################################################
aligned_regions = find_aligned_regions(maf_parse, tmp_dir, gff_file, maf_file, lincRNA_maf_outdir, aligned_regions)

aligned_regions.each_pair do |taxon, maf_item_objs|
  outfile = File.join([blast_outdir, taxon+'.blast8'])
  out_fh = File.open(outfile, 'w')
  maf_item_objs = get_new_maf_item_objs(maf_item_objs)
  maf_item_objs.each do |maf_item_obj|
    lincRNA = maf_item_obj.lincRNA
    chr = maf_item_obj.chr
    start = maf_item_obj.start
    stop = maf_item_obj.stop
    length = maf_item_obj.stop - maf_item_obj.start + 1
    out_fh.puts [lincRNA, chr, 100, length, 0, 0, 0, 0, start, stop, 1e-100, 1000].map{|i|i.to_s}.join("\t")
  end
  out_fh.close
end


