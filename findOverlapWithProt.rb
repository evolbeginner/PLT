#! /bin/env ruby


require 'getoptlong'

require 'Hash'


####################################################################
linc_gff_indir = nil
gene_gff_indir = nil
attr = 'ID'
maf_indir = nil

linc_rela = Hash.new{|h1,k1|h1[k1]=Hash.new}
taxa_ordered = %w[Aarabicum Sirio tsa tpa Lalabamica CRU ALY]
linc_aligned_info = Hash.new{|h1,k1|h1[k1]=Hash.new{|h2,k2|h2[k2]=[]}}
lincs_overlapping_with_prot = Hash.new

coors = multi_D_Hash(3)
indices = Array.new

ref_taxon = nil


####################################################################
class Intersect_item
  attr_accessor :linc, :chr, :start, :stop, :intersect_start, :intersect_stop, :prot
end


####################################################################
def get_overlap_info(bedtools_intersect_output, attr, corename, linc_rela, linc_aligned_info, type)
  bedtools_intersect_output.each_line do |line|
    #scaffold_1                      5087943 5087987         +               ID=scaffold_1:5087943-5088216:At1NC017980       scaffold_1      phytozomev10    gene    5087772 5087987 .       -       0       ID=Carubv10012680m;pacid=20891911
    line.chomp!
    line_arr = line.split("\t")
    intersect_item_obj = Intersect_item.new

    intersect_start, intersect_stop = line_arr.values_at(3,4).map{|i|i.to_i}
    chr = line_arr[0]
    linc_attr_str = line_arr[8]
    linc_attr_str =~ /#{attr}=([^;]+)/
    linc = $1
    linc =~ /[^:]+$/
    linc_corename = $&
    start, stop = linc.split(':')[1].split('-').map{|i|i.to_i}

    intersect_item_obj = Intersect_item.new
    intersect_item_obj.linc = linc
    intersect_item_obj.chr = chr
    intersect_item_obj.start = start
    intersect_item_obj.stop = stop
    intersect_item_obj.intersect_start = intersect_start
    intersect_item_obj.intersect_stop = intersect_stop

    if type == 'P'
      prot_attr_str = line_arr[17]
      prot_attr_str =~ /#{attr}=([^;]+)/
      prot = $1
      intersect_item_obj.prot = prot
    end

    linc_rela[linc_corename][corename] = type
    linc_aligned_info[linc_corename][corename] << intersect_item_obj

  end
  return([linc_rela, linc_aligned_info])
end


####################################################################
opts = GetoptLong.new(
  ['--linc_gff_indir', '--linc_gff_dir', GetoptLong::REQUIRED_ARGUMENT],
  ['--gene_gff_indir', '--gene_gff_dir', GetoptLong::REQUIRED_ARGUMENT],
  ['--attr', GetoptLong::REQUIRED_ARGUMENT],
  ['--maf_dir', '--maf_indir', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '--linc_gff_indir', '--linc_gff_dir'
      linc_gff_indir = value
    when '--gene_gff_indir', '--gene_gff_dir'
      gene_gff_indir = value
    when '--attr'
      attr = value
    when '--maf_indir', '--maf_dir'
      maf_indir = value
  end
end


####################################################################
Dir.foreach(linc_gff_indir) do |b|
  next if b =~ /^\./
  next if b =~ /\.saf$/
  b =~ /(^[^.]+)/
  corename = $1
  linc_gff_file = File.join([linc_gff_indir, b])
  gene_gff_file = File.join([gene_gff_indir, corename+'.gene'+'.gff'])

  bedtools_intersect_output = `bedtools intersect -a #{linc_gff_file} -b #{gene_gff_file} -wb`
  linc_rela, linc_aligned_info = get_overlap_info(bedtools_intersect_output, attr, corename, linc_rela, linc_aligned_info, 'P')

  bedtools_intersect_output = `bedtools intersect -a #{linc_gff_file} -b #{gene_gff_file} -v`
  linc_rela, linc_aligned_info = get_overlap_info(bedtools_intersect_output, attr, corename, linc_rela, linc_aligned_info, 'L')
end


####################################################################
puts ['', taxa_ordered].flatten.join("\t")
linc_rela.each_pair do |linc, v|
  output_arr = Array.new
  output_arr << linc
  taxa_ordered.each do |taxon|
    if v.include?(taxon)
      output_arr << v[taxon]
      #p '!' if linc_aligned_info[linc][taxon][0].linc =~ /\w/
      #output_arr << linc_aligned_info[linc][taxon].map{|i|i.prot}.join("-")
      if linc_aligned_info[linc][taxon].map{|i|i.prot}.any?{|i|i=~/\w/}
        lincs_overlapping_with_prot[linc] = ''
      end
    else
      output_arr << 'N'
    end
  end
  #puts output_arr.join("\t")
end


####################################################################
lincs_overlapping_with_prot.each_pair do |linc, v|
  taxa = linc_rela[linc].keys
  maf_infile_basename = linc+'.maf'
  maf_infile = `find #{maf_indir} -name #{maf_infile_basename}`.chomp

  count = 0
  in_fh = File.open(maf_infile, 'r')
  sub_indices = Hash.new{|h,k|h[k]=[]}
  in_fh.each_line do |line|
    next if line =~ /^#/
    if count == 0
      sub_indices = Hash.new{|h,k|h[k]=[]}
    end

    if line =~ /^$/
      count = 0
      indices << sub_indices
    else
      count+=1
      if line =~ /^a score=/ #count==1
        ;
      else
        line_arr = line.split(/\s+/)
        chr_full_name, start, length, strand, chr_length, seq = line_arr.values_at(1,2,3,4,5,6)
        chr_full_name =~ /(^[^.]+)\.(.+)/
        taxon = $1
        chr = $2
        length = length.to_i
        chr_length = chr_length.to_i
        line.chomp!
        if strand == '+'
          start = start.to_i + 1
          stop = start + length - 1
        else
          stop = chr_length - start.to_i
          start = stop - length + 1
        end

        if count == 2
          ref_taxon = taxon
        end

        if count >= 2
          seq_arr = seq.split("")
          ref_seq_arr = seq_arr
          (start..stop).each_with_index do |coor, index|
            if seq_arr[index] =~ /[a-zA-Z]/
              coors[taxon][chr][coor] = coor
              sub_indices[taxon] << coor
            else
              if coors[taxon][chr].include?(coor-1)
                coors[taxon][chr][coor] = coors[taxon][chr][coor-1]
              else
                coors[taxon][chr][coor] = coor-1
              end
              if sub_indices[taxon].include?(index-1)
                sub_indices[taxon] << sub_indices[taxon][index-1]
              else
                sub_indices[taxon] << coor-1
              end
            end
          end
        end

      end
    end
  end

  linc_aligned_info[linc].each_pair do |taxon, list|
    p taxon
    list.each do |i|
      indices.each_with_index do |sub_indices, ind1|
        intersect_start_ind = sub_indices[taxon].find_index(i.intersect_start)
        intersect_stop_ind = sub_indices[taxon].find_index(i.intersect_stop)
        if not intersect_stop_ind.nil?
          ref_intersect_start_coor = sub_indices[ref_taxon][intersect_start_ind]
          ref_intersect_stop_coor = sub_indices[ref_taxon][intersect_stop_ind]
          puts [linc, i.intersect_start, i.intersect_stop].map{|i|i.to_s}.join("\t")
          puts [linc, intersect_start_ind, intersect_stop_ind].map{|i|i.to_s}.join("\t")
          puts [linc, ref_intersect_start_coor, ref_intersect_stop_coor].map{|i|i.to_s}.join("\t")
        end
      end
    end
  end
  indices = Array.new

end


