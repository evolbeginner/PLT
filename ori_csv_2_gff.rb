#! /bin/env ruby


require 'getoptlong'


############################################################
infile = nil


############################################################
opts = GetoptLong.new(
  ['-i', GetoptLong::REQUIRED_ARGUMENT],
)


opts.each do |opt, value|
  case opt
    when '-i'
      infile = value
  end
end


############################################################
in_fh = File.open(infile, 'r')
in_fh.each_line do |line|
  next if $. == 1
  line.chomp!
  line_arr = line.split('",')
  line_arr.each do |i|
    i.gsub!('"', '')
  end
  region1 = line_arr[1]
  region1 =~ /([^,]+),(.+)/
  chr = $1
  start = $2
  start.gsub!(',', '')
  start = start.to_i
  region2 = line_arr[2]
  stop = region2
  stop.gsub!(',', '')
  stop = stop.to_i

  lincRNA = 'lincRNA' + ($.-1).to_s
  puts [chr, '.', 'gene', start, stop, '.', '+', '.', "ID=#{lincRNA}"].join("\t")
end
in_fh.close



