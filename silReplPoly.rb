#!/usr/bin/env ruby
#
# this outputs polymorphism stats

load "#{File.join(__dir__, 'SequenceMatrix.rb')}"

if ARGV.empty?
    print "Outputs estimates of pi and the number of segregating sites at silent sites, replacement sites, and the entire coding region.\n\nusage: ./silReplPoly.rb fastaFile\n\nfastaFile: a path to a fasta-formatted file containing an alignment of coding sequences from a population sample\n"
    exit
end
print "segSites\tsil_s\trepl_s\tpi\tsilPi\trPi\n"

seq = SequenceMatrix.new.initializeFromFasta(ARGV[0])
seq.asCodingSequence(nil,nil)
muts = seq.silReplMutations
new = seq.silReplPi2
# print "sil_s\trepl_s\tsegSites\tsilentPi_old\treplacementPi_old\tpi\tsilPi_new\trPi_new\n"

print seq.segSites,"\t",muts["silents"].size,"\t",muts["replacements"].size,"\t",seq.piClean.round(4),"\t",new[1].round(4),"\t",new[0].round(4),"\n"

