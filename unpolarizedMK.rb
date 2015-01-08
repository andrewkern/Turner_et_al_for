#!/usr/bin/env ruby
#
# unpolarizedMK.rb
#
# takes a pair of directories of multiple alignments as it's argument: one for outgroup seqs and one for ingroup seqs

load "#{File.join(__dir__, 'SequenceMatrix.rb')}"

if ARGV.size < 2
        print "Performs an unpolarized MK test for each coding sequence for which both ingroup and outgroup alignments are provided.\n\nusage: ./unpolarizedMK.rb outGroupDir inGroupDir\n\noutGroupDir: a directory of cross-species alignments, one file for each transcript. The first sequence in each alignment is the reference sequence for the ingroup (which is ignored). A sequence from a related species must follow.\ninGroupDir: a directory of within-species alignments, one file for each transcript. The file names within both directories must be of the format <transcriptId>.fa\n"
	exit
end


genes = Dir[ARGV[0]+"/*"]
print "locus\tmel aaFix\tmel aaPoly\tmel silFix\tmel silPoly\n"
genes.each{ | outGroupPath |
        geneName = outGroupPath.split("/")[-1]
        inGroupPath = ARGV[1]+"/"+geneName
        if File.exists?(inGroupPath)
		outGroupSeqs = SequenceMatrix.new.initializeFromFasta(outGroupPath)
                #first one is mel, so we only want the second two (sim and yak)
		sim = SequenceMatrix.new.initializeWith([outGroupSeqs.nameVector[1]],[outGroupSeqs.matrix[1]])
		sim.asCodingSequence(nil,nil)
                #now to get all of the mel seqs out of the ingroup file
                inGroupSeqs = SequenceMatrix.new.initializeFromFasta(inGroupPath)
                melNamesArray = Array.new
                inGroupSeqs.nameVector.each { | seqName |
			melNamesArray << seqName
		}
		melSeqsArray = Array.new
		inGroupSeqs.matrix.each { | seq |
			#melSeqsArray << seq
			melSeqsArray << seq
		}
		mel = SequenceMatrix.new.initializeWith(melNamesArray,melSeqsArray)
		mel.asCodingSequence(nil,nil)
		align = Alignment.new(mel,sim,nil)
		#print align.sequenceMatrix1.matrix[0],"\n"
		array = align.mkTestGreedy
		print geneName.chomp.split(".fa")[0],"\t",array[0],"\t",array[1],"\t",array[2],"\t",array[3],"\n"
		#print array[1][0],"\t",array[1][1],"\t",array[1][2],"\t",array[1][3],"\n"
	end
	}
