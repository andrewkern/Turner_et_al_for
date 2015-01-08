# GFF class 
# 
# 10/11/04 Andrew Kern
#
#
# meant to be compatible with GFF 3
#
#

load "#{File.join(__dir__, 'SequenceMatrix.rb')}"

class Array
	include Enumerable
end

class GFF
	attr_reader :records, :sequenceDict, :sequenceNames
	attr_writer	:records, :sequenceDict, :sequenceNames
	
	def initialize(aFilename)
		self.records = Array.new
		aFile = File.new(aFilename)
		aFile.each_line{ | line |
			if ! line.strip.empty?
				if  !(line[0,1] == "#")
					self.records.push(Record.new(line))
				end
			end
			}
		aFile.close
		
	end
	

	def initializeFromString(aString)
		aString.each_line{ | line |
			if ! line.strip.empty?
				if  !(line[0,1] == "#")
					self.records.push(Record.new(line))
				end
			end
			}
	end
		
		
	def inheritSequence(aString)
		self.sequence = aString.upcase
	end
	
	def inheritSequenceFromFasta(aFilename)			#note that this passes the sequence/names as hash
		seqMat = SequenceMatrix.new.initializeFromFasta(aFilename)
		self.sequenceDict = Hash.new
		self.sequenceNames = seqMat.nameVector
		self.sequenceNames.each_index{ | i | 
			self.sequenceDict[self.sequenceNames[i].split.first] = seqMat.matrix[i]
			}
	end
		
	def getRecordSequence(aRecord)
		coords = Array.new
		coords << aRecord.start.to_i 
		coords << aRecord.end.to_i
		coords.sort!
		size = self.sequenceDict[aRecord.seqid].size
		max = coords.last
		if size < max
			max = size
		end
		l = max - coords.first + 1
		seq = self.sequenceDict[aRecord.seqid][coords.first - 1, l]		
			
		if aRecord.strand == "+"
			return seq
		else
			return seq.reverse.tr('ATGC','TACG')
		end
	end
	
	def getSequenceBetweenRecords(record1, record2)
		coords = Array.new
		if record1.strand == "+"
			coords << record2.start.to_i - 1
			coords << record1.end.to_i + 1
		else
			coords << record2.end.to_i + 1
			coords << record1.start.to_i - 1
		end
		coords.sort!
		l = coords.last - coords.first + 1
		seq = self.sequenceDict[record1.seqid][coords.first - 1, l]
		if record1.strand == "+"
			return seq
		else
			return seq.reverse.tr('ATGC','TACG')
		end
	end
	
	def returnRecordType(aString)
		return records.find_all{ | each | each.type == aString }
	end
	
	def returnRecordSeqID(aString)
		return records.find_all{ | each | each.seqid == aString }
	end
	
	def returnTypeSet
		set = Array.new
		self.records.each{ | each |
			set << each.type
			}
		return set.uniq
	end
	
	def returnSeqIDSet
		set = Array.new
		self.records.each{ | each |
			set << each.seqid
			}
		return set.uniq
	end
	
	def returnSourceSet
		set = Array.new
		self.records.each{ | each |
			set << each.source
			}
		return set.uniq
	end
	
	def outputRecordTypeFasta(aString)
		records = self.returnRecordType(aString)
		records.sort!{ | x, y |	x.start.to_i <=> y.start.to_i }
		records.each{ | each |
			print ">"+each.attributes["ID"]+"\n"
			print self.getRecordSequence(each),"\n"
			}
	end
	def returnIntergenicRegions		#returns a hash based on each seqid
		genes = self.returnRecordType("gene") | self.returnRecordType("transposable_element")
		genes.sort!{ | x, y | x.start.to_i <=> y.start.to_i}
		geneHash = Hash.new
		resultHash = Hash.new
		self.returnSeqIDSet.each{ | x | geneHash[x] = Array.new
										resultHash[x] = Array.new}
		genes.each{ | each | geneHash[each.seqid] << each }
		geneHash.each_pair{ | key, value |
			start = 0
			value.each{ | each |
				resultHash[key] << [start, each.start.to_i - 1]
				start = each.end.to_i + 1
				}
			resultHash[key] << [start, self.sequenceDict[key].size - 1]
		}
		return resultHash
	end
	
	def outputIntergenicFasta
		
		interH = self.returnIntergenicRegions
		interH.each_pair{ | key, value |
				i = 0
				value.each{ | x |
					x.sort!
					print ">"+key+".intergenic."+i.to_s+"\n"
					l = x.last - x.first + 1
					seq = self.sequenceDict[key][x.first - 1, l] 
					print seq,"\n"
					i += 1
					}
				}
	end
		
	def returnTranscriptSet      # note this is probably not general. designed specifically for flybase
		set = Array.new
		self.returnRecordType("CDS").each{ | x | x.attributes["Parent"].each{|x| set << x }}
		return set.uniq
	end
	
	#note there have been substantial changes to the GFF format necessitating the bizarre changes below
	# I wish flybase would stay consistent....
	def arrangeTranscriptHash    # returns hash of transcripts in record form
            hash = Hash.new
            allTrans = self.returnTranscriptSet
            cdsArray = self.returnRecordType("CDS")
            cdsHash = Hash.new
            cdsArray.each{ | eachCDS |
                eachCDS.attributes["Parent"].each{ |aPar | 
                  if cdsHash[aPar].nil?
                    cdsHash[aPar] = Array.new
                  end
                  cdsHash[aPar] << eachCDS
                }
            }
            cdsKeys = cdsHash.keys
            #allTrans.each{ | transName | hash[transName] = Array.new }
            #keys = hash.keys
            #exons = self.returnRecordType("exon")
            #exons.each{ | eachExon | 
              #  eachExon.attributes["Parent"].each{ | eachMRNA |
               #     if keys.include?(eachMRNA.strip)
                #        hash[eachMRNA.strip] << eachExon.clone 
                 #   end
                #}
            #}
            #trying this
            hash = cdsHash.clone
            allTrans.each{ | transName |
                if cdsKeys.include?(transName)
                 # print transName,"\n"
                  if hash[transName].first.strand == "+"
                          hash[transName].sort!{ | x, y | (x.start.to_i) <=> (y.start.to_i) }
                          #hash[transName].reject!{ | x | x.end.to_i < cdsHash[transName].start.to_i or x.start.to_i > cdsHash[transName].end.to_i}
                          #hash[transName].first.start = cdsHash[transName].start
                          #hash[transName].last.end = cdsHash[transName].end.to_i + 3
                  else
                          hash[transName].sort!{ | x, y | (y.start.to_i) <=> (x.start.to_i) }
                          #hash[transName].reject!{ | x | x.end.to_i < cdsHash[transName].start.to_i or x.start.to_i > cdsHash[transName].end.to_i}
                          #hash[transName].first.end = cdsHash[transName].end
                          #hash[transName].last.start = cdsHash[transName].start.to_i - 3
                  end
                end
            }
            return hash
	end
	
	def transcriptGenomicLocations			#returns a Hash of arrays indicating the coding sites for each transcript
		locHash = Hash.new
		transHash = arrangeTranscriptHash
		transHash.keys.each{ | key |
			tempArray = Array.new
			transHash[key].each{ | record |
				if record.strand == "+"
					(record.start.to_i).upto(record.end.to_i){ | i | tempArray << i }
				else
					(record.end.to_i).downto(record.start.to_i){ | i | tempArray << i }
				end
			}
		locHash[key] = tempArray
		}
		return locHash
	end
	
	def transcriptStartStop					#returns a hash of chromosome oriented start and stops
		transHash = arrangeTranscriptHash
		transHash.keys.each{ | key |
			tempArray = Array.new
			transHash[key].each{ | record |
					print record.seqid,"\t",record.start,"\t",record.end,"\n"
			}
		}
	end
	
	
	def arrangeTranscriptHashStrings  #returns hash of transcript strings
		hash = self.arrangeTranscriptHash
		stringHash = Hash.new
		hash.each_pair{ | key, value | 
                #print value.first.seqid,"\n"
                if self.sequenceDict[value.first.seqid] != nil
                  l = self.sequenceDict[value.first.seqid].size
                  stringHash[key] = String.new("")
                  value.each{ | eachRecord |
                      if eachRecord.start.to_i < l and eachRecord.end.to_i < l
                          stringHash[key] << self.getRecordSequence(eachRecord) 
                      end
                  }
                end
            }
            stringHash.delete_if{ |key,value| value.empty? }
            return stringHash
	end
	
  def arrangeTranscriptCodons  #returns hash of transcript strings
    hash = self.arrangeTranscriptHash
    codonHash = Hash.new
    posHash = Hash.new
    hash.each_pair{ | key, value | 
      #print value.first.seqid,"\n"
      if self.sequenceDict[value.first.seqid] != nil
        l = self.sequenceDict[value.first.seqid].size
        posHash[key] = Array.new
          value.each{ | eachRecord |
            if eachRecord.start.to_i < l and eachRecord.end.to_i < l
              eachRecord.start.upto(eachRecord.end){| i | posHash[key] << i}
            end
          }
      end
    }
    posHash.delete_if{ |key,value| value.empty? }
    posHash.each_pair{ |key,value|
      value.sort!
      codonHash[key] = Array.new
      tempStrand = hash[key].first.strand
      while value.length > 0
        pos = [value.shift,value.shift,value.shift]
        seq = self.sequenceDict[hash[key].first.seqid][pos[0].to_i - 1, 1] + self.sequenceDict[hash[key].first.seqid][pos[1].to_i - 1, 1] + self.sequenceDict[hash[key].first.seqid][pos[2].to_i - 1, 1]
        if tempStrand == "-"
          seq = seq.reverse.tr('ATGC','TACG')
          pos = pos.reverse  #maybe?
        end
        codonHash[key] << Codon.new(seq,pos,tempStrand).translate
      end
      }
    return codonHash
  end
	
	def spliceSiteHash			#returns hash of arrays of splice site strings
		hash = self.arrangeTranscriptHash
		result = Hash.new
		hash.each_pair{ | key, value |
			if value.size > 1
				l = self.sequenceDict[value.first.seqid].size
				spliceStringArray = Array.new
				0.upto(value.size - 2){ | i |
					if (value[i].start.to_i < l and value[i].end.to_i < l) and (value[i + 1].start.to_i < l and value[i + 1].end.to_i < l)
						intronString = self.getSequenceBetweenRecords(value[i], value[i + 1])
						string = String.new("")
						string << intronString[0,2]
						string  << intronString[intronString.size - 2,2]
						spliceStringArray << [i,string]
					end
					}
				if spliceStringArray.size > 0	
					result[key] = spliceStringArray
				end
			end
			}
		return result
	end
	
	def relativePosition(sequenceName, anIndex)			#returns the position of anIndex relative to the gaps inserted during alignment in sequence, useful for annotations but note it is zero based
		tempIndex = relIndex = 0
		string = self.sequenceDict[sequenceName]
		while tempIndex <= anIndex.to_i
			if string[relIndex,1] == "-"
				relIndex += 1
			else
				relIndex += 1
				tempIndex += 1
			end
			
		end
		return relIndex - 1
	end
	
	def hashPositions			#goes through sequence makes hash(orig) -> aligned
		posHash = Hash.new
		self.sequenceDict.keys.each{ | aKey |
			posHash[aKey] = Hash.new
			}
		posHash.keys.each{ | aKey |
		
			tempIndex = relIndex = 1
			last = self.sequenceDict[aKey].size
			while tempIndex < last
				if self.sequenceDict[aKey][relIndex - 1,1] == "-"
					relIndex += 1
				else
					posHash[aKey][tempIndex] = relIndex
					relIndex += 1
					tempIndex += 1
				end
				
			end
		}
	
		return posHash
	end
	
	def adjustPositions			#adjusts the coordinate space given gaps from alignment
		posHash = self.hashPositions
		self.records.each{ | eachRecord |
			eachRecord.start = posHash[eachRecord.seqid][eachRecord.start.to_i]
			eachRecord.end = posHash[eachRecord.seqid][eachRecord.end.to_i]
			}
	end	
	
end

$thegeneticcode = { "GCT" => "A",
	 "GCC" => "A",
	 "GCA" => "A",
	 "GCG" => "A",
	 "TGT" => "C",
	 "TGC" => "C",
	 "GAT" => "D",
	 "GAC" => "D",
	 "GAA" => "E",
	 "GAG" => "E",
	 "TTT" => "F",
	 "TTC" => "F",
	 "GGT" => "G",
	 "GGC" => "G",
	 "GGA" => "G",
	 "GGG" => "G",
	 "CAT" => "H",
	 "CAC" => "H",
	 "ATT" => "I",
	 "ATC" => "I",
	 "ATA" => "I",
	 "AAA" => "K",
	 "AAG" => "K",
	 "TTG" => "L",
	 "TTA" => "L",
	 "CTT" => "L",
	 "CTC" => "L",
	 "CTA" => "L",
	 "CTG" => "L",
	 "ATG" => "M",
	 "AAT" => "N",
	 "AAC" => "N",
	 "CCT" => "P",
	 "CCC" => "P",
	 "CCA" => "P",
	 "CCG" => "P",
	 "CAA" => "Q",
	 "CAG" => "Q",
	 "CGT" => "R",
	 "CGC" => "R",
	 "CGA" => "R",
	 "CGG" => "R",
	 "AGA" => "R",
	 "AGG" => "R",
	 "TCT" => "S",
	 "TCC" => "S",
	 "TCA" => "S",
	 "TCG" => "S",
	 "AGT" => "S",
	 "AGC" => "S",
	 "ACT" => "T",
	 "ACC" => "T",
	 "ACA" => "T",
	 "ACG" => "T",
	 "GTT" => "V",
	 "GTC" => "V",
	 "GTA" => "V",
	 "GTG" => "V",
	 "TGG" => "W",
	 "TAT" => "Y",
	 "TAC" => "Y",
	 "TAA" => "*",
	 "TAG" => "*",
	 "TGA" => "*"}			


			
class Record
	attr_reader :seqid, :source, :type, :start, :end, :score, :strand, :phase, :attributes 
	attr_writer  :seqid, :source, :type, :start, :end, :score, :strand, :phase, :attributes
	
	def initialize(aString)
		tokens = aString.chomp.split
		self.seqid = tokens[0]
		self.source = tokens[1]
		self.type = tokens[2]
		self.start = tokens[3]
		self.end = tokens[4]
		self.score = tokens[5]
		self.strand = tokens[6]
		self.phase = tokens[7]
		atr = tokens[8]
		if atr
			self.attributes = parseAttributes(atr)
		end
		return self
	end
	
	def printRecord
		print self.seqid,"\t",self.source,"\t",self.type,"\t",self.start,"\t",self.end,"\t",self.score,"\t",self.strand,"\t",self.phase,"\t"
		if self.attributes
			self.attributes.each_pair{ | key, value |
					print key,"=",value,"; "
					}
		end
		print "\n"
	end
	
	private
		
	def parseAttributes(aString)
		aHash = Hash.new
		aString.split(";").each{ | atr |
				key, value = atr.split("=", 2)
        if key == "Parent"
				  aHash[key] = value.split(",")
        else
				  aHash[key] = value
        end
		}
		return aHash
	end
	

	
end
		
class Codon
	attr_reader :trans, :seq, :strand, :positions		
	attr_writer	:trans, :seq, :strand, :positions	
	
	def initialize(aSeq, aPos, aStrand)
	  self.seq = aSeq
	  self.positions = aPos
	  self.strand = aStrand
	end
	
	def translate
			self.trans = $thegeneticcode[self.seq]
			return self
	end
	
end
