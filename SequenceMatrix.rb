# SequenceMatrix class - This is the begining of a Ruby polymorphism data implementation
# The Goal here is to make a lightweight, scriptable, population genetic analysis environment
#
#
#

include Math
load "#{File.join(__dir__, 'arrayADK.rb')}"
load "#{File.join(__dir__, 'mathADK.rb')}"
load "#{File.join(__dir__, 'CodingSequence.rb')}"
load "#{File.join(__dir__, 'Alignment.rb')}"

class Array
  include ArrayADK
  include Enumerable
end
class SequenceMatrix
  include CodingSequence
end

class MultipleMatrix
  attr_reader :samples
  attr_writer :samples

  def initialize(aFile)
    fileLines = File.readlines(aFile)
    self.samples = Array.new
    while ! fileLines.empty?
      tempArray = fileLines.slice!(0,fileLines.index("//\n") + 1)
      tempArray.pop
      self.samples << SequenceMatrix.new.initializeFromString(tempArray)
    end
  end
end

class SequenceMatrix
  attr_reader :matrix, :nameVector, :codingRegions, :readingFrame, :geneticCode, :codons, :features,:filenameID, :codonDists, :puDict, :silentSiteDict
  attr_writer :matrix, :nameVector, :codingRegions, :readingFrame, :geneticCode, :codons, :features, :filenameID, :codonDists, :puDict, :silentSiteDict

  def initializeFromFasta(aFilename)
    aFile = File.new(aFilename)
    array = Array.new
    string = String.new("")
    nameVector = Array.new
    aFile.each_line{ | line |
      if ! line.empty?
        if line =~ /^>.+/
          line.chomp!.slice!(/>/)
          nameVector.push(line)
          if string != ""
            array.push(string)
          end
          string = String.new("")
        else 
          string << line.chomp.gsub(/ /,"").upcase
        end
      end }
      array.push(string.upcase)
      self.matrix = array
      self.nameVector = nameVector
      self.filenameID= aFilename
      aFile.close
      return self
  end

    def initializeFromString(aString)
      if aString.class == String
        lines = aString.split("\n")
      else
        lines = aString
      end
      array = Array.new
      string = String.new("")
      nameVector = Array.new
      lines.each{ | line |
        if ! line.empty?
          if line =~ /^>.+/
            line.slice!(/>/)
            nameVector.push(line.chomp)
            if string != ""
              array.push(string)
            end
            string = String.new("")
          else 
            string << line.chomp.gsub(/ /,"").upcase
          end
        end }
        array.push(string.upcase)
        self.matrix = array
        self.nameVector = nameVector
        return self
    end

      def asCodingSequence(anArray,frame)
        if frame == nil
          self.readingFrame = 1
        else
          self.readingFrame = frame
        end
        if anArray != nil
          self.codingRegions = anArray
        else
          self.codingRegions = [1,self.matrix[0].length]
        end
        self.setGeneticCode("standard")
        self.setCodonDists
        self.setPUDict
        self.setSilentSiteDict
        self.setCodons
      end

      #	def initializeFromGenbank(aFilename)
      #		
      #		aFile = File.new(aFilename)
      #		array = Array.new
      #		string = String.new("")
      #		nameVector = Array.new
      #		self.codingRegions = Array.new
      #		featureDict = Hash.new
      #		aFile.each_line{ | line |
      #				if ! line.empty?
      #					if line.strip =~ /^CDS/
      #						line.chomp!.gsub(/d+/){ | match | self.codingRegions.push(match)
      #															print self.codingRegions}
      #					end	
      #					if line.strip =~ /^ORIGIN/
      #							array.push(string.upcase)
      #						end
      #						string = String.new("")
      #					else 
      #						string << line.chomp.gsub(" ","")
      #					end
      #				 
      #				}
      #		array.push(string.upcase)
      #		self.matrix = array
      #		self.nameVector = nameVector

      #		self.filenameID= aFilename
      #		self.readingFrame= 1
      #		aFile.close
      #		return self
      #	end


      def initializeWith(nameVector, matrix)
        self.matrix = matrix
        self.nameVector = nameVector
        self.codingRegions = Array.new
        return self
      end

      # manipulations
      def concatenate			#returns a new instance of a SequenceMatrix with all of the sequences of the original concatenated into one long sequence
        temp = String.new("")
        name = self.nameVector[0]
        self.matrix.each{ | eachString |
          temp << eachString
        }
        return SequenceMatrix.new.initializeWith([name],[temp])
      end

      def sequenceSubSet(anArray)			#takes an array of sites, returns a new instance of a SequenceMatrix with only those sites
        newMat = Array.new
        0.upto(self.sampleSize - 1){ | r |
          newSeq = String.new("")
          anArray.each{ | c |
            if c < self.length - 1
              newSeq << self.matrix[r][c,1]
            else
              newSeq << "N"
            end
          }
          newMat << newSeq
        }
        return SequenceMatrix.new.initializeWith(self.nameVector,newMat)
      end

      def sequenceMatrixSplit(anArray) #returns a new instance of a SequenceMatrix with only those seqs specified in array
        newMat = Array.new
        newNames = Array.new
        anArray.each{ | anIndex |
          if self.nameVector[anIndex] != nil
            newNames << self.nameVector[anIndex].clone
            newMat << self.matrix[anIndex].clone
          end
        }
        return(SequenceMatrix.new.initializeWith(newNames,newMat.flatten))  
      end

      def	sequenceSlice(anArray)			#meant to take a start and stop array
        newMat = Array.new
        0.upto(self.sampleSize - 1){ | r |
          l = anArray.last - anArray.first 
          newSeq = self.matrix[r][anArray.first, l]
          newMat << newSeq
        }
        return SequenceMatrix.new.initializeWith(self.nameVector,newMat)
      end

      def returnSequences(anArray)  #returns a new instance of sequenceMatrix with the sequences indicated in the array
        newMat = Array.new
        newNames = Array.new
        anArray.each{ | x | 
          newMat << self.matrix[x]
          newNames << self.nameVector[x]
        }
        return SequenceMatrix.new.initializeWith(newNames,newMat)
      end

      def sampleSize
        return @matrix.size
      end

      def length
        return self.matrix[0].size
      end

      def siteSet(anIndex)
        set = Array.new
        @matrix.each{ | aSeq |
          set.push(aSeq[anIndex,1])
        }
        set.uniq!
        return set
      end
      def siteSetClean(anIndex)		#cleans out ambiguous bases
        set = Array.new
        @matrix.each{ | aSeq |
          set.push(aSeq[anIndex,1])
        }
        set.uniq!
        set.delete("N")
        set.delete("-")
        return set
      end

      def siteArray(anIndex)
        array = Array.new
        @matrix.each{ | aSeq |
          array.push(aSeq[anIndex,1])
        }
        return array
      end

      def siteArrayClean(anIndex)
        array = Array.new
        @matrix.each{ | aSeq |
          array.push(aSeq[anIndex,1])
        }
        array.delete("N")
        array.delete("-")
        return array
      end

      def aminoSet(anIndex)
        set = Array.new
        self.codons.each{ | anAllele |
          if anAllele[anIndex] =~ /-/ 
            set.push("-")
          else 
            if anAllele[anIndex] =~ /N/
              set.push("N")
            else set.push(self.geneticCode[anAllele[anIndex]])
            end
          end
        }
        return set.uniq
      end
# end

    def aminoSetClean(anIndex)
      set = Array.new
      self.codons.each{ | anAllele |
        if ! (anAllele[anIndex] =~ /-/ or anAllele[anIndex] =~ /N/)
          set.push(self.geneticCode[anAllele[anIndex]])
        end
      }
      return set.uniq
    end

    def translate  # returns a new instance of sequence matrix
      self.setCodons
      trans = SequenceMatrix.new
      oc = Array.new
      self.codons.each { | anAllele |	
        temp = String.new("")
        anAllele.each{ | aCodon |
          if aCodon.include?("-") or aCodon.include?("N")
            temp << "-"
          else 
            if self.geneticCode[aCodon.upcase]
              temp << (self.geneticCode[aCodon.upcase])
            else 
              print aCodon
            end
          end }
          oc.push(temp)
        }
        oc.each{ | row |
          last = row.length - 1
          if row[last,1] == "-" or row[last,1]  == "*"
            row.slice!(last)
          end
        }
        trans.matrix = oc
        trans.nameVector = self.nameVector
        return trans
    end

      def has_stops?    #data check. translates and looks for stop codons
        trans = self.translate
        count = 0
        trans.matrix.each{ | r |
          count += r.count("*") }
          return count
      end

        def percentGaps			#this is for data checking, returns an array where each element is for an allele

          array = Array.new
          matrix.each{ | r |
            gaps = r.count("-")
            count = r.size.to_f
            array.push(gaps/count)
          }
          return array
        end

        def percentAlignment   #returns a float with the percentage of nucleotides aligned in all alleles
          count = 0
          length = matrix[0].length
          0.upto(length){ | i |
            if self.siteArray(i).include?("-")
              count += 1
            end
          }
          return 1.0 - (count.to_f / length.to_f)
        end

        def averageSampleSize   #returns a float with the average sample size - N's thrown away
          count = 0.0
          n = self.sampleSize
          length = self.length
          0.upto(length - 1){ | i |
            count += self.siteArrayClean(i).size
          }
          return count.to_f / length
        end	

        def distancesBetweenSegSites  	#returns an array of distances, caluculated right to left
          dists = Array.new
          oc = self.segSitesLocations
          oc.reverse!
          top = oc.size
          count = 0
          while count < top - 1
            dists.push((oc[count].to_f - oc[count + 1].to_f))
            count += 1
          end
          return dists
        end

        def getHaplotypes(aFlag) 		#returns an instance of sequenceMatrix with only the segregating sites included
          if aFlag.nil?
            locs = self.segSitesLocations
          else 
            locs = self.segSitesLocationsWithGaps
          end
          mat = Array.new
          0.upto(self.sampleSize - 1){ | i | mat.push(String.new("")) }
          locs.each{ | c |
            0.upto(self.sampleSize - 1){ | r |
              mat[r]<<(self.matrix[r])[c,1]
            }
          }
          newSeqMat = SequenceMatrix.new
          newSeqMat.nameVector = self.nameVector
          newSeqMat.matrix = mat
          newSeqMat.readingFrame = 1
          return newSeqMat
        end

        def getHaplotypesPretty(aFlag)	# returns an instance of sequenceMatrix with only the segregating sites but with .'s for the same as ref

          haps = self.getHaplotypes(aFlag)
          ref = haps.matrix[0]
          mat = Array.new			
          0.upto(haps.sampleSize - 1){ | i | mat.push(String.new("")) }
          mat[0] = ref
          0.upto(ref.length - 1){ | c |
            1.upto(haps.sampleSize - 1){ | r |
              if (haps.matrix[r])[c,1] == ref[c,1]
                mat[r]<< "."
              else
                mat[r]<<(haps.matrix[r])[c,1]
              end
            }
          }
          newSeqMat = SequenceMatrix.new
          newSeqMat.nameVector = haps.nameVector
          newSeqMat.matrix = mat
          newSeqMat.readingFrame = 1
          return newSeqMat
        end

        def majorAlleleFreqSite(aSite)
          tempArray = self.siteArrayClean(aSite)
          set = self.siteSetClean(aSite)
          n = tempArray.size.to_f
          major = 0.0
          set.each{ | each | 
            if tempArray.occurrencesOf(each) > major 
              major = tempArray.occurrencesOf(each)
            end
          }
          return major/n
        end


        def majorAlleleFreqSiteDirty(aSite)
          tempArray = self.siteArray(aSite)
          set = self.siteSet(aSite)
          n = tempArray.size.to_f
          major = 0.0
          set.each{ | each | 
            if tempArray.occurrencesOf(each) > major 
              major = tempArray.occurrencesOf(each)
            end
          }
          return major/n
        end

        def minorAlleleFreqSite(aSite)
          tempArray = self.siteArrayClean(aSite)
          set = self.siteSetClean(aSite)
          n = tempArray.size.to_f
          minor = self.sampleSize
          set.each{ | each | 
            if tempArray.occurrencesOf(each) < minor 
              minor = tempArray.occurrencesOf(each)
            end
          }
          return minor
        end
        def minorAlleleFreqSiteDirty(aSite)
          tempArray = self.siteArray(aSite)
          set = self.siteSet(aSite)
          n = tempArray.size.to_f
          minor = self.sampleSize
          set.each{ | each | 
            if tempArray.occurrencesOf(each) < minor 
              minor = tempArray.occurrencesOf(each)
            end
          }
          return minor/n
        end

        def majorAlleleNumberSite(aSite)
          tempArray = self.siteArray(aSite)
          set = self.siteSet(aSite)
          edit = tempArray.delete_if { |x | x =="N" }
          major = 0.0
          set.each{ | each | 
            if edit.occurrencesOf(each) > major 
              major = edit.occurrencesOf(each)
            end
          }
          return major
        end


        def majorAlleleStateSite(aSite)
          tempArray = self.siteArray(aSite)
          set = self.siteSet(aSite)
          edit = tempArray.delete_if { |x | x =="N" }
          n = edit.size.to_f
          major = 0
          state = nil
          set.each{ | each | 
            if edit.occurrencesOf(each) > major 
              major = edit.occurrencesOf(each)
              state = each
            end
          }
          return state
        end

        def stateNumberSite(aSite, aState)
          tempArray = self.siteArray(aSite)
          return tempArray.occurrencesOf(aState)
        end

        def rSquaredSites(aSite1, aSite2)

          pA = self.majorAlleleFreqSite(aSite1)
          stateA = self.majorAlleleStateSite(aSite1)
          pB = self.majorAlleleFreqSite(aSite2)
          stateB = self.majorAlleleStateSite(aSite2)

          array1 = self.siteArrayClean(aSite1)
          array2 = self.siteArrayClean(aSite2)
          count = 0.0
          siteCount = [array1.length,array2.length].min

          0.upto(self.sampleSize - 1) { | i |
            if (array1[i] == stateA) and (array2[i] == stateB)
              count += 1.0
            end
          }
          pAB = count/siteCount
          d = pAB - (pA * pB)
          denom = Math.sqrt(pA * (1 - pA) * pB * (1 - pB))
          r = d / denom
          return (r * r)
        end

        def reverseComplement   #returns sequence matrix

          anArray = Array.new
          self.matrix.each{ | string |
            temp = string.reverse.tr('ATGC','TACG')
            anArray.push(temp)
          }
          return SequenceMatrix.new.initializeWith(self.nameVector, anArray)
        end

        def relativePosition(aSequenceNumber, anIndex)			#returns the position of anIndex relative to the gaps inserted during alignment in sequence aSequenceNumber, useful for annotations but note it is zero based
          tempIndex = relIndex = 0
          string = self.matrix[aSequenceNumber]
          while tempIndex <= anIndex
            if string[relIndex,1] == "-"
              relIndex += 1
            else
              relIndex += 1
              tempIndex += 1
            end

          end
          return relIndex - 1
        end

        def padAlignment		#evens sequence length with N's from the back
          max = 0
          self.matrix.each{ | eachString |
            if eachString.length > max
              max = eachString.length
            end
          }
          self.matrix.each{ | eachString |
            if eachString.length < max
              l = max - eachString.length
              l.times{ | i | eachString << "N" }
            end
          }
        end

        def padAlignmentLength(anInt)		#adds sequence length with N's from the back to length anInt
          self.matrix.each{ | eachString |
            if eachString.length < anInt
              l = anInt - eachString.length
              l.times{ | i | eachString << "N" }
            end
          }
        end

        # conversions
        def codonsToPaml		#this attempts to cut off termination codons- still needs debugging
          fn = self.filenameID, ".paml"
          ws = File.new(fn.to_s, "w")
          n = self.sampleSize 
          size = self.codons[0].size - 1
          l = (size * 3) + 3
          if (self.aminoSet((size - 1))).include?("*")
            l = l - 3
            ws.print(n.to_s,"\t",l.to_s,"\n")
            self.codons.each_index{ | i |
              name = nameVector[i].to_s
              ws.print(name.slice(0,10),"  ")
              0.upto(self.codons[0].size - 2){ | c |
                ws.print(self.codons[i][c])
              }
              ws.print("\n") }
          else
              ws.print(n.to_s,"\t",l.to_s,"\n")
              self.codons.each_index{ | i |
                name = nameVector[i].to_s
                ws.print(name.slice(0,10),"  ")
                0.upto(self.codons[0].size - 1){ | c |
                  ws.print(self.codons[i][c])
                }
                ws.print("\n") }
          end
              ws.close
        end

            def sequenceNameHash		#returns a hash with the sequence names as keys and the sequence strings as values
              i = 0
              aHash = Hash.new
              self.nameVector.each{ | key |
                aHash[key] = self.matrix[i]
                i += 1
              }
              return aHash
            end

            def returnFasta    #spits out a new fasta file reflecting any changes that might have occured
              fn = self.filenameID, ".new"
              ws = File.new(fn.to_s, "w")
              self.nameVector.each_index{ | i |
                ws.print(">",nameVector[i],"\n",matrix[i],"\n")
              }
              ws.close
            end

            def returnFastaNamed(aName)    #spits out a new fasta file reflecting any changes that might have occured
              ws = File.new(aName, "w")
              self.nameVector.each_index{ | i |
                ws.print(">",nameVector[i],"\n",matrix[i],"\n")
              }
              ws.close
            end

            def outputFasta		
              i = 0
              self.nameVector.each{ | name |
                print ">",name,"\n"
                count = 60
                0.upto(self.matrix[i].length - 1){ | j |
                  print self.matrix[i][j,1]
                  count -= 1
                  if count == 0
                    print "\n"
                    count = 60
                  end
                }
                i += 1
                print "\n\n"
              }
            end


            # polymorphism stats

            def segSites   #this ignores N's but keeps -'s. will be different from estimator S's
              s = 0
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSet(count)
                temp.delete("N")
                if 1 < temp.size
                  s += 1
                end
                count += 1
              end
              return s
            end

            def segSitesKillN
              s = 0
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSet(count)
                if temp.include?("N") or temp.include?("-")
                  count += 1
                else
                  if 1 < temp.size
                    s += 1
                  end
                  count += 1
                end
              end
              return s
            end

            def sitesNArray   #returns array length n of number of sites at each sampleSize
              count = 0
              n = self.sampleSize
              oc = Array.new
              n.times{ | i | oc << 0 }
              while count < self.matrix[0].size  
                ni = self.siteArrayClean(count).size
                oc[ni - 1] += 1
                count += 1
              end
              return oc
            end    

            def segSitesNArray   #returns array length N of number of segSites at each sampleSize
              count = 0
              n = self.sampleSize
              oc = Array.new
              n.times{ | i | oc << 0 }
              while count < self.matrix[0].size  
                temp = self.siteSetClean(count)
                if 1 < temp.size
                  ni = self.siteArrayClean(count).size
                  oc[ni - 1] += 1
                end
                count += 1
              end
              return oc
            end

            def segSitesNArrayWindow(start, fin)   #returns array length N of number of segSites at each sampleSize, note zero indexing- doesn't look at seq at position "fin"
              count = start
              n = self.sampleSize
              oc = Array.new
              n.times{ | i | oc << 0 }
              while count < fin  
                temp = self.siteSetClean(count)
                if 1 < temp.size
                  ni = self.siteArrayClean(count).size
                  oc[ni - 1] += 1
                end
                count += 1
              end
              return oc
            end

            def segSitesLocations 			#returns array of locations from segSitesKillN zero indexed
              oc = Array.new				# note that this is for sites counted in estimators like theta 
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSet(count)
                if temp.include?("N") or temp.include?("-")
                  count += 1
                else
                  if 1 < temp.size
                    oc.push(count)
                  end
                  count += 1
                end
              end
              return oc
            end

            def segSitesLocationsWithGaps 			#returns array of locations of segSites with gaps includes (still no N's)
              oc = Array.new
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSet(count)
                if temp.include?("N")
                  count += 1
                else
                  if 1 < temp.size
                    oc.push(count)
                  end
                  count += 1
                end
              end
              return oc
            end

            def segSitesLocationsAllSites 			#returns array of locations from segSitesKillN zero indexed
              oc = Array.new				# note that this is for all columns (even those that contain -'s or N's)
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSetClean(count)
                if 1 < temp.size
                  oc.push(count)
                end
                count += 1			
              end
              return oc
            end

            def segSitesLocationsAllSitesNoSingletons 	#returns array of locations from segSitesKillN zero indexed
              oc = Array.new				# note that this is for all columns (even those that contain -'s or N's)
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSetClean(count)
                if 1 < temp.size and minorAlleleFreqSite(count) > 1
                  oc.push(count)
                end
                count += 1			
              end
              return oc
            end

            def segSitesLocationsAllSitesWithGaps 			# same as above but count gaps. returns array of locations from segSitesKillN zero indexed
              oc = Array.new								# note that this is for all columns (even those that contain -'s or N's)
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSet(count)
                temp.delete("N")
                if 1 < temp.size
                  oc.push(count)
                end
                count += 1

              end
              return oc
            end

            def fractionSegSitesKillN		#used in calculated Watterson's theta
              s = 0
              count = 0
              siteCount = 0
              while count < self.matrix[0].size  
                temp = self.siteSet(count)
                if temp.include?("N") or temp.include?("-")
                  count += 1
                else
                  siteCount += 1
                  if 1 < temp.size
                    s += 1
                  end
                  count += 1
                end
              end
              return s.to_f/siteCount
            end

            def SequenceMatrix.a1(n)			#this is the harmonic sum in Watterson's estimator
              sum = 0
              1.upto(n - 1){ | i | sum += 1.0/i }
              return sum
            end

            def thetaKillN 		#this is Watterson's estimator where sites which contain indels and N's are excluded
              s = self.fractionSegSitesKillN
              n = self.sampleSize
              if n < 2
                return 0.0
              else
                sum = 0.0
                count = 1
                while count < n
                  sum += (1.0/count)
                  count += 1
                end
              end
              return s/sum
            end

            def thetaMissingData		#adapted version of Watterson's estimator
              sArray = self.segSitesNArray
              l = self.length
              sum = 0
              1.upto(sArray.size - 1){ | i | sum += sArray[i].to_f / SequenceMatrix.a1(i + 1) }
              return sum / l	
            end	

            def thetaMissingDataWindow(start, fin)
              sArray = self.segSitesNArrayWindow(start, fin)
              l = fin - start
              sum = 0
              1.upto(sArray.size - 1){ | i | sum += sArray[i].to_f / SequenceMatrix.a1(i + 1) }
              return sum / l	
            end	

            def piKillN  #this is Tajima's (1983) estimator 
              n = self.sampleSize
              oc = Array.new
              l = self.matrix[0].size
              index = 0
              while index < self.matrix[0].size
                siteArray = self.siteArray(index)
                if siteArray.include?("-") or siteArray.include?("N")
                  l -= 1
                  index += 1
                else
                  siteSet = siteSet(index)
                  ni = siteArray.size.to_f
                  if siteSet.size > 1
                    pSum = 0.0
                    siteSet.each{ | state |
                      p = (siteArray.occurrencesOf(state)).to_f/ ni
                      pSum += (p * p)
                    }
                    oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                  end
                  index += 1
                end
              end
              return (oc.sum/l)
            end

            def piClean				 #pi where sample size varies per column, also cleans gaps
              oc = Array.new
              l = self.matrix[0].size
              index = 0
              while index < self.matrix[0].size
                siteArray = self.siteArrayClean(index)
                siteSet = siteArray.uniq
                ni = siteArray.size.to_f
                if siteSet.size > 1
                  pSum = 0.0
                  siteSet.each{ | state |
                    p = (siteArray.occurrencesOf(state)).to_f/ ni
                    pSum += (p * p)
                  }
                  oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                end
                index += 1

              end
              if oc.empty?
                return "0.0"
              else
                return (oc.sum/l)
              end
            end

            def piCleanNoSingletons				 #pi where sample size varies per column, also cleans gaps
              oc = Array.new
              l = self.matrix[0].size
              self.segSitesLocationsAllSitesNoSingletons.each{ | index |
                siteArray = self.siteArrayClean(index)
                siteSet = siteArray.uniq
                ni = siteArray.size.to_f
                if siteSet.size > 1
                  pSum = 0.0
                  siteSet.each{ | state |
                    p = (siteArray.occurrencesOf(state)).to_f/ ni
                    pSum += (p * p)
                  }
                  oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                end
              }
              if oc.empty?
                return "0.0"
              else
                return (oc.sum/l)
              end
            end

            def piCleanWindow(start, fin)				 #pi where sample size varies per column, also cleans gaps, zero indexing- doesn't look at site "fin"
              oc = Array.new
              l = fin - start 
              index = start
              nCount = 0
              while index < fin
                siteArray = self.siteArrayClean(index)
                if siteArray.empty?
                  nCount += 1
                else
                  siteSet = siteArray.uniq
                  ni = siteArray.size.to_f
                  if siteSet.size > 1
                    pSum = 0.0
                    siteSet.each{ | state |
                      p = (siteArray.occurrencesOf(state)).to_f/ ni
                      pSum += (p * p)
                    }
                    oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                  end
                end
                index += 1
              end
              if nCount > l * 0.5
                return "NA"
              else
                return (oc.sum/l)
              end
            end
            
            def tajimasK 	#this is pi total for Tajima's D
              n = self.sampleSize
              oc = Array.new
              l = self.matrix[0].size
              index = 0
              while index < self.matrix[0].size
                siteArray = self.siteArray(index)
                if siteArray.include?("-") or siteArray.include?("N")
                  index += 1
                else
                  siteSet = siteSet(index)
                  ni = siteArray.size.to_f
                  if siteSet.size > 1
                    pSum = 0.0
                    siteSet.each{ | state |
                      p = (siteArray.occurrencesOf(state)).to_f/ ni
                      pSum += (p * p)
                    }
                    oc.push((1.0 - pSum) * (ni/(ni - 1))) 
                  end
                  index += 1
                end
              end
              return oc.sum
            end

            def thetaH(state) 	#Fay and Wu's H estimator, variable is the derived state
              n = self.sampleSize
              l = self.matrix[0].size
              index = 0
              pSum = 0.0
              while index < self.matrix[0].size
                siteArray = self.siteArray(index)
                if siteArray.include?("-") or siteArray.include?("N")
                  index += 1
                else
                  siteSet = siteSet(index)
                  ni = siteArray.size.to_f
                  if siteSet.size > 1
                    p = (siteArray.occurrencesOf(state)).to_f
                    pSum += (p * p)
                  end
                  index += 1
                end
              end
              return (pSum * 2.0) / (n * (n-1.0))
            end

            def hFay(state) 	#Fay and Wu's H 
              n = self.sampleSize
              nnm1 = n / (n-1.0)
              l = self.matrix[0].size
              index = 0
              pSum = 0.0
              while index < self.matrix[0].size
                siteArray = self.siteArray(index)
                if siteArray.include?("-") or siteArray.include?("N")
                  index += 1
                else
                  siteSet = siteSet(index)
                  ni = siteArray.size.to_f
                  if siteSet.size > 1
                    p = (siteArray.occurrencesOf(state)).to_f / n
                    pSum += 2.0 * p * (2.0 * p - 1) * nnm1
                  end
                  index += 1
                end
              end
              return (-pSum)
            end

            def numberHaplotypes
              haps = self.getHaplotypes(nil)
              return haps.matrix.uniq.size
            end

            def haplotypeDiversity
              haps = self.getHaplotypes(nil)
              set = haps.matrix.uniq
              sum = 0.0
              set.each{ | string | 
                p = haps.matrix.select{ | x | x == string }.size.to_f / haps.sampleSize.to_f
                sum += p * p
              }
              return 1.0 - sum
            end

            #foldedSFS- returns array of site frequency spectrum assuming complete sample size
            def sfs
              s = self.segSitesLocationsAllSites
              sfs = Array.new
              0.upto(self.sampleSize - 1){ | i | sfs[i] = 0}
              s.each{ | i |
                sfs[self.minorAlleleFreqSite(i)] += 1
              }
              return sfs
            end

            def tajimasDKillN
              s = self.segSitesKillN
              if s < 3
                return "na"
              else
                n = self.sampleSize.to_f
                if n < 4
                  return "na"
                else
                  a1 = 0.0
                  a2 = 0.0
                  1.upto(n - 1){ | i |
                    a1 +=  (1.0/i)
                    a2 += (1.0/ (i*i))
                  }
                  b1 = (n + 1.0)/(3.0 * ( n - 1.0))
                  b2 = (2*( (n * n) + n + 3.0))/ (9.0 * n * (n - 1.0))
                  c1 = b1 - (1.0 / a1)
                  c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
                  e1 = c1 / a1
                  e2 = c2 / ((a1 * a1) + a2)
                  d = self.tajimasK - (s.to_f / a1)
                  stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
                  return d/stdev
                end
              end
            end

            def tajimasDHack
              s = self.segSitesLocationsAllSites.size
              if s < 3
                return "na"
              else
                n = self.sampleSize.to_f
                if n < 4
                  return "na"
                else
                  a1 = 0.0
                  a2 = 0.0
                  1.upto(n - 1){ | i |
                    a1 +=  (1.0/i)
                    a2 += (1.0/ (i*i))
                  }
                  b1 = (n + 1.0)/(3.0 * ( n - 1.0))
                  b2 = (2*( (n * n) + n + 3.0))/ (9.0 * n * (n - 1.0))
                  c1 = b1 - (1.0 / a1)
                  c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
                  e1 = c1 / a1
                  e2 = c2 / ((a1 * a1) + a2)
                  d = (self.piClean * self.length) - (s.to_f / a1)
                  stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
                  return d/stdev
                end
              end
            end

            def tajimasDHackNoSingletons
              s = self.segSitesLocationsAllSitesNoSingletons.size
              if s < 3
                return "na"
              else
                n = self.sampleSize.to_f
                if n < 4
                  return "na"
                else
                  a1 = 0.0
                  a2 = 0.0
                  1.upto(n - 1){ | i |
                    a1 +=  (1.0/i)
                    a2 += (1.0/ (i*i))
                  }
                  b1 = (n + 1.0)/(3.0 * ( n - 1.0))
                  b2 = (2*( (n * n) + n + 3.0))/ (9.0 * n * (n - 1.0))
                  c1 = b1 - (1.0 / a1)
                  c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
                  e1 = c1 / a1
                  e2 = c2 / ((a1 * a1) + a2)
                  d = (self.piCleanNoSingletons * self.length) - (s.to_f / a1)
                  stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
                  return d/stdev
                end
              end
            end

            def slidingWindowPi(windowSize, offset)    #note this uses piClean, prints array with [site, pi]
              i = 0
              l = self.length - windowSize
              oc = Array.new
              while i <= l
                oc << [i, self.piCleanWindow(i, i + windowSize)]
                i += offset
              end
              oc.reject!{ | each | each[1] == "NA" }
              oc.each{ | each | print each[0],"\t",each[1],"\n"}
            end

            def slidingWindowTheta(windowSize, offset)    #note this uses thetaMissingData, prints array with [site, theta]
              i = 0
              l = self.length - windowSize
              oc = Array.new
              while i <= l
                oc << [i, self.thetaMissingDataWindow(i, i + windowSize)]
                i += offset
              end
              oc.each{ | each | print each[0],"\t",each[1],"\n"}
            end

            def polySummary    #returns an array with the name of the file, n, length, S, Stotal, H, thetaKillN, piKillN, tajimasDKillN
              array = Array.new
              array.push(self.filenameID)
              array.push(self.sampleSize)
              array.push(self.matrix[0].length)
              array.push(self.segSitesKillN)
              array.push(self.segSites)
              array.push(self.numberHaplotypes)
              array.push(self.haplotypeDiversity)
              array.push(self.thetaKillN)
              array.push(self.piKillN)
              array.push(self.tajimasDKillN)
              return array
            end

            def printPolySummary(anInt)
              if anInt.nil?
                sum = self.polySummary
                print "file ID \tN \tbp \tS \tStotal \tH \thaploDiveristy \tthetaKillN \tpiKillN \ttajimasDKillN\n"
                sum.each{| x | print x,"\t" }
                print "\n"
              else 
                if anInt == "-p"
                sum = self.polySummary
                string = "ms ",sum[1].to_s," 10000 -s ",sum[3].to_s," | sample_stats > ~/rubyStuff/simTemp"
                #	print "now running simulations: ",string,"\n"
                system(string.to_s)
                aFile = File.new("/Users/adk/rubyStuff/simTemp","r")
                lines = aFile.readlines
                array = Array.new
                lines.each{| x | tokens = x.chomp.split 
                  array.push(tokens[5].to_f)}
                  array.sort!
                  sum.push(array[499])
                  sum.push(array[9499])
                  print "file ID \tN \tbp \tS \tStotal \tH \thaploDiveristy  \tthetaKillN \tpiKillN \ttajimasDKillN\tlower 95% d\tupper 95% d\n"
                  sum.each{| x | print x,"\t" }
                  print "\n"
                  `rm ~/rubyStuff/simTemp`
                end

              end
            end	

            #end	


            #
            # divergence (still not sure if this should be moved to separate class)
            #

            def averagePairwiseDifferencesPerSite			#assumes all seqs are the same length, kills "-"s and "N"s

              difPerSite = Array.new
              n = self.sampleSize - 1 
              width = self.matrix[0].size - 1
              0.upto(n - 1){ | i |
                (i+1).upto(n){ | j |
                  siteCount = 0
                  diffs = 0
                  0.upto(width){ | c |
                    testArray = [self.matrix[i][c,1],self.matrix[j][c,1]]
                    if ! testArray.include?("-") or ! testArray.include?("N")
                      siteCount += 1
                      if testArray.uniq.size > 1
                        diffs += 1
                      end
                    end
                  }
                  difPerSite << diffs.to_f / siteCount
                }
              }
              return difPerSite.sampleMean
            end

            def averagePairwiseDifferencesArray			#assumes all seqs are the same length, kills "-"s and "N"s

              difPerSite = Array.new
              n = self.sampleSize - 1
              width = self.matrix[0].size - 1
              0.upto(n - 1){ | i |
                (i+1).upto(n){ | j |
                  siteCount = 0
                  diffs = 0
                  0.upto(width){ | c |
                    testArray = [self.matrix[i][c,1],self.matrix[j][c,1]]
                    if ! testArray.include?("-") and ! testArray.include?("N")
                      siteCount += 1
                      if testArray.uniq.size > 1
                        diffs += 1
                      end
                    end
                  }
                  difPerSite << diffs.to_f 
                }
              }
              return difPerSite
            end

            def averagePairwiseDifferencesPerSiteSites(anArray)			#assumes all seqs are the same length, doesn't kill N's but kills "-"s

              difPerSite = Array.new
              n = self.sampleSize - 1 
              0.upto(n - 1){ | i |
                (i+1).upto(n){ | j |
                  siteCount = 0
                  diffs = 0
                  anArray.each{ | c |
                    testArray = [self.matrix[i][c,1],self.matrix[j][c,1]]
                    if ! testArray.include?("-")
                      siteCount += 1
                      if testArray.uniq.size > 1
                        diffs += 1
                      end
                    end
                  }
                  if siteCount == 0
                    difPerSite << "NA"
                  else
                    difPerSite << diffs.to_f/siteCount
                  end
                }
              }
              temp = difPerSite.checkArray
              if temp
                return temp.sampleMean
              else
                return "NA"
              end
            end

            def averagePairwiseDifferencesSites(anArray)			#assumes all seqs are the same length, doesn't kill N's but kills "-"s

              difPerSite = Array.new
              n = self.sampleSize - 1 
              0.upto(n - 1){ | i |
                (i+1).upto(n){ | j |
                  diffs = 0
                  anArray.each{ | c |
                    testArray = [self.matrix[i][c,1],self.matrix[j][c,1]]
                    if ! testArray.include?("-") or ! testArray.include?("N")
                      if testArray.uniq.size > 1
                        diffs += 1
                      end
                    end
                  }
                  difPerSite << diffs.to_f
                }
              }
              return difPerSite.sampleMean
            end


            #
            # stats 
            #
            #

            def pairwiseRSquaredByDistance		#returns an Array

              temp = self.segSitesLocationsAllSitesWithGaps

              sites = temp.delete_if{ | x | self.majorAlleleFreqSite(x) > 0.95 }
              n = sites.size - 1
              counter = 0
              anArray = Array.new
              i = 0
              while i < n
                j = i + 1
                while j <= n
                  print "comparison ",counter.to_s,"\n"
                  counter+=1
                  temp = Array.new
                  temp.push(self.rSquaredSites(sites[i],sites[j]))
                  temp.push(sites[j] - sites[i])
                  anArray.push(temp)
                  j += 1
                end
                i += 1
              end
              return anArray
            end

            def pairwiseRSquared		#returns lower diag. matrix
              temp = self.segSitesLocationsAllSitesWithGaps
              sites = temp.delete_if{ | x | self.majorAlleleFreqSite(x) > 0.95 }
              n = sites.size - 1
              counter = 0
              anArray = Array.new
              i = 0
              while i < n
                j = i + 1
                temp = Array.new(sites.length,"NA")
                while j <= n    
                  temp[j] = self.rSquaredSites(sites[i],sites[j])
                  #temp.push(sites[j] - sites[i])
                  j += 1
                end
                anArray.push(temp)
                i += 1
              end
              return anArray
            end


            ###############################################	
            # Protein Sequence Stuff	
            def segSitesLocationsAA 			#returns array of locations from segSitesKillN zero indexed
              oc = Array.new				# note that this is for sites counted in estimators like theta 
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSet(count)
                if temp.include?("X") or temp.include?("-")
                  count += 1
                else
                  if 1 < temp.size
                    oc.push(count)
                  end
                  count += 1
                end
              end
              return oc
            end

            def segSitesLocationsWithGapsAA 			#returns array of locations of segSites with gaps includes (still no N's)
              oc = Array.new
              count = 0
              while count < self.matrix[0].size  
                temp = self.siteSet(count)
                if temp.include?("X")
                  count += 1
                else
                  if 1 < temp.size
                    oc.push(count)
                  end
                  count += 1
                end
              end
              return oc
            end	

            def getHaplotypesAA(aFlag) 		#returns an instance of sequenceMatrix with only the segregating sites included
              if aFlag.nil?
                locs = self.segSitesLocationsAA
              else 
                locs = self.segSitesLocationsWithGapsAA
              end
              mat = Array.new
              0.upto(self.sampleSize - 1){ | i | mat.push(String.new("")) }
              locs.each{ | c |
                0.upto(self.sampleSize - 1){ | r |
                  mat[r]<<(self.matrix[r])[c,1]
                }
              }
              newSeqMat = SequenceMatrix.new
              newSeqMat.nameVector = self.nameVector
              newSeqMat.matrix = mat
              newSeqMat.readingFrame = 1
              return newSeqMat
            end	

###########################################
###   SFS stuff
##
            ##unfoldedSFS-- for ms output
            def unfoldedSFSHudson
              s = self.segSitesLocationsAllSites
              sfs = Array.new
              0.upto(self.sampleSize - 1){ | i | sfs[i] = 0}
              s.each{ | i |
                sfs[self.stateNumberSite(i,"1")] += 1
              }
              return sfs
            end
            
            ##Uses Achaz's "system" to calculate theta_w based on sfs
            def sfs2ThetaW(sfsArray)
              weightSum = 0.0
              n = sfsArray.length
              1.upto(n-1){|i|
                weightSum += (1.0 / i)
                 }
              thetaW = 0.0
              1.upto(n-1){|i|
                thetaW += sfsArray[i] * i * (1.0 / i)
                 }
              return(thetaW * (1.0 / weightSum))
            end
            
            ##Uses Achaz's "system" to calculate theta_w based on sfs
            # doesn't include singletons
            def sfs2ThetaWNoSingletons(sfsArray)
              weightSum = 0.0
              n = sfsArray.length
              2.upto(n-1){|i|
                weightSum += (1.0 / i)
                 }
              thetaW = 0.0
              2.upto(n-1){|i|
                thetaW += sfsArray[i] * i * (1.0 / i)
                 }
              return(thetaW * (1.0 / weightSum))
            end
            
            ##Uses Achaz's "system" to calculate theta_pi based on sfs
            # doesn't include singletons
            def sfs2ThetaPiNoSingletons(sfsArray)
              weightSum = 0.0
              n = sfsArray.length
              2.upto(n-1){|i|
                weightSum += (n - i)
                 }
              thetaW = 0.0
              2.upto(n-1){|i|
                thetaW += sfsArray[i] * i * (n- i)
                 }
              return(thetaW * (1.0 / weightSum))
            end
            ##Uses Achaz's "system" to calculate theta_h based on sfs
            # doesn't include singletons
            def sfs2ThetaHNoSingletons(sfsArray)
              weightSum = 0.0
              n = sfsArray.length
              2.upto(n-1){|i|
                weightSum +=  i
                 }
              thetaW = 0.0
              2.upto(n-1){|i|
                thetaW += sfsArray[i] * i *  i
                 }
              return(thetaW * (1.0 / weightSum))
            end
            
            ##Uses Achaz's "system" to calculate Fay and Wu's H based on sfs
            # doesn't include singletons
            def sfs2FayWuHNoSingletons(sfsArray)
              return(self.sfs2ThetaPiNoSingletons(sfsArray) - self.sfs2ThetaHNoSingletons(sfsArray) )
            end
            
            ##Uses Achaz's "system" to calculate theta_h based on sfs
            # doesn't include singletons
            def sfs2TajDNoSingletons(sfsArray)
              return(self.sfs2ThetaPiNoSingletons(sfsArray) - self.sfs2ThetaWNoSingletons(sfsArray) )
            end
end				



