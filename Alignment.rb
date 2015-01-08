#Alignment.rb- a class for divergence statistics
#
#
class Alignment
  attr_reader :sequenceMatrix1, :sequenceMatrix2, :outgroup, :codingRegions, :readingFrame, :geneticCode
  attr_writer :sequenceMatrix1, :sequenceMatrix2, :outgroup, :codingRegions, :readingFrame, :geneticCode

  def initialize(seqMat1, seqMat2, outG)
    self.sequenceMatrix1 = seqMat1
    self.sequenceMatrix2 = seqMat2
    self.outgroup = outG
    self.readingFrame = seqMat1.readingFrame		#inherits annotations from seqMat1
    self.codingRegions = seqMat1.codingRegions
    self.geneticCode = seqMat1.geneticCode
  end


  # manipulations
  def length
    return self.sequenceMatrix1.matrix[0].length
  end

  def alignSet(aSite)
    a = Array.new
    self.sequenceMatrix1.siteSet(aSite).each{ | x | a << x }
    self.sequenceMatrix2.siteSet(aSite).each{ | x | a << x }
    return a.uniq
  end

  def alignCodonSet(aCodonIndex)
    a = Array.new
    self.sequenceMatrix1.codonSet(aCodonIndex).each{ | x | a << x }
    self.sequenceMatrix2.codonSet(aCodonIndex).each{ | x | a << x }
    return a.uniq
  end


  def alignCodonSetPlus(aCodonIndex)
    a = alignCodonSet(aCodonIndex)
    a << self.outgroup.codonSet(aCodonIndex).first
    return a.uniq
  end

  def alignAminoSet(aCodonIndex)
    c = self.alignCodonSet(aCodonIndex)
    c.reject!{ | x | x.include?("N")}
    temp = Array.new
    c.each{ | eachCodon |
      temp << self.geneticCode[eachCodon]
    }
    return temp
  end

  def alignAminoSetPlus(aCodonIndex)
    a = self.alignAminoSet(aCodonIndex)
    a << self.outgroup.aminoSet(aCodonIndex).first
    return a.uniq
  end

  def alignAminoSetPlusClean(aCodonIndex)
    a = self.alignAminoSetPlus(aCodonIndex)
    a.reject!{ | x | x == ("-") }
    return a.uniq
  end




  #
  # divergence
  #
  def countTransitions					#gets counts of differences according to basepair changes, only uses first sequence in each matrix
    aHash = Hash.new
    params = ["AT","AC","AG","CG","CT","GT","AA","CC","TT","GG"]
    params.each{ | each | aHash[each] = 0 }
    l = self.sequenceMatrix1.length - 1
    0.upto(l){ | i | 
      temp = Array.new
      temp << self.sequenceMatrix1.matrix[0][i,1]
      temp << self.sequenceMatrix2.matrix[0][i,1]
      string = temp.sort.to_s
      if params.include?(string)
        aHash[string] += 1
      end
    }
    params.each{ | each |
      print each,"\t",aHash[each],"\n"
    }
  end


  def fixedDiffSites(sampleSizeFilter)			#returns a dict with fixes and gaps/ambiguous
    fixes = Array.new
    gaps = Array.new
    width = self.sequenceMatrix1.length
    height1 = self.sequenceMatrix1.sampleSize
    height2 = self.sequenceMatrix2.sampleSize
    0.upto(width - 1){ | c |
      sampSize = self.sequenceMatrix1.siteArrayClean(c).size
      if sampSize > sampleSizeFilter
        #use cleaned up data
        testSet1 = self.sequenceMatrix1.siteSet(c).reject{| x | x == "N"}
        testSet2 = self.sequenceMatrix2.siteSet(c).reject{| x | x == "N"}
        if testSet1.empty? or testSet2.empty?
          gaps << c
        else
          if testSet1.size == 1 and testSet2.size == 1 						#both sites fixed?
            if  (testSet1.include?("-") or testSet2.include?("-"))			#include gaps? if not...
              gaps << c	
            else
              if testSet1 != testSet2	  #are the sets equal? if not...
                fixes << c
              end
            end
          end
        end
      end
    }
    dict = Hash.new
    dict["fixes"] = fixes
    dict["gaps"] = gaps
    return dict
  end

  def fixedDiffSitesWindow(start, fin)			#returns a dict with fixes and gaps/ambiguous
    fixes = Array.new
    gaps = Array.new
    height1 = self.sequenceMatrix1.sampleSize
    height2 = self.sequenceMatrix2.sampleSize
    start.upto(fin - 1){ | c |
      testSet1 = self.sequenceMatrix1.siteSet(c)
      testSet2 = self.sequenceMatrix2.siteSet(c)
      #clean up ambiguous sites
      testSet1.delete("N")
      testSet2.delete("N")
      if testSet1.empty? or testSet2.empty?
        gaps << c
      else
        if testSet1.size == 1 and testSet2.size == 1 						#both sites fixed?
          if  (testSet1.include?("-") or testSet2.include?("-"))			#include gaps? if not...
            gaps << c	
          else
            if testSet1 != testSet2										#are the sets equal? if not...
              fixes << c
            end
          end
        end
      end
    }
    dict = Hash.new
    dict["fixes"] = fixes
    dict["gaps"] = gaps
    return dict
  end

  def fixedDiffsWindow(start, fin)										#returns number of fixed differences
    fixes = 0
    nCount = 0
    start.upto(fin-1){ | c |
      testSet1 = self.sequenceMatrix1.siteSet(c)
      testSet2 = self.sequenceMatrix2.siteSet(c)
      #clean up ambiguous sites
      testSet1.delete("N")
      testSet2.delete("N")
      testSet1.delete("-")
      testSet2.delete("-")
      if testSet1.empty? or testSet2.empty?
        nCount += 1
      else
        if testSet1.size == 1 and testSet2.size == 1 						#both sites fixed?
         # print testSet1[0],"\t",testSet2[0],"\n"
            if testSet1[0] != testSet2[0]										#are the sets equal? if not...
              fixes += 1
            end
        end
      end
    }
    if nCount > ((fin - start) * 0.9)
      return "NA"
    else
      return fixes
    end
  end

  def fixedDiffsWindowPer(start, fin)										#returns number of fixed differences
    fixes = 0
    nCount = 0
    bi = 0
    start.upto(fin - 1){ | c |
      testSet1 = self.sequenceMatrix1.siteSet(c)
      testSet2 = self.sequenceMatrix2.siteSet(c)
      #clean up ambiguous sites
      testSet1.delete("N")
      testSet2.delete("N")
      if testSet1.empty? or testSet2.empty?
        nCount += 1
      else
        if testSet1.size == 1 and testSet2.size == 1 						#both sites fixed?
          if  (testSet1.include?("-") or testSet2.include?("-"))			#include gaps? if not...
          else
            if testSet1 != testSet2										#are the sets equal? if not...
              fixes += 1
            end
            bi += 1
          end
        end
      end
    }
    if nCount > ((fin - start) * 0.5)
      return "NA"
    else
      return fixes.to_f / bi
    end
  end

  def fixedDiffCodingDict
    pos = Hash.new
    fixDict = self.fixedDiffSites(1)   #currently only finds fixations with a sample size of 2 or greater in seqmat 1
    fixSites = fixDict["fixes"]
    codingSites = self.sequenceMatrix1.codingSites
    fixSites = fixSites.select{ | each | codingSites.include?(each) }
    codonPos = Array.new
    fixSites.each{ | eachSite | codonPos << (codingSites.index(eachSite) / 3) }
    pos["codonFixes"] = codonPos
    pos["locusFixes"] = fixSites
    pos["gaps"] = fixDict["gaps"]
    return pos
  end

  def silReplFixations					#returns a hash of silent and replacement fixations
    fixes = Hash.new
    ["replacements","silents","missingData"].each{ | x | fixes[x] = Array.new }
    dict = self.fixedDiffCodingDict
    codingSites = self.sequenceMatrix1.codingSites
    lastCodon = nil
    dict["codonFixes"].each_with_index{ | site, index |
      if site != lastCodon			#fixes the problem of multiple substitutions per codon
        lastCodon = site
        testSet = self.alignCodonSet(site)
        if testSet.any?{ | x | x.include?("N") or x.include?("-") }		#clean up ambiguous data
          fixes["missingData"] << dict["locusFixes"][index]
        else
          if self.sequenceMatrix1.codonSegSitesPositions(testSet).size == 1		#one fixed diff in the codon?
            amino1 = self.sequenceMatrix1.aminoSet(site)
            amino2 = self.sequenceMatrix2.aminoSet(site)
            if amino1.size < amino2.size
              if amino1.any?{ | aa | amino2.include?(aa) }
                fixes["silents"] << dict["locusFixes"].uniq[index]	
              else
                fixes["replacements"] << dict["locusFixes"].uniq[index]
              end
            else
              if amino2.any?{ | aa | amino1.include?(aa) }
                fixes["silents"] << dict["locusFixes"].uniq[index]
              else
                fixes["replacements"] << dict["locusFixes"].uniq[index]
              end
            end
          else				#more than one fixed difference in codon, use paths
            paths = Array.new
            set1 = self.sequenceMatrix1.codonSet(site)
            set2 = self.sequenceMatrix2.codonSet(site)
            set1.each{ | codon1 |
              set2.each{ | codon2 |
                print codon1,"\t",codon2,"\n"
                paths << self.sequenceMatrix1.lookupPath2Codons([codon1,codon2])
              }
            }
            #set up the path length variable
            short = 100
            shortPath = nil
            paths = paths.reject{ | x | x.nil? or x.include?("*") }

            #go through paths, find shortest with respect to Replacements
            paths.each{ | aPath |
              rLength = 0
              aPath.each{ | eachArray | if eachArray.include?("R")
                rLength += 1 
              end
            }
            if short > rLength
              short = rLength
              shortPath = aPath
            else
              if shortPath.size > aPath.size		# more replacements but less changes overall
                shortPath = aPath
              end
            end
          }

          shortPath.each{ | eachArray |
            if eachArray.first == "R"
              fixes["replacements"] << codingSites[ (site * 3) + eachArray.last.to_i]
            else
              if eachArray.first == "S"
                fixes["silents"] << codingSites[ (site * 3) + eachArray.last.to_i]
              end
            end
          }
        end
      end
    end
  }
  fixes["missingData"] = fixes["missingData"] | dict["gaps"]
  return fixes
end

def silReplFixationsGreedy	#returns a hash of silent and replacement fixations - greedy cleans out N's and Gaps
  fixes = Hash.new
  ["replacements","silents","missingData"].each{ | x | fixes[x] = Array.new }
  dict = self.fixedDiffCodingDict
  codingSites = self.sequenceMatrix1.codingSites
  lastCodon = nil
  dict["codonFixes"].each_with_index{ | site, index |
    if site != lastCodon		#fixes the problem of multiple substitutions per codon
      lastCodon = site
      testSet = self.alignCodonSet(site)
      testSet.reject!{ | x | x.include?("N") or x.include?("-") }		#clean up ambiguous data
      if ! testSet.empty?
        if self.sequenceMatrix1.codonSegSitesPositions(testSet).size == 1		#one fixed diff in the codon?
          amino1 = self.sequenceMatrix1.aminoSetClean(site)
          amino2 = self.sequenceMatrix2.aminoSetClean(site)
          if amino1.size < amino2.size
            if amino1.any?{ | aa | amino2.include?(aa) }
              fixes["silents"] << dict["locusFixes"].uniq[index]	
            else
              fixes["replacements"] << dict["locusFixes"].uniq[index]
            end
          else
            if amino2.any?{ | aa | amino1.include?(aa) }
              fixes["silents"] << dict["locusFixes"].uniq[index]
            else
              fixes["replacements"] << dict["locusFixes"].uniq[index]
            end
          end
        else			#more than one fixed difference in codon, use paths
          paths = Array.new
          set1 = self.sequenceMatrix1.codonSetClean(site)
          set2 = self.sequenceMatrix2.codonSetClean(site)
          set1.each{ | codon1 |
            set2.each{ | codon2 |
              paths << self.sequenceMatrix1.lookupPath2Codons([codon1,codon2])
            }
          }
          #set up the path length variable
          short = 100
          shortPath = nil
          paths.reject!{ | x | x.nil? }
          paths.reject!{  | x | x.include?("*") }

          #go through paths, find shortest with respect to Replacements
          paths.each{ | aPath |
            rLength = 0
            aPath.each{ | eachArray | if eachArray.include?("R")
              rLength += 1 
            end
          }
          if short > rLength
            short = rLength
            shortPath = aPath
          else
            if shortPath.size > aPath.size		# more replacements but less changes overall
              shortPath = aPath
            end
          end
        }
        if shortPath != nil							#bug here!
          shortPath.each{ | eachArray |
            if eachArray.first == "R"
              fixes["replacements"] << codingSites[ (site * 3) + eachArray.last.to_i]
            else
              if eachArray.first == "S"
                fixes["silents"] << codingSites[ (site * 3) + eachArray.last.to_i]
              end
            end

          }
        end
      end
    end
  end
}
fixes["missingData"] = fixes["missingData"] | dict["gaps"]
return fixes
end
def polarizeAllFixations			#returns a hash of polarized fixations, designed for one outgroup
  polarized = Hash.new
  temp1 = Array.new
  temp2 = Array.new
  notPol = Array.new
  sites = self.fixedDiffSites(1)									#use fixed differences between ingroups as starting material
  sites["fixes"].each{ | eachSite |
    os = self.outgroup.siteSet(eachSite)
    if os != ["-"] 
      set1 = self.sequenceMatrix1.siteSet(eachSite)			
      set2 = self.sequenceMatrix2.siteSet(eachSite)			#skip those sites which have a gap in the outgroup
      if os.any?{ | each | set1.include?(each) }				#outgroup includes set1?
        if os.any?{ | each | set2.include?(each) }			#outgroup includes set2?
          notPol << eachSite
        else
          temp2 << eachSite								# fixation along lineage 2
        end
      else
        if os.any?{ | each | set2.include?(each) }								
          temp1 << eachSite
        end
      end
    else
      notPol << eachSite
    end
  }
  polarized["species1"] = temp1
  polarized["species2"] = temp2
  polarized["notPol"] = notPol
  return polarized
end

def polarizeAllFixationsGreedy			#clean 'em up! ## returns a hash of polarized fixations, designed for one outgroup
  polarized = Hash.new
  temp1 = Array.new
  temp2 = Array.new
  notPol = Array.new
  sites = self.fixedDiffSites(1)									#use fixed differences between ingroups as starting material
  sites["fixes"].each{ | eachSite |
    os = self.outgroup.siteSetClean(eachSite)
    if os != ["-"] 
      set1 = self.sequenceMatrix1.siteSetClean(eachSite)			
      set2 = self.sequenceMatrix2.siteSetClean(eachSite)			#skip those sites which have a gap in the outgroup
      if os.any?{ | each | set1.include?(each) }				#outgroup includes set1?
        if os.any?{ | each | set2.include?(each) }			#outgroup includes set2?
          notPol << eachSite
        else
          temp2 << eachSite								# fixation along lineage 2
        end
      else
        if os.any?{ | each | set2.include?(each) }								
          temp1 << eachSite
        end
      end
    else
      notPol << eachSite
    end
  }
  polarized["species1"] = temp1
  polarized["species2"] = temp2
  polarized["notPol"] = notPol
  return polarized
end
def polarizeCodingFixations
  pol = Hash.new
  polFixes = self.polarizeAllFixations											#use all polarized fixations
  srHash = self.silReplFixations
  codingSites = self.sequenceMatrix1.codingSites
  pol["notPol"] = Array.new
  ["species1","species2"].each{ | name |											#set up pol hash with an internal hash
    aHash = Hash.new
    ["silents","replacements","preferreds","unpreferreds","noChange"].each{ | x | aHash[x] = Array.new }
    pol[name] = aHash
  }

  polFixes["species1"].each{ | eachFix |
    if srHash["replacements"].include?(eachFix)									#silreplHash has site as repl
      if self.alignAminoSetPlus(codingSites.index(eachFix) / 3).size < 3	#aminoSet can be polarized?
        pol["species1"]["replacements"] << eachFix									
      else
        pol["notPol"] << eachFix
      end
    else
      if srHash["silents"].include?(eachFix)
        pol["species1"]["silents"] << eachFix								
        #move to pref/unpref analysis
        puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachFix) / 3)	
        puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachFix) / 3)
        puSetO = self.outgroup.puSet(codingSites.index(eachFix) / 3)
        if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
          if puSet2 == puSet1
            pol["species1"]["noChange"] << eachFix
          else
            if ! puSetO.include?(puSet1.first)						#is it derived?
              if puSet1.first == 1
                pol["species1"]["preferreds"] << eachFix
              else
                pol["species1"]["unpreferreds"] << eachFix
              end
            end
          end
        end
      end
    end
  }
  polFixes["species2"].each{ | eachFix |
    if srHash["replacements"].include?(eachFix)									#silreplHash has site as repl
      if self.alignAminoSetPlus( codingSites.index(eachFix) / 3 ).size < 3	#aminoSet can be polarized?
        pol["species2"]["replacements"] << eachFix		
      else
        pol["notPol"] << eachFix
      end
    else
      if srHash["silents"].include?(eachFix)
        pol["species2"]["silents"] << eachFix								
        #move to pref/unpref analysis
        puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachFix) / 3)	
        puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachFix) / 3)
        puSetO = self.outgroup.puSet(codingSites.index(eachFix) / 3)
        if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
          if puSet2 == puSet1
            pol["species2"]["noChange"] << eachFix
          else
            if ! puSetO.include?(puSet2.first)						#is it derived?
              if puSet2.first == 1
                pol["species2"]["preferreds"] << eachFix
              else
                pol["species2"]["unpreferreds"] << eachFix
              end
            end
          end
        end
      end
    end
  }
  pol["notPol"] = pol["notPol"] | polFixes["notPol"]
  return pol
end

def polarizeCodingFixationsGreedy			#clean 'em!
  pol = Hash.new
  polFixes = self.polarizeAllFixationsGreedy				#use all polarized fixations
  srHash = self.silReplFixationsGreedy
  codingSites = self.sequenceMatrix1.codingSites
  pol["notPol"] = Array.new
  ["species1","species2"].each{ | name |				#set up pol hash with an internal hash
    aHash = Hash.new
    ["silents","replacements","preferreds","unpreferreds","noChange"].each{ | x | aHash[x] = Array.new }
    pol[name] = aHash
  }

  polFixes["species1"].each{ | eachFix |
    if srHash["replacements"].include?(eachFix)									#silreplHash has site as repl
      if self.alignAminoSetPlusClean(codingSites.index(eachFix) / 3).size < 3	#aminoSet can be polarized?
        pol["species1"]["replacements"] << eachFix									
      else
        pol["notPol"] << eachFix
      end
    else
      if srHash["silents"].include?(eachFix)
        pol["species1"]["silents"] << eachFix								
        #move to pref/unpref analysis
        puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachFix) / 3)	
        puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachFix) / 3)
        puSetO = self.outgroup.puSet(codingSites.index(eachFix) / 3)
        if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
          if puSet2 == puSet1
            pol["species1"]["noChange"] << eachFix
          else
            if ! puSetO.include?(puSet1.first)						#is it derived?
              if puSet1.first == 1
                pol["species1"]["preferreds"] << eachFix
              else
                pol["species1"]["unpreferreds"] << eachFix
              end
            end
          end
        end
      end
    end
  }
  polFixes["species2"].each{ | eachFix |
    if srHash["replacements"].include?(eachFix)									#silreplHash has site as repl
      if self.alignAminoSetPlusClean( codingSites.index(eachFix) / 3 ).size < 3	#aminoSet can be polarized?
        pol["species2"]["replacements"] << eachFix		
      else
        pol["notPol"] << eachFix
      end
    else
      if srHash["silents"].include?(eachFix)
        pol["species2"]["silents"] << eachFix								
        #move to pref/unpref analysis
        puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachFix) / 3)	
        puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachFix) / 3)
        puSetO = self.outgroup.puSet(codingSites.index(eachFix) / 3)
        if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
          if puSet2 == puSet1
            pol["species2"]["noChange"] << eachFix
          else
            if ! puSetO.include?(puSet2.first)						#is it derived?
              if puSet2.first == 1
                pol["species2"]["preferreds"] << eachFix
              else
                pol["species2"]["unpreferreds"] << eachFix
              end
            end
          end
        end
      end
    end
  }
  pol["notPol"] = pol["notPol"] | polFixes["notPol"]
  return pol
end

#
#polymorphism
#

def silReplPolymorphismWithin		#returns a hash with silent and replacement polymorphisms that are specific to one of the seq mats
  polys = Hash.new
  ["replacements","silents","missingData"].each{|each | polys[each] = Array.new }
  codingSites = self.sequenceMatrix1.codingSites
  checkArray = Array.new
  all = Array.new
  hash1 = self.sequenceMatrix1.silReplMutations
  hash2 = self.sequenceMatrix2.silReplMutations
  hash1["silents"].each{ | each | 
    all << each
    if hash2["silents"].include?(each)
      checkArray << each
    end
  }
  hash2["silents"].each{ | each | all << each }
  checkArray.uniq.each{ | x | all.delete(x) }
  checkArray.uniq.each{ | site | 
    codon1 = self.sequenceMatrix1.codonSet( codingSites.index(site) / 3 )
    codon2 = self.sequenceMatrix2.codonSet( codingSites.index(site)/ 3 )
    if codon1.size == codon2.size
      if codon1.all?{ | each | codon2.include?(each) }
        hash1["silents"].occurrencesOf(site).times{ all << site }
      else
        hash1["silents"].occurrencesOf(site).times{ all << site }
        hash2["silents"].occurrencesOf(site).times{ all << site }
      end
    else
      siteSet1 = self.sequenceMatrix1.siteSet(site)
      siteSet2 = self.sequenceMatrix2.siteSet(site)
      bothSet = Array.new
      siteSet1.each{ | eachSite | bothSet << eachSite }
      siteSet2.each{ | eachSite | bothSet << eachSite }
      (bothSet.uniq.size - 1).times{ all << site }
    end
  }
  all.each{ | eachSite | polys["silents"] << eachSite }

  #add replacements

  checkArray = Array.new
  all = Array.new
  hash1["replacements"].each{ | each | 
    all << each
    if hash2["replacements"].include?(each)
      checkArray << each
    end
  }
  hash2["replacements"].each{ | each | all << each }
  checkArray.uniq.each{ | x | all.delete(x) }
  checkArray.uniq.each{ | site | 
    codon1 = self.sequenceMatrix1.codonSet( codingSites.index(site) / 3 )
    codon2 = self.sequenceMatrix2.codonSet( codingSites.index(site) / 3 )
    if codon1.size == codon2.size
      if codon1.all?{ | each | codon2.include?(each) }
        hash1["replacements"].occurrencesOf(site).times{ all << site }
      else
        hash1["replacements"].occurrencesOf(site).times{ all << site }
        hash2["replacements"].occurrencesOf(site).times{ all << site }
      end
    else
      siteSet1 = self.sequenceMatrix1.siteSet(site)
      siteSet2 = self.sequenceMatrix2.siteSet(site)
      bothSet = siteSet1 | siteSet2 
      (bothSet.size - 1).times{ all << site }
    end
  }
  all.each{ | eachSite | polys["replacements"] << eachSite }
  polys["complex"] = hash1["complex"] | hash2["complex"]
  return polys
end


def polarizeAllPolymorphisms
  polys = Hash.new
  ["species1","species2","notPol"].each{ | x | polys[x] = Array.new }
  both = self.sequenceMatrix1.segSitesLocationsAllSites | self.sequenceMatrix2.segSitesLocationsAllSites
  both.each{ | eachSite |
    os = self.outgroup.siteSetClean(eachSite)
    set1 = self.sequenceMatrix1.siteSetClean(eachSite)
    set2 = self.sequenceMatrix2.siteSetClean(eachSite)
    if set1.size == 1 or set2.size == 1				#one ingroup fixed?
      if set1.size == 1					#set1 fixed? if true
        if  set2.all?{ | each | os.include?(each) }	#outgroup has all set2 states? 
          polys["notPol"]  << eachSite
        else	
          if set1.any?{ | each | os.include?(each) }		#outgroup has any set1 states? if true...
            polys["species2"] << eachSite 			#species2 polymorphism
          else
            polys["notPol"]  << eachSite
          end
        end
      else							#set2 fixed
        if  set1.all?{ | each | os.include?(each) }
          poly["notPol"]  << eachSite
        else
          if set2.any?{ | each | os.include?(each) }
            polys["species1"] << eachSite 
          else
            polys["notPol"]  << eachSite
          end
        end
      end
    else								#both sets larger than 1 state
      set1.each{ | x | 					# delete overlapping states
        if set2.include?(x)
          set1.delete(x)
          set2.delete(x)
        end
      }
      if os.any?{ | each | set1.include?(each) }			#outgroup state in set1? if true
        if ! os.any?{ | each | set2.include?(each) }		#outgroup state in set2? if false
          polys["species2"] << eachSite 			#species2 polymorphism
        else 
          poly["notPol"]  << eachSite			#outgroup has state in set1 and set2 -> nonPol
        end
      else
        if  os.any?{ | each | set2.include?(each) }
          polys["species1"] << eachSite 
        end
      end
    end
  }
  return polys
end

def polarizeCodingPolymorphisms
  pol = Hash.new
  pol["notPol"] = Array.new
  polPolys = self.polarizeAllPolymorphisms					#use all polarized polys
  srHash = self.silReplPolymorphismWithin
  codingSites = self.sequenceMatrix1.codingSites
  ["species1","species2"].each{ | name |								#set up pol hash with an internal hash
    aHash = Hash.new
    ["silents","replacements","preferreds","unpreferreds","noChange"].each{ | x | aHash[x] = Array.new }
    pol[name] = aHash
  }
  polPolys["species1"].each{ | eachPoly |
    if srHash["replacements"].include?(eachPoly)									#silreplHash has site as repl
      oAA = self.outgroup.aminoSet(codingSites.index(eachPoly) / 3)
      inAA = self.sequenceMatrix2.aminoSet(codingSites.index(eachPoly) / 3)
      if oAA.any?{ | each | inAA.include?(each) } 
        pol["species1"]["replacements"] << eachPoly
      else 
        pol["notPol"] << eachPoly
      end
    else
      if srHash["silents"].include?(eachPoly)
        pol["species1"]["silents"] << eachPoly								
        #move to pref/unpref analysis
        puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachPoly) / 3)	
        puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachPoly) / 3)
        puSetO = self.outgroup.puSet(codingSites.index(eachPoly) / 3)
        if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
          if puSet1.size == puSet2.size
            pol["species1"]["noChange"] << eachPoly
          else
            puSet1 = puSet1.reject{ | x | puSet2.include?(x) } 
            if ! puSet1.empty?
              if puSet1.first == 1
                pol["species1"]["preferreds"] << eachPoly
              else
                pol["species1"]["unpreferreds"] << eachPoly
              end
            end
          end
        end
      end
    end
  }
  polPolys["species2"].each{ | eachPoly |
    if srHash["replacements"].include?(eachPoly)						#silreplHash has site as repl
      oAA = self.outgroup.aminoSet(codingSites.index(eachPoly) / 3)
      inAA = self.sequenceMatrix1.aminoSet(codingSites.index(eachPoly) / 3)
      if oAA.any?{ | each | inAA.include?(each) } 
        pol["species2"]["replacements"] << eachPoly
      else 
        pol["notPol"] << eachPoly
      end								
    else
      if srHash["silents"].include?(eachPoly)
        pol["species2"]["silents"] << eachPoly								
        #move to pref/unpref analysis
        puSet1 = self.sequenceMatrix1.puSet(codingSites.index(eachPoly) / 3)	
        puSet2 = self.sequenceMatrix2.puSet(codingSites.index(eachPoly) / 3)
        puSetO = self.outgroup.puSet(codingSites.index(eachPoly) / 3)
        if ! (puSet1.empty? or puSet2.empty? or puSetO.empty?)
          if puSet1.size == puSet2.size
            pol["species2"]["noChange"] << eachPoly
          else
            puSet2 = puSet2.reject{ | x | puSet1.include?(x) } 
            if ! puSet1.empty?
              if puSet2.first == 1
                pol["species2"]["preferreds"] << eachPoly
              else
                pol["species2"]["unpreferreds"] << eachPoly
              end
            end
          end
        end
      end
    end
  }
  pol["notPol"] = pol["notPol"] | polPolys["notPol"]
  return pol
end

def fst
  n = self.sequenceMatrix1.sampleSize + self.sequenceMatrix2.sampleSize
  hb = self.averagePairwiseDist 
  print hb,"\t"
  hw = self.sequenceMatrix1.averagePairwiseDifferencesArray  + self.sequenceMatrix2.averagePairwiseDifferencesArray  

  print hw.sampleMean,"\n"
  if hw == "NA" or hb == "NA"
    return("NA")
  else
    return(1.0 - (hw.sampleMean/hb))
  end
end


def averagePairwiseDist
  ds = Array.new      
  d = 0
  #go through all pairwise combos
  self.sequenceMatrix1.matrix.each_index{ | i |
    self.sequenceMatrix2.matrix.each_index{ | j |
      #print i,"\t",j,"\n"
      tempSeq1 = self.sequenceMatrix1.matrix[i]
      tempSeq2 = self.sequenceMatrix2.matrix[j]
      ds << Alignment.stringComp(tempSeq1,tempSeq2)
    }
  }
  return ds.sampleMean
end

def averagePairwiseJCDist
  ds = Array.new

  #go through all pairwise combos
  self.sequenceMatrix1.matrix.each_index{ | i |
    self.sequenceMatrix2.matrix.each_index{ | j |
      tempSeq1 = self.sequenceMatrix1.returnSequences([i])
      tempSeq2 = self.sequenceMatrix2.returnSequences([j])
      tempAlign = Alignment.new(tempSeq1,tempSeq2,nil)
      d = tempAlign.jcDiv
      ds << d

    }
  }
  return ds.sampleMean
end

def averagePairwiseNGDist     #returns an array (dn,ds), this is algorithm 1
  dns = Array.new
  dss = Array.new

  #go through all the pairwise combos
  self.sequenceMatrix1.matrix.each_index{ | i |
    self.sequenceMatrix2.matrix.each_index{ | j |
      dists = silReplDifsAllPaths(i,j)
      if ! dists.include?(nil)
        dns << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[0]))
        dss << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[1]))
      end
    }
  }
  return [dns.sampleMean, dss.sampleMean]
end
def averageWeightedPairwiseNGDist     #same as above but weights final average by sil/repl sites
  dns = Array.new
  dss = Array.new

  #make arrays to store average sites for comparisons
  ssComps = Array.new
  rsComps = Array.new
  #go through all the pairwise combos
  self.sequenceMatrix1.matrix.each_index{ | i |
    self.sequenceMatrix2.matrix.each_index{ | j |
      dists = silReplDifsAllPaths(i,j)
      if ! dists.include?(nil)
        sites = silReplSitesInComp(i,j)
        rsComps << sites[0]
        ssComps << sites[1]
        dns << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[0]))
        dss << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[1]))
      end
    }
  }
  #calc weights
  ssSum = ssComps.sum
  rsSum = rsComps.sum
  ssComps.each_index{ | i |
    dss[i] = dss[i] * (ssComps[i] / ssSum)
    dns[i] = dns[i] * (rsComps[i] / rsSum)
  }

  return [dns.sum, dss.sum]
end
def averagePairwiseNGDist2     #returns an array (dn,ds), this is algorithm 2
  dns = Array.new
  dss = Array.new

  #go through all the pairwise combos
  self.sequenceMatrix1.matrix.each_index{ | i |
    self.sequenceMatrix2.matrix.each_index{ | j |
      dists = silReplDifsNG(i,j)
      dns << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[0]))
      dss << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*dists[1]))
    }
  }
  return [dns.sampleMean, dss.sampleMean]
end
def averageWeightedPairwiseNGDist2     #same as above but weights final average by sil/repl sites
  dns = Array.new
  dss = Array.new
  #get silent and replacement sites for each species
  ss1 = self.sequenceMatrix1.silentSites
  rs1 = self.sequenceMatrix1.replacementSites
  ss2 = self.sequenceMatrix2.silentSites
  rs2 = self.sequenceMatrix2.replacementSites
  #make arrays to store average sites for comparisons
  ssComps = Array.new
  rsComps = Array.new
  #go through all the pairwise combos
  self.sequenceMatrix1.matrix.each_index{ | i |
    self.sequenceMatrix2.matrix.each_index{ | j |
      numerator = silReplDifsNG(i,j)
      denom = [rs1[i],rs2[j]].sampleMean
      rsComps << denom
      dns << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*(numerator[0].to_f / denom)))
      denom = [ss1[i],ss2[j]].sampleMean
      ssComps << denom
      dss << (-3.0/4.0) * Math.log(1.0 - ((4.0/3.0)*(numerator[1].to_f / denom)))

    }
  }
  #calc weights
  ssSum = ssComps.sum
  rsSum = rsComps.sum
  ssComps.each_index{ | i |
    dss[i] = dss[i] * (ssComps[i] / ssSum)
    dns[i] = dns[i] * (rsComps[i] / rsSum)
  }

  return [dns.sum, dss.sum]
end


def silReplDifsNG(specAIndex, specBIndex)  #this returns an array (nonSyn, Syn) as estimated by Nei & Gojobori alg. 2

  codons1 = self.sequenceMatrix1.codons[specAIndex]
  codons2 = self.sequenceMatrix2.codons[specBIndex]

  #initialize counts
  ms = [0,0,0]
  mm = [0,0,0]
  s = [0,0,0]

  #go through codons, count up number of diffs in each codon and their positions
  codons1.each_index{ | i |
    set = Array.new 
    set << codons1[i]
    set << codons2[i]
    #get rid of any ambiguous codons
    set.reject!{ | x | x.include?("N") or x.include?("X") }
    ssPos = self.sequenceMatrix1.codonSegSitesPositions(set)
    #one change?
    if ssPos.size == 1
      ms[ssPos.first] += 1
      #is it silent?
      if self.sequenceMatrix1.aminoCodonSet(set).size == 1
        s[ssPos.first] += 1
      end
    else
      #two or three changes in codon?
      if ssPos.size > 1
        ssPos.each{ | j| 
          mm[j] += 1
        }
      end
    end
  }
  #get prob of synonymous per position
  if s[0] == 0
    probSyn1 = 0.0
  else
    probSyn1 = s[0].to_f / ms[0]
  end

  if s[2] == 0
    probSyn3 = 0.0
  else
    probSyn3 = s[2].to_f / ms[2]
  end
  probNonSyn1 = 1.0 - probSyn1
  probNonSyn3 = 1.0 - probSyn3

  #estimate numbers
  sd = (probSyn1 * (ms[0] + mm[0])) + (probSyn3 * (ms[2] + mm[2]))
  nd = (probNonSyn1 * (ms[0] + mm[0])) +  (ms[1] + mm[1]) + (probNonSyn3 * (ms[2] + mm[2])) #assumes all second positions are nonsyn
  return [nd,sd]
end

def silReplDifsAllPaths(specAIndex, specBIndex)  #this returns an array (nonSyn, Syn) as estimated by Nei & Gojobori alg. 1

  codons1 = self.sequenceMatrix1.codons[specAIndex]
  codons2 = self.sequenceMatrix2.codons[specBIndex]

  #initialize arrays
  nsCount = 0
  sCount = 0
  nsSites = 0
  sSites = 0
  #go through codons, count up number of diffs in each codon and their positions
  codons1.each_index{ | i |
    set = Array.new 
    set << codons1[i]
    set << codons2[i]
    #get rid of any ambiguous codons
    set.reject!{ | x | x.include?("N") or x.include?("X") or self.geneticCode[x] == "*" }
    if set.size > 1
      #tally up sites along the way
      nsSites += self.sequenceMatrix1.silentSiteDict[set[0]][1]
      nsSites += self.sequenceMatrix1.silentSiteDict[set[1]][1]
      sSites += self.sequenceMatrix1.silentSiteDict[set[0]][0]
      sSites += self.sequenceMatrix1.silentSiteDict[set[1]][0]
    end
    #more than one codon state?
    if set.uniq.size > 1
      #get allPaths between codons
      paths = self.sequenceMatrix1.allPaths2Codons(set)
      tempRCount = 0
      tempSCount = 0
      #clean out paths which go through stop codons
      paths.reject!{ | x | x[0] == "*"}
      #go through paths, count stuff
      paths.each{ | aPath |
        #go through path steps, count em
        aPath.each{ | aStep |
          if aStep[0] == "R"
            tempRCount += 1
          else
            tempSCount += 1
          end
        }
      }
      nsCount += tempRCount.to_f / paths.size
      sCount += tempSCount.to_f/ paths.size
    end
  }

  if (nsSites == 0 and sSites == 0)
    return [nil,nil]
  else
    return [nsCount / (nsSites.to_f / 2), sCount / (sSites.to_f / 2)]    
  end
end
def silReplSitesInComp(specAIndex, specBIndex)  #this returns an array (nonSyn, Syn) of sites used in comparison

  codons1 = self.sequenceMatrix1.codons[specAIndex]
  codons2 = self.sequenceMatrix2.codons[specBIndex]

  #initialize arrays
  nsCount = 0
  sCount = 0
  nsSites = 0
  sSites = 0
  #go through codons, count up number of diffs in each codon and their positions
  codons1.each_index{ | i |
    set = Array.new 
    set << codons1[i]
    set << codons2[i]
    #get rid of any ambiguous codons
    set.reject!{ | x | x.include?("N") or x.include?("X") or self.sequenceMatrix1.geneticCode[x] == "*"}

    #more than one codon state?
    if set.size > 1
      #tally up sites along the way
      nsSites += self.sequenceMatrix1.silentSiteDict[set[0]][1]
      nsSites += self.sequenceMatrix1.silentSiteDict[set[1]][1]
      sSites += self.sequenceMatrix1.silentSiteDict[set[0]][0]
      sSites += self.sequenceMatrix1.silentSiteDict[set[1]][0]
    end
  }
  if (nsSites == 0 and sSites == 0)
    return [nil,nil]
  else
    return [nsSites.to_f / 2, sSites.to_f / 2]    
  end
end                        


#statistics


def mkTest				#returns an array in the order aaFix, aaPoly, silFix, silPoly
  fix = self.silReplFixations
  poly = self.silReplPolymorphismWithin
  return [fix["replacements"].size,poly["replacements"].size,fix["silents"].size,poly["silents"].size]
end

def mkTestGreedy				#returns an array in the order aaFix, aaPoly, silFix, silPoly
  fix = self.silReplFixationsGreedy
  poly = self.silReplPolymorphismWithin
  return [fix["replacements"].size,poly["replacements"].size,fix["silents"].size,poly["silents"].size]
end

def polMKTest			#returns hash
  fix = self.polarizeCodingFixations
  poly = self.polarizeCodingPolymorphisms
  dict = Hash.new
  ["species1","species2"].each{ | x |
    temp = Hash.new
    temp["aaFix"] = fix[x]["replacements"].size
    temp["silFix"] = fix[x]["silents"].size
    temp["prefFix"] = fix[x]["preferreds"].size
    temp["unFix"] = fix[x]["unpreferreds"].size
    temp["ncFix"] = fix[x]["noChange"].size
    temp["aaPoly"] = poly[x]["replacements"].size
    temp["silPoly"] = poly[x]["silents"].size
    temp["prefPoly"] = poly[x]["preferreds"].size
    temp["unPoly"] = poly[x]["unpreferreds"].size
    temp["ncPoly"] = poly[x]["noChange"].size
    dict[x] = temp
  }
  return dict
end

def polMKTestGreedy			#returns hash
  fix = self.polarizeCodingFixationsGreedy
  poly = self.polarizeCodingPolymorphisms
  dict = Hash.new
  ["species1","species2"].each{ | x |
    temp = Hash.new
    temp["aaFix"] = fix[x]["replacements"].size
    temp["silFix"] = fix[x]["silents"].size
    temp["prefFix"] = fix[x]["preferreds"].size
    temp["unFix"] = fix[x]["unpreferreds"].size
    temp["ncFix"] = fix[x]["noChange"].size
    temp["aaPoly"] = poly[x]["replacements"].size
    temp["silPoly"] = poly[x]["silents"].size
    temp["prefPoly"] = poly[x]["preferreds"].size
    temp["unPoly"] = poly[x]["unpreferreds"].size
    temp["ncPoly"] = poly[x]["noChange"].size
    dict[x] = temp
  }
  return dict
end
def polMKTestArray 		#returns 2 dimensional mkArray: aaFix, aaPoly, silFix, silPoly
  polMK = polMKTest
  out = Array.new
  out << [polMK["species1"]["aaFix"],polMK["species1"]["aaPoly"],polMK["species1"]["silFix"],polMK["species1"]["silPoly"]]
  out << [polMK["species2"]["aaFix"],polMK["species2"]["aaPoly"],polMK["species2"]["silFix"],polMK["species2"]["silPoly"]]
  return out
end
def polMKTestArrayGreedy 		#returns 2 dimensional mkArray: aaFix, aaPoly, silFix, silPoly
  polMK = polMKTestGreedy
  out = Array.new
  out << [polMK["species1"]["aaFix"],polMK["species1"]["aaPoly"],polMK["species1"]["silFix"],polMK["species1"]["silPoly"]]
  out << [polMK["species2"]["aaFix"],polMK["species2"]["aaPoly"],polMK["species2"]["silFix"],polMK["species2"]["silPoly"]]
  return out
end
def polPUTestArray 		#returns 2 dimensional mkArray: aaFix, aaPoly, silFix, silPoly
  polMK = polMKTest
  out = Array.new
  out << [polMK["species1"]["prefFix"],polMK["species1"]["prefPoly"],polMK["species1"]["unFix"],polMK["species1"]["unPoly"],polMK["species1"]["ncFix"],polMK["species1"]["ncPoly"]]
  out << [polMK["species2"]["prefFix"],polMK["species2"]["prefPoly"],polMK["species2"]["unFix"],polMK["species2"]["unPoly"],polMK["species2"]["ncFix"],polMK["species2"]["ncPoly"]]
  return out
end

#
#	
#misc
#
#
def padAlignment		#evens sequence length with N's from the back
  max = 0
  self.sequenceMatrix1.matrix.each{ | eachString |
    if eachString.length > max
      max = eachString.length
    end
  }
  self.sequenceMatrix2.matrix.each{ | eachString |
    if eachString.length > max
      max = eachString.length
    end
  }
  self.sequenceMatrix1.matrix.each{ | eachString |
    if eachString.length < max
      l = max - eachString.length
      l.times{ | i | eachString << "N" }
    end
  }
  self.sequenceMatrix2.matrix.each{ | eachString |
    if eachString.length < max
      l = max - eachString.length
      l.times{ | i | eachString << "N" }
    end
  }	
end
def slidingWindowFixedDiffs(windowSize, offset)    # [site, fd]
  i = 0
  l = self.length - windowSize
  oc = Array.new
  while i <= l
    oc << [i, self.fixedDiffsWindow(i, i + windowSize)]
    i += offset
  end
  oc.reject!{ | each | each[1] == "NA"}
  oc.each{ | each | print each[0],"\t",each[1].to_f / windowSize,"\n"}
end

def slidingWindowJCDiv(windowSize, offset)    # only use with two sequences!!! [site, fd]
  i = 0
  l = self.length - windowSize
  oc = Array.new
  while i <= l
    d =  self.fixedDiffsWindowPer(i, i + windowSize)
    if d != "NA"
      oc << [i, (-3.0/4) * Math.log(1 - ((4.0 *d )/3))]
    end
    i += offset
  end
  oc.each{ | each | print each[0],"\t",each[1],"\n"}
end

def jcDiv
  d =  self.fixedDiffsWindowPer(0, self.length - 1)
  if d != "NA"
    return  (-3.0/4) * Math.log(1 - ((4.0 *d )/3))
  else
    return "NA"
  end
end

#class methods
def Alignment.stringComp(string1,string2)
  count = 0
  0.upto(string1.length-1){|i|
    if string1[i,1] != "N" and string2[i,1] != "N" and string1[i,1] != "-" and string2[i,1] != "-"
      if string1[i,1] != string2[i,1]
        count += 1
      end
    end
  }
  return(count)
end

def Alignment.gTest(anArray)		#takes an array of length 4, [a,c,b,d]
  cellSum = (anArray[0] * Math.log(anArray[0])) + (anArray[1] * Math.log(anArray[1])) + (anArray[2] * Math.log(anArray[2])) + (anArray[3] * Math.log(anArray[3]))   
  r1 = anArray[0] + anArray[2]
  r2 = anArray[1] + anArray[3]
  c1 = anArray[0] + anArray[1]
  c2 = anArray[2] + anArray[3]
  n = anArray.sum.to_f
  nln = n * Math.log(n)
  rcSum = (r1 * Math.log(r1)) + (r2 * Math.log(r2)) + (c1 * Math.log(c1)) + (c2 * Math.log(c2))
  g = 2 * (cellSum - rcSum + nln)
  q = 1 + ((((n/r1) + (n/r2) - 1) * ((n/c1)+(n/c2) - 1)) / (6 * n))
  return g/q
end

def Alignment.fishersExactTest(anArray)		#takes an array of length 4, [a,c,b,d]
  tmpArray= Array.new
  anArray.each{ | x | tmpArray << x}
  n = anArray.sum.to_i
  q1 = (anArray[0] + anArray[2]).logFactorial + (anArray[1] + anArray[3]).logFactorial + (anArray[0] + anArray[1]).logFactorial + (anArray[2] + anArray[3]).logFactorial - n.logFactorial
  q2 = anArray[0].logFactorial + anArray[1].logFactorial + anArray[2].logFactorial + anArray[3].logFactorial
  pTail1 = 10 ** (q1 - q2)
  origP1 = 10 ** (q1 - q2)
  while (! anArray.include?(0))
    if (anArray[0] * anArray[3]) - (anArray[1] * anArray[2]) < 0
      anArray[0] -= 1
      anArray[3] -= 1
      anArray[1] += 1
      anArray[2] += 1
    else
      anArray[0] += 1
      anArray[3] += 1
      anArray[1] -= 1
      anArray[2] -= 1
    end
    q2 = anArray[0].logFactorial + anArray[1].logFactorial + anArray[2].logFactorial + anArray[3].logFactorial
    pTail1 += 10 ** (q1 - q2)

  end
  #now add up second tail
  if (tmpArray[0] * tmpArray[3]) - (tmpArray[1] * tmpArray[2]) < 0
    adj = [anArray[1],anArray[2]].min
    tmpArray[1] -= adj
    tmpArray[2] -= adj
    tmpArray[0] += adj
    tmpArray[3] += adj
  else
    adj = [anArray[0],anArray[3]].min
    tmpArray[1] += adj
    tmpArray[2] += adj
    tmpArray[0] -= adj
    tmpArray[3] -= adj
  end
  q1 = (tmpArray[0] + tmpArray[2]).logFactorial + (tmpArray[1] + tmpArray[3]).logFactorial + (tmpArray[0] + tmpArray[1]).logFactorial + (tmpArray[2] + tmpArray[3]).logFactorial - n.logFactorial
  q2 = tmpArray[0].logFactorial + tmpArray[1].logFactorial + tmpArray[2].logFactorial + tmpArray[3].logFactorial
  pTail2 = 10 ** (q1 - q2)
  origP2 = 10 ** (q1 - q2)
  while (origP2 < origP1)
    if (tmpArray[0] * tmpArray[3]) - (tmpArray[1] * tmpArray[2]) < 0
      tmpArray[0] += 1
      tmpArray[3] += 1
      tmpArray[1] -= 1
      tmpArray[2] -= 1
    else
      tmpArray[0] -= 1
      tmpArray[3] -= 1
      tmpArray[1] += 1
      tmpArray[2] += 1
    end
    q2 = tmpArray[0].logFactorial + tmpArray[1].logFactorial + tmpArray[2].logFactorial + tmpArray[3].logFactorial
    pTail2 += 10 ** (q1 - q2)
    origP2 = 10 ** (q1 - q2)
  end
  pTail2 -= origP2
  return pTail1 + pTail2
end
end	

