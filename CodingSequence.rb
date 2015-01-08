# CodingSequence Module - This gives methods to the SequenceMatrix class which describe coding sequence
##
##

module CodingSequence

  #coding sequences
  def codingSites
    cr = Array.new
    self.codingRegions.each{ | x | cr.push(x - 1) }
    oc = Array.new
    exonNumber = cr.size / 2
    exonNumber.times{ | i |
      b = cr.shift
      e = cr.shift
      b.upto(e){ | site |
      oc.push(site) }
    }
    return oc
  end

  def intronSites           # this assumes only those sites surrounded by coding sequence are introns
    is = Array.new
    cs = self.codingSites
    if cs.empty?
      0.upto(self.matrix[1].size - 1){ | i | is.push(i) }
    end
    (self.codingRegions.first - 1).upto(self.codingRegions.last - 1){ | i | is.push(i) }
    cs.each{ | x | is.delete(x)}
    return is
  end

  def setCodons
    trips = Array.new
    oc = Array.new
    pad = self.readingFrame - 1

    self.matrix.each{ | anAllele |
      temp = String.new("")
      self.codingSites.each{ | c |
        temp << anAllele[c,1].to_s
      }
    oc.push(temp)}
    oc.each{ | cds |
      pad.times{
      cds.insert(0,"-") }
    }
    codonNumber = (oc[0].size / 3.0 ).floor
    top = oc[0].size
    lastCompleteCodonBase = codonNumber * 3
    if (top - lastCompleteCodonBase) != 0
      oc.each{ | cds |
        (3 - (top - lastCompleteCodonBase)).times{ cds << ("-") }
      }
    end

    oc.each{ | cds |
      temp = cds.scan(/.../)
      trips.push(temp)
    }
    self.codons = trips
    return self
  end
  def codonSet(aCodonIndex)
    cSet = Array.new
    self.codons.each{ | each | cSet.push(each[aCodonIndex])}
    return cSet.uniq
  end

  def codonSetClean(aCodonIndex)
    cSet = codonSet(aCodonIndex)
    cSet.reject!{ | x | x.include?("N") or x.include?("-") }
    return cSet
  end

  def aminoSet(aCodonIndex)
    aSet = Array.new
    self.codons.each{ | each | aSet.push(self.geneticCode[each[aCodonIndex]]) }
    return aSet.uniq
  end

  def puSet(aCodonIndex)
    c = self.codonSet(aCodonIndex)
    aSet = Array.new
    edit = c.reject{ | each | each.include?("N") or each.include?("-") }
    edit.each{ | each | aSet.push(self.puDict[each]) }
    return aSet.uniq
  end

  def aminoCodonSet(aSet)
    edit = aSet.reject{ | each | each.include?("N") or each.include?("-") }
    anArray = Array.new
    edit.each{ | each |
      anArray << self.geneticCode[each.upcase]
    }
    return anArray.uniq
  end

  def codonSegSites(aCodonIndex)
    cSet = Array.new
    self.codons.each{ | anArray |
      cSet.push(anArray[aCodonIndex])
    }
    cSet.uniq!
    setHold = Array.new(3,Array.new)

    sum = 0
    0.upto(2){ | i |
      test = Array.new
      cSet.each{ | codon | test.push(codon[i,1])}
      test.uniq!
      if test.size > 1
        sum += 1
    end }

    return sum
  end

  def codonSegSitesSet(aCodonSet)                       #returns number of sites segregating in a codon set
    sum = 0
    0.upto(2){ | c |
      temp = Array.new
      aCodonSet.each{ | eachCodon |
        temp << eachCodon[c,1]
      }
      if temp.include?("N") or temp.include?("-")
        return 0
      else
        if temp.uniq.size > 1
          sum += 1
        end
      end
    }
    return sum
  end

  def codonSegSitesPositions(aCodonSet)             # returns an Array of positions
    pos = Array.new
    sum = 0
    0.upto(2){ | c |
      temp = Array.new
      aCodonSet.each{ | eachCodon |
        temp << eachCodon[c,1]
      }
      if temp.include?("N") or temp.include?("-")
        return Array.new
      else
        if temp.uniq.size > 1
          pos << c
        end
      end
    }
    return pos
  end

  def codonSegSitesPositionDict                         #returns array of segSites in both codon and nuc. indexes
    pos = Hash.new
    segSites = self.segSitesLocationsAllSites
    codingSites = self.codingSites
    segs = segSites.find_all{ | x | codingSites.include?(x) }
    if segs.empty?
      return nil
    else
      codonPos = Array.new
      segs.each{ | each | codonPos << (((codingSites.index(each)).to_f) / 3) }
      pos["codonSegSites"] = codonPos
      pos["locusSegSites"] = segs
      return pos
    end
  end

  def codonIndexToCodingSequenceIndex(anIndex)                          #returns first nucleotide index of codon
    codingSites = self.codingSites
    nucIndex = anIndex * 3
    return codingSites[nucIndex]
  end

  def shortestPath2Codons(aCodonSet)                    # takes a set of codons and returns an array indicating a path

    aminoSet = self.aminoCodonSet(aCodonSet)
    pos = self.codonSegSitesPositions(aCodonSet)
    segSites = pos.size
    paths = Array.new

    if segSites == 1
      path = Array.new
      if aminoSet.size > 1
        path << ["R",pos[0].to_i]
      else
        path << ["S",pos[0].to_i]
      end
      paths << path
    else
      if segSites == 2
        perms = [[pos[0],pos[1]],[pos[1],pos[0]]]    #permutations for swaps
        paths = Array.new
        tempString = String.new("")
        tempString << aCodonSet[0]
        perms.each{ | eachArray |

          path = Array.new
          tempString = String.new("")
          tempString << aCodonSet[0]
          eachArray.each{ | i |
            tempAmino1 = String.new("")
            tempAmino1 << self.geneticCode[tempString]
            tempString[i,1] = aCodonSet[1][i,1]
            tempAmino2 = String.new("")
            tempAmino2 << self.geneticCode[tempString]
            if tempAmino2 == "*"
              path << ["*",i]
            end
            if tempAmino1 == tempAmino2
              path << ["S", i]
            else
              path << ["R", i]
            end
          }
        paths << path }
      else  #3 segSites!!!
        perms = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,1,0],[2,0,1]]    #permutations for swaps
        paths = Array.new
        perms.each{ | eachArray |
          path = Array.new
          tempString = String.new("")
          tempString << aCodonSet[0]
          eachArray.each{ | i |
            tempAmino1 = String.new("")
            tempAmino1 << self.geneticCode[tempString]
            tempString[i,1] = aCodonSet[1][i,1]
            tempAmino2 = String.new("")
            tempAmino2 << self.geneticCode[tempString]
            if tempAmino2 == "*"
              path << ["*",i]
            end
            if tempAmino1 == tempAmino2
              path << ["S", i]
            else
              path << ["R", i]
            end
          }
        paths << path }
      end
    end
    #set up the path length variable
    short = 100
    shortPath = nil
    noStopsPaths = paths.delete_if{ | x | x.first.first == "*" }
    #go through paths, find shortest with respect to Replacements
    noStopsPaths.each{ | aPath |
      rLength = 0
      aPath.each{ | eachArray | if eachArray.include?("R")
        rLength += 1
        end
      }
      if short > rLength
        short = rLength
        shortPath = aPath
      else
        if shortPath.size > aPath.size      # more replacements but less changes overall
          shortPath = aPath
        end
      end
    }
    return shortPath

  end
  def lookupPath2Codons(aCodonSet)
    if aCodonSet.include?("TAA") or aCodonSet.include?("TAG") or aCodonSet.include?("TGA")
      return nil
    end
    string = self.codonDists[aCodonSet.first][aCodonSet.last]
    path = Array.new
    steps = string.length / 2
    0.upto(steps - 1){ | i | path << [string[(2*i),1],string[(2*i)+1,1]] }
    return path
  end

  def shortestPath3Codons(aCodonSet)
    segSites = self.codonSegSitesSet(aCodonSet)
    aminoSet = self.aminoCodonSet(aCodonSet)
    pos = self.codonSegSitesPositions(aCodonSet)
    paths = Array.new
    shortPath = nil

    if segSites == 1
      if aminoSet.size > 1
        path = Array.new
        r = aminoSet.size - 1
        s = aCodonSet.size - 1 - r
        r.times{ path << ["R",pos[0]] }
        s.times{ path << ["S",pos[0]] }
        paths << path
      else
        path = Array.new
        s = aCodonSet.size - 1
        s.times{ path << ["S",pos[0]] }
        paths << path
      end
    else
      if segSites == 2
        if aCodonSet.size == 3
          #set up all possible paths, use 2 codon paths

          codonPaths = Array.new
          codonPaths << [[aCodonSet[0],aCodonSet[1]],[aCodonSet[1],aCodonSet[2]]]
          codonPaths << [[aCodonSet[0],aCodonSet[2]],[aCodonSet[1],aCodonSet[2]]]
          codonPaths << [[aCodonSet[0],aCodonSet[2]],[aCodonSet[0],aCodonSet[1]]]
          codonPaths.each{ | eachPath |
            temp = Array.new
            eachPath.each{ | testArray |
              if self.lookupPath2Codons(testArray).nil?
                print testArray,"\n"
              end
              self.lookupPath2Codons(testArray).each{ | x | temp << x }
            }
            paths << temp
          }
        else
          if aCodonSet.size == 4
            #set up all possible paths, use 2 codon paths
            codonPaths = Array.new
            codonPaths << [[aCodonSet[0],aCodonSet[1]],[aCodonSet[1],aCodonSet[2]],[aCodonSet[2],aCodonSet[3]]]
            codonPaths << [[aCodonSet[0],aCodonSet[1]],[aCodonSet[1],aCodonSet[3]],[aCodonSet[3],aCodonSet[2]]]
            codonPaths << [[aCodonSet[0],aCodonSet[2]],[aCodonSet[2],aCodonSet[1]],[aCodonSet[1],aCodonSet[3]]]
            codonPaths << [[aCodonSet[0],aCodonSet[2]],[aCodonSet[2],aCodonSet[3]],[aCodonSet[3],aCodonSet[1]]]
            codonPaths << [[aCodonSet[0],aCodonSet[3]],[aCodonSet[3],aCodonSet[2]],[aCodonSet[2],aCodonSet[1]]]
            codonPaths << [[aCodonSet[0],aCodonSet[3]],[aCodonSet[3],aCodonSet[1]],[aCodonSet[1],aCodonSet[2]]]
            codonPaths << [[aCodonSet[1],aCodonSet[0]],[aCodonSet[0],aCodonSet[2]],[aCodonSet[2],aCodonSet[3]]]
            codonPaths << [[aCodonSet[1],aCodonSet[3]],[aCodonSet[3],aCodonSet[0]],[aCodonSet[0],aCodonSet[2]]]
            codonPaths << [[aCodonSet[1],aCodonSet[3]],[aCodonSet[1],aCodonSet[0]],[aCodonSet[0],aCodonSet[2]]]
            codonPaths << [[aCodonSet[1],aCodonSet[2]],[aCodonSet[2],aCodonSet[0]],[aCodonSet[0],aCodonSet[3]]]
            codonPaths << [[aCodonSet[1],aCodonSet[0]],[aCodonSet[0],aCodonSet[3]],[aCodonSet[3],aCodonSet[2]]]
            codonPaths << [[aCodonSet[3],aCodonSet[0]],[aCodonSet[0],aCodonSet[1]],[aCodonSet[1],aCodonSet[2]]]
            codonPaths << [[aCodonSet[0],aCodonSet[1]],[aCodonSet[0],aCodonSet[3]],[aCodonSet[0],aCodonSet[2]]]
            codonPaths << [[aCodonSet[1],aCodonSet[0]],[aCodonSet[1],aCodonSet[2]],[aCodonSet[1],aCodonSet[3]]]
            codonPaths << [[aCodonSet[3],aCodonSet[2]],[aCodonSet[3],aCodonSet[1]],[aCodonSet[3],aCodonSet[0]]]
            codonPaths << [[aCodonSet[2],aCodonSet[0]],[aCodonSet[2],aCodonSet[1]],[aCodonSet[2],aCodonSet[3]]]
            codonPaths.each{ | eachPath |
              temp = Array.new
              eachPath.each{ | testArray |
                self.lookupPath2Codons(testArray).each{ | x | temp << x }
              }
              paths << temp
            }
          else
            return ["C",pos[0]]
          end
        end
      else
        return ["C",pos[0]]
      end
    end

    #set up the path length variable
    short = 100
    shortPath = nil
    paths = paths.reject{ | x | x.include?("*") }

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
        if shortPath.size > aPath.size      # more replacements but less changes overall
          shortPath = aPath
        end
      end
    }
    return shortPath
  end

  def allPaths2Codons(aCodonSet)        # takes a set of codons and returns an array of paths
    #initialize arrays
    aminoSet = self.aminoCodonSet(aCodonSet)
    pos = self.codonSegSitesPositions(aCodonSet)
    segSites = pos.size
    paths = Array.new

    #1 segSite?
    if segSites == 1
      path = Array.new
      if aminoSet.size > 1  #is it a aa replacement?
        path << ["R",pos[0].to_i]
      else
        path << ["S",pos[0].to_i]
      end
      paths << path
    else
      #2 segsites?
      if segSites == 2
        perms = [[pos[0],pos[1]],[pos[1],pos[0]]]    #permutations for swaps
        tempString = String.new("")
        tempString << aCodonSet[0]
        perms.each{ | eachArray |
          path = Array.new
          tempString = String.new("")
          tempString << aCodonSet[0]
          eachArray.each{ | i |
            tempAmino1 = String.new("")
            tempAmino1 << self.geneticCode[tempString]
            tempString[i,1] = aCodonSet[1][i,1]
            tempAmino2 = String.new("")
            tempAmino2 << self.geneticCode[tempString]
            if tempAmino2 == "*"
              path << ["*",i]
            end
            if tempAmino1 == tempAmino2
              path << ["S", i]
            else
              path << ["R", i]
            end
          }
        paths << path }
      else  #3 segSites!!!
        perms = [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,1,0],[2,0,1]]    #permutations for swaps
        paths = Array.new
        perms.each{ | eachArray |
          path = Array.new
          tempString = String.new("")
          tempString << aCodonSet[0]
          eachArray.each{ | i |
            tempAmino1 = String.new("")
            tempAmino1 << self.geneticCode[tempString]
            tempString[i,1] = aCodonSet[1][i,1]
            tempAmino2 = String.new("")
            tempAmino2 << self.geneticCode[tempString]
            if tempAmino2 == "*"
              path << ["*",i]
            end
            if tempAmino1 == tempAmino2
              path << ["S", i]
            else
              path << ["R", i]
            end
          }
          paths << path
        }
      end
    end
    return paths
  end

  def silReplMutations      #returns a dictionary
    subs = Hash.new
    last = self.codons[0].length - 1
    ["replacements","silents","missingData","complex"].each{ | x | subs[x] = Array.new}
    codingSites = self.codingSites
    dict = self.codonSegSitesPositionDict
    if dict
      segLocus = dict["locusSegSites"]
      segCodons = dict["codonSegSites"]
      #reduce segCodons
      segCodons.collect!{| x | x.floor }.uniq!
      #go through reduced set and do stuff
      segCodons.each_with_index{ | site, index |
        if site.floor != last
          #  print site,"\n"
          set = self.codonSet(site)
          mis = 0

          #is there ambiguous data? if yes cleanup
          # p site
          set.each{ | each | if each.include?("N") or each.include?("-")
            mis += 1
            end
          }
          if mis > 0
            subs["missingData"] << segLocus[index]
            set = set.reject{ | x | x.include?("N") or x.include?("-") }
          end
    	  set.reject!{| x | x == "TAG"}
	  set.reject!{| x | x == "TAA"}
	  set.reject!{| x | x == "TGA"}
          #how many states?

          if set.uniq.size == 2
            shortPath = lookupPath2Codons(set.uniq)
            if shortPath
              shortPath.each{ | each |
                if each.include?("R")
                  subs["replacements"] << (site * 3) + each[1].to_i #segLocus[index]
                else
                  subs["silents"] <<  (site * 3) + each[1].to_i # segLocus[index]
                end
              }
            end
          else
            shortPath = self.shortestPath3Codons(set.uniq)
            if shortPath
              if shortPath.include?("C")
                subs["complex"] << (site * 3) #segLocus[index]
              else
                shortPath.each{ | each |
                  if each.include?("R")
                    subs["replacements"] << (site * 3) + each[1].to_i #segLocus[index]
                  else
                    subs["silents"] << (site * 3) + each[1].to_i #segLocus[index]
                  end
                }
              end
            end
          end
        end
      }
    end
    return subs
  end

  def silReplDifsAllPaths2Seqs(seqIndex1, seqIndex2)  #this returns an array (nonSyn/nsSites, Syn/sSites) as estimated by Nei & Gojobori alg. 1

    codons1 = self.codons[seqIndex1]
    codons2 = self.codons[seqIndex2]

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
      if set.size > 0
        #tally up sites along the way
        nsSites += self.silentSiteDict[set[0]][1]
        sSites += self.silentSiteDict[set[0]][0]
      end
      if set.size > 1
	nsSites += self.silentSiteDict[set[1]][1]
	sSites += self.silentSiteDict[set[1]][0]
      end
      #more than one codon state?
      if set.uniq.size > 1
        #get allPaths between codons
        paths = self.allPaths2Codons(set)
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


  #accessing info

  def setGeneticCode(aCodeTable)

    if aCodeTable == "standard"
      gc = { "GCT" => "A",
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
      self.geneticCode = gc
    end
    return self
  end

  def setPUDict
    pu = { "GCT" => 0,
           "GCC" => 1,
           "GCA" => 0,
           "GCG" => 0,
           "TGT" => 0,
           "TGC" => 1,
           "GAT" => 0,
           "GAC" => 1,
           "GAA" => 0,
           "GAG" => 1,
           "TTT" => 0,
           "TTC" => 1,
           "GGT" => 0,
           "GGC" => 1,
           "GGA" => 0,
           "GGG" => 0,
           "CAT" => 0,
           "CAC" => 1,
           "ATT" => 0,
           "ATC" => 1,
           "ATA" => 0,
           "AAA" => 0,
           "AAG" => 1,
           "TTG" => 0,
           "TTA" => 0,
           "CTT" => 0,
           "CTC" => 1,
           "CTA" => 0,
           "CTG" => 1,
           "ATG" => 0,
           "AAT" => 0,
           "AAC" => 1,
           "CCT" => 0,
           "CCC" => 1,
           "CCA" => 0,
           "CCG" => 0,
           "CAA" => 0,
           "CAG" => 1,
           "CGT" => 1,
           "CGC" => 1,
           "CGA" => 0,
           "CGG" => 0,
           "AGA" => 0,
           "AGG" => 0,
           "TCT" => 0,
           "TCC" => 1,
           "TCA" => 0,
           "TCG" => 1,
           "AGT" => 0,
           "AGC" => 0,
           "ACT" => 0,
           "ACC" => 1,
           "ACA" => 0,
           "ACG" => 0,
           "GTT" => 0,
           "GTC" => 1,
           "GTA" => 0,
           "GTG" => 1,
           "TGG" => 0,
           "TAT" => 0,
           "TAC" => 1,
           "TAA" => 0,
           "TAG" => 0,
           "TGA" => 0 }
    self.puDict = pu
    return self
  end

  def degenCodonSite(aCodon, aSite)         #sets the synonymous "degeneracy" of each position in each codon; Li 1993
    h ={
      "TAC" => [0,0,2],
      "GTC" => [0,0,4],
      "CTG" => [2,0,4],
      "CAT" => [0,0,2],
      "GCG" => [0,0,4],
      "ACC" => [0,0,4],
      "AGG" => [2,0,2],
      "CCA" => [0,0,4],
      "TTA" => [2,0,2],
      "AAA" => [0,0,2],
      "ATT" => [0,0,3],
      "GGA" => [0,0,4],
      "TGT" => [0,0,2],
      "TCG" => [0,0,4],
      "GAG" => [0,0,2],
      "GCT" => [0,0,4],
      "TGA" => [0,2,0],
      "AGT" => [0,0,2],
      "CGG" => [2,0,4],
      "CAA" => [0,0,2],
      "CCC" => [0,0,4],
      "AAC" => [0,0,2],
      "CTT" => [0,0,4],
      "ATA" => [0,0,3],
      "GGC" => [0,0,4],
      "TTC" => [0,0,2],
      "TAG" => [0,0,2],
      "GTG" => [0,0,4],
      "GAT" => [0,0,2],
      "ACG" => [0,0,4],
      "TCT" => [0,0,4],
      "AGA" => [2,0,2],
      "CGT" => [0,0,4],
      "CTA" => [2,0,4],
      "ATC" => [0,0,3],
      "CAC" => [0,0,2],
      "TGC" => [0,0,2],
      "GCA" => [0,0,4],
      "TAT" => [0,0,2],
      "GTT" => [0,0,4],
      "ACT" => [0,0,4],
      "AGC" => [0,0,2],
      "TCA" => [0,0,4],
      "CGA" => [2,0,4],
      "CCG" => [0,0,4],
      "CTC" => [0,0,4],
      "TTG" => [2,0,2],
      "AAG" => [0,0,2],
      "GGG" => [0,0,4],
      "GAA" => [0,0,2],
      "GCC" => [0,0,4],
      "TAA" => [0,2,2],
      "TGG" => [0,0,0],
      "GTA" => [0,0,4],
      "ACA" => [0,0,4],
      "TCC" => [0,0,4],
      "CGC" => [0,0,4],
      "CAG" => [0,0,2],
      "CCT" => [0,0,4],
      "AAT" => [0,0,2],
      "ATG" => [0,0,0],
      "GGT" => [0,0,4],
      "TTT" => [0,0,2],
    "GAC" => [0,0,2]}

    return h[aCodon][aSite]
  end


  def setSilentSiteDict      #again Nei and Gojorbori style creation, altered a little to accomodate stops; returns an array [sSites,nSites]
    dict = Hash.new
    codons = self.geneticCode.keys.sort
    codons.delete("TAA")
    codons.delete("TGA")
    codons.delete("TAG")
    dict["TAA"] = [nil,nil]
    dict["TAG"] = [nil,nil]
    dict["TGA"] = [nil,nil]

    #go through codons
    codons.each{ | aCodon |
      s = 0
      n = 9
      origAA = self.geneticCode[aCodon]
      0.upto(2){ | i |
        temp = String.new
        0.upto(2){  | j | temp[j,1] = aCodon[j,1]}
        #set up possible states, delete actual
        states = ["A","C","T","G"]
        states.delete(temp[i,1])
        states.each{ | aState |
          temp[i,1] = aState
          tempAA = self.geneticCode[temp]
          #if the test aa is the same add a silent site
          if tempAA == origAA
            s += 1
            n -= 1
          else
            if tempAA == '*'
              n -= 1
            end
          end
        }

      }
      dict[aCodon] = [s.to_f / 3, n.to_f / 3]
    }
    self.silentSiteDict =  dict
  end


  def setCodonDists                     #reads in and stores paths as Hash, this creates a dependancy on the file codonMatrix.txt
    keys = self.geneticCode.keys.sort
    keys.delete("TAA")
    keys.delete("TGA")
    keys.delete("TAG")
    codonDictionary = Hash.new
    keys.each{ | key | codonDictionary[key] = Hash.new }
    # aFile = File.new("/Users/adk/rubyStuff/codonMatrix.txt")  #/Network/Servers/i-dpgp.ucdavis.edu/Users/adk/rubyStuff/codonMatrix.txt")
    aFile = File.readlines(File.join(__dir__, 'codonMatrix.txt'))
    # aFile.each_line{ | line |
    aFile.each{ | line |
      array = line.chomp.split("\t")
      if  array.first != ""
        from = array.first
        array.each_with_index{ | item, index |
          if index > 0
            if item != "nil"
              codonDictionary[from][keys[index - 1]] = item
              codonDictionary[keys[index - 1]][from] = item
            end
          end
        }
      end
    }
    self.codonDists = codonDictionary
  end

  def degeneracyMatrix          #returns an array of arrays corresponding to the 'degeneracy' of each sequence. [nondegen., 2fold, 4fold, codonCount]
    a = Array.new
    self.codons.each{ | eachAlleleArray |
      nd = 0
      f2 = 0
      f4 = 0
      cc = 0
      eachAlleleArray.each{ | codon |
        if ! (codon.include?("N") or codon.include?("-"))
          cc += 1
          0.upto(2){ | i |
            x = self.degenCodonSite(codon,i)
            case x
            when 0
              nd += 1
            when 2
              f2 += 1
            when 3
              f2 += 1
            else
              f4 += 1
            end
          }
        end
      }
      a << [nd,f2,f4,cc]
    }
    return a
  end

  #silent site stuff, like sil poly
  #

  def gcThree    #average gc 3 values for seqs
    freqs = Array.new
    self.codons.each{ | eachRow |
      gcCount = 0
      codonCount = 0
      eachRow.each{ |eachCodon |
        if ! (eachCodon.include?("N") or eachCodon.include?("-"))
          codonCount += 1
          if eachCodon[2,1] =~ /[GC]/
            gcCount += 1
          end
        end
      }

      if codonCount != 0
        freqs << gcCount.to_f / codonCount
      end
    }


    return freqs.sampleMean
  end

  def silentSites           #returns array  of silent sites
    silSites = Array.new
    self.codons.each{ | eachRow |
      temp = Array.new
      eachRow.each{ | eachCodon |
        if ! ( eachCodon.include?("N") or eachCodon.include?("-") )
          temp << self.silentSiteDict[eachCodon].first
        end
      }
      silSites << temp.reject{ | x | x == nil }.sum
    }
    return silSites
  end

  def replacementSites
    replSites = Array.new
    self.codons.each{ | eachRow |
      temp = Array.new
      eachRow.each{ | eachCodon |
        if ! ( eachCodon.include?("N") or eachCodon.include?("-") )
          temp << self.silentSiteDict[eachCodon].last
        end
      }
      replSites << temp.reject{ | x | x == nil}.sum
    }
    return replSites
  end

  def averageNumberSilentSites      #returns a float
    return self.silentSites.sampleMean
  end

  def averageNumberReplacementSites
    return self.replacementSites.sampleMean
  end

  def silentPi      #cleans ambiguites, uses parsimony counts and uses n&g silent site counts
    sites = self.silReplMutations["silents"]
    l = self.averageNumberSilentSites
    oc = Array.new
    sites.each{ | c |
      siteArray = self.siteArrayClean(c)
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

  def silReplPi2   #this routine uses the NG unweighted estimates of silent changes
    comps = 0
    silDifs = 0.0
    replDifs = 0.0
    #do all pairwise comparisons
      0.upto(self.sampleSize - 2){ | i |
        (i+1).upto(self.sampleSize - 1){ | j |
          difs = self.silReplDifsAllPaths2Seqs(i, j)
          if ! difs.include?(nil)
            replDifs += difs[0]
            silDifs += difs[1]
            comps += 1
          end
        }
      }
      silDifs = silDifs / comps
      if silDifs.nan?
        silDifs = 0
      end
      replDifs = replDifs / comps
      if replDifs.nan?
        replDifs = 0
      end
      return [replDifs, silDifs]
    end


    def replacementPi                   #cleans ambiguites and uses n&g silent site counts
      sites = self.silReplMutations["replacements"]
      l = self.averageNumberReplacementSites
      oc = Array.new
      sites.each{ | c |
        siteArray = self.siteArrayClean(c)
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

    def replacementTheta
      sites = self.silReplMutations["replacements"]
      l = self.averageNumberReplacementSites
      s = sites.size / l
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

    def silentTheta
      sites = self.silReplMutations["silents"]
      l = self.averageNumberSilentSites
      s = sites.size / l
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

    def silentTajD

      s = self.silReplMutations["silents"].size
      if s < 3
        return -666
      else
        n = self.sampleSize.to_f
        if n < 4
          return -666
        else
          k = self.silReplPi2[1] * self.averageNumberSilentSites
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
          d = k - (s.to_f / a1)
          stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
          return d/stdev
        end
      end
    end
    def replTajD

      s = self.silReplMutations["replacements"].size
      if s < 3
        return -666
      else
        n = self.sampleSize.to_f
        if n < 4
          return -666
        else
          k = self.silReplPi2[0] * self.averageNumberReplacementSites
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
          d = k - (s.to_f / a1)
          stdev = Math.sqrt((e1 * s) + ((e2 * s) * (s - 1)))
          return d/stdev
        end
      end
    end






  end
