# Array add ons from A.D. Kern

include Math

load "#{File.join(__dir__, 'mathADK.rb')}"

module ArrayADK

	def checkArray		#returns edited self or nil 
		temp = self.delete_if{ | x | x == "NA" }
		if temp.empty?
			return nil
		else
			return temp
		end
	end
	def inject(n)
		each do | value |
			n = yield(n, value)
		end
		n
	end
	def sum(initial = 0)
		inject(initial){ | n, value | n + value }
	end
	def max
		temp = self.checkArray
		if temp
			max = 0.0
			temp.each{ | x | 
					if x > max
						max = x
					end
					}
			return max
		else
			return "NA"
		end
	end
	def product(initial = 1)
		inject(initial){ | n, value | n * value }
	end
	def sampleMean
		temp = self.delete_if{ | x | x.class == String }
		if temp.empty?
			return "NA"
		else
			return (temp.sum).to_f/(temp.size).to_f
		end
	end
        def weightedMean(weightArray)
            numSum = denSum = 0.0
            0.upto(self.size - 1){ | i |
                numSum += self[i].to_f * weightArray[i].to_f
                denSum += weightArray[i].to_f
            }
            return(numSum / denSum)
        end
                
	def sampleVariance
		temp = self.delete_if{ | x | x.class == String }
		mean = temp.sampleMean
		count = 0
		temp.each{ | x |
				count += (x - mean) * (x - mean)
				}
		count/(temp.size - 1)
	end
	def stdev
		temp = self.delete_if{ | x | x.class == String }
		mean = temp.sampleMean
		count = 0
		temp.each{ | x |
				count += (x - mean) * (x - mean)
				}
		return Math.sqrt(count/temp.size)
	end
	
	def occurrencesOf(anObject)
		count = 0
		self.each{ | x | 
			if x == anObject
				count += 1
			end }
		return count
	end
	def pearsonCorrelation(secondArray)
		if self.size == secondArray.size
			naArray = Array.new
			if self.include?("NA") or secondArray.include?("NA")
				self.each_index{ | i | if self[i] == "NA" or secondArray[i] == "NA"
											naArray<<i
										end
									}
			end
			naArray.reverse.each{ | x | self.delete_at(x)
										secondArray.delete_at(x)
									}
			
			
			temp1 = self.checkArray
			temp2 = secondArray.checkArray
			
			mean1 = temp1.sampleMean
			mean2 = temp2.sampleMean
			
			stdev1 = temp1.stdev
			stdev2 = temp2.stdev
			sum = 0
			
			0.upto(temp1.size - 1){ | i |
				sum += ((temp1[i] - mean1) / stdev1) * ((temp2[i] - mean2) / stdev2)
				}
			return sum / temp1.size.to_f
			
			
			
		else
			print "Error: Arrays different lengths\n"
			exit()
		end
	end
	
	
	def distancesBetween(anObject)   #returns an array of index distances between anObject in current array
		dists = Array.new
		positions = Array.new
		self.each_index{ | i |
				if (self[i] == anObject)
					positions.push(i)
				end
				}
		positions.reverse!
		positions.each_index { | i |
					if i < (positions.size - 1)
						dists.push( positions[i] - positions[i + 1] + 1)
					end
					}
		return dists
	end
	
	def distancesBetweenTwo(anObject, anotherObject)   #returns an array of index distances between anObject in current array
		dists = Array.new
		positions = Array.new
		current = anObject
		self.each_index{ | i |
				if (self[i].to_s == current)
					positions.push(i)
					if current == anObject
						current = anotherObject
					else
						current = anObject
					end
				end
				}
		if positions.empty?
			print "distancesBetween(anObject): no match in to object!\n"
			exit
		end
		positions.reverse!
		positions.each_index { | i |
					if i < (positions.size - 1)
						dists.push( positions[i] - positions[i + 1])
					end
					}
		return dists
	end
	
	def chisquare				#takes an array with 4 entries
		row1 = self[0] + self[1]
		row2 = self[2] + self[3]
		column1 = self[0] + self[2]
		column2 = self[1] + self[3]
		tot = column1 + column2
		top = (((self[0] * self[3]) - (self[1] * self[2])) * ((self[0] * self[3]) - (self[1] * self[2]))) * tot
		bottom = row1 * row2 * column1 * column2
		return top.to_f/bottom
	end
	
	def randomArray
		rSum = self[0] + self[1]  #returns a random array of length 2
		cSum = self[0] + self[2]
		min = [rSum, cSum].min
		new = Array.new(2,0)
		new[0] = rand(min)
		new[1] = rSum - new[0] 
		return new
	end
	
	def randomTable   #takes an array with 4 entries, returns random array conditional on marginals
		new = Array.new
		col1 = self[0] + self[2]
		col2 = self[1] + self[3]
		randRow = self.randomArray
		new.push(randRow[0])
		new.push(randRow[1])
		new.push(col1 - randRow[0])
		new.push(col2 - randRow[1])
		return new
	end
	
	def monteCarloContingency		#take an array with 4 entries
		chiArray = Array.new
		10000.times { 
			chiArray.push(self.randomTable.chisquare)
			}
		chiArray.sort!
		
		testStat = self.chisquare
		chiArray.each_index{ | i |
					if chiArray[i] >= testStat
						return (1 - (i.to_f/10000.0))
					end
					}
	
	end
	
		
	def fishersExactTest  #takes an array with 4 entries order: A11 A12 A21 A22#
		
	if self.include?(0)
		return 1
	end
		
		
		pArray = Array.new
		while ! self.include?(-1) 

			row1 = self[0] + self[1]
			
			row2 = self[2] + self[3]
			column1 = self[0] + self[2]
			column2 = self[1] + self[3]
			tot = column1 + column2
			top = row1.factorial * row2.factorial * column2.factorial * column1.factorial
			bottom = tot.factorial * (self[0].factorial * self[1].factorial * self[2].factorial * self[3].factorial)
			pArray.push(top.to_f/bottom.to_f)
			
			
			
			if (self[0] * self[3]) - (self[1] * self[2]) > 0
				self[0] += 1
				self[3] += 1
				self[1] -= 1
				self[2] -= 1
			else
				self[0] -= 1
				self[3] -= 1
				self[1] += 1
				self[2] += 1
			end
			
		end				
		
		return pArray.sum * 2.0
	end
end
