# Math add ons from A. D. Kern
#
#

include Math

class Fixnum
  def factorial
  	if self == 0
  		return 1
  	end
  	sum = self
  	a = self - 1
    while a > 0
     sum *= a
     a -= 1
    end
    return sum
  end
  
  def logFactorial		#this returns log factorial
  	sum = 0
  	1.upto(self){ | i |
  		sum += Math.log10(i)
  		}
  	return sum
  end
end

