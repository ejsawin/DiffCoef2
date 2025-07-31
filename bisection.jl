function bisection(func,x1,x2,precision,max_iters)
    
    #Set iter count and initialize values for beginning/end points 
    iters=0
    val1=func(x1)
    val2=func(x2)
    
    while iters < max_iters

        #Calculate midpoint
        midpoint=(x1+x2)/2
        valmid=func(midpoint)

        #Check if root 
        if abs(valmid) < precision || abs(x2-x1)/2 < precision
            return midpoint
        end

        #Shift midpoint depending on sign
        if val1*valmid < 0
            x2=midpoint
            val2=valmid
        else
            x1=midpoint
            val1=valmid
        end
    end
end