        integer nhistsmax, nbinsmax
        parameter (nhistsmax=100, nbinsmax=1000)
        character*(100) histlist(nhistsmax)
        integer nhists, unitlist(nhistsmax)
        
        common /c_hist/ nhists,unitlist,histlist