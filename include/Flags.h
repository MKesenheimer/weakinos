c############### Flags.h ###############################################
c last modified by MK, 01.12.2015

c store here the flags which are necessary for this process
        
        ! whether to use fake virtuals or not
        logical flg_fakevirtuals

        ! flg_btilde: if flg_split!=0: 
        ! distinguish real contributions for btilde and remnants...
        logical flg_btilde
        
        ! flg_rad: when subtracting: for radiation-generation and 
        ! modified reals: R might turn negative 
        ! -> cut around the resonant gluino-mass 
        logical flg_rad
        
        common/flags/ flg_fakevirtuals, flg_btilde, flg_rad
        save /flags/

c############### end Flags.h ###########################################