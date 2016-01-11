c############### Analysis-dummy.f ######################################
c last modified by MK, 07.01.2015
c adapted from dislepton
c open some histograms and fill them with data
c + several auxiliary subroutines

c############### subroutine init_hist ##################################
c book all histograms, will be filled later

      subroutine init_hist
        implicit none
      end

c############### end subroutine init_hist ##############################

c############### subroutine analysis ###################################
c extract all data required for the histograms
c calculate quantities that shall be plotted
c fill histograms

      subroutine analysis(dsig)
        implicit none
        double precision dsig ! total cross section
      end

c############### end subroutine analysis ###############################
