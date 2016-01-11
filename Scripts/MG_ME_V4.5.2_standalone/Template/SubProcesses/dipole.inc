c *************************************************************                                                         
c  Variables used for the integrated dipoles                                                                            
c *************************************************************                                                         
      logical intdipole                                                                                                 
      common /intd/ intdipole                                                                                           
                                                                                                                        
                                                                                                                        
c *************************************************************                                                         
c Mass and flavour parameters                                                                                           
c *************************************************************                                                         
                                                                                                                        
      real*8 mu              !Renormalization scale                                                                     
      real*8 muf             !Factorization scale                                                                       
      common /c_scales/mu,muf                                                                                           
                                                                                                                        
c Number of flavours in the integrated dipoles for initial                                                              
c state splittings of gluons. Should be set equal to                                                                    
c the number of quark flavors in the initial state.                                                                     
c Note that for final state splittings, this factor is not used,                                                        
c and included gluon to quark splittings are included via                                                               
c the explicit real-emission subprocesses. For the                                                                      
c initial state this cannot be done because of the PDFs.                                                                
c (Running of alpha_s is independent of this factor)                                                                    
      integer Nf                                                                                                        
      parameter (Nf=4)                                                                                                  
                                                                                                                        
                                                                                                                        
c **************************************************************                                                        
c Alpha parameters for phase space restriction                                                                          
c **************************************************************                                                        
                                                                                                                        
      real*8 alpha_ff,alpha_fi,alpha_if,alpha_ii                                                                        
      parameter (alpha_ff=1d0)                                                                                          
      parameter (alpha_fi=1d0)                                                                                          
      parameter (alpha_if=1d0)                                                                                          
      parameter (alpha_ii=1d0)                                                                                          
                                                                                                                        
                                                                                                                        
c ***************************************************************                                                       
c Renormalisation scheme                                                                                                
c Two renormalisation schemes are implemented:                                                                          
c 1. The 't Hooft-Veltman scheme (scheme=HV)                                                                            
c 2. dimensional reduction      (scheme=DR)                                                                             
c ***************************************************************                                                       
      character*2 scheme                                                                                                
      parameter (scheme='HV')                                                                                           
                                                                                                                        
                                                                                                                        
                                                                                                                        
c ***************************************************************                                                       
c Variables only needed for Photons                                                                                     
c ***************************************************************                                                       
                                                                                                                        
c ***************************************************************                                                       
c Regularisation scheme                                                                                                 
c Two regularisation schemes are implemented:                                                                           
c 1. dimensional regularisation (rscheme= DREG)                                                                         
c 2. mass regularisation        (rscheme= MREG)                                                                         
c ***************************************************************                                                       
      character*4 rscheme                                                                                               
      parameter (rscheme='MREG')                                                                                        
                                                                                                                        
c ***************************************************************                                                       
c Photons in the initial state can come out of a hadron. In that                                                        
c case collinear singularities are cancelled by redefinition of                                                         
c the PDF.                                                                                                              
c ***************************************************************                                                       
      logical photonpdf                                                                                                 
      parameter (photonpdf=.FALSE.)                                                                                     
                                                                                                                        
                                                                                                                        
c ***************************************************************                                                       
c Masses of particles, used for regularisation in integrated                                                            
c subtraction terms when using mass regularisation. The masses                                                          
c of the n+1 process should be defined here, a photon must be declared                                                  
c as a massless particle. Units are GeV.                                                                                
c ***************************************************************                                                       
      real*8 massint( 4)                                                                                                
      data massint/ 1d0,1d0,1d0,1d0/                                                                                    
                                                                                                                        
                                                                                                                        
                                                                                                                        
                                                                                                                        
c ***************************************************************                                                       
c Photon mass                                                                                                           
c ***************************************************************                                                       
      real*8 mass_a                                                                                                     
      parameter (mass_a=1d0)                                                                                            
                                                                                                                        
                                                                                                                        
                                                                                                                        
c ***************************************************************                                                       
c Non-collinear safe observables.                                                                                       
c ***************************************************************                                                       
      logical ncs                                                                                                       
      parameter (ncs=.false.)                                                                                           
                                                                                                                        
                                                                                                                        
c ***************************************************************                                                       
c Photon fragmenation function                                                                                          
c ***************************************************************                                                       
      logical photonfrag                                                                                                
      parameter (photonfrag=.true.)                                                                                     
