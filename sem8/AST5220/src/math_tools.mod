
  n*  Z   k820309    P          16.0         «EW                                                                                                           
       math_tools.f90 MATH_TOOLS                                                     
                                                              u #INVERT_MATRIX_DPC    #INVERT_MATRIX_DP    #INVERT_MATRIX_SP                                                           u #INVERT_MATRIX_WITH_MASK_DPC    #INVERT_MATRIX_WITH_MASK_DP                                                           u #SOLVE_LINEAR_SYSTEM_DP    #SOLVE_LINEAR_SYSTEM_EIGEN    #SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK 	                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                   
                  ÿÿÿÿÿÿï                                                          HUGE                                                  
                
          )       -DTû!	@        3.141592653589793238462643383279502884197                                                 
                
          *       mBP×Ò?        0.2820947917738781434740397257803862929220#         @      X                                                 #MATRIX              
D@                                                  "              & p                  & p                                          #         @      X                                                 #MATRIX              
D@                                                 
 %              & p                  & p                                          #         @      X                                                 #MATRIX              
D@                                                 	 (              & p                  & p                                          #         @      X                                                 #MATRIX              
D@                                                  +              & p                  & p                                          #         @      X                                                 #MATRIX              
D@                                                 
 /              & p                  & p                                          #         @      X                                                 #A    #X    #B              
                                                    
              & p                  & p                                                    D@                                                 
               & p                                                    
                                                    
              & p                                          #         @      X                                                 #MYID    #A    #X    #B    #CUTOFF     #EIGENVALS !             
  @                                                   
  @                                                 
              & p                  & p                                                    D@                                                 
 	              & p                                                    
                                                    
              & p                                                    
                                       
                F @                              !                   
 
              & p                                          #         @      X                             	                    #MYID "   #A #   #X $   #B %   #MASK &   #CUTOFF '   #EIGENVALS (             
  @                              "                     
                                 #                   
              & p                  & p                                                    D@                              $                   
               & p                                                    
                                 %                   
              & p                                                    
  @                              &                                 & p                                                    
                                 '     
                F @                              (                   
               & p                                          #         @                                  )                    #MYID *   #MATRIX +   #EIGENVALS ,   #EIGENVECTORS -             
                                 *                     
                                 +                   
              & p                  & p                                                    D@                              ,                   
               & p                                                    D @                              -                   
               & p                  & p                                          #         @                                   .                    #MATRIX /   #THRESHOLD 0             
D@                              /                   
               & p                  & p                                                    
                                 0     
      #         @                                   1                    #MATRIX 2             
D@                              2                   	 3              & p                  & p                                          #         @                                   3                    #A 4   #X 5   #B 6                                             4                    7              &                   &                                                     D                                5                    8              &                                                  0  @                              6                    9              &                                           #         @                                  7                    #A 8   #X 9   #B :                                             8                   
 =              &                   &                                                     D                                9                   
 >              &                                                  0  @                              :                   
 ?              &                                           #         @                                   ;                    #M <   #LEN =             D @                              <                   
 C              &                   &                                                                                     =            #         @                                   >                    #A ?   #L @             
                                ?                   
 F             &                   &                                                     D @                              @                   
 G              &                   &                                           #         @                                   A                    #A B   #L C             
 @                              B                   
 I             & p                  & p                                                    D @                              C                   
 J              & p                  & p                                          #         @                                   D                    #L E   #B F   #X G             
                                E                   
 L             &                   &                                                     
                                 F                   
 M             &                                                     D                                G                   
 N              &                                           #         @                                   H                    #NLMAX I   #M J   #THETA K   #PLM L             
                                 I                     
  @                              J                     
  @                              K     
               D                                L                    
 Q    p           & p         5  p        r I         5  p        r I   p         p                          #         @                                   M                    #FRAC N   #SIGMA O             D                                N     	                 
  @                              O     	      %         @                               P                    
       #X Q             
  @                              Q     
      #         @                                   R                    #SIGMA S   #FRACT T             D                                S     	                 
                                 T     	      %         @                                U                    
       #XX V             
                                 V     
             "      fn#fn    Â   @   J   HEALPIX_TYPES "            gen@INVERT_MATRIX ,            gen@INVERT_MATRIX_WITH_MASK (     ¤       gen@SOLVE_LINEAR_SYSTEM !   ª  p       DP+HEALPIX_TYPES "     p       I4B+HEALPIX_TYPES "     p       I8B+HEALPIX_TYPES "   ú  p       LGT+HEALPIX_TYPES !   j  p       SP+HEALPIX_TYPES "   Ú  p       DPC+HEALPIX_TYPES %   J  p       MAX_DP+HEALPIX_TYPES #   º  =       HUGE+HEALPIX_TYPES !   ÷         PI+HEALPIX_TYPES (            SQ4PI_INV+HEALPIX_TYPES "   *  T       INVERT_MATRIX_DPC )   ~  ¬   a   INVERT_MATRIX_DPC%MATRIX !   *  T       INVERT_MATRIX_DP (   ~  ¬   a   INVERT_MATRIX_DP%MATRIX !   *	  T       INVERT_MATRIX_SP (   ~	  ¬   a   INVERT_MATRIX_SP%MATRIX ,   *
  T       INVERT_MATRIX_WITH_MASK_DPC 3   ~
  ¬   a   INVERT_MATRIX_WITH_MASK_DPC%MATRIX +   *  T       INVERT_MATRIX_WITH_MASK_DP 2   ~  ¬   a   INVERT_MATRIX_WITH_MASK_DP%MATRIX '   *  ]       SOLVE_LINEAR_SYSTEM_DP )     ¬   a   SOLVE_LINEAR_SYSTEM_DP%A )   3     a   SOLVE_LINEAR_SYSTEM_DP%X )   Ã     a   SOLVE_LINEAR_SYSTEM_DP%B *   S         SOLVE_LINEAR_SYSTEM_EIGEN /   Õ  @   a   SOLVE_LINEAR_SYSTEM_EIGEN%MYID ,     ¬   a   SOLVE_LINEAR_SYSTEM_EIGEN%A ,   Á     a   SOLVE_LINEAR_SYSTEM_EIGEN%X ,   Q     a   SOLVE_LINEAR_SYSTEM_EIGEN%B 1   á  @   a   SOLVE_LINEAR_SYSTEM_EIGEN%CUTOFF 4   !     a   SOLVE_LINEAR_SYSTEM_EIGEN%EIGENVALS 4   ±         SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK 9   =  @   a   SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK%MYID 6   }  ¬   a   SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK%A 6   )     a   SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK%X 6   ¹     a   SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK%B 9   I     a   SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK%MASK ;   Ù  @   a   SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK%CUTOFF >        a   SOLVE_LINEAR_SYSTEM_EIGEN_WITH_MASK%EIGENVALS (   ©         GET_EIGEN_DECOMPOSITION -   (  @   a   GET_EIGEN_DECOMPOSITION%MYID /   h  ¬   a   GET_EIGEN_DECOMPOSITION%MATRIX 2        a   GET_EIGEN_DECOMPOSITION%EIGENVALS 5   ¤  ¬   a   GET_EIGEN_DECOMPOSITION%EIGENVECTORS '   P  c       INVERT_SINGULAR_MATRIX .   ³  ¬   a   INVERT_SINGULAR_MATRIX%MATRIX 1   _  @   a   INVERT_SINGULAR_MATRIX%THRESHOLD +     T       INVERT_MATRIX_WITH_MASK_SP 2   ó  ¬   a   INVERT_MATRIX_WITH_MASK_SP%MATRIX      ]       SOLVE_SYSTEM    ü  ¤   a   SOLVE_SYSTEM%A          a   SOLVE_SYSTEM%X    ,     a   SOLVE_SYSTEM%B "   ¸  ]       SOLVE_SYSTEM_REAL $     ¤   a   SOLVE_SYSTEM_REAL%A $   ¹     a   SOLVE_SYSTEM_REAL%X $   E     a   SOLVE_SYSTEM_REAL%B #   Ñ  X       INVERT_MATRIX_REAL %   )  ¤   a   INVERT_MATRIX_REAL%M '   Í  @   a   INVERT_MATRIX_REAL%LEN #      V       CHOLESKY_DECOMPOSE %   c   ¤   a   CHOLESKY_DECOMPOSE%A %   !  ¤   a   CHOLESKY_DECOMPOSE%L -   «!  V       CHOLESKY_DECOMPOSE_WITH_MASK /   "  ¬   a   CHOLESKY_DECOMPOSE_WITH_MASK%A /   ­"  ¬   a   CHOLESKY_DECOMPOSE_WITH_MASK%L    Y#  ]       CHOLESKY_SOLVE !   ¶#  ¤   a   CHOLESKY_SOLVE%L !   Z$     a   CHOLESKY_SOLVE%B !   æ$     a   CHOLESKY_SOLVE%X $   r%  n       COMP_NORMALISED_PLM *   à%  @   a   COMP_NORMALISED_PLM%NLMAX &    &  @   a   COMP_NORMALISED_PLM%M *   `&  @   a   COMP_NORMALISED_PLM%THETA (    &  ä   a   COMP_NORMALISED_PLM%PLM #   '  ]       CONVERT_SIGMA2FRAC (   á'  @   a   CONVERT_SIGMA2FRAC%FRAC )   !(  @   a   CONVERT_SIGMA2FRAC%SIGMA    a(  W       CORR_ERF    ¸(  @   a   CORR_ERF%X $   ø(  ^       CONVERT_FRACT2SIGMA *   V)  @   a   CONVERT_FRACT2SIGMA%SIGMA *   )  @   a   CONVERT_FRACT2SIGMA%FRACT    Ö)  X       GAMMLN    .*  @   a   GAMMLN%XX 