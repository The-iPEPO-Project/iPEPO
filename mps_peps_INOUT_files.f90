MODULE mps_peps_INOUT_files

 USE utility
 USE definitions_mps_mpo
 USE definitions_peps_pepo

 IMPLICIT NONE

 !!! MPS ROUTINES !!!
 interface mps_print_function
    module procedure mps_print_function_GamLam
    module procedure mps_print_function_PLAIN
 end interface mps_print_function

 interface output_mps
    module procedure output_mps_GamLam
    module procedure output_mps_PLAIN
 end interface output_mps

 interface reload_mps
    module procedure reload_mps_GamLam_OLDFORMATTING
    module procedure reload_mps_GamLam
    module procedure reload_mps_PLAIN
 end interface reload_mps

 !!! MPO ROUTINES !!!
 interface mpo_print_function
    module procedure mpo_print_function_GamLam
    module procedure mpo_print_function_PLAIN
 end interface mpo_print_function

 interface output_mpo
    module procedure output_mpo_GamLam
    module procedure output_mpo_PLAIN
 end interface output_mpo

 interface reload_mpo
    module procedure reload_mpo_GamLam
    module procedure reload_mpo_PLAIN
 end interface reload_mpo

 !!! PEPS ROUTINES !!!
 interface peps_print_function
    module procedure peps_print_function_GamLam
    module procedure peps_print_function_PLAIN
 end interface peps_print_function

 interface output_peps
    module procedure output_peps_GamLam
    module procedure output_peps_PLAIN
 end interface output_peps

 interface reload_peps
    module procedure reload_peps_GamLam
    module procedure reload_peps_PLAIN
 end interface reload_peps

 private
 public mps_print_function,  reload_mps
 public mpo_print_function,  reload_mpo
 public peps_print_function, reload_peps

CONTAINS

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPS BLOCK (PLAIN, no lambdas) OUTPUT/RELOAD ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! This routine sets up the file for printing MPS, 
 !!! then prints it using output_mps routine
 SUBROUTINE mps_print_function_PLAIN(mps_block, mps_name)

  TYPE(block_mps),  INTENT(IN) :: mps_block(:)
  CHARACTER(LEN=*), INTENT(IN) :: mps_name

  INTEGER   :: N_sites
  CHARACTER :: dump_mps_filename*256, char_N*16

  !! Size of imps
  N_sites = SIZE(mps_block)

  !! Determine file name for imps dump
  WRITE(char_N,'(I6)') N_sites
  dump_mps_filename=TRIM(ADJUSTL(mps_name))//"_N="//TRIM(ADJUSTL(char_N))//".dat"
  
  !! Output MPS
  CALL output_mps(dump_mps_filename, mps_block)

  !! Write to an index file where we dumped the state.
  WRITE(DUMP_INDEX,*) "Output to:", dump_mps_filename
  WRITE(DUMP_INDEX,*)

 END SUBROUTINE mps_print_function_PLAIN



 !!! Output MPS to file !!!
 SUBROUTINE output_mps_PLAIN(filename, mps_block)

  CHARACTER,       INTENT(IN) :: filename*256
  TYPE(block_mps), INTENT(IN) :: mps_block(:)

  !! output file handle
  INTEGER, PARAMETER :: OUTPUT=90

  !! dims & indices
  INTEGER :: N_sites, site

  !! size of imps
  N_sites = SIZE(mps_block)

  OPEN(OUTPUT, FILE=filename,form="unformatted",status="replace")

  !! Num of sites
  WRITE(OUTPUT) N_sites

  !! Size of each matrix
  DO site=1,N_sites
      !! Dimensions:
      WRITE(OUTPUT) mps_block(site)%SNdim
      WRITE(OUTPUT) mps_block(site)%Wdim, mps_block(site)%Edim
  end DO

  !! Output each site of mps_block
  DO site=1,N_sites
      WRITE(OUTPUT) mps_block(site)%m
  end DO

  CLOSE(OUTPUT)
   
 END SUBROUTINE output_mps_PLAIN



 !!! Reload MPS from file !!!
 SUBROUTINE reload_mps_PLAIN(filename, mpo_block, mps_block, compactMpoFlag)

  CHARACTER,                    INTENT(IN)           :: filename*256
  TYPE(block_mpo),              INTENT(IN)           :: mpo_block(:)
  TYPE(block_mps), ALLOCATABLE, INTENT(INOUT)        :: mps_block(:)
  CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: compactMpoFlag

  INTEGER :: N_sites, site
  INTEGER :: INPUT=90
  INTEGER :: temp, SNdim, Wdim, Edim

  OPEN(INPUT,FILE=filename,form="unformatted",status="old")

  !! Number of t-steps = length of MPS
  READ(INPUT) N_sites

  !! Allocate empty mps block
  CALL allocate_empty_mps_block(mps_block, N_sites)

  !! Sizes of matrices - allocate each matrix
  DO site=1,N_sites

      !! Read in SNdim
      READ(INPUT) temp

      !! Get the expected SNdim from input MPO
      IF(PRESENT(compactMpoFlag)) THEN 
         SNdim = mpo_block(getSiteMpo(site, N_sites)) % Ndim
      ELSE
         SNdim = mpo_block(site) % Ndim
      end IF

      !! Check South-North dim to ensure consistency
      IF(temp .NE. SNdim)  THEN
         WRITE(*,*) "Attempt to read in mps_block from file with wrong South-North dim", filename
         WRITE(*,*) "temp & SNdim: ", temp, SNdim
         WRITE(*,*) "at site: ", site
         STOP
      end IF

      !! Allocate a new site of mps_block (also sets the dims of each site)
      READ(INPUT) Wdim, Edim 
      CALL allocate_mps_site(mps_block(site), SNdim, Wdim, Edim)

  end DO

  !! Read mps block from file
  DO site=1,N_sites
     READ(INPUT) mps_block(site)%m
  end DO
    
  CLOSE(INPUT)

  WRITE(*,*) "Succesfully reloaded nolambda mps_block from file: ", filename

 END SUBROUTINE reload_mps_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! iMPS (Gamma-Lambda form) OUTPUT/RELOAD ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! This routine sets up the file for printing iMPS, 
 !!! then prints it using output_mps routine
 SUBROUTINE mps_print_function_GamLam(iGamma, iLambda, imps_name)

  TYPE(block_mps),    INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)
  CHARACTER(LEN=*),   INTENT(IN) :: imps_name

  INTEGER :: N_sites
  CHARACTER :: dump_imps_filename*256, char_sites*16

  !! Size of imps
  N_sites = SIZE(iGamma)

  !! Determine file name for imps dump
  WRITE(char_sites,'(I6)') N_sites
  dump_imps_filename=TRIM(ADJUSTL(imps_name))//".dat"
  
  !! Output iMPS
  CALL output_mps(dump_imps_filename, iGamma, iLambda)

  !! Write to an index file where we dumped the state.
  WRITE(*,*) "Output to:", dump_imps_filename
  WRITE(*,*)

 END SUBROUTINE mps_print_function_GamLam




 !!! Output iMPS to file !!!
 SUBROUTINE output_mps_GamLam(filename, iGamma, iLambda)

  CHARACTER,          INTENT(IN) :: filename*256
  TYPE(block_mps),    INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)

  !! Output file handle
  INTEGER, PARAMETER :: OUTPUT=90

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Size of imps
  N_sites = SIZE(iGamma)

  OPEN(OUTPUT, FILE=filename, form="unformatted", status="replace")

  !! Num of sites
  WRITE(OUTPUT) N_sites

  !! Size of each matrix
  DO site=1,N_sites
     !! Dimensions:
     WRITE(OUTPUT) iGamma(site)%SNdim
     WRITE(OUTPUT) iGamma(site)%Wdim, iGamma(site)%Edim
  end DO

  !! Output each iGamma
  DO site=1,N_sites
     WRITE(OUTPUT) iGamma(site)%m
  end DO

  !! Output each iLambda
  DO site=1,N_sites
     WRITE(OUTPUT) iLambda(site)%m
  end DO

  CLOSE(OUTPUT)
   
 END SUBROUTINE output_mps_GamLam




 !!! Reload iMPS from file !!!
 SUBROUTINE reload_mps_GamLam(filename, mpo, iGamma, iLambda)

  CHARACTER,                       INTENT(IN)    :: filename*256
  TYPE(block_mpo),                 INTENT(IN)    :: mpo(:)
  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)

  INTEGER :: N_sites, site
  INTEGER :: INPUT=90
  INTEGER :: temp, SNdim, Wdim, Edim

  OPEN(INPUT, FILE=filename, form="unformatted", status="old")

  !! Number of t-steps = length of imps 
  READ(INPUT) N_sites

  !! Allocate empty mps block for imps
  CALL allocate_empty_mps_block(iGamma, N_sites)
  CALL allocate_empty_lambda_block(iLambda, N_sites)

  !! Sizes of matrices - allocate each matrix
  DO site=1,N_sites

      !! Read in SNdim
      READ(INPUT) temp

      !! Get the expected SNdim from input MPO
      SNdim = mpo(site)%Ndim

      !! Check South-North dim to ensure consistency
      IF(temp .NE. SNdim)  THEN
         WRITE(*,*) "Attempt to read in iMPS from file with wrong South-North dim", filename
         WRITE(*,*) "temp & SNdim: ", temp, SNdim
         WRITE(*,*) "at site: ", site
         STOP
      end IF

      !! Allocate a new site of iGamma, iLambda (also sets the dims of each site)
      READ(INPUT) Wdim, Edim 
      CALL allocate_mps_site(iGamma(site), SNdim, Wdim, Edim)
      CALL allocate_lambda_site(iLambda(site), Edim)
  end DO

  !! Read iGamma from file
  DO site=1,N_sites
     READ(INPUT) iGamma(site)%m
  end DO

  !! Read iLambda from file
  DO site=1,N_sites
     READ(INPUT) iLambda(site)%m
  end DO
    
  CLOSE(INPUT)

  WRITE(*,*) "Succesfully reloaded iMPS from file: ", filename

 END SUBROUTINE reload_mps_GamLam





 !!! Reload iMPS from files with old formatting, produced by previous TEBD code !!!
 SUBROUTINE reload_mps_GamLam_OLDFORMATTING(filename, mpo, iGamma, iLambda, OLDFORMATTING)

  CHARACTER,                       INTENT(IN)    :: filename*256
  TYPE(block_mpo),                 INTENT(IN)    :: mpo(:)
  TYPE(block_mps),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  CHARACTER(LEN=*),                INTENT(IN)    :: OLDFORMATTING !! Pass an extra flag
                                                                  !! to call the old formatting reload function 
                                                                  !! instead of the new one
  INTEGER :: N_sites, N_bonds, site
  INTEGER :: INPUT=90
  INTEGER :: temp, SNdim, Wdim, Edim

  OPEN(INPUT, FILE=filename, form="unformatted", status="old")

  !! Number of t-steps = length of imps 
  READ(INPUT) N_sites, N_bonds

  !! Allocate empty mps block for imps
  CALL allocate_empty_mps_block(iGamma, N_sites)
  CALL allocate_empty_lambda_block(iLambda, N_sites)

  !! Read in SNdim
  READ(INPUT) temp

  !! Get the expected SNdim from input MPO
  SNdim = mpo(1)%Ndim

  !! Check South-North dim to ensure consistency
  IF(temp .NE. SNdim)  THEN
     WRITE(*,*) "Attempt to read in iMPS from file with wrong South-North dim", filename
     WRITE(*,*) "temp & SNdim: ", temp, SNdim
     STOP
  end IF

  !! Sizes of matrices - allocate each matrix
  DO site=1,N_sites

      !! Allocate a new site of iGamma, iLambda (also sets the dims of each site)
      READ(INPUT) Wdim, Edim 
      CALL allocate_mps_site(iGamma(site), SNdim, Wdim, Edim)
      CALL allocate_lambda_site(iLambda(site), Edim)
  end DO

  !! Read iGamma from file
  DO site=1,N_sites
     READ(INPUT) iGamma(site)%m
  end DO

  !! Read iLambda from file
  DO site=1,N_sites
     READ(INPUT) iLambda(site)%m
  end DO
    
  CLOSE(INPUT)

  WRITE(*,*) "Succesfully reloaded iMPS from file: ", filename

 END SUBROUTINE reload_mps_GamLam_OLDFORMATTING

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MPO BLOCK (PLAIN, no lambdas) OUTPUT/RELOAD ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! This routine sets up the file for printing MPO, 
 !!! then prints it using output_mpo routine
 SUBROUTINE mpo_print_function_PLAIN(mpo_block, mpo_name)

  TYPE(block_mpo),  INTENT(IN) :: mpo_block(:)
  CHARACTER(LEN=*), INTENT(IN) :: mpo_name

  INTEGER   :: N_sites
  CHARACTER :: dump_mpo_filename*256, char_N*16

  !! Size of impo
  N_sites = SIZE(mpo_block)

  !! Determine file name for impo dump
  WRITE(char_N,'(I6)') N_sites
  dump_mpo_filename=TRIM(ADJUSTL(mpo_name))//"_N="//TRIM(ADJUSTL(char_N))//".dat"
  
  !! Output MPO
  CALL output_mpo(dump_mpo_filename, mpo_block)

  !! Write to an index file where we dumped the state.
  WRITE(*,*) "Output to:", dump_mpo_filename
  WRITE(*,*)

 END SUBROUTINE mpo_print_function_PLAIN






 !!! Output MPO to file !!!
 SUBROUTINE output_mpo_PLAIN(filename, mpo_block)

  CHARACTER,       INTENT(IN) :: filename*256
  TYPE(block_mpo), INTENT(IN) :: mpo_block(:)

  !! output file handle
  INTEGER, PARAMETER :: OUTPUT=90

  !! dims & indices
  INTEGER :: N_sites, site

  !! size of impo
  N_sites = SIZE(mpo_block)

  OPEN(OUTPUT, FILE=filename, form="unformatted", status="replace")

  !! Num of sites
  WRITE(OUTPUT) N_sites

  !! Size of each matrix
  DO site=1,N_sites
      !! Dimensions:
      WRITE(OUTPUT) mpo_block(site)%Sdim, mpo_block(site)%Ndim
      WRITE(OUTPUT) mpo_block(site)%Wdim, mpo_block(site)%Edim
  end DO

  !! Output each site of mpo_block
  DO site=1,N_sites
      WRITE(OUTPUT) mpo_block(site)%m
  end DO

  CLOSE(OUTPUT)
   
 END SUBROUTINE output_mpo_PLAIN





 !!! Reload MPO from file !!!
 SUBROUTINE reload_mpo_PLAIN(filename, mpo_block, Sdims, Ndims, compactMpoFlag)

  CHARACTER,                    INTENT(IN)           :: filename*256
  TYPE(block_mpo), ALLOCATABLE, INTENT(INOUT)        :: mpo_block(:)
  INTEGER,                      INTENT(IN)           :: Sdims(:), Ndims(:)
  CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: compactMpoFlag

  INTEGER :: N_sites, site
  INTEGER :: INPUT=90
  INTEGER :: tempS, tempN
  INTEGER :: Sdim, Ndim, Wdim, Edim

  OPEN(INPUT,FILE=filename,form="unformatted",status="old")

  !! Number of t-steps = length of MPO
  READ(INPUT) N_sites

  !! Allocate empty mpo block
  CALL allocate_empty_mpo_block(mpo_block, N_sites)

  !! Sizes of matrices - allocate each matrix
  DO site=1,N_sites

      !! Read in localdim
      READ(INPUT) tempS, tempN

      !! Get the expected SNdim from input MPO
      IF(PRESENT(compactMpoFlag)) THEN 
         Sdim = Sdims(getSiteMpo(site, N_sites))
         Ndim = Ndims(getSiteMpo(site, N_sites))
      ELSE
         Sdim = Sdims(site)
         Ndim = Ndims(site)
      end IF

      !! Check South && North dims to ensure consistency
      IF((tempS .NE. Sdim) .OR. (tempN .NE. Ndim)) THEN
         WRITE(*,*) "Attempt to read in mpo from file with wrong South-North dims ", filename
         WRITE(*,*) "tempS & Sdim: ", tempS, Sdim
         WRITE(*,*) "tempN & Ndim: ", tempN, Ndim
         WRITE(*,*) "at site: ", site
         STOP
      end IF

      !! Allocate a new site of mpo_block (also sets the dims of each site)
      READ(INPUT) Wdim, Edim 
      CALL allocate_mpo_site(mpo_block(site), (/Sdim, Ndim, Wdim, Edim/))
  end DO

  !! Read mpo block from file
  DO site=1,N_sites
     READ(INPUT) mpo_block(site)%m
  end DO
    
  CLOSE(INPUT)
  WRITE(*,*) "Succesfully reloaded nolambda mpo_block from file: ", filename

 END SUBROUTINE reload_mpo_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! iMPO (Gamma-Lambda form) OUTPUT/RELOAD ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! This routine sets up the file for printing iMPO, 
 !!! then prints it using output_mpo routine
 SUBROUTINE mpo_print_function_GamLam(iGamma, iLambda, impo_name)

  TYPE(block_mpo),    INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)
  CHARACTER(LEN=*),   INTENT(IN) :: impo_name

  INTEGER :: N_sites
  CHARACTER :: dump_impo_filename*256, char_sites*16

  !! Size of impo
  N_sites = SIZE(iGamma)

  !! Determine file name for impo dump
  WRITE(char_sites,'(I6)') N_sites
  dump_impo_filename=TRIM(ADJUSTL(impo_name))//".dat"
  
  !! Output iMPO
  CALL output_mpo(dump_impo_filename, iGamma, iLambda)

  !! Write to an index file where we dumped the state.
  WRITE(*,*) "Output to:", dump_impo_filename
  WRITE(*,*)

 END SUBROUTINE mpo_print_function_GamLam




 !!! Output iMPO to file !!!
 SUBROUTINE output_mpo_GamLam(filename, iGamma, iLambda)

  CHARACTER,          INTENT(IN) :: filename*256
  TYPE(block_mpo),    INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)

  !! Output file handle
  INTEGER, PARAMETER :: OUTPUT=90

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Size of impo
  N_sites = SIZE(iGamma)

  !! Open file
  OPEN(OUTPUT, FILE=filename, form="unformatted", status="replace")

  !! Num of sites
  WRITE(OUTPUT) N_sites

  !! Size of each matrix
  DO site=1,N_sites
     !! Dimensions:
     WRITE(OUTPUT) iGamma(site)%Sdim, iGamma(site)%Ndim
     WRITE(OUTPUT) iGamma(site)%Wdim, iGamma(site)%Edim
  end DO

  !! Output each iGamma
  DO site=1,N_sites
     WRITE(OUTPUT) iGamma(site)%m
  end DO

  !! Output each iLambda
  DO site=1,N_sites
     WRITE(OUTPUT) iLambda(site)%m
  end DO

  CLOSE(OUTPUT)
   
 END SUBROUTINE output_mpo_GamLam




 !!! Reload iMPO from file !!!
 SUBROUTINE reload_mpo_GamLam(filename, iGamma, iLambda, Sdims, Ndims)

  CHARACTER,                       INTENT(IN)    :: filename*256
  TYPE(block_mpo),    ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)
  INTEGER,                         INTENT(IN)    :: Sdims(:), Ndims(:)

  INTEGER :: N_sites, site
  INTEGER :: INPUT=90
  INTEGER :: tempS, tempN
  INTEGER :: Sdim, Ndim, Wdim, Edim

  !! Open file
  OPEN(INPUT, FILE=filename, form="unformatted", status="old")

  !! Number of t-steps = length of impo
  READ(INPUT) N_sites

  !! Allocate empty mpo block for impo
  CALL allocate_empty_mpo_block(iGamma, N_sites)
  CALL allocate_empty_lambda_block(iLambda, N_sites)

  !! Sizes of matrices - allocate each matrix
  DO site=1,N_sites

      !! Read in localdim
      READ(INPUT) tempS, tempN

      !! Get the expected Sdim, Ndim from input
      Sdim = Sdims(site)
      Ndim = Ndims(site)

      !! Check South && North dims to ensure consistency
      IF((tempS .NE. Sdim) .OR. (tempN .NE. Ndim)) THEN
         WRITE(*,*) "Attempt to read in mpo from file with wrong South-North dims ", filename
         WRITE(*,*) "tempS & Sdim: ", tempS, Sdim
         WRITE(*,*) "tempN & Ndim: ", tempN, Ndim
         WRITE(*,*) "at site: ", site
         STOP
      end IF

      !! Allocate a new site of iGamma, iLambda (also sets the dims of each site)
      READ(INPUT) Wdim, Edim 
      CALL allocate_mpo_site(iGamma(site), (/Sdim, Ndim, Wdim, Edim/))
      CALL allocate_lambda_site(iLambda(site), Edim)
  end DO

  !! Read iGamma from file
  DO site=1,N_sites
     READ(INPUT) iGamma(site)%m
  end DO

  !! Read iLambda from file
  DO site=1,N_sites
     READ(INPUT) iLambda(site)%m
  end DO
    
  CLOSE(INPUT)

  WRITE(*,*) "Succesfully reloaded iMPO from file: ", filename

 END SUBROUTINE reload_mpo_GamLam

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PEPS BLOCK (PLAIN, no lambdas) OUTPUT/RELOAD ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! This routine sets up the file for printing iPEPS, 
 !!! then prints it using output_peps routine
 SUBROUTINE peps_print_function_PLAIN(peps, peps_name)

  TYPE(block_peps),   INTENT(IN) :: peps(:)
  CHARACTER(LEN=*),   INTENT(IN) :: peps_name

  CHARACTER :: dump_peps_filename*256

  !! Determine file name for ipeps dump
  dump_peps_filename=TRIM(ADJUSTL(peps_name))//".dat"
  
  !! Output PEPS
  CALL output_peps(dump_peps_filename, peps)

  !! Write to an index file where we dumped the state.
  WRITE(*,*) "Output to:", dump_peps_filename
  WRITE(*,*)

 END SUBROUTINE peps_print_function_PLAIN




 !!! Output iPEPS to file !!!
 SUBROUTINE output_peps_PLAIN(filename, peps)

  CHARACTER,          INTENT(IN) :: filename*256
  TYPE(block_peps),   INTENT(IN) :: peps(:)

  !! Output file handle
  INTEGER, PARAMETER :: OUTPUT=90

  !! Dims & indices
  INTEGER :: N_sites, site

  !! Size of ipeps
  N_sites = SIZE(peps)

  OPEN(OUTPUT, FILE=filename, form="unformatted", status="replace")

  !! Num of sites && bonds
  WRITE(OUTPUT) N_sites

  !! Size of each peps matrix
  DO site=1,N_sites
     WRITE(OUTPUT) peps(site) % LocalDim
     WRITE(OUTPUT) peps(site) % Sdim, peps(site) % Ndim
     WRITE(OUTPUT) peps(site) % Wdim, peps(site) % Edim
  end DO

  !! Output each peps site
  DO site=1,N_sites
     WRITE(OUTPUT) peps(site) % m
  end DO

  CLOSE(OUTPUT)
   
 END SUBROUTINE output_peps_PLAIN




 !!! Reload iPEPS from file !!!
 SUBROUTINE reload_peps_PLAIN(filename, LocalDim, peps)

  CHARACTER,                       INTENT(IN)    :: filename*256
  INTEGER,                         INTENT(IN)    :: LocalDim
  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: peps(:)

  INTEGER :: N_sites, site
  INTEGER :: INPUT=90
  INTEGER :: temp
  INTEGER :: Sdim, Ndim, Wdim, Edim

  OPEN(INPUT,FILE=filename,form="unformatted",status="old")

  !! Number of iPEPS sites
  READ(INPUT) N_sites

  !! Allocate empty peps block
  CALL allocate_empty_peps_block(peps, N_sites)

  !! Sizes of peps matrices - allocate each matrix
  DO site=1,N_sites

     !! Read in SNdim
     READ(INPUT) temp

     !! Check LocalDim to ensure consistency
     IF(temp .NE. LocalDim)  THEN
        WRITE(*,*) "Attempt to read in iPEPS from file with wrong LocalDim", filename
        WRITE(*,*) "temp & LocalDim: ", temp, LocalDim
        WRITE(*,*) "at site: ", site
        STOP
     end IF

     !! Allocate a new site of peps (also sets the dims of each site)
     READ(INPUT) Sdim, Ndim 
     READ(INPUT) Wdim, Edim 
     CALL allocate_peps_site(peps(site),  (/ LocalDim, Sdim, Ndim, Wdim, Edim /))
  end DO

  !! Read peps from file
  DO site=1,N_sites
     READ(INPUT) peps(site) % m
  end DO
    
  CLOSE(INPUT)

  WRITE(*,*) "Succesfully reloaded PEPS from file: ", filename

 END SUBROUTINE reload_peps_PLAIN

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! iPEPS (Gamma-Lambda form) OUTPUT/RELOAD ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!! This routine sets up the file for printing iPEPS, 
 !!! then prints it using output_peps routine
 SUBROUTINE peps_print_function_GamLam(iGamma, iLambda, ipeps_name)

  TYPE(block_peps),   INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)
  CHARACTER(LEN=*),   INTENT(IN) :: ipeps_name

  CHARACTER :: dump_ipeps_filename*256

  !! Determine file name for ipeps dump
  dump_ipeps_filename=TRIM(ADJUSTL(ipeps_name))//".dat"
  
  !! Output iPEPS
  CALL output_peps(dump_ipeps_filename, iGamma, iLambda)

  !! Write to an index file where we dumped the state.
  WRITE(*,*) "Output to:", dump_ipeps_filename
  WRITE(*,*)

 END SUBROUTINE peps_print_function_GamLam




 !!! Output iPEPS to file !!!
 SUBROUTINE output_peps_GamLam(filename, iGamma, iLambda)

  CHARACTER,          INTENT(IN) :: filename*256
  TYPE(block_peps),   INTENT(IN) :: iGamma(:)
  TYPE(block_lambda), INTENT(IN) :: iLambda(:)

  !! Output file handle
  INTEGER, PARAMETER :: OUTPUT=90

  !! Dims & indices
  INTEGER :: N_sites, site
  INTEGER :: N_bonds, bond

  !! Size of ipeps
  N_sites = SIZE(iGamma)
  N_bonds = SIZE(iLambda)

  OPEN(OUTPUT, FILE=filename, form="unformatted", status="replace")

  !! Num of sites && bonds
  WRITE(OUTPUT) N_sites
  WRITE(OUTPUT) N_bonds

  !! Size of each gamma matrix
  DO site=1,N_sites
     WRITE(OUTPUT) iGamma(site)%LocalDim
     WRITE(OUTPUT) iGamma(site)%Sdim, iGamma(site)%Ndim
     WRITE(OUTPUT) iGamma(site)%Wdim, iGamma(site)%Edim
  end DO

  !! Size of each lambda matrix
  DO bond=1,N_bonds
     WRITE(OUTPUT) iLambda(bond)%ChiDim
  end DO

  !! Output each iGamma
  DO site=1,N_sites
     WRITE(OUTPUT) iGamma(site)%m
  end DO

  !! Output each iLambda
  DO bond=1,N_bonds
     WRITE(OUTPUT) iLambda(bond)%m
  end DO

  CLOSE(OUTPUT)
   
 END SUBROUTINE output_peps_GamLam




 !!! Reload iPEPS from file !!!
 SUBROUTINE reload_peps_GamLam(filename, LocalDim, iGamma, iLambda)

  CHARACTER,                       INTENT(IN)    :: filename*256
  INTEGER,                         INTENT(IN)    :: LocalDim
  TYPE(block_peps),   ALLOCATABLE, INTENT(INOUT) :: iGamma(:)
  TYPE(block_lambda), ALLOCATABLE, INTENT(INOUT) :: iLambda(:)

  INTEGER :: N_sites, N_bonds, site, bond
  INTEGER :: INPUT=90
  INTEGER :: temp
  INTEGER :: Sdim, Ndim, Wdim, Edim, EWdim

  OPEN(INPUT,FILE=filename,form="unformatted",status="old")

  !! Number of iPEPS sites = num of gammas, Num of bonds = num of lambdas
  READ(INPUT) N_sites
  READ(INPUT) N_bonds

  !! Allocate empty pesps block for ipeps
  CALL allocate_empty_peps_block(iGamma, N_sites)
  CALL allocate_empty_lambda_block(iLambda, N_bonds)

  !! Sizes of gamma matrices - allocate each matrix
  DO site=1,N_sites

     !! Read in SNdim
     READ(INPUT) temp

     !! Check LocalDim to ensure consistency
     IF(temp .NE. LocalDim)  THEN
        WRITE(*,*) "Attempt to read in iPEPS from file with wrong LocalDim", filename
        WRITE(*,*) "temp & LocalDim: ", temp, LocalDim
        WRITE(*,*) "at site: ", site
        STOP
     end IF

     !! Allocate a new site of iGamma (also sets the dims of each site)
     READ(INPUT) Sdim, Ndim 
     READ(INPUT) Wdim, Edim 
     CALL allocate_peps_site(iGamma(site),  (/ LocalDim, Sdim, Ndim, Wdim, Edim /))
  end DO

  !! Sizes of lambda matrices - allocate each matrix
  DO bond=1,N_bonds
     READ(INPUT) EWdim
     CALL allocate_lambda_site(iLambda(bond), EWdim)
  end DO

  !! Read iGamma from file
  DO site=1,N_sites
     READ(INPUT) iGamma(site)%m
  end DO

  !! Read iLambda from file
  DO bond=1,N_bonds
     READ(INPUT) iLambda(bond)%m
  end DO
    
  CLOSE(INPUT)

  WRITE(*,*) "Succesfully reloaded iPEPS from file: ", filename

 END SUBROUTINE reload_peps_GamLam

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE mps_peps_INOUT_files
