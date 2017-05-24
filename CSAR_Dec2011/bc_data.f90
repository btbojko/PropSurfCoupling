MODULE bc_data

  IMPLICIT NONE

  SAVE

  REAL*8,ALLOCATABLE :: capd(:,:),capdd(:,:),ffnew(:,:),cape(:,:),capee(:,:)
  REAL*8,ALLOCATABLE ::  capa(:,:,:), capb(:,:,:), capc(:,:,:), capr(:,:,:)

END MODULE bc_data
