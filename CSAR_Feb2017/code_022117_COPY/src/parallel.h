!This header file contains the interfaces necessary to use pointers
! as dummy argument for functions
!****************************************
interface
   subroutine parallel_swap(px,col, current_mesh)
     USE data_types
     TYPE(mesh_pointers),INTENT(IN) :: current_mesh
     REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px
     INTEGER, INTENT(IN) :: col(6)
   end subroutine parallel_swap
end interface

interface
   subroutine parallel_swap4(col,current_mesh)
     USE GLOBAL_DATA
     USE MYMPI
     INTEGER, INTENT(IN) :: col(7)
     TYPE(mesh_pointers),INTENT(IN) :: current_mesh
   end subroutine parallel_swap4
end interface

interface
   subroutine twolayer_swap(col,current_mesh)
     USE GLOBAL_DATA
     USE MYMPI
     INTEGER, INTENT(IN) :: col(4)
     TYPE(mesh_pointers),INTENT(IN) :: current_mesh
   end subroutine twolayer_swap
end interface
