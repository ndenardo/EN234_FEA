
subroutine user_print(n_steps)
  use Types
  use ParamIO
  use Globals, only : TIME, DTIME
  use Printparameters
  use Staticstepparameters, only: current_step_number

  implicit none
  
  integer, intent(in) :: n_steps                                 ! Current step number
  
  integer ::  left_node,right_node

  real (prec) :: vol_averaged_strain(6)               ! Contribution to strainrate from crystal plasticity
  real (prec) :: vol_averaged_stress(6)                           ! Average stress in polycrystal
  real (prec) :: uniaxial_strain                                  ! Uniaxial strain calculated from boundary conditions
!
!
! Compute the total extensional strain and strainrate from the boundary conditions
!
   call compute_hypoelastic_averages_3D(vol_averaged_strain,vol_averaged_stress)

   if (TIME<1.d-12) then
     write(user_print_units(1),'(A)') 'VARIABLES = e11,s11'
     write(user_print_units(1),'(2(1x,D12.5))') 0.0,0.0
   endif

   write(user_print_units(1),'(2(1x,D12.5))') vol_averaged_strain(1),vol_averaged_stress(1)

end subroutine user_print

subroutine compute_hypoelastic_averages_3D(vol_averaged_strain,vol_averaged_stress)
    use Types
    use ParamIO
    use Globals
    use Mesh, only : n_elements
    use Mesh, only : extract_element_data
    use Mesh, only : extract_node_data
    use User_Subroutine_Storage
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none


    real (prec), intent( out )  ::  vol_averaged_strain(6)
    real (prec), intent( out )  ::  vol_averaged_stress(6)

    ! Local variables

    integer    :: node_identifier                              ! Flag identifying node type
    integer    :: element_identifier                           ! Flag identifying element type (specified in .in file)
    integer    :: n_nodes                                      ! # nodes on the element
    integer    :: n_properties                                 ! # properties for the element
    integer    :: n_state_variables                            ! # state variables for the element
    integer    :: n_state_vars_per_intpt                       ! # state variables for one integration point

    integer, allocatable    :: node_list(:)                                ! List of nodes on the element (connectivity)

    real( prec ), allocatable   :: element_properties(:)                  ! Element or material properties, stored in order listed in input file
    real( prec ), allocatable   :: initial_state_variables(:)             ! Element state variables at start of step.
    real( prec ), allocatable   :: updated_state_variables(:)             ! State variables at end of time step

    real( prec ), allocatable   :: x(:,:)                                  ! Nodal coords x(i,a) is ith coord of ath node
    real( prec ), allocatable   :: dof_increment(:)                        ! DOF increment, using usual element storage convention
    real( prec ), allocatable   :: dof_total(:)                            ! accumulated DOF, using usual element storage convention

    integer      :: n_points,kint,i
    integer      :: n_coords, n_dof
    integer      :: iof
    integer      :: status
    integer      :: lmn

    real (prec)  ::  stress(6)                         ! Stress vector [s11, s22, e33, s12]
    real (prec)  ::  strain(6)                          ! Plastic strain rate vector
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  total_vol                         ! Total vol (actually area) of the polycrystal (in undeformed config)
    !
    !  Allocate memory to store element data.
    !  The variables specifying the size of the arrays are stored in the module user_subroutine_storage
    !  They are initialized when the input file is read, and specify the size of the arrays required to store data
    !  for any element in the mesh.  Some elements may require less storage.

    allocate(node_list(length_node_array), stat=status)
    allocate(element_properties(length_property_array), stat=status)
    allocate(initial_state_variables(length_state_variable_array), stat=status)
    allocate(updated_state_variables(length_state_variable_array), stat=status)
    allocate(x(3,length_coord_array/3), stat=status)
    allocate(dof_increment(length_dof_array), stat=status)
    allocate(dof_total(length_dof_array), stat=status)

    if (status/=0) then
       write(IOW,*) ' Error in subroutine compute_crystal_averages_2D'
       write(IOW,*) ' Unable to allocate memory for element variables '
       stop
    endif


    vol_averaged_stress = 0.d0
    vol_averaged_strain = 0.d0

    total_vol = 0.d0

    do lmn = 1,n_elements

        !
        ! Extract element and node data from global storage (see module Mesh.f90 for the source code for these subroutines)

        call extract_element_data(lmn,element_identifier,n_nodes,node_list,n_properties,element_properties, &
            n_state_variables,initial_state_variables,updated_state_variables)

 
            do i = 1, n_nodes
                iof = 3*(i-1)+1     ! Points to first DOF for the node in the dof_increment and dof_total arrays
                call extract_node_data(node_list(i),node_identifier,n_coords,x(1:3,i),n_dof, &
                    dof_increment(iof:iof+2),dof_total(iof:iof+2))
            end do

            if (n_nodes == 4) n_points = 1
            if (n_nodes == 10) n_points = 4
            if (n_nodes == 8) n_points = 8
            if (n_nodes == 20) n_points = 27

            call initialize_integration_points(n_points, n_nodes, xi, w)

            n_state_vars_per_intpt = n_state_variables/n_points

            !     --  Loop over integration points
            do kint = 1, n_points
                call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
                dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
                call invert_small(dxdxi,dxidx,determinant)

                iof = n_state_vars_per_intpt*(kint-1)+1
                stress(1:6) = updated_state_variables(iof:iof+5)
                strain(1:6) = updated_state_variables(iof+6:iof+11)

                vol_averaged_stress(1:6) = vol_averaged_stress(1:6) + stress(1:6)*w(kint)*determinant
                vol_averaged_strain(1:6) = vol_averaged_strain(1:6) + strain(1:6)*w(kint)*determinant

                total_vol = total_vol + w(kint)*determinant

            end do


    end do


    vol_averaged_stress = vol_averaged_stress/total_vol
    vol_averaged_strain = vol_averaged_strain/total_vol

    deallocate(node_list)
    deallocate(element_properties)
    deallocate(initial_state_variables)
    deallocate(updated_state_variables)
    deallocate(x)
    deallocate(dof_increment)
    deallocate(dof_total)

    return




end subroutine compute_hypoelastic_averages_3D