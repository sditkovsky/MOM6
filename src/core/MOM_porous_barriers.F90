!> Function for calculating curve fit for porous barriers.
!written by sjd
module MOM_porous_barriers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs, porous_barrier_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use MOM_interface_heights, only : find_eta

implicit none ; private

#include <MOM_memory.h>

public porous_widths

!> Calculates curve fit from D_min, D_max, D_avg
interface porous_widths
  module procedure por_widths, calc_por_layer
end interface porous_widths

contains


subroutine por_widths(h, tv, G, GV, US, eta, &
    uv, por_i, por_j, Dmin, Dmax, Davg, pbv, eta_bt, halo_size, eta_to_m)
    !eta_bt, halo_size, eta_to_m not currently used
    
  !variables needed to call find_eta
  type(ocean_grid_type),                      intent(in)  :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                    intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(unit_scale_type),                      intent(in)  :: US     !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),   intent(in)  :: h      !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                      intent(in)  :: tv     !< A structure pointing to various
                                                                    !! thermodynamic variables.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(out) :: eta    !< layer interface heights
                                                                    !! [Z ~> m] or 1/eta_to_m m).
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in)  :: eta_bt !< optional barotropic
             !! variable that gives the "correct" free surface height (Boussinesq) or total water
             !! column mass per unit area (non-Boussinesq).  This is used to dilate the layer.
             !! thicknesses when calculating interfaceheights [H ~> m or kg m-2].
  integer,                          optional, intent(in)  :: halo_size !< width of halo points on
                                                                       !! which to calculate eta.

  real,                             optional, intent(in)  :: eta_to_m  !< The conversion factor from
             !! the units of eta to m; by default this is US%Z_to_m.


  !variables used to modify porous interface
  !Currently using len 100, but need a more robust solution
  integer, dimension(100),              intent(in)  :: uv    !< indicator of whether stored 
                                                             !! coordinates belong to U or V grid
  integer, dimension(100),              intent(in)  :: por_i !< zonal coordinate of modified cell
  integer, dimension(100),              intent(in)  :: por_j !< meridional coordinate of modified cell
  real, dimension(100),                 intent(in)  :: Dmin  !< parameterized topography mimimum depth
                                                             !![Z ~> m]
  real, dimension(100),                 intent(in)  :: Dmax  !< parameterized topography maximum depth
                                                             !! [Z ~> m]
  real, dimension(100),                 intent(in)  :: Davg  !< parameterized topography average depth
                                                             !! [Z ~> m]
  type(porous_barrier_ptrs),           intent(inout) :: pbv
  !should the following be intent(out) or inout?
  !real, dimension(G%IsdB:G%IedB, G%jsd:G%jed, SZK_(G)+1), & 
  !                                      intent(inout) :: por_layer_widthU !< fractional porous open 
                                                         !!width of interface on zonal flux grid [nondim] 
  !real, dimension(G%isd:G%ied, G%JsdB:G%JedB, SZK_(G)+1), & 
  !                                      intent(inout) :: por_layer_widthV !< fractional porous open
                                                         !!width of interface on merid. flux grid [nondim]
  !real, dimension(G%IsdB:G%IedB, G%jsd:G%jed, SZK_(G)), &
  !                                      intent(inout) :: por_face_areaU   !< fractional porous open
                                                         !!area on zonal flux grid [nondim]
  !real, dimension(G%isd:G%ied, G%JsdB:G%JedB, SZK_(G)), &
  !                                      intent(inout) :: por_face_areaV   !< fractional porous open
                                                         !!area on merid. flux grid [nondim]

  !local variables
  integer ii, i, j, k, nk, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real w_layer, & !< frac. open interface width of current layer [nondim]
       A_layer, & !< integral of fractional open width for current layer [nondim]
       A_layer_cum, & !< cumulative integral of fractional open width  [nondim]
       eta_s !< shallower of upwind/downwind layer [Z ~> m]
  
  !eta is zero at surface and decreases downward


  nk = SZK_(G)

  !currently no treatment for using optional find_eta arguments if present
  call find_eta(h, tv, G, GV, US, eta)


  do ii = 1, size(uv)
     if (uv(ii) == 1) then !zonal face
        I = por_i(ii) !in [IsdB, IedB]
        j = por_j(ii) !in [jsd, jed]
        do K = nk+1,1,-1
           eta_s = max(eta(I,j,K), eta(I+1,j,K)) !take shallower layer height
           if (eta_s <= Dmin(ii)) then
              pbv%por_layer_widthU(I,j,K) = 0.0
              A_layer_cum = 0.0
              if (K < nk+1) then
                 pbv%por_face_areaU(I,j,k) = 0.0; endif
           else
              call calc_por_layer(Dmin(ii), Dmax(ii), Davg(ii), eta_s, w_layer, A_layer)
              pbv%por_layer_widthU(I,j,K) = w_layer
              if (k <= nk) then
                 pbv%por_face_areaU(I,j,k) = (A_layer - A_layer_cum)/(eta(I,j,K)-eta(I,j,K+1))
              endif
               A_layer_cum = A_layer
           endif; enddo
     elseif (uv(ii) == 2) then !meridional face
        i = por_i(ii) !in [isd, ied]
        J = por_j(ii) !in [JsdB, JedB] 
        do K = nk+1,1,-1
           eta_s = max(eta(i,J,K), eta(i,J+1,K)) !take shallower layer height
           if (eta_s <= Dmin(ii)) then
              pbv%por_layer_widthV(i,J,K) = 0.0
              A_layer_cum = 0.0
              if (K < nk+1) then
                 pbv%por_face_areaV(i,J,k) = 0.0; endif
           else
              call calc_por_layer(Dmin(ii), Dmax(ii), Davg(ii), eta_s, w_layer, A_layer)
              pbv%por_layer_widthV(i,J,K) = w_layer
              if (k <= nk) then
                 pbv%por_face_areaV(i,J,k) = (A_layer - A_layer_cum)/(eta(i,J,K)-eta(i,J,K+1))
              endif
               A_layer_cum = A_layer
           endif; enddo     
     endif
  enddo

end subroutine por_widths

subroutine calc_por_layer(D_min, D_max, D_avg, h_layer, w_layer, A_layer)
!subroutine to calculate the profile fit for a layer
 
  real,            intent(in) :: D_min !< minimum topographic height [Z ~> m]
  real,            intent(in) :: D_max !< maximum topographic height [Z ~> m]
  real,            intent(in) :: D_avg !< mean topographic height [Z ~> m]
  real,            intent(in) :: h_layer !< height of interface [Z ~> m] 
  real,            intent(out) :: w_layer !< frac. open interface width of current layer [nondim]
  real,            intent(out) :: A_layer !< frac. open face area of current layer [nondim]
  !local variables
  real m, a, zeta, psi, psi_int 
  !psi_int is the integral of psi from 0 to zeta

  !three parameter fit from Adcroft 2013
  m = (D_avg - D_min)/(D_max - D_min)
  a = (1. - m)/m

  zeta = (h_layer - D_min)/(D_max - D_min)
  
  if (h_layer <= D_min) then
     w_layer = 0.0
     A_layer = 0.0
  elseif (h_layer >= D_max) then
     w_layer = 1.0
     A_layer = h_layer - D_avg
  else
     if (m < 0.5) then
        psi = zeta**(1./a)
        psi_int = (1.-m)*zeta**(1./(1.-m))
     elseif (m == 0.5) then
        psi = zeta
        psi_int = 0.5*zeta*zeta
     else 
        psi = 1. - (1. - zeta)**a
        psi_int = zeta - m + m*((1-zeta)**(1/m))
     endif
     w_layer = psi
     A_layer = (D_max - D_min)*psi_int
  endif


end subroutine calc_por_layer

end module MOM_porous_barriers
