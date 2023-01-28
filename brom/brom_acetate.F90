!-----------------------------------------------------------------------
! brom_acetate is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Elizaveta Protsenko, Alisa Ilinskaya
!-----------------------------------------------------------------------
#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE:
!
! !INTERFACE:
   module fabm_niva_brom_acetate
!
! !DESCRIPTION:
!
! brom_acetate parameterize  the Chemical Oxygen Demand, COD (https://en.wikipedia.org/wiki/Chemical_oxygen_demand) 
! OM mineralization, nitrification, and oxidation of reduced specied of S, Mn, Fe, present in suboxic conditions.
! brom_acetate consists of 1 state variables ( in oxygen-units):
! - acetate - is an organic compound CnHaObNc that can be can be fully oxidized to inorganic nutrients 
!   with a strong oxidizing agent under acidic condition (oxygen).
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev, Elizaveta Protsenko
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_acetate
!     Variable identifiers
      type (type_state_variable_id)        :: id_oxy, id_dom, id_acetate
      type (type_dependency_id)            :: id_temp, id_salt
      type (type_diagnostic_variable_id)   :: id_acetate_miner_rate, id_mg, id_magn
!     Model parameters
      !---Organic matter mineralization---- !
       real(rk) :: r_acetate_miner, Tda, beta_da, Wacetate, Bu, mg_s_ratio
   contains
      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the OXYDEP-COD model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_acetate), intent(inout), target :: self
   integer,                      intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
       real(rk),parameter :: d_per_s = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
  ! acetate
   call self%get_parameter(self%mg_s_ratio, 'mg_s_ratio', '-', 'mg to s ratio',               default=0.10_rk)
   call self%get_parameter(self%r_acetate_miner, 'r_acetate_miner', '1/d', 'Specific rate of acetate mineralization',           default=0.10_rk,scale_factor=d_per_s)
   call self%get_parameter(self%beta_da,        'beta_da',         'nd', 'Coefficient for dependence of mineralization on t ', default=20._rk)
   call self%get_parameter(self%Tda,            'Tda',             'nd', 'Coefficient for dependence of mineralization on t ', default=13._rk)
   call self%get_parameter(self%Wacetate,      'Wacetate',          'm/s', 'vertical velocity of acetate (<0 for sinking)',         default=-1.0_rk,scale_factor=d_per_s)
   call self%get_parameter(self%Bu,             'Bu',             'nd',  'Burial coeficient for lower boundary',               default=0.25_rk)
   ! Register state variables
   call self%register_state_variable(self%id_acetate,'acetate','mmol/m**3','Acetate', 0.0_rk, minimum=0.0_rk, vertical_movement=self%Wacetate)
   ! Register link to external variables
   call self%register_state_dependency(self%id_oxy,'Oxy','mmol/m**3','OXY')
!   call self%register_state_dependency(self%id_oxy,'DOM','mmol/m**3','DOM')
   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_acetate_miner_rate,'acetate_miner_rate','mmol/m**3/d',  'acetate_miner_rate,  Mineralization of acetate with oxygen',           &
                    output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_mg,'Mg','mmol/m**3',  'Mg, concentration', output=output_time_step_integrated)

   ! Register environmental dependencies
   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_acetate),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk)                   :: acetate, oxy, dom, t, s

   real(rk) :: doxy,dacetate
 ! Rates of biogeochemical processes
   real(rk) :: acetate_miner_rate    ! oxic mineralization of acetate (1/d)
   real(rk) :: Mg            ! Mg (1/d)
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_oxy,oxy)
   _GET_(self%id_acetate,acetate)
!   _GET_(self%id_dom,dom)

   ! Retrieve current environmental conditions.
   _GET_(self%id_temp,t)              ! temperature
   _GET_(self%id_salt,s)              ! salinity

   Mg=self%mg_s_ratio*s

!--------------------------------------------------------------
! Oxic mineralization of acetate depends on T
   acetate_miner_rate   = self%r_acetate_miner*(1.+self%beta_da*yy(self%tda,t))*acetate

! Now we can summarize processes and write state variables sink/sources:
!--------------------------------------------------------------
! OXY
!--------------------------------------------------------------
! Changes of OXY due to OM production and decay!
   dacetate = -acetate_miner_rate
   doxy     = -2.0_rk*acetate_miner_rate

!-------------------------------------------------------------


!derivatives for FABM
   _SET_ODE_(self%id_oxy, doxy)
   _SET_ODE_(self%id_acetate,dacetate)

   ! Export diagnostic variables
   
   _SET_DIAGNOSTIC_(self%id_acetate_miner_rate,acetate_miner_rate)   
   _SET_DIAGNOSTIC_(self%id_Mg,Mg)
   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
!

! !INPUT PARAMETERS:
   class (type_niva_brom_acetate),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
   real(rk)                   :: acetate

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_acetate,acetate)

   ! BURYING into the sediments, mmol/m2/s (sinking rates "Wxxx" are in m/s and positive upward)
   _SET_BOTTOM_EXCHANGE_(self%id_acetate,self%Bu*self%Wacetate*acetate)

   _HORIZONTAL_LOOP_END_

   end subroutine
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: Saturation function squared
!
! !INTERFACE:
   real(rk) function yy(a,x)
!
! !DESCRIPTION:
! This is a squared Michaelis-Menten type of limiter:
! \begin{equation}\label{Y}
! Y(x_w,x) = \frac{x^2}{x_w^2+x^2}.
! \end{equation}
!
! !IN2PUT PARAMETERS:
   real(rk), intent(in)                :: a,x
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard, Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   yy=x**2/(a**2+x**2)

   end function yy
!EOC

   end module fabm_niva_brom_acetate

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
