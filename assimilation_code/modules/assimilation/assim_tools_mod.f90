! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!>  A variety of operations required by assimilation.
module assim_tools_mod

!> \defgroup assim_tools assim_tools_mod
!>
!> @{
use      types_mod,       only : r8, i8, digits12, PI, missing_r8

use    options_mod,       only : get_missing_ok_status

use  utilities_mod,       only : file_exist, get_unit, check_namelist_read, do_output,    &
                                 find_namelist_in_file, error_handler,   &
                                 E_ERR, E_MSG, nmlfileunit, do_nml_file, do_nml_term,     &
                                 open_file, close_file, timestamp
use       sort_mod,       only : index_sort 
use random_seq_mod,       only : random_seq_type, random_gaussian, init_random_seq,       &
                                 random_uniform

use obs_sequence_mod,     only : obs_sequence_type, obs_type, get_num_copies, get_num_qc, &
                                 init_obs, get_obs_from_key, get_obs_def, get_obs_values, &
                                 destroy_obs

use          obs_def_mod, only : obs_def_type, get_obs_def_location, get_obs_def_time,    &
                                 get_obs_def_error_variance, get_obs_def_type_of_obs

use         obs_kind_mod, only : get_num_types_of_obs, get_index_for_type_of_obs,                   &
                                 get_quantity_for_type_of_obs, assimilate_this_type_of_obs

use       cov_cutoff_mod, only : comp_cov_factor

use       reg_factor_mod, only : comp_reg_factor

use       obs_impact_mod, only : allocate_impact_table, read_impact_table, free_impact_table

use sampling_error_correction_mod, only : get_sampling_error_table_size, &
                                          read_sampling_error_correction

use         location_mod, only : location_type, get_close_type, query_location,           &
                                 operator(==), set_location_missing, write_location,      &
                                 LocationDims, is_vertical, vertical_localization_on,     &
                                 set_vertical, has_vertical_choice, get_close_init,       &
                                 get_vertical_localization_coord, get_close_destroy,      &
                                 set_vertical_localization_coord, get_dist

use ensemble_manager_mod, only : ensemble_type, get_my_num_vars, get_my_vars,             &
                                 compute_copy_mean_var, get_var_owner_index,              &
                                 map_pe_to_task

use mpi_utilities_mod,    only : my_task_id, broadcast_send, broadcast_recv,              &
                                 sum_across_tasks, task_count, start_mpi_timer,           &
                                 read_mpi_timer, task_sync

use adaptive_inflate_mod, only : do_obs_inflate,  do_single_ss_inflate, do_ss_inflate,    &
                                 do_varying_ss_inflate,                                   &
                                 update_inflation, update_single_state_space_inflation,   &
                                 update_varying_state_space_inflation,                    &
                                 inflate_ens, adaptive_inflate_type,                      &
                                 deterministic_inflate, solve_quadratic

use time_manager_mod,     only : time_type, get_time

use assim_model_mod,      only : get_state_meta_data,                                     &
                                 get_close_obs,         get_close_state,                  &
                                 convert_vertical_obs,  convert_vertical_state

use distributed_state_mod, only : create_mean_window, free_mean_window

use quality_control_mod, only : good_dart_qc, DARTQC_FAILED_VERT_CONVERT

use inner_domain_mod, only : get_num_vars_inner_domain, get_var_index_inner_domain,       &
                             get_var_ens_inner_domain

implicit none
private

public :: filter_assim, &
          set_assim_tools_trace, &
          test_state_copies, &
          update_ens_from_weights

! Indicates if module initialization subroutine has been called yet
logical :: module_initialized = .false.


integer :: print_timestamps    = 0
integer :: print_trace_details = 0

! True if random sequence needs to be initialized
logical                :: first_inc_ran_call = .true.
type (random_seq_type) :: inc_ran_seq

integer                :: num_types = 0
real(r8), allocatable  :: cutoff_list(:)
logical                :: has_special_cutoffs
logical                :: close_obs_caching = .true.

! true if we have multiple vert choices and we're doing vertical localization
! (make it a local variable so we don't keep making subroutine calls)
logical                :: is_doing_vertical_conversion = .false.

character(len=512)     :: msgstring, msgstring2, msgstring3

! Need to read in table for off-line based sampling correction and store it
integer                :: sec_table_size
real(r8), allocatable  :: exp_true_correl(:), alpha(:)

! if adjust_obs_impact is true, read in triplets from the ascii file
! and fill this 2d impact table.
real(r8), allocatable  :: obs_impact_table(:,:)

character(len=*), parameter :: source = 'assim_tools_mod.f90'

!============================================================================

!---- namelist with default values

! Filter kind selects type of observation space filter
!      1 = EAKF filter
!      2 = ENKF
!      3 = Kernel filter
!      4 = particle filter
!      5 = random draw from posterior
!      6 = deterministic draw from posterior with fixed kurtosis
!      8 = Rank Histogram Filter (see Anderson 2011)
!
!  special_localization_obs_types -> Special treatment for the specified observation types
!  special_localization_cutoffs   -> Different cutoff value for each specified obs type
!
integer  :: filter_kind                     = 1
real(r8) :: cutoff                          = 0.2_r8
logical  :: sort_obs_inc                    = .false.
logical  :: spread_restoration              = .false.
logical  :: sampling_error_correction       = .false.
integer  :: adaptive_localization_threshold = -1
real(r8) :: adaptive_cutoff_floor           = 0.0_r8
integer  :: print_every_nth_obs             = 0

! since this is in the namelist, it has to have a fixed size.
integer, parameter   :: MAX_ITEMS = 300
character(len = 129) :: special_localization_obs_types(MAX_ITEMS)
real(r8)             :: special_localization_cutoffs(MAX_ITEMS)

logical              :: output_localization_diagnostics = .false.
character(len = 129) :: localization_diagnostics_file = "localization_diagnostics"

! Following only relevant for filter_kind = 8
logical  :: rectangular_quadrature          = .true.
logical  :: gaussian_likelihood_tails       = .false.

! False by default; if true, expect to read in an ascii table
! to adjust the impact of obs on other state vector and obs values.
logical            :: adjust_obs_impact  = .false.
character(len=256) :: obs_impact_filename = ''
logical            :: allow_any_impact_values = .false.

! These next two only affect models with multiple options
! for vertical localization:
!
! "convert_state" is false by default; it depends on the model
! what is faster - do the entire state up front and possibly
! do unneeded work, or do the conversion during the assimilation
! loop. we think this depends heavily on how much of the state
! is going to be adjusted by the obs.  for a global model
! we think false may be better; for a regional model with
! a lot of obs and full coverage true may be better.
!
! "convert_obs" is true by default; in general it seems to
! be better for each task to convert the obs vertical before
! going into the loop but again this depends on how many
! obs per task and whether the mean is distributed or
! replicated on each task.
logical :: convert_all_state_verticals_first = .false.
logical :: convert_all_obs_verticals_first   = .true.

! Not in the namelist; this var disables the experimental
! linear and spherical case code in the adaptive localization
! sections.  to try out the alternatives, set this to .false.
logical  :: only_area_adapt  = .true.

! Option to distribute the mean.  If 'false' each task will have a full
! copy of the ensemble mean, which speeds models doing vertical conversion.
! If 'true' the mean will be spread across all tasks which reduces the
! memory needed per task but requires communication if the mean is used
! for vertical conversion.  We have changed the default to be .false.
! compared to previous versions of this namelist item.
logical  :: distribute_mean  = .false.

! CCHU:
! parameters used for the PFF: define the maximum number of inner domain
! variables within all the obs
integer, parameter :: max_ni = 50

! CCHU: namelist for PFF-DART:
real(r8) :: min_kernel_value, learning_rate_fac, max_learning_rate, min_eig_ratio
real(r8) :: eakffg_inf, fixed_ker_alpha
integer  :: obs_adj_kind
logical  :: eakffg_io, adaptive_ker_io

namelist / assim_tools_nml / filter_kind, cutoff, sort_obs_inc, &
   spread_restoration, sampling_error_correction,                          &
   adaptive_localization_threshold, adaptive_cutoff_floor,                 &
   print_every_nth_obs, rectangular_quadrature, gaussian_likelihood_tails, &
   output_localization_diagnostics, localization_diagnostics_file,         &
   special_localization_obs_types, special_localization_cutoffs,           &
   distribute_mean, close_obs_caching,                                     &
   adjust_obs_impact, obs_impact_filename, allow_any_impact_values,        &
   convert_all_state_verticals_first, convert_all_obs_verticals_first,     &
   min_eig_ratio, adaptive_ker_io, min_kernel_value, fixed_ker_alpha,      &
   obs_adj_kind, learning_rate_fac, max_learning_rate, eakffg_io, eakffg_inf

!============================================================================

contains

!-------------------------------------------------------------

subroutine assim_tools_init()

integer :: iunit, io, i, j
integer :: num_special_cutoff, type_index
logical :: cache_override = .false.


! do this up front
module_initialized = .true.

! give these guys initial values at runtime *before* we read
! in the namelist.  this is to help detect how many items are
! actually given in the namelist.
special_localization_obs_types(:)  = 'null'
special_localization_cutoffs(:)    =  missing_r8

! Read the namelist entry
call find_namelist_in_file("input.nml", "assim_tools_nml", iunit)
read(iunit, nml = assim_tools_nml, iostat = io)
call check_namelist_read(iunit, io, "assim_tools_nml")

! Write the namelist values to the log file
if (do_nml_file()) write(nmlfileunit, nml=assim_tools_nml)
if (do_nml_term()) write(     *     , nml=assim_tools_nml)

! Forcing distributed_mean for single processor.
! Note null_win_mod.f90 ignores distibute_mean.
if (task_count() == 1) distribute_mean = .true.

! FOR NOW, can only do spread restoration with filter option 1 (need to extend this)
if(spread_restoration .and. .not. filter_kind == 1) then
   write(msgstring, *) 'cannot combine spread_restoration and filter_kind ', filter_kind
   call error_handler(E_ERR,'assim_tools_init:', msgstring, source)
endif

! allocate a list in all cases - even the ones where there is only
! a single cutoff value.  note that in spite of the name these
! are specific types (e.g. RADIOSONDE_TEMPERATURE, AIRCRAFT_TEMPERATURE)
! because that's what get_close() is passed.   and because i've confused
! myself several times -- we define generic kinds starting at 0, but
! the specific types are autogenerated and always start at 1.  so the
! cutoff list is never (0:num_types); it is always (num_types).
num_types = get_num_types_of_obs()
allocate(cutoff_list(num_types))
cutoff_list(:) = cutoff
has_special_cutoffs = .false.

! Go through special-treatment observation kinds, if any.
num_special_cutoff = 0
j = 0
do i = 1, MAX_ITEMS
   if(special_localization_obs_types(i) == 'null') exit
   if(special_localization_cutoffs(i) == MISSING_R8) then
      write(msgstring, *) 'cutoff value', i, ' is uninitialized.'
      call error_handler(E_ERR,'assim_tools_init:', &
                         'special cutoff namelist for types and distances do not match', &
                         source, &
                         text2='kind = '//trim(special_localization_obs_types(i)), &
                         text3=trim(msgstring))
   endif
   j = j + 1
enddo
num_special_cutoff = j

if (num_special_cutoff > 0) has_special_cutoffs = .true.

do i = 1, num_special_cutoff
   type_index = get_index_for_type_of_obs(special_localization_obs_types(i))
   if (type_index < 0) then
      write(msgstring, *) 'unrecognized TYPE_ in the special localization namelist:'
      call error_handler(E_ERR,'assim_tools_init:', msgstring, source, &
                         text2=trim(special_localization_obs_types(i)))
   endif
   cutoff_list(type_index) = special_localization_cutoffs(i)
end do

! cannot cache previous obs location if different obs types have different
! localization radii.  change it to false, and warn user why.
if (has_special_cutoffs .and. close_obs_caching) then
   cache_override = .true.
   close_obs_caching = .false.
endif

if(sampling_error_correction) then
   sec_table_size = get_sampling_error_table_size()
   allocate(exp_true_correl(sec_table_size), alpha(sec_table_size))
   ! we can't read the table here because we don't have access to the ens_size
endif

is_doing_vertical_conversion = (has_vertical_choice() .and. vertical_localization_on())

call log_namelist_selections(num_special_cutoff, cache_override)

end subroutine assim_tools_init

!-------------------------------------------------------------

subroutine filter_assim(ens_handle, obs_ens_handle, obs_seq, keys,           &
   ens_size, num_groups, obs_val_index, inflate, ENS_MEAN_COPY, ENS_SD_COPY, &
   ENS_INF_COPY, ENS_INF_SD_COPY, OBS_KEY_COPY, OBS_GLOBAL_QC_COPY,          &
   OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END, OBS_PRIOR_VAR_START,            &
   OBS_PRIOR_VAR_END, inflate_only,                                          &
   iter, max_iter, &
   initial_ker_alpha, outer_update)

type(ensemble_type),         intent(inout) :: ens_handle, obs_ens_handle
type(obs_sequence_type),     intent(in)    :: obs_seq
integer,                     intent(in)    :: keys(:)
integer,                     intent(in)    :: ens_size, num_groups, obs_val_index
! JLA: At present, this only needs to be inout because of the possible use of
! non-determinstic obs_space adaptive inflation that is not currently supported.
! Implementing that would require communication of the info about the inflation
! values as each observation updated them.
type(adaptive_inflate_type), intent(inout) :: inflate
integer,                     intent(in)    :: ENS_MEAN_COPY, ENS_SD_COPY, ENS_INF_COPY
integer,                     intent(in)    :: ENS_INF_SD_COPY
integer,                     intent(in)    :: OBS_KEY_COPY, OBS_GLOBAL_QC_COPY
integer,                     intent(in)    :: OBS_PRIOR_MEAN_START, OBS_PRIOR_MEAN_END
integer,                     intent(in)    :: OBS_PRIOR_VAR_START, OBS_PRIOR_VAR_END
logical,                     intent(in)    :: inflate_only

! changed the ensemble sized things here to allocatable

real(r8) :: obs_prior(ens_size), obs_inc(ens_size)
real(r8) :: final_factor
real(r8) :: net_a(num_groups), correl(num_groups)
real(r8) :: obs(1), obs_err_var, my_inflate, my_inflate_sd
real(r8) :: obs_qc, cutoff_rev, cutoff_orig
real(r8) :: orig_obs_prior_mean(num_groups), orig_obs_prior_var(num_groups)
real(r8) :: obs_prior_mean(num_groups), obs_prior_var(num_groups)
real(r8) :: vertvalue_obs_in_localization_coord, whichvert_real
real(r8), allocatable :: close_obs_dist(:)
real(r8), allocatable :: close_state_dist(:)

integer(i8) :: state_index
integer(i8), allocatable :: my_state_indx(:)
integer(i8), allocatable :: my_obs_indx(:)

integer :: my_num_obs, i, j, owner, owners_index, my_num_state
integer :: obs_mean_index, obs_var_index
integer :: grp_beg(num_groups), grp_end(num_groups), grp_size, grp_bot, grp_top, group
integer :: num_close_obs, obs_index, num_close_states
integer :: last_num_close_obs, last_num_close_states
integer :: base_obs_kind, base_obs_type, nth_obs
integer :: num_close_obs_cached, num_close_states_cached
integer :: num_close_obs_calls_made, num_close_states_calls_made
integer :: whichvert_obs_in_localization_coord
integer :: istatus, localization_unit
integer, allocatable :: close_obs_ind(:)
integer, allocatable :: close_state_ind(:)
integer, allocatable :: my_obs_kind(:)
integer, allocatable :: my_obs_type(:)
integer, allocatable :: my_state_kind(:)
integer, allocatable :: vstatus(:)

type(location_type)  :: base_obs_loc, last_base_obs_loc, last_base_states_loc
type(location_type)  :: dummyloc
type(location_type), allocatable :: my_obs_loc(:)
type(location_type), allocatable :: my_state_loc(:)

type(get_close_type) :: gc_obs, gc_state
type(obs_type)       :: observation
type(obs_def_type)   :: obs_def
type(time_type)      :: obs_time

logical :: allow_missing_in_state
logical :: local_single_ss_inflate
logical :: local_varying_ss_inflate
logical :: local_ss_inflate
logical :: local_obs_inflate

! CCHU for PFF-DART new arguments:
integer,  intent(in),    optional :: iter, max_iter     ! the PFF iteration
real(r8), intent(inout), optional :: initial_ker_alpha  ! for adpative kernel width
logical,  intent(in),    optional :: outer_update       ! if true, update the outer domain

! variables for broadcast send/recv: 
integer(i8) :: inner_index(max_ni)
real(r8)    :: inner_index_r8(max_ni) ! temporay storage for index in real number
! the ensemble information for inner doamin variables temporarily stored for
! broadcasting. There are two blocks: one for current pseudo time step, the
! other is for the prior values.  In each block, we start with the first inner
! domain variable, and go through all the ensemble member, then the next inner
! domain variable, ...
real(r8) :: inner_p(max_ni*ens_size), inner_c(max_ni*ens_size)
real(r8) :: hx_c(ens_size), hx_p(ens_size)
integer  :: Ni    ! the size of inner domain
real(r8) :: Ni_r8 ! the real(Ni)

! variables for calculating increments
real(r8) :: one_inner_to_state_inc(ens_size)
! output from obs_increment_pff
real(r8), allocatable :: inner_increment(:,:)      ! inner domain increment
real(r8), allocatable :: Binv_inner_increment(:,:) ! Binv* inner domain increment
real(r8), allocatable :: norm_inner_increment(:)   ! norm for increment
real(r8), allocatable :: inner_p_mean(:)           ! inner domain prior mean
real(r8), allocatable :: prior_cov(:,:), prior_cov_inv(:,:), input_inverse(:,:)
real(r8)              :: initial_alpha             ! temporary storage for initial alpha

! for adapative learning rate
real(r8) :: scalar_norm_tmp

! dummy vars:
integer :: ii, jj, inner,cc, ss
integer :: n_total_obs ! total number of global obs
character :: output_name*50
real(r8) :: ccwu, ccwu_tmp, ccwu_add
integer  :: cwui, cchu
integer :: ccwu1_ct, ccwu2_ct, n_my_state
real(r8), dimension(3) :: ccwu_array


! allocate rather than dump all this on the stack
allocate(close_obs_dist(     obs_ens_handle%my_num_vars), &
         close_obs_ind(      obs_ens_handle%my_num_vars), &
         vstatus(            obs_ens_handle%my_num_vars), &
         my_obs_indx(        obs_ens_handle%my_num_vars), &
         my_obs_kind(        obs_ens_handle%my_num_vars), &
         my_obs_type(        obs_ens_handle%my_num_vars), &
         my_obs_loc(         obs_ens_handle%my_num_vars))

allocate(close_state_dist(     ens_handle%my_num_vars), &
         close_state_ind(      ens_handle%my_num_vars), &
         my_state_indx(        ens_handle%my_num_vars), &
         my_state_kind(        ens_handle%my_num_vars), &
         my_state_loc(         ens_handle%my_num_vars))
! end alloc

! Initialize assim_tools_module if needed
if (.not. module_initialized) call assim_tools_init()

!HK make window for mpi one-sided communication
! used for vertical conversion in get_close_obs
! Need to give create_mean_window the mean copy
call create_mean_window(ens_handle, ENS_MEAN_COPY, distribute_mean)

! filter kinds 1 and 8 return sorted increments, however non-deterministic
! inflation can scramble these. the sort is expensive, so help users get better
! performance by rejecting namelist combinations that do unneeded work.
if (sort_obs_inc) then
   if(deterministic_inflate(inflate) .and. ((filter_kind == 1) .or. (filter_kind == 8))) then
      write(msgstring,  *) 'With a deterministic filter [assim_tools_nml:filter_kind = ',filter_kind,']'
      write(msgstring2, *) 'and deterministic inflation [filter_nml:inf_deterministic = .TRUE.]'
      write(msgstring3, *) 'assim_tools_nml:sort_obs_inc = .TRUE. is not needed and is expensive.'
      call error_handler(E_MSG,'', '')  ! whitespace
      call error_handler(E_MSG,'WARNING filter_assim:', msgstring, source, &
                         text2=msgstring2,text3=msgstring3)
      call error_handler(E_MSG,'', '')  ! whitespace
      sort_obs_inc = .FALSE.
   endif
endif

! Open the localization diagnostics file
if(output_localization_diagnostics .and. my_task_id() == 0) &
  localization_unit = open_file(localization_diagnostics_file, action = 'append')

! For performance, make local copies of these settings which
! are really in the inflate derived type.
local_single_ss_inflate  = do_single_ss_inflate(inflate)
local_varying_ss_inflate = do_varying_ss_inflate(inflate)
local_ss_inflate         = do_ss_inflate(inflate)
local_obs_inflate        = do_obs_inflate(inflate)

! Default to printing nothing
nth_obs = -1

! Divide ensemble into num_groups groups.
! make sure the number of groups and ensemble size result in
! at least 2 members in each group (to avoid divide by 0) and
! that the groups all have the same number of members.
grp_size = ens_size / num_groups
if ((grp_size * num_groups) /= ens_size) then
   write(msgstring,  *) 'The number of ensemble members must divide into the number of groups evenly.'
   write(msgstring2, *) 'Ensemble size = ', ens_size, '  Number of groups = ', num_groups
   write(msgstring3, *) 'Change number of groups or ensemble size to avoid remainders.'
   call error_handler(E_ERR,'filter_assim:', msgstring, source, &
                         text2=msgstring2,text3=msgstring3)
endif
if (grp_size < 2) then
   write(msgstring,  *) 'There must be at least 2 ensemble members in each group.'
   write(msgstring2, *) 'Ensemble size = ', ens_size, '  Number of groups = ', num_groups
   write(msgstring3, *) 'results in < 2 members/group.  Decrease number of groups or increase ensemble size'
   call error_handler(E_ERR,'filter_assim:', msgstring, source, &
                         text2=msgstring2,text3=msgstring3)
endif
do group = 1, num_groups
   grp_beg(group) = (group - 1) * grp_size + 1
   grp_end(group) = grp_beg(group) + grp_size - 1
enddo

! Put initial value of state space inflation in copy normally used for SD
! This is to avoid weird storage footprint in filter
ens_handle%copies(ENS_SD_COPY, :) = ens_handle%copies(ENS_INF_COPY, :)

! For single state or obs space inflation, the inflation is like a token
! Gets passed from the processor with a given obs on to the next
if(local_single_ss_inflate) then
   my_inflate    = ens_handle%copies(ENS_INF_COPY,    1)
   my_inflate_sd = ens_handle%copies(ENS_INF_SD_COPY, 1)
end if

! Get info on my number and indices for obs
my_num_obs = get_my_num_vars(obs_ens_handle)
call get_my_vars(obs_ens_handle, my_obs_indx)

! Construct an observation temporary
call init_obs(observation, get_num_copies(obs_seq), get_num_qc(obs_seq))

! Get the locations for all of my observations
! HK I would like to move this to before the calculation of the forward operator so you could
! overwrite the vertical location with the required localization vertical coordinate when you
! do the forward operator calculation
call get_my_obs_loc(obs_ens_handle, obs_seq, keys, my_obs_loc, my_obs_kind, my_obs_type, obs_time)

if (convert_all_obs_verticals_first .and. is_doing_vertical_conversion) then
   ! convert the vertical of all my observations to the localization coordinate
   if (obs_ens_handle%my_num_vars > 0) then
      call convert_vertical_obs(ens_handle, obs_ens_handle%my_num_vars, my_obs_loc, &
                                my_obs_kind, my_obs_type, get_vertical_localization_coord(), vstatus)
      do i = 1, obs_ens_handle%my_num_vars
         if (good_dart_qc(nint(obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i)))) then
            !> @todo Can I just use the OBS_GLOBAL_QC_COPY? Is it ok to skip the loop?
            if (vstatus(i) /= 0) obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, i) = DARTQC_FAILED_VERT_CONVERT
         endif
      enddo
   endif 
endif

! Get info on my number and indices for state
my_num_state = get_my_num_vars(ens_handle)
call get_my_vars(ens_handle, my_state_indx)

! Get the location and kind of all my state variables
do i = 1, ens_handle%my_num_vars
   call get_state_meta_data(my_state_indx(i), my_state_loc(i), my_state_kind(i))
end do

!> optionally convert all state location verticals
if (convert_all_state_verticals_first .and. is_doing_vertical_conversion) then
   if (ens_handle%my_num_vars > 0) then
      call convert_vertical_state(ens_handle, ens_handle%my_num_vars, my_state_loc, my_state_kind,  &
                                  my_state_indx, get_vertical_localization_coord(), istatus)
   endif
endif

! Get mean and variance of each group's observation priors for adaptive inflation
! Important that these be from before any observations have been used
if(local_ss_inflate) then
   do group = 1, num_groups
      obs_mean_index = OBS_PRIOR_MEAN_START + group - 1
      obs_var_index  = OBS_PRIOR_VAR_START  + group - 1
         call compute_copy_mean_var(obs_ens_handle, grp_beg(group), grp_end(group), &
           obs_mean_index, obs_var_index)
   end do
endif

! Initialize the method for getting state variables close to a given ob on my process
if (has_special_cutoffs) then
   call get_close_init(gc_state, my_num_state, 2.0_r8*cutoff, my_state_loc, 2.0_r8*cutoff_list)
else
   call get_close_init(gc_state, my_num_state, 2.0_r8*cutoff, my_state_loc)
endif

! Initialize the method for getting obs close to a given ob on my process
if (has_special_cutoffs) then
   call get_close_init(gc_obs, my_num_obs, 2.0_r8*cutoff, my_obs_loc, 2.0_r8*cutoff_list)
else
   call get_close_init(gc_obs, my_num_obs, 2.0_r8*cutoff, my_obs_loc)
endif

if (close_obs_caching) then
   ! Initialize last obs and state get_close lookups, to take advantage below
   ! of sequential observations at the same location (e.g. U,V, possibly T,Q)
   ! (this is getting long enough it probably should go into a subroutine. nsc.)
   last_base_obs_loc           = set_location_missing()
   last_base_states_loc        = set_location_missing()
   last_num_close_obs          = -1
   last_num_close_states       = -1
   num_close_obs_cached        = 0
   num_close_states_cached     = 0
   num_close_obs_calls_made    = 0
   num_close_states_calls_made = 0
endif

allow_missing_in_state = get_missing_ok_status()


! CCHU: initialize the state increment for PFF-DART
ens_handle%state_inc(:,:) = 0 ! initialzation

! Loop through all the (global) observations sequentially

SEQUENTIAL_OBS: do i = 1, obs_ens_handle%num_vars

   ! Some compilers do not like mod by 0, so test first.
   if (print_every_nth_obs > 0) nth_obs = mod(i, print_every_nth_obs)

   ! If requested, print out a message every Nth observation
   ! to indicate progress is being made and to allow estimates
   ! of how long the assim will take.
   if (nth_obs == 0) then
      write(msgstring, '(2(A,I8))') 'Processing observation ', i, &
                                         ' of ', obs_ens_handle%num_vars
      if (print_timestamps == 0) then
         call error_handler(E_MSG,'filter_assim',msgstring)
      else
         call timestamp(trim(msgstring), pos="brief")
      endif
   endif

   ! CCHU: 2022/01/18:
   ! However, inner domain info is ONLY stored in the pe which owns the obs, so need to broadcast

   ! Every pe has information about the global obs sequence
   call get_obs_from_key(obs_seq, keys(i), observation)
   call get_obs_def(observation, obs_def)
   base_obs_loc = get_obs_def_location(obs_def)

   ! CCHU 2023/03/24
   obs_err_var = get_obs_def_error_variance(obs_def)
   base_obs_type = get_obs_def_type_of_obs(obs_def)
   if (base_obs_type > 0) then
      base_obs_kind = get_quantity_for_type_of_obs(base_obs_type)
   else
      call get_state_meta_data(-1 * int(base_obs_type,i8), dummyloc, base_obs_kind)  ! identity obs
   endif
   ! Get the value of the observation
   call get_obs_values(observation, obs, obs_val_index)

   ! should be different for different obs operator
   !obs_err_var = (0.25*obs(1))**2

   ! Find out who has this observation and where it is
   call get_var_owner_index(ens_handle, int(i,i8), owner, owners_index)


   ! Following block is done only by the owner of this observation
   !-----------------------------------------------------------------------
   if(ens_handle%my_pe == owner) then
      ! each task has its own subset of all obs.  if they were converted in the
      ! vertical up above, then we need to broadcast the new values to all the other
      ! tasks so they're computing the right distances when applying the increments.
      if (is_doing_vertical_conversion) then
         vertvalue_obs_in_localization_coord = query_location(my_obs_loc(owners_index), "VLOC")
         whichvert_obs_in_localization_coord = query_location(my_obs_loc(owners_index), "WHICH_VERT")
      else
         vertvalue_obs_in_localization_coord = 0.0_r8
         whichvert_obs_in_localization_coord = 0
      endif

      obs_qc = obs_ens_handle%copies(OBS_GLOBAL_QC_COPY, owners_index)

      ! Only value of 0 for DART QC field should be assimilated
      IF_QC_IS_OKAY: if(nint(obs_qc) ==0) then

         ! CCHU: This part is revised for PFF-DART
         ! The info needs to be broadcast (ie, local to the owner of obs)
         ! includes two (1) ensemble in obs space (2) ensemble in inner domain

         ! === information of ensemble in observation space H(x) ===

         ! CCHU: 2022/01/18: for PFF-DART obs_prior is confusing because of the
         ! iterations. To differentiate:
         ! hx_p = prior H(x)
         ! hx_c = current H(x)
         ! Here obs_ens_handle%copies refers to the ensemble of H(x) at the current iterations. 
         ! So, rename obs_prior as hx_c (H(x) current)

         ! The following is the original DART code: 
         ! obs_prior = obs_ens_handle%copies(1:ens_size, owners_index)

         hx_p = obs_ens_handle%hx_prior(1:ens_size)
         hx_c = obs_ens_handle%copies(1:ens_size, owners_index)

         ! Note that these are before DA starts, so can be different from current obs_prior
         orig_obs_prior_mean = obs_ens_handle%copies(OBS_PRIOR_MEAN_START: &
            OBS_PRIOR_MEAN_END, owners_index)
         orig_obs_prior_var  = obs_ens_handle%copies(OBS_PRIOR_VAR_START:  &
            OBS_PRIOR_VAR_END, owners_index)

         ! === information of ensemble in inner domain ===
         ! we need: prior values, current values, inner domain size, inner domain (global) indices

         ! Ni = # of variables in the inner domain for this obs
         Ni = get_num_vars_inner_domain(owners_index)

         ! Get the indices of state variables for inner domain:
         inner_index = 0
         call get_var_index_inner_domain(owners_index, inner_index(1:Ni), Ni)

         ! It is easier to broadcast real numbers than integer, so change the
         ! format of the integer to real number (r8):
         inner_index_r8 = 1.0_r8*inner_index
         Ni_r8          = 1.0_r8*Ni

         ! The following info has to be "column vector", so need to "reshape" matrices
         ! (1) inner_p: prior of inner domain variables "of this obs"
         ! (2) inner_c: current values of inner domain variables "of this obs"

         do jj = 1, Ni
            inner_p(ens_size*(jj-1)+1: ens_size*jj) = obs_ens_handle%inner_prior(:,jj)
            ! this is obtained from inner_domain module:
            call get_var_ens_inner_domain(owners_index, jj, inner_c(ens_size*(jj-1)+1: ens_size*jj))
         enddo

         ! CCHU: save inner domain variables:
         !if (my_task_id().eq.0) then

            !n_total_obs = obs_ens_handle%num_vars

            !if ((iter.eq.1).and.(i.eq.1)) then
               ! open a new file:
               !output_name='./output/PFF_linear_two_obs_x.dat'
               !output_name='./output/PFF_linear_two_obs_unobs_x_test.dat'
               !open(12,file=output_name,status='new',form='unformatted',access='direct',recl=4*Ni*ens_size)
            !endif

            !write(12, rec=n_total_obs*(iter-1)+i ) inner_current(1:Ni*ens_size)
            !write(12, rec=n_total_obs*(iter-1)+i ) ens_handle%copies(2,1:ens_size)
            !write(*,*) sum(ens_handle%copies(2,1:ens_size))/ens_size
         !end if

      endif IF_QC_IS_OKAY

      !Broadcast the info from this obs to all other processes
      ! orig_obs_prior_mean and orig_obs_prior_var only used with adaptive inflation
      ! my_inflate and my_inflate_sd only used with single state space inflation
      ! vertvalue_obs_in_localization_coord and whichvert_real only used for
      ! vertical coordinate transformation
      whichvert_real = real(whichvert_obs_in_localization_coord, r8)

      ! broadcast in the following:

      call broadcast_send(map_pe_to_task(ens_handle, owner), hx_p, hx_c,  &
         inner_p, inner_c, inner_index_r8, orig_obs_prior_mean, orig_obs_prior_var,   &
         scalar1=obs_qc, scalar2=vertvalue_obs_in_localization_coord,      &
         scalar3=whichvert_real, scalar4=my_inflate, scalar5=my_inflate_sd, scalar6=Ni_r8)

      !if(my_task_id()==0) then 
          !print*, 'scalar1 = ',obs_qc
          !print*, 'scalar2 = ',vertvalue_obs_in_localization_coord
          !print*, 'scalar3 = ',whichvert_real
          !print*, 'scalar4 = ',my_inflate
          !print*, 'scalar5 = ',my_inflate_sd
          !print*, 'scalar6 = ',Ni_r8
          !print*, 'norm_obs = ', norm_obs(1)
      !endif

!      call broadcast_send(map_pe_to_task(ens_handle, owner), obs_prior,    &
!         orig_obs_prior_mean, orig_obs_prior_var,                          &
!         scalar1=obs_qc, scalar2=vertvalue_obs_in_localization_coord,      &
!         scalar3=whichvert_real, scalar4=my_inflate, scalar5=my_inflate_sd)


   ! Next block is done by processes that do NOT own this observation
   !-----------------------------------------------------------------------
   else

      call broadcast_recv(map_pe_to_task(ens_handle, owner), hx_p, hx_c,  &
         inner_p, inner_c, inner_index_r8, orig_obs_prior_mean, orig_obs_prior_var,   &
         scalar1=obs_qc, scalar2=vertvalue_obs_in_localization_coord,      &
         scalar3=whichvert_real, scalar4=my_inflate, scalar5=my_inflate_sd, scalar6=Ni_r8)

      ! CCHU: turns the r8 real number back to integer:
      inner_index = nint(inner_index_r8)
      Ni          = nint(Ni_r8)

!      call broadcast_recv(map_pe_to_task(ens_handle, owner), obs_prior,    &
!         orig_obs_prior_mean, orig_obs_prior_var,                          & 
!         scalar1=obs_qc, scalar2=vertvalue_obs_in_localization_coord,      &
!         scalar3=whichvert_real, scalar4=my_inflate, scalar5=my_inflate_sd)

      whichvert_obs_in_localization_coord = nint(whichvert_real)

   endif

      !if(my_task_id()==1) then
          !print*, 'scalar1 = ',obs_qc
          !print*, 'scalar2 = ',vertvalue_obs_in_localization_coord
          !print*, 'scalar3 = ',whichvert_real
          !print*, 'scalar4 = ',my_inflate
          !print*, 'scalar5 = ',my_inflate_sd
          !print*, 'scalar6 = ',Ni_r8
          !print*, 'norm_obs = ', norm_obs(1)
      !endif

   !-----------------------------------------------------------------------

   ! Everybody is doing this section, cycle if qc is bad
   if(nint(obs_qc) /= 0) cycle SEQUENTIAL_OBS

   !> all tasks must set the converted vertical values into the 'base' version of this loc
   !> because that's what we pass into the get_close_xxx() routines below.
   if (is_doing_vertical_conversion) &
      call set_vertical(base_obs_loc, vertvalue_obs_in_localization_coord, whichvert_obs_in_localization_coord)

   ! CCHU
   ! once every pe knows Ni "of this obs" (each obs has different Ni),
   !  allocate the dimension of inner domain increment:
   allocate(inner_increment     (ens_size, Ni) )
   allocate(norm_inner_increment(Ni)           )

   if (.not. outer_update) then

      ! Compute observation space increments for each group
      do group = 1, num_groups
         grp_bot = grp_beg(group); grp_top = grp_end(group)

         ! CCHU:: call the filter (PFF) here:

         ! initial value of kernel alpha (a random number, doesn't matter)
         initial_alpha = initial_ker_alpha
 
         call obs_increment(hx_c(grp_bot:grp_top), grp_size,                              &
              obs(1), obs_err_var, obs_inc(grp_bot:grp_top), inflate, my_inflate,         &
              my_inflate_sd, net_a(group),                                                &
              iter=iter, max_iter=max_iter, initial_alpha=initial_alpha,                  &
              Ni=Ni, inner_index=inner_index(1:Ni),                                       &
              inner_p=inner_p(1:ens_size*Ni), inner_c=inner_c(1:ens_size*Ni), hx_p=hx_p,  & 
              inner_increment=inner_increment, &
              norm_inner_increment=norm_inner_increment, &
              base_obs_type=base_obs_type, base_obs_loc=base_obs_loc )

         ! For the 1st iteration PFF, will determine an kernel alpha value
         ! save this number and this will be used (not changed anymore) for the remaining iterations
         if ( iter.eq.1 ) initial_ker_alpha = initial_alpha

         ! save inner_increment into the owner pe for verfication:
         if ( ens_handle%my_pe == owner) then
            obs_ens_handle%inner_inc(1:ens_size,1:Ni,iter) = inner_increment
         endif    

         ! save the vector norm for all pe:
         obs_ens_handle%vector_norm(1:Ni, iter) = norm_inner_increment

         ! calculate the "scalar norm" by combining the vector norm for all pe:
         scalar_norm_tmp = 0.0_r8
         do jj=1,Ni
            scalar_norm_tmp = scalar_norm_tmp + norm_inner_increment(jj)/ &
                                                obs_ens_handle%vector_norm(jj,1)
         enddo
         obs_ens_handle%scalar_norm(iter) = scalar_norm_tmp/Ni

         ! Also compute prior mean and variance of obs for efficiency here
         obs_prior_mean(group) = sum(hx_c(grp_bot:grp_top)) / grp_size
         obs_prior_var(group) = sum((hx_c(grp_bot:grp_top) - obs_prior_mean(group))**2) / &
            (grp_size - 1)
         if (obs_prior_var(group) < 0.0_r8) obs_prior_var(group) = 0.0_r8
      end do
   
      ! CCHU 2023/02/14:
      ! update the inner domain 

      do inner = 1, Ni
         call get_var_owner_index(ens_handle, inner_index(inner), owner, owners_index)
         if (my_task_id()==owner) then
           !print*, 'inner (global) id  =', inner_index(inner), 'owner = pe ',owner
           !print*, 'owner (local) index =', owners_index, '(gloabl) index =',ens_handle%my_vars(owners_index)
            ens_handle%state_inc(1:ens_size, owners_index) = inner_increment(1:ens_size, inner)
         endif
      enddo

      cycle SEQUENTIAL_OBS

   endif ! if .not. outer_update


   ! Compute updated values for single state space inflation
   if(local_single_ss_inflate) then
      ! Update for each group separately
      do group = 1, num_groups
         call update_single_state_space_inflation(inflate, my_inflate, my_inflate_sd, &
            ens_handle%copies(ENS_SD_COPY, 1), orig_obs_prior_mean(group), &
            orig_obs_prior_var(group), obs(1), obs_err_var, grp_size, inflate_only)
      end do
   endif
  
   ! Adaptive localization needs number of other observations within localization radius.
   ! Do get_close_obs first, even though state space increments are computed before obs increments.
   call  get_close_obs_cached(gc_obs, base_obs_loc, base_obs_type,      &
      my_obs_loc, my_obs_kind, my_obs_type, num_close_obs, close_obs_ind, close_obs_dist,  &
      ens_handle, last_base_obs_loc, last_num_close_obs, num_close_obs_cached,             &
      num_close_obs_calls_made)

   ! set the cutoff default, keep a copy of the original value, and avoid
   ! looking up the cutoff in a list if the incoming obs is an identity ob
   ! (and therefore has a negative kind).  specific types can never be 0;
   ! generic kinds (not used here) start their numbering at 0 instead of 1.
   if (base_obs_type > 0) then
      cutoff_orig = cutoff_list(base_obs_type)
   else
      cutoff_orig = cutoff
   endif

   ! JLA, could also cache for adaptive_localization which may be expensive?
   call adaptive_localization_and_diags(cutoff_orig, cutoff_rev, adaptive_localization_threshold, &
      adaptive_cutoff_floor, num_close_obs, close_obs_ind, close_obs_dist, my_obs_type, &
      i, base_obs_loc, obs_def, localization_unit)

   ! Find state variables on my process that are close to observation being assimilated
   call  get_close_state_cached(gc_state, base_obs_loc, base_obs_type,      &
      my_state_loc, my_state_kind, my_state_indx, num_close_states, close_state_ind, close_state_dist,  &
      ens_handle, last_base_states_loc, last_num_close_states, num_close_states_cached,              &
      num_close_states_calls_made)

   !call test_close_obs_dist(close_state_dist, num_close_states, i)

   ! CCHU: test the radius of influence
   ! write(*,*) 'close state ind = ', close_state_ind(1:num_close_states)

   ! Loop through to update each of my state variables that are close to obs location (current)
   ! inner domain location"s" (need to write the code for this in the future)


   ! === Binv * inner domain increment ===
   allocate(prior_cov           (Ni,Ni)       )
   allocate(prior_cov_inv       (Ni,Ni)       )
   allocate(input_inverse       (Ni,Ni)       )
   allocate(Binv_inner_increment(ens_size, Ni)) 
   allocate(inner_p_mean        (Ni)          ) 
 
   ! calculate the prior covariance (for inner domain)
   
   do inner = 1, Ni
      inner_p_mean(inner) = sum(inner_p(ens_size*(inner-1)+1: ens_size*inner))/ens_size
   enddo

   do ii=1,Ni
      do jj=1,Ni
         if (jj.ge.ii) then
            prior_cov(jj,ii) =  &
               dot_product( inner_p(ens_size*(ii-1)+1: ens_size*ii) - inner_p_mean(ii), &
                            inner_p(ens_size*(jj-1)+1: ens_size*jj) - inner_p_mean(jj) )/(ens_size-1)
         else
            prior_cov(jj,ii) = prior_cov(ii,jj)
         endif
      enddo
   enddo

   ! calculate the inverse of prior covariance matrix:
   input_inverse = prior_cov ! input_inverse is a dummy variable that passes prior_cov to the inverse

   ! use SVD to calculate the inverse of prior covariance (of inner domain)
   call svd_pseudo_inverse(input_inverse,prior_cov_inv,Ni,Ni,min_eig_ratio)
 
   ! Binv_increment:
   do inner = 1, Ni
      inner_increment(:,inner) = inner_c(ens_size*(inner-1)+1: ens_size*inner) - &
                                 inner_p(ens_size*(inner-1)+1: ens_size*inner)
   enddo

   ! print out the inner domain increment!
   if (my_task_id()==0) print*, ' '
   if (my_task_id()==0) print*, '       obs value =', obs
   if (my_task_id()==0) print*, '       H(x) prior mean =', sum(hx_p)/(ens_size*1.0_r8)
   if (my_task_id()==0) print*, '       H(x) post  mean =', sum(hx_c)/(ens_size*1.0_r8)
   if (my_task_id()==0) print*, '       |inner inc| =', sum(abs(inner_increment),dim=1)/(ens_size*1.0_r8)
!   if (my_task_id()==0) print*, 'inner inc =',  sum(inner_increment(1:3,:),dim=2)/4

   Binv_inner_increment = matmul(inner_increment, prior_cov_inv)


   INNER_DOMAIN_INC: do inner = 1, Ni

      ! For PFF (the following should be the inner domain update):
      ! overwrite "obs_prior", "obs_inc", "obs_prior_mean" and "obs_prior_var"
      ! In PFF, obs_prior = inner domain prior
      !         obs_inc   = inner domain increment
      !         and so on...
      obs_prior = inner_p(ens_size*(inner-1)+1: ens_size*inner)  ! overwrite obs_prior for PFF
      obs_inc   = Binv_inner_increment(:,inner)                  ! overwrite obs_inc for PFF
      obs_prior_mean = sum( obs_prior )/(1.0_r8*ens_size)        ! overwrite obs_prior_mean for PFF
      obs_prior_var = 0.0_r8 ! overwrite obs_prior_var for PFF

      do cc = 1, ens_size
         obs_prior_var = obs_prior_var + (obs_prior(cc)-obs_prior_mean)**2 / (1.0_r8*ens_size-1)
      enddo
      
      ! not sure why the below does not work...
!     ! obs_prior_var  = sum( (obs_prior - obs_prior_mean)**2 )/(1.0_r8*ens_size-1)
!      write(*,*) 'obs_prior - mean  = ', obs_prior - obs_prior_mean
!      write(*,*) 'mean =', obs_prior_mean, 'variance = ', obs_prior_var

      ! Need to modify below espeically for non-local observations
      ! uncomment the below two lines to get "B localization"
      ! comment out the below two lines to do the original DART 'increment localization'
   
      !call get_state_meta_data(inner_index(inner),base_obs_loc) 
      !base_obs_type  = -inner_index(inner)

      !if(my_task_id()==0) print*, 'base_obs_type = ', base_obs_type

      
      STATE_UPDATE: do j = 1, num_close_states
         state_index = close_state_ind(j)

         ! CCHU: originally below line does not exist!
         ! base_obs_loc should be the location of the inner domain variable
         ! my_state_loc(state_index) should be the location of the j-th state
         ! needs to be updated (j=1,num_close_states)

         !close_state_dist(j) = get_dist(base_obs_loc, my_state_loc(state_index))

         !print*, 'dist = ',close_state_dist(j)
         !write(*,*) 'dist between ',state_index,'and inner domain',base_obs_type, 'is ', close_state_dist(j)

         if ( allow_missing_in_state ) then
            ! Don't allow update of state ensemble with any missing values
            if (any(ens_handle%copies(1:ens_size, state_index) == MISSING_R8)) cycle STATE_UPDATE
         endif

         ! Compute the covariance localization and adjust_obs_impact factors
         ! (module storage)
         final_factor = cov_and_impact_factors(base_obs_loc, base_obs_type, my_state_loc(state_index), &
            my_state_kind(state_index), close_state_dist(j), cutoff_rev)
         !final_factor = 1.0_r8


         !if (my_task_id()==0) then
         !   if (any(inner_index(1:Ni)==state_index)) then
         !      print*,'state_index = ', state_index
         !      print*,'inner_index = ', -base_obs_type
         !      print*,'distance    = ', close_state_dist(j)
         !      print*,'final_factor= ', final_factor
         !   endif
         !endif

         if ( final_factor<=0.0_r8 ) cycle STATE_UPDATE

         ! OLD obs_updates_ens:
         !call obs_updates_ens(ens_size, num_groups, ens_handle%copies(1:ens_size, state_index), &
         !   updated_ens, my_state_loc(state_index), my_state_kind(state_index), obs_prior, obs_inc, &
         !   obs_prior_mean, obs_prior_var, base_obs_loc, base_obs_type, obs_time, &
         !   close_state_dist(j), cutoff_rev, net_a, adjust_obs_impact, obs_impact_table, &
         !   grp_size, grp_beg, grp_end, i, my_state_indx(state_index), final_factor, correl, &
         !   state_prior=pstate(1:ens_size, state_index))

         ! CCWU: Just for PFF test (NEED TO REMOVE THIS in the future)
         !local_varying_ss_inflate = .false.
         !inflate_only             = .false.     
        
         !print*,'inner_ind = ', -base_obs_type, 'state index =', state_index
 
         call obs_updates_ens(ens_size, num_groups, ens_handle%copies(1:ens_size, state_index), &
            my_state_loc(state_index), my_state_kind(state_index), obs_prior, obs_inc, &
            obs_prior_mean, obs_prior_var, base_obs_loc, base_obs_type, obs_time, &
            net_a, grp_size, grp_beg, grp_end, i, &
            my_state_indx(state_index), final_factor, correl, .false., .false., &
            state_prior=ens_handle%state_prior(1:ens_size, state_index), &
            one_inner_to_state_inc=one_inner_to_state_inc)

         ! Compute spatially-varying state space inflation
         if(local_varying_ss_inflate) then
            do group = 1, num_groups
               call update_varying_state_space_inflation(inflate,                     &
                  ens_handle%copies(ENS_INF_COPY, state_index),                       &
                  ens_handle%copies(ENS_INF_SD_COPY, state_index),                    &
                  ens_handle%copies(ENS_SD_COPY, state_index),                        &
                  orig_obs_prior_mean(group), orig_obs_prior_var(group), obs(1),      &
                  obs_err_var, grp_size, final_factor, correl(group), inflate_only)
            end do
         endif

         ! CCHU: important change here!
         ! Do NOT update the state variable immediately, ONLY store the increment for now, and will
         ! be added to the state variable later

         !if(.not. inflate_only) ens_handle%copies(1:ens_size, state_index) = updated_ens

         if (.not. inflate_only) then
            ens_handle%state_inc(1:ens_size,state_index) = &
               ens_handle%state_inc(1:ens_size,state_index) + one_inner_to_state_inc
         endif

         ! Compute spatially-varying state space inflation
         if(local_varying_ss_inflate .and. final_factor > 0.0_r8) then
            do group = 1, num_groups
               call update_varying_state_space_inflation(inflate,                     &
                  ens_handle%copies(ENS_INF_COPY, state_index),                       &
                  ens_handle%copies(ENS_INF_SD_COPY, state_index),                    &
                  ens_handle%copies(ENS_SD_COPY, state_index),                        &
                  orig_obs_prior_mean(group), orig_obs_prior_var(group), obs(1),      &
                  obs_err_var, grp_size, final_factor, correl(group), inflate_only)
            end do
         endif

      end do STATE_UPDATE
   end do INNER_DOMAIN_INC

   if (allocated(inner_increment) )      deallocate( inner_increment      )
   if (allocated(Binv_inner_increment) ) deallocate( Binv_inner_increment )
   if (allocated(norm_inner_increment) ) deallocate( norm_inner_increment )
   if (allocated(prior_cov) )            deallocate( prior_cov            )
   if (allocated(prior_cov_inv) )        deallocate( prior_cov_inv        )
   if (allocated(input_inverse) )        deallocate( input_inverse        )
   if (allocated(inner_p_mean) )         deallocate( inner_p_mean         )

! BELOW is the original DART code (use obs_inc to update states)
!         ! Compute the covariance localization and adjust_obs_impact factors (module storage)
!            final_factor = cov_and_impact_factors(base_obs_loc, base_obs_type, my_obs_loc(obs_index), &
!            my_obs_kind(obs_index), close_obs_dist(j), cutoff_rev)
!
!            if(final_factor <= 0.0_r8) cycle OBS_UPDATE
!
!            call obs_updates_ens(ens_size, num_groups, obs_ens_handle%copies(1:ens_size, obs_index), &
!               my_obs_loc(obs_index), my_obs_kind(obs_index), obs_prior, obs_inc, &
!               obs_prior_mean, obs_prior_var, base_obs_loc, base_obs_type, obs_time, &
!               net_a, grp_size, grp_beg, grp_end, i, &
!               -1*my_obs_indx(obs_index), final_factor, correl, .false., inflate_only)
!         endif
!      end do OBS_UPDATE
!   endif

end do SEQUENTIAL_OBS

! Every pe needs to get the current my_inflate and my_inflate_sd back
if(local_single_ss_inflate) then
   ens_handle%copies(ENS_INF_COPY, :) = my_inflate
   ens_handle%copies(ENS_INF_SD_COPY, :) = my_inflate_sd
end if

! Free up the storage
call destroy_obs(observation)
call get_close_destroy(gc_state)
call get_close_destroy(gc_obs)

! do some stats - being aware that unless we do a reduce() operation
! this is going to be per-task.  so only print if something interesting
! shows up in the stats?  maybe it would be worth a reduce() call here?

! Assure user we have done something
if (print_trace_details >= 0) then

if (iter.eq.99999) then
   write(msgstring, '(A,I8,A)') 'Processed', obs_ens_handle%num_vars, ' total observations'
   call error_handler(E_MSG,'filter_assim:',msgstring)
endif
endif

! diagnostics for stats on saving calls by remembering obs at the same location.
! change .true. to .false. in the line below to remove the output completely.

! CCWU: do not want the output!
close_obs_caching = .false.

if (close_obs_caching) then
   if (num_close_obs_cached > 0 .and. do_output()) then
      print *, "Total number of calls made    to get_close_obs for obs/states:    ", &
                num_close_obs_calls_made + num_close_states_calls_made
      print *, "Total number of calls avoided to get_close_obs for obs/states:    ", &
                num_close_obs_cached + num_close_states_cached
      if (num_close_obs_cached+num_close_obs_calls_made+ &
          num_close_states_cached+num_close_states_calls_made > 0) then
         print *, "Percent saved: ", 100.0_r8 * &
                   (real(num_close_obs_cached+num_close_states_cached, r8) /  &
                   (num_close_obs_calls_made+num_close_obs_cached +           &
                    num_close_states_calls_made+num_close_states_cached))
      endif
   endif
endif

!call test_state_copies(ens_handle, 'end')

! Close the localization diagnostics file
if(output_localization_diagnostics .and. my_task_id() == 0) call close_file(localization_unit)

! get rid of mpi window
call free_mean_window()

! deallocate space
deallocate(close_obs_dist,      &
           my_obs_indx,         &
           my_obs_kind,         &
           my_obs_type,         &
           close_obs_ind,       &
           vstatus,             &
           my_obs_loc)

deallocate(close_state_dist,      &
           my_state_indx,         &
           close_state_ind,       &
           my_state_kind,         &
           my_state_loc)

! end dealloc

end subroutine filter_assim

!-------------------------------------------------------------

subroutine obs_increment(hx_c, ens_size, obs, obs_var, obs_inc,            &
   inflate, my_cov_inflate, my_cov_inflate_sd,  net_a,                     &
   iter, max_iter, initial_alpha, Ni, inner_index, inner_p, inner_c, hx_p, &
   inner_increment, norm_inner_increment, base_obs_type, base_obs_loc )

! Given the ensemble prior for an observation, the observation, and
! the observation error variance, computes increments and adjusts
! observation space inflation values

integer,                     intent(in)    :: ens_size
real(r8),                    intent(in)    :: hx_c(ens_size), obs, obs_var
real(r8),                    intent(out)   :: obs_inc(ens_size)
type(adaptive_inflate_type), intent(inout) :: inflate
real(r8),                    intent(inout) :: my_cov_inflate, my_cov_inflate_sd
real(r8),                    intent(out)   :: net_a

real(r8) :: ens(ens_size), inflate_inc(ens_size)
real(r8) :: prior_mean, prior_var, new_val(ens_size)
integer  :: i, ens_index(ens_size), new_index(ens_size)

real(r8) :: rel_weights(ens_size)

! CCHU
! PFF-DART new input and output
! the size of inner_p (prior) and inner_c (current) should be (# ens_size) * (# inner domain var):
integer,     intent(in), optional :: iter, max_iter
real(r8), intent(inout), optional :: initial_alpha  ! adaptive kernel width
integer,     intent(in)           :: Ni
integer(i8), intent(in), optional :: inner_index(Ni)
real(r8),    intent(in), optional :: inner_p(ens_size*Ni), inner_c(ens_size*Ni)
real(r8),    intent(in), optional :: hx_p(ens_size) ! prior H(x); not at this iteration
real(r8),   intent(out), optional :: inner_increment(ens_size, Ni)
real(r8),   intent(out), optional :: norm_inner_increment(Ni)

integer,     intent(in), optional :: base_obs_type ! the obs type
type(location_type), intent(in), optional :: base_obs_loc

! new dummy variables for PFF-DART
real(r8) :: inner_pmatrix(ens_size, Ni),inner_cmatrix(ens_size, Ni)
integer  :: j


! CCHU
! first reshape the input vector to the matrix form
do i=1,Ni
    do j=1,ens_size
        inner_pmatrix(j,i) = inner_p(ens_size*(i-1)+j)
        inner_cmatrix(j,i) = inner_c(ens_size*(i-1)+j)
    enddo
enddo

! Copy the input ensemble to something that can be modified
ens = hx_c

! Null value of net spread change factor is 1.0
net_a = 0.0_r8

! Compute prior variance and mean from sample
prior_mean = sum(ens) / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

! If observation space inflation is being done, compute the initial
! increments and update the inflation factor and its standard deviation
! as needed. my_cov_inflate < 0 means don't do any of this.
if(do_obs_inflate(inflate)) then
   ! If my_cov_inflate_sd is <= 0, just retain current my_cov_inflate setting
   if(my_cov_inflate_sd > 0.0_r8) &
      ! Gamma set to 1.0 because no distance for observation space
      call update_inflation(inflate, my_cov_inflate, my_cov_inflate_sd, prior_mean, &
         prior_var, ens_size, obs, obs_var, gamma_corr = 1.0_r8)

   ! Now inflate the ensemble and compute a preliminary inflation increment
   call inflate_ens(inflate, ens, prior_mean, my_cov_inflate, prior_var)
   ! Keep the increment due to inflation alone
   inflate_inc = ens - hx_c

   ! Need to recompute variance if non-deterministic inflation (mean is unchanged)
   if(.not. deterministic_inflate(inflate)) &
      prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)
endif

! If obs_var == 0, delta function.  The mean becomes obs value with no spread.
! If prior_var == 0, obs has no effect.  The increments are 0.
! If both obs_var and prior_var == 0 there is no right thing to do, so Stop.
if ((obs_var == 0.0_r8) .and. (prior_var == 0.0_r8)) then

   ! fail if both obs variance and prior spreads are 0.
   write(msgstring,  *) 'Observation value is ', obs, ' ensemble mean value is ', prior_mean
   write(msgstring2, *) 'The observation has 0.0 error variance, and the ensemble members have 0.0 spread.'
   write(msgstring3, *) 'These require inconsistent actions and the algorithm cannot continue.'
   call error_handler(E_ERR, 'obs_increment', msgstring, &
           source, text2=msgstring2, text3=msgstring3)

else if (obs_var == 0.0_r8) then

   ! new mean is obs value, so increments are differences between obs
   ! value and current value.  after applying obs, all state will equal obs.
   obs_inc(:) = obs - ens

else if (prior_var == 0.0_r8) then

   ! if all state values are the same, nothing changes.
   obs_inc(:) = 0.0_r8
   inner_increment = 0.0_r8
   norm_inner_increment  = 0.0_r8

else
! CCWU (test PFF) SO...
filter_kind = 9
   ! Call the appropriate filter option to compute increments for ensemble
   ! note that at this point we've taken care of the cases where either the
   ! obs_var or the prior_var is 0, so the individual routines no longer need
   ! to have code to test for those cases.
   if(filter_kind == 1) then
      call obs_increment_eakf(ens, ens_size, prior_mean, prior_var, &
         obs, obs_var, obs_inc, net_a)
   else if(filter_kind == 2) then
      call obs_increment_enkf(ens, ens_size, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == 3) then
      call obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
   else if(filter_kind == 4) then
      call obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
   else if(filter_kind == 5) then
      call obs_increment_ran_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == 6) then
      call obs_increment_det_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == 7) then
      call obs_increment_boxcar(ens, ens_size, obs, obs_var, obs_inc, rel_weights)
   else if(filter_kind == 8) then
      call obs_increment_rank_histogram(ens, ens_size, prior_var, obs, obs_var, obs_inc)
   else if(filter_kind == 9) then
      call obs_increment_pff(iter, max_iter, ens_size, Ni,              &
                             hx_p, hx_c, inner_pmatrix, inner_cmatrix,  &
                             inner_index,                               &
                             obs, obs_var, base_obs_type, base_obs_loc, &
                             initial_alpha,                             &
                             inner_increment, norm_inner_increment )

   else
      call error_handler(E_ERR,'obs_increment', &
              'Illegal value of filter_kind in assim_tools namelist [1-8 OK]', source)
   endif
endif


! Add in the extra increments if doing observation space covariance inflation
if(do_obs_inflate(inflate)) obs_inc = obs_inc + inflate_inc

! To minimize regression errors, may want to sort to minimize increments
! This makes sense for any of the non-deterministic algorithms
! By doing it here, can take care of both standard non-deterministic updates
! plus non-deterministic obs space covariance inflation. This is expensive, so
! don't use it if it's not needed.
if (sort_obs_inc) then
   new_val = hx_c + obs_inc
   ! Sorting to make increments as small as possible
   call index_sort(hx_c, ens_index, ens_size)
   call index_sort(new_val, new_index, ens_size)
   do i = 1, ens_size
      obs_inc(ens_index(i)) = new_val(new_index(i)) - hx_c(ens_index(i))
   end do
endif

! Get the net change in spread if obs space inflation was used
if(do_obs_inflate(inflate)) net_a = net_a * sqrt(my_cov_inflate)
end subroutine obs_increment



! CCHU: subroutine PFF here:
subroutine obs_increment_pff(iter, max_iter, ens_size, Ni,              &
                             hx_p, hx_c, inner_pmatrix, inner_cmatrix,  &
                             inner_index,                               &
                             obs, obs_var, base_obs_type, base_obs_loc, &
                             initial_alpha,                             &
                             inner_increment, norm_inner_increment )
!========================================================================
!
! PFF version of "inner domain" increment

!Use lapack_interfaces, Only: dgesvd


integer,      intent(in)  :: iter, max_iter, ens_size, Ni 
real(r8),     intent(in)  :: hx_p(ens_size), hx_c(ens_size)
real(r8),     intent(in)  :: inner_pmatrix(ens_size, Ni), inner_cmatrix(ens_size, Ni)
integer(i8),  intent(in)  :: inner_index(Ni)      ! the iteration index

real(r8),           intent(in)  :: obs, obs_var
integer,             intent(in) :: base_obs_type  ! the observation type
type(location_type), intent(in) :: base_obs_loc

real(r8),   intent(inout) :: initial_alpha
real(r8),     intent(out) :: inner_increment(ens_size, Ni)
real(r8),     intent(out) :: norm_inner_increment(Ni)

! the following are temporary variables
real(r8) :: norm_nodim_inner_increment(Ni)
real(r8) :: prior_cov(Ni,Ni), prior_cov_inv(Ni,Ni), input_inverse(Ni,Ni)
!real(r8) :: kernel(ens_size, ens_size, Ni) ! matrix-valued kernel
real(r8) :: kernel(ens_size, ens_size) ! scalar-valued kernel
real(r8) :: argument
real(r8) :: dx(Ni)
real(r8) :: kernel_width(Ni), prior_mean(Ni)
real(r8) :: obs_mean, inner_mean(Ni), HT(ens_size, Ni), BHT(ens_size, Ni), HBHT(ens_size)
real(r8) :: post_pdf(ens_size, Ni), like_pdf(ens_size, Ni), prir_pdf(ens_size, Ni)
real(r8) :: ker_width(Ni)
real(r8) :: ker_alpha
real(r8) :: grad_ker(ens_size, ens_size, Ni)
real(r8) :: eps_assim, eps_err
!real(r8) :: eps_type
real(r8) :: hx_var, hx_mean, innov_std ! variance and mean of ens
real(r8) :: dd(Ni,Ni), cov_factor(Ni,Ni)
real(r8) :: dd2(Ni), cov_factor2(Ni)
real(r8) :: input_x(ens_size, Ni)
real(r8) :: inner_inc_T(Ni, ens_size) ! temporarily used for Binv multiplication
real(r8) :: particle_dis(ens_size, ens_size)
!real(r8) :: min_kernel_value
integer  :: i,j,k

type(location_type):: inner_loc(Ni)

! adaptive kernel width
integer  :: inner_var_type
real(r8) :: prior_hx_mean, prior_hx_var, post_hx_var, post_hx_mean, current_hx_mean, current_hx_var
real(r8) :: inf_obs_var, inf_post_hx_mean
real(r8) :: var_ratio, obs_space_inc(ens_size)
real(r8) :: cov_xy(Ni)


! Although not necessarily used, the following info can be handy
! if assuming all Gaussian solution, the posterior variance/mean will be:
prior_hx_mean = sum(hx_p)/ens_size
prior_hx_var  = sum((hx_p-prior_hx_mean)**2)/(ens_size-1)
post_hx_var   = obs_var*prior_hx_var/(obs_var + prior_hx_var)
post_hx_mean  = prior_hx_var/(prior_hx_var + obs_var)*obs + &
                 obs_var    /(prior_hx_var + obs_var)*prior_hx_mean

current_hx_mean = sum(hx_c)/ens_size
current_hx_var  = sum((hx_c-current_hx_mean)**2)/(ens_size-1)

! special case: only for the second iteration
! examine if the learning rate at the 1st iteration is too large:

!if (iter .eq. 2) then
!   if ( (current_hx_mean-prior_hx_mean)*(current_hx_mean-obs) .gt. 0 ) then
!      inner_increment = 0.0_r8
!      norm_inner_increment = 0.0_r8
!      if (my_task_id()==0) print*, '1st iteration learning rate too large, return...'
!      return

!   endif
!endif



! Get the location informaiton for inner domain:
do i=1,Ni
   call get_state_meta_data(inner_index(i), inner_loc(i),inner_var_type)
   !IF (my_task_id()==0) print*, 'inner index = ', inner_index
   !if (my_task_id()==0) print*, 'inner domain var ',i,' var type =', inner_var_type
   !if (my_task_id()==0 .and. iter ==1 ) print*, '  location lat  = ', inner_loc(i)%lat/3.1415*180., ' lon =',inner_loc(i)%lon/3.1415*180.
   !if (my_task_id()==0) print*, '          height = ', inner_loc(i)%vloc,'( vert = ', inner_loc(i)%which_vert,' )'
enddo

! caluclate the prior mean, prior covariance matrix
prior_mean = sum(inner_pmatrix,dim=1)/ens_size

! the "unlocalized" prior covariance matrix
! Note that in this version the localization is done directly on the final increment
do i=1,Ni
   do j=1,Ni
      if (j.ge.i) then
         prior_cov(j,i) = dot_product(inner_pmatrix(:,j)-prior_mean(j), &
                                      inner_pmatrix(:,i)-prior_mean(i))/(ens_size-1)
      else
         prior_cov(j,i) = prior_cov(i,j)
      endif
   enddo
enddo

! calculate the inverse of prior covariance matrix:
! CAUTIOUS!! by calling the subroutine "inverse" might actually change the value of
! prior_cov. Need to check the code for why?

input_inverse = prior_cov ! input_inverse is a dummy variable that passes prior_cov to the inverse

! use SVD to calculate the inverse of prior covariance (of inner domain)
call svd_pseudo_inverse(input_inverse,prior_cov_inv,Ni,Ni,min_eig_ratio)

! ====== estimation of the adjoint of the observation operator ======

input_x = inner_cmatrix

if ( obs_adj_kind == 0 ) then
   ! ----- METHOD 1: linear regression: -----
   call HT_regress(HT, input_x, hx_c, ens_size, Ni, min_eig_ratio)

elseif ( obs_adj_kind == 1 ) then
   ! the analytical adjoint for exponential H(x)
   ! CAUTION: need to manually check if the adjoint is correct!!
   do i=1,ens_size
      do j=1,4
   !      HT(i,j) = 0.25* hx_c(i)/400
         HT(i,j) = 0.25* hx_c(i)/100
      enddo
   enddo

   !   do j=5,8
   !      HT(i,j) = 0.25*(inner_cmatrix(i,5)+inner_cmatrix(i,6)+ &
   !                      inner_cmatrix(i,7)+inner_cmatrix(i,8))/sqrt(ens(i))
   !      HT(i,j) =  0.5*(inner_cmatrix(i,5)+inner_cmatrix(i,6)+ &
   !                      inner_cmatrix(i,7)+inner_cmatrix(i,8))
   !   enddo
   !enddo

elseif ( obs_adj_kind == 2 ) then
   ! ----- METHOD 2: kernel approx: -----
   call HT_kernel(HT, input_x, hx_c, ens_size, Ni, 0.05*1.0_r8)

endif

if ((iter==1).and.(my_task_id()==0)) then
   !print*, 'x (1st)=', inner_cmatrix(:,1)
   !print*, 'x (2nd)=', inner_cmatrix(:,2)
   !print*, 'x (3rd)=', inner_cmatrix(:,3)
   !print*, 'x (4th)=', inner_cmatrix(:,4)
   !print*, 'H(x) =', hx_c(1:ens_size)
   print*,'HT=', HT(1,:)
endif

! ====== END OF estimate the adjoint of the observation operator ======


! in the following, calculate the gradient of the posterior pdf
! Decompose the B* (log posterior) = B*log likelihood (like_pdf) + B*log prior (prir_pdf)
do i=1,ens_size
   like_pdf(i,:) = matmul( prior_cov, HT(i,:)*(obs-hx_c(i))/obs_var)
   prir_pdf(i,:) = -(inner_cmatrix(i,:)-prior_mean(:))
enddo

post_pdf = like_pdf + prir_pdf


! ====== Adaptive kernel width algorithm ======
! 2022/11/11 adaptive kernel width:
do i=1,ens_size
   do j=1,ens_size
      dx                = inner_cmatrix(i,:) - inner_cmatrix(j,:)
      particle_dis(i,j) = dot_product(dx,matmul(prior_cov_inv, dx))
   enddo
enddo

if ( adaptive_ker_io ) then
   if ( iter.eq.1 ) then
   !   min_kernel_value = 0.1 ! set the minimum kernel value
      ker_alpha = -maxval(particle_dis)/log(min_kernel_value)
      initial_alpha = ker_alpha
   elseif ( current_hx_var .le. post_hx_var ) then
      if (my_task_id()==0) print*, 'caution: kernel values are set to 1'
      ker_alpha = -maxval(particle_dis)/log(0.999999)
      !ker_alpha = initial_alpha
   else
      ker_alpha = initial_alpha
   endif

else
   ker_alpha = fixed_ker_alpha

endif

if ((my_task_id().eq.0).and.(iter.eq.1)) then
   print*, 'kernel alpha =',ker_alpha
endif

! ====== END OF Adaptive kernel width algorithm ======


! evaluate the pairwise kernels and their gradients:
do i=1,ens_size
   do j=1,ens_size
      dx              = inner_cmatrix(i,:) - inner_cmatrix(j,:)
      kernel(i,j)     = exp(-particle_dis(i,j)/ker_alpha)
      grad_ker(i,j,:) = 2*dx/ker_alpha*kernel(i,j)

      !grad_ker(i,j,:) = 2*matmul(prior_cov_inv,dx)/ker_alpha*kernel(i,j)
   enddo
enddo

if ((my_task_id().eq.0).and.(iter.eq.1)) then
   print*, 'minimum kernel = ', minval(kernel)
endif


! ====== Adaptive learning rate for eps_assim ======

! eps_assim = eps_type * eps_err (does NOT depend on iteration anymore)

!if (base_obs_type.le.0) then
!   eps_type = 0.1 ! for identity obs
!    eps_type = 0.1
!endif

! 2022/11/21: fix eps_type for now (the value is via input.nml)
!             in the future, better fix an optimal eps_type for each obs type
!select case (base_obs_type)
!  case(4)
!    !print*,'obs type = surface pressure'
!    eps_type = 0.1
!  case default
!    !print*,'unknown obs type'
!    eps_type = 0.2
!end select

!hx_mean = sum(ens)/(1.0_r8*ens_size)
!hx_var  = sum((ens-hx_mean)**2)/(1.0_r8*(ens_size-1))
!innov_std = sqrt(sum( (hx_c-obs)**2 )/(1.0_r8*(ens_size-1)))

! maybe no need eps_err?

! min is to avoid too large learning rate for small prior variance
!eps_err = min(obs_var/hx_var,1.0_r8)
!eps_err = min(obs_var/innov_std, 5.0_r8)

!eps_assim = eps_type

!eps_assim = min(eps_type/ (sum(abs(HT(1,:)))/Ni)*1e-3, eps_type)

! ====== END OF Adaptive learning rate for eps_assim ======


! Calculate the increment (and its norm) for scalar-valued kernel:
do i=1,ens_size
   inner_increment(i,:) = 0.0_r8

   do j=1,ens_size
!      inner_increment(i,:) = inner_increment(i,:) + &
!                                eps_assim*( kernel(i,j)*post_pdf(j,:) + grad_ker(i,j,:) )/(1.0_r8*ens_size)
      inner_increment(i,:) = inner_increment(i,:) + &
                                           ( kernel(i,j)*post_pdf(j,:) + grad_ker(i,j,:) )/(1.0_r8*ens_size)
!      inner_increment(i,:) = inner_increment(i,:) + &
!                                           ( kernel(i,j)*post_pdf(j,:) +  grad_ker(i,j,:) )/sum(kernel(i,:))
   enddo
end do

!Binv_inner_increment = matmul(inner_increment, prior_cov_inv)

! Calculate the norm of inner domain increment for each component:
do i=1,Ni
!   norm_inner_increment(i) = norm2(inner_increment(:,i)/eps_assim)
    norm_inner_increment(i) = norm2(inner_increment(:,i)/sqrt(ens_size*1.0_r8))
!    norm_Binv_inner_increment(i) = norm2(Binv_inner_increment(:,i)/sqrt(ens_size*1.0_r8))
    norm_nodim_inner_increment(i) = &
       norm2(inner_increment(:,i)/sqrt(ens_size*1.0_r8))/sqrt(prior_cov(i,i))
enddo

! Try to make learning rate depend on the actual magnitude of particle flow
!eps_assim = min( 1.0_r8/(sum(norm_inner_increment)/Ni)*10, eps_type)
!eps_assim = min(1.0_r8/(sum(norm_Binv_inner_increment)/Ni)*learning_rate_fac, max_learning_rate)
eps_assim = min(1.0_r8/(norm2(norm_nodim_inner_increment/sqrt(Ni*1.0_r8)))*learning_rate_fac, max_learning_rate)

if (my_task_id()==0) print*, 1.0_r8/(norm2(norm_nodim_inner_increment/sqrt(Ni*1.0_r8)))*learning_rate_fac

inner_increment = eps_assim*inner_increment


! ======== To speed up the iteration, try EAKF update in the first iteration ========

! first step EAKF/EnKF update
if ( eakffg_io .and. (iter ==1) ) then

   inf_obs_var = eakffg_inf*obs_var ! first-guess EAKF uses inflated obs error var

   ! re-evaluate EAKF solution if use inflated EAKF first-guess:
   inf_post_hx_mean  = prior_hx_var/(prior_hx_var + inf_obs_var)*obs + &
                        inf_obs_var/(prior_hx_var + inf_obs_var)*prior_hx_mean

   ! Compute the new mean
   var_ratio = inf_obs_var / (prior_hx_var + inf_obs_var)

   ! Compute sd ratio and shift ensemble
   obs_space_inc = sqrt(var_ratio) * (hx_p - prior_hx_mean) + inf_post_hx_mean - hx_p

   ! update the inner domain
   ! Note: the exact EAKF for the inner domain:

   ! compute the localization factor for state-obs pair:

   ! original EAKF =========
   do i=1,Ni
      !dd2(i)         = get_dist(inner_loc(i), base_obs_loc)
      !cov_factor2(i) = comp_cov_factor(dd2(i), cutoff)
      !cov_xy(i)      = cov_factor2(i)*dot_product(inner_pmatrix(:,i)-prior_mean(i), &
      !                                                          hx_p-prior_hx_mean)/(ens_size-1)
      cov_xy(i)      = dot_product(inner_pmatrix(:,i)-prior_mean(i), &
                                                                hx_p-prior_hx_mean)/(ens_size-1)
   enddo

   ! linear regression:
   do i=1,ens_size
      inner_increment(i,:) = cov_xy(:)/prior_hx_var*obs_space_inc(i)
   enddo


   ! exact-adjoint EAKF ========
   !do i=1,ens_size
   !   BHT(i,:) = matmul(prior_cov, HT(i,:))
   !   HBHT(i)  = dot_product(HT(i,:), BHT(i,:))
   !   inner_increment(i,:) = BHT(i,:)/HBHT(i)*obs_space_inc(i)
   !enddo

   !IF (my_task_id()==0) print*, inner_increment(:,1)

endif

! print all the diagnostics

!if (iter .eq. max_iter) then
if (my_task_id()==0) print*, ' '
if (my_task_id()==0) print*, '      obs value =', obs
if (my_task_id()==0) print*, '      eps_assim = ',eps_assim
!if (my_task_id()==0) print*, '      current std = ', sqrt(current_hx_var), 'posterior std =', sqrt(post_hx_var)
if (my_task_id()==0) print*, '      current mean = ', current_hx_mean, 'posterior mean =', post_hx_mean
if (my_task_id()==0) print*, ' '

!endif

end subroutine obs_increment_pff





subroutine HT_kernel(HT, X, Hx, Np, Ni, alpha)
! estimate the adjoint of the observation operator based on kernel method
! (specifically, diagonal B matrix for kernel, and linear samples added)
! input:  X   (ensemble of inner domain variables, size: #particle, #Ni)
!         Hx  (ensemble in obs space, size: #particle)
!         alpha (a parameter for the kernel width)
!         Np = #particle
!         Ni = #Ni
! output: HT (the ensemble mean adjoint, size: #particle, #Ni)
! latest update: cchu 2022/07/06

   implicit none
   real(r8),    intent(in) :: X(Np, Ni), Hx(Np), alpha
   integer,     intent(in) :: Np, Ni
   real(r8),   intent(out) :: HT(Np,Ni)
   real(r8) :: varB(Ni), prior_mean(Ni)
   real(r8) :: ker_width, Np_r, Ni_r
   real(r8) :: kernel(Np,Np)
   real(r8) :: grad_kernel(Ni,Np,Np)
   integer  :: i,j,k
   real(r8) :: sum_K
   real(r8) :: sum_H_grad_K(Ni)

Np_r = Np*1.0_r8
Ni_r = Ni*1.0_r8
ker_width = alpha*Np/Ni
!write(*,*) ker_width

!write(*,*) 'X=',X(1,1)
!write(*,*) 'X=',X(1,2)
!write(*,*) 'hx=',Hx(1)

! variance of each state variable:
! caluclate the prior mean, prior covariance matrix
prior_mean = 0.0_r8
prior_mean = sum(X,dim=1)/Np_r

! prior covariance matrix (need to see how to localize!)
do i=1,Ni
   varB(i) = sum( (X(:,i)-prior_mean(i))**2 )/(Np_r-1)
enddo

! calculate the pairwise kernel matrix
do i=1,Np
   do j=1,Np
      if (j.ge.i) then
         kernel(i,j) = exp( - ker_width*sum( (X(i,:)-X(j,:))**2/varB ) )
      else
         kernel(i,j) = kernel(j,i)
      endif
   enddo
enddo

!write(*,*) 'X(1,:) =',X(1,:)
!write(*,*) 'X(3,:) =',X(3,:)
!write(*,*) 'varB = ',varB
!write(*,*) 'kernel =',kernel(1,:)
!write(*,*) 'kernel =',kernel(1,2)

! calculate the gradient of the kernel:
do k=1,Ni
   do i=1,Np
      do j=1,Np
         if (j.ge.i) then
            grad_kernel(k,i,j) = -2*ker_width*(X(i,k)-X(j,k))/varB(k)*kernel(i,j)
         else
            grad_kernel(k,i,j) = - grad_kernel(k,j,i)
         endif
      enddo
   enddo
enddo

!write(*,*) 'grad ker (1,2) = ', grad_kernel(:,1,2)
!write(*,*) 'grad ker (2,1) = ', grad_kernel(:,2,1)
!write(*,*) 'grad ker (3,3) = ', grad_kernel(:,3,3)

! estimation of the gradient:
do i=1,Np
   sum_K = 0
   sum_H_grad_K = 0

   do j=1,Np
      sum_K        = sum_K + kernel(i,j)
      sum_H_grad_K = sum_H_grad_K + (Hx(j)-Hx(i))*grad_kernel(:,i,j)
   enddo
   HT(i,:) = sum_H_grad_K/sum_K
   !write(*,*) 'sum K =',sum_K
   !write(*,*) 'sum H grad K =', sum_H_grad_K

enddo
!write(*,*) HT(1,:) ! print the adjoint of the first ensemble member

end subroutine HT_kernel

subroutine HT_regress(HT, X, Hx, Np, Ni, inv_max_cond_num)
! estimate the adjoint of the observation operator, based on the ensemble
! input:  X   (ensemble of inner domain variables, size: #particle, #Ni)
!         Hx  (ensemble in obs space, size: #particle)
!         inv_cond_num (the inverse of the maximum condition number for the matrix). 
!                      Note that the condition number of a matrix A, cond(A) is the 
!                      ratio of the max eigenvalue of A to the min eigenvalue of A
!                      By setting a "maximum" condition number, we enforce those 
!                      "too small" eigenvalues to 0 when inverting the matrix
!         Np = #particle
!         Ni = #Ni
! output: HT (the ensemble mean adjoint, size: #particle, #Ni)

   implicit none
   real(r8),    intent(in) :: X(Np, Ni), Hx(Np), inv_max_cond_num
   integer,     intent(in) :: Np, Ni
   real(r8),   intent(out) :: HT(Np,Ni)
   real(r8) :: pX(Ni,Np), HT_ensmean(Ni), Xp(Np,Ni), Hxp(Np)
   integer  :: i

! first calculate the perturbation matrices:
do i=1,Ni
   Xp(:,i) = X(:,i) - sum(X(:,i))/Np
enddo
Hxp = Hx - sum(Hx)/Np

! then obtain the pseudo inverse of Xp:
call svd_pseudo_inverse(Xp,pX,Np,Ni,inv_max_cond_num)
HT_ensmean = matmul(pX, Hxp)

do i=1,Np
   HT(i,:) = HT_ensmean
enddo

end subroutine HT_regress


subroutine svd_pseudo_inverse(A,pA,m,n,inv_max_cond_num)
      ! A is the input matrix, which has dimension mxn
      ! cond_num is the largest condition number that A can be
      ! if cond(A)>cond_num, this subroutine makes cond(A)=cond_num
      ! before calculating the (pseudo)inverse
      ! pA is the pseudo-inverse of A (dimension: nxm)
      implicit none

      real(r8), intent(in)   :: A(m,n)
      real(r8), intent(in)   :: inv_max_cond_num
      integer, intent(in)   :: m,n
      real(r8), intent(out)  :: pA(n,m)
      real(r8) :: U(m,m), UT(m,m), V(n,n), VT(n,n), S_inv(n,m)
      real(r8), allocatable:: S(:), work(:)
      integer  :: info_svd, LWORK, len_s, i, j
      character:: stringA*50, stringS*50
!      real(8) :: test1(m,n), test2(m,n), SS(m,n)

!      WRITE(*,*) 'A ='
!      WRITE(stringA, '( "(" ,I4, "f10.4)" )' ) n ! # of column A
      !write(*,*) stringA
!      WRITE(*,stringA) ((A(I,J),J=1,n),I=1,m)

      len_s = min(m,n)
      allocate(S(len_s))

      ! calculate the svd of A:

      LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
      allocate(work(LWORK))

      call dgesvd('A','A',m,n,A,m,S,U,m,VT,n,WORK,LWORK,info_svd)
     ! calculate the inverse of A (pA) below:
      V  = transpose(VT)
      UT = transpose(U)

      S_inv = 0
      do i=1,len_s
         if ( S(i).ge.inv_max_cond_num*S(1) ) then
             S_inv(i,i) = 1./S(i)
         endif
      enddo

!      write(*,*) 'pinv(S) = '
!      WRITE(stringS, '( "(" ,I4, "f10.4)" )' ) m ! # of column S
!      write(*,stringS) ((S_inv(i,j),j=1,m),i=1,n)

      pA = matmul(V,matmul(S_inv,UT))

!      write(*,*) 'pinv(A) = '
!      write(*,stringS) ((pA(i,j),j=1,m),i=1,n)

end subroutine svd_pseudo_inverse


subroutine inverse(a_in,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer , intent(in) :: n
real(r8), intent(in) :: a_in(n,n)
real(r8), intent(out):: c(n,n)
real(r8):: a(n,n)
real(r8):: L(n,n), U(n,n), b(n), d(n), x(n)
real(r8):: coeff
integer:: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0_r8
U=0.0_r8
b=0.0_r8

a = a_in

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0_r8
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0_r8
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0_r8
end do

end subroutine inverse





subroutine obs_increment_eakf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc, a)
!========================================================================
!
! EAKF version of obs increment

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_mean, prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(out) :: a

real(r8) :: new_mean, var_ratio

! Compute the new mean
var_ratio = obs_var / (prior_var + obs_var)
new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)

! Compute sd ratio and shift ensemble
a = sqrt(var_ratio)
obs_inc = a * (ens - prior_mean) + new_mean - ens

end subroutine obs_increment_eakf


subroutine obs_increment_ran_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
!========================================================================
!
! Forms a random sample of the Gaussian from the update equations.
! This is very close to what a true 'ENSEMBLE' Kalman Filter would
! look like. Note that outliers, multimodality, etc., get tossed.

integer,   intent(in)  :: ens_size
real(r8),  intent(in)  :: prior_mean, prior_var
real(r8),  intent(in)  :: ens(ens_size), obs, obs_var
real(r8),  intent(out) :: obs_inc(ens_size)

real(r8) :: new_mean, var_ratio
real(r8) :: temp_mean, temp_var, new_ens(ens_size), new_var
integer  :: i

var_ratio = obs_var / (prior_var + obs_var)
new_var = var_ratio * prior_var
new_mean  = var_ratio * (prior_mean  + prior_var*obs / obs_var)

! This will reproduce exactly for multiple runs with the same task count,
! but WILL NOT reproduce for a different number of MPI tasks.
! To make it independent of the number of MPI tasks, it would need to
! use the global ensemble number or something else that remains constant
! as the processor count changes.  this is not currently an argument to
! this function and so we are not trying to make it task-count invariant.

! Form a random sample from the updated distribution
! Then adjust the mean (what about adjusting the variance?)!
! Definitely need to sort with this; sort is done in main obs_increment
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq, my_task_id() + 1)
   first_inc_ran_call = .false.
endif

do i = 1, ens_size
   new_ens(i) = random_gaussian(inc_ran_seq, new_mean, sqrt(prior_var*var_ratio))
end do

! Adjust the mean of the new ensemble
temp_mean = sum(new_ens) / ens_size
new_ens(:) = new_ens(:) - temp_mean + new_mean

! Compute prior variance and mean from sample
temp_var  = sum((new_ens - new_mean)**2) / (ens_size - 1)
! Adjust the variance, also
new_ens = (new_ens - new_mean) * sqrt(new_var / temp_var) + new_mean

! Get the increments
obs_inc = new_ens - ens

end subroutine obs_increment_ran_kf



subroutine obs_increment_det_kf(ens, ens_size, prior_mean, prior_var, obs, obs_var, obs_inc)
!========================================================================
!
! Does a deterministic ensemble layout for the updated Gaussian.
! Note that all outliers, multimodal behavior, etc. get tossed.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: prior_mean, prior_var
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: new_mean, var_ratio, temp_var, new_ens(ens_size), new_var
integer :: i

var_ratio = obs_var / (prior_var + obs_var)
new_var = var_ratio * prior_var
new_mean = var_ratio * (prior_mean  + prior_var*obs / obs_var)

! Want a symmetric distribution with kurtosis 3 and variance new_var and mean new_mean
if(ens_size /= 20) then
   write(*, *) 'EXPERIMENTAL version obs_increment_det_kf only works for ens_size 20 now'
   stop
endif

! This has kurtosis of 3.0, verify again from initial uniform
!new_ens(1) = -2.146750_r8
!new_ens(2) = -1.601447_r8
!new_ens(3) = -1.151582_r8
!new_ens(4) = -0.7898650_r8
!new_ens(5) = -0.5086292_r8
!new_ens(6) = -0.2997678_r8
!new_ens(7) = -0.1546035_r8
!new_ens(8) = -6.371084E-02_r8
!new_ens(9) = -1.658448E-02_r8
!new_ens(10) = -9.175255E-04_r8

! This has kurtosis of 3.0, verify again from initial inverse gaussian
!new_ens(1) = -2.188401_r8
!new_ens(2) = -1.502174_r8
!new_ens(3) = -1.094422_r8
!new_ens(4) = -0.8052422_r8
!new_ens(5) = -0.5840152_r8
!new_ens(6) = -0.4084518_r8
!new_ens(7) = -0.2672727_r8
!new_ens(8) = -0.1547534_r8
!new_ens(9) = -6.894587E-02_r8
!new_ens(10) = -1.243549E-02_r8

! This has kurtosis of 2.0, verify again
new_ens(1) = -1.789296_r8
new_ens(2) = -1.523611_r8
new_ens(3) = -1.271505_r8
new_ens(4) = -1.033960_r8
new_ens(5) = -0.8121864_r8
new_ens(6) = -0.6077276_r8
new_ens(7) = -0.4226459_r8
new_ens(8) = -0.2598947_r8
new_ens(9) = -0.1242189_r8
new_ens(10) = -2.539018E-02_r8

! This has kurtosis of 1.7, verify again
!new_ens(1) = -1.648638_r8
!new_ens(2) = -1.459415_r8
!new_ens(3) = -1.272322_r8
!new_ens(4) = -1.087619_r8
!new_ens(5) = -0.9056374_r8
!new_ens(6) = -0.7268229_r8
!new_ens(7) = -0.5518176_r8
!new_ens(8) = -0.3816142_r8
!new_ens(9) = -0.2179997_r8
!new_ens(10) = -6.538583E-02_r8
do i = 11, 20
   new_ens(i) = -1.0_r8 * new_ens(20 + 1 - i)
end do

! Right now, this ensemble has mean 0 and some variance
! Compute prior variance and mean from sample
temp_var  = sum((new_ens)**2) / (ens_size - 1)

! Adjust the variance of this ensemble to match requirements and add in the mean
new_ens = new_ens * sqrt(new_var / temp_var) + new_mean

! Get the increments
obs_inc = new_ens - ens

end subroutine obs_increment_det_kf




subroutine obs_increment_particle(ens, ens_size, obs, obs_var, obs_inc)
!------------------------------------------------------------------------
!
! A observation space only particle filter implementation for a
! two step sequential update filter. Second version, 2 October, 2003.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: weight(ens_size), rel_weight(ens_size), cum_weight(0:ens_size)
real(r8) :: base, frac, new_val(ens_size), weight_sum
integer  :: i, j, indx(ens_size)

! Begin by computing a weight for each of the prior ensemble members
do i = 1, ens_size
   weight(i) = exp(-1.0_r8 * (ens(i) - obs)**2 / (2.0_r8 * obs_var))
end do

! Compute relative weight for each ensemble member
weight_sum = sum(weight)
do i = 1, ens_size
   rel_weight(i) = weight(i) / weight_sum
end do

! Compute cumulative weights at boundaries
cum_weight(0) = 0.0_r8
do i = 1, ens_size
   cum_weight(i) = cum_weight(i - 1) + rel_weight(i)
!   write(*,'(1x,i3,3(e10.4,1x))') i, weight(i), rel_weight(i), cum_weight(i)
end do
! Fix up for round-off error if any
cum_weight(ens_size) = 1.0_r8

! Do a deterministic implementation: just divide interval into ens_size parts and see
! which interval this is in (careful to offset; not start at 0)
base = 1.0_r8 / (ens_size * 2.0_r8)

do i = 1, ens_size

   frac = base + (i - 1.0_r8) / ens_size

   ! Now search in the cumulative range to see where this frac falls
   ! Can make this search more efficient by limiting base
   do j = 1, ens_size
      if(cum_weight(j - 1) < frac .and. frac < cum_weight(j)) then
         indx(i) = j
!         write(*, *) i, frac, 'gets index ', j
         goto 111
      end if
   end do

111 continue

end do

! Set the new values for the ensemble members
do i = 1, ens_size
   new_val(i) = ens(indx(i))
!   write(*, *) 'new_val ', i, new_val(i)
end do

! Generate increments
obs_inc = new_val - ens

end subroutine obs_increment_particle



subroutine obs_increment_enkf(ens, ens_size, prior_var, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment_enkf(ens, ens_size, obs, obs_var, obs_inc)
!

! ENKF version of obs increment

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

real(r8) :: obs_var_inv, prior_var_inv, new_var, new_mean(ens_size)
! real(r8) :: sx, s_x2
real(r8) :: temp_mean, temp_obs(ens_size)
integer  :: i

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var
prior_var_inv = 1.0_r8 / prior_var

new_var       = 1.0_r8 / (prior_var_inv + obs_var_inv)

! If this is first time through, need to initialize the random sequence.
! This will reproduce exactly for multiple runs with the same task count,
! but WILL NOT reproduce for a different number of MPI tasks.
! To make it independent of the number of MPI tasks, it would need to
! use the global ensemble number or something else that remains constant
! as the processor count changes.  this is not currently an argument to
! this function and so we are not trying to make it task-count invariant.
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq, my_task_id() + 1)
   first_inc_ran_call = .false.
endif

! Generate perturbed obs
do i = 1, ens_size
    temp_obs(i) = random_gaussian(inc_ran_seq, obs, sqrt(obs_var))
end do

! Move this so that it has original obs mean
temp_mean = sum(temp_obs) / ens_size
temp_obs(:) = temp_obs(:) - temp_mean + obs

! Loop through pairs of priors and obs and compute new mean
do i = 1, ens_size
   new_mean(i) = new_var * (prior_var_inv * ens(i) + temp_obs(i) / obs_var)
   obs_inc(i)  = new_mean(i) - ens(i)
end do

! Can also adjust mean (and) variance of final sample; works fine
!sx         = sum(new_mean)
!s_x2       = sum(new_mean * new_mean)
!temp_mean = sx / ens_size
!temp_var  = (s_x2 - sx**2 / ens_size) / (ens_size - 1)
!new_mean = (new_mean - temp_mean) * sqrt(new_var / temp_var) + updated_mean
!obs_inc = new_mean - ens


end subroutine obs_increment_enkf



subroutine obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
!========================================================================
! subroutine obs_increment_kernel(ens, ens_size, obs, obs_var, obs_inc)
!

! Kernel version of obs increment

integer, intent(in)             :: ens_size
real(r8), intent(in)            :: ens(ens_size), obs, obs_var
real(r8), intent(out)           :: obs_inc(ens_size)

real(r8) :: obs_var_inv
real(r8) :: prior_mean, prior_cov_inv, new_cov, prior_cov
real(r8) :: sx
real(r8) :: weight(ens_size), new_mean(ens_size)
real(r8) :: cum_weight, total_weight, cum_frac(ens_size)
real(r8) :: unif, norm, new_member(ens_size)

integer :: i, j, kernel

! Compute mt_rinv_y (obs error normalized by variance)
obs_var_inv = 1.0_r8 / obs_var

! Compute prior mean and covariance
sx         = sum(ens)
prior_mean = sx / ens_size
prior_cov  = sum((ens - prior_mean)**2) / (ens_size - 1)

prior_cov     = prior_cov / 10.0_r8     ! For kernels, scale the prior covariance
prior_cov_inv = 1.0_r8 / prior_cov

! Compute new covariance once for these kernels
new_cov = 1.0_r8 / (prior_cov_inv + obs_var_inv)

! New mean is computed ens_size times as is weight
do i = 1, ens_size
   new_mean(i) = new_cov*(prior_cov_inv * ens(i) + obs / obs_var)
   weight(i) =  2.71828_r8 ** (-0.5_r8 * (ens(i)**2 * prior_cov_inv + &
      obs**2 * obs_var_inv - new_mean(i)**2 / new_cov))
end do

! Compute total weight
total_weight = sum(weight)
cum_weight   = 0.0_r8
do i = 1, ens_size
   cum_weight  = cum_weight + weight(i)
   cum_frac(i) = cum_weight / total_weight
end do

! If this is first time through, need to initialize the random sequence.
! This will reproduce exactly for multiple runs with the same task count,
! but WILL NOT reproduce for a different number of MPI tasks.
! To make it independent of the number of MPI tasks, it would need to
! use the global ensemble number or something else that remains constant
! as the processor count changes.  this is not currently an argument to
! this function and so we are not trying to make it task-count invariant.
if(first_inc_ran_call) then
   call init_random_seq(inc_ran_seq, my_task_id() + 1)
   first_inc_ran_call = .false.
endif

! Generate a uniform random number and a Gaussian for each new member
do i = 1, ens_size
   unif = random_uniform(inc_ran_seq)
   ! Figure out which kernel it's in
   whichk: do j = 1, ens_size
      if(unif < cum_frac(j)) then
         kernel = j
         exit whichk
      end if
   end do whichk

   ! Next calculate a unit normal in this kernel
   norm = random_gaussian(inc_ran_seq, 0.0_r8, sqrt(new_cov))
   ! Now generate the new ensemble member
   new_member(i) = new_mean(kernel) + norm
end do

! Generate the increments
obs_inc = new_member - ens

end subroutine obs_increment_kernel



subroutine update_from_obs_inc(obs, obs_prior_mean, obs_prior_var, obs_inc, &
               state, ens_size, state_inc, reg_coef, net_a_in, correl_out)
!========================================================================

! Does linear regression of a state variable onto an observation and
! computes state variable increments from observation increments

integer,            intent(in)    :: ens_size
real(r8),           intent(in)    :: obs(ens_size), obs_inc(ens_size)
real(r8),           intent(in)    :: obs_prior_mean, obs_prior_var
real(r8),           intent(in)    :: state(ens_size)
real(r8),           intent(out)   :: state_inc(ens_size), reg_coef
real(r8),           intent(in) :: net_a_in
real(r8), optional, intent(inout) :: correl_out

real(r8) :: obs_state_cov, intermed
real(r8) :: restoration_inc(ens_size), state_mean, state_var, correl
real(r8) :: factor, exp_true_correl, mean_factor, net_a


! For efficiency, just compute regression coefficient here unless correl is needed

state_mean = sum(state) / ens_size
obs_state_cov = sum( (state - state_mean) * (obs - obs_prior_mean) ) / (ens_size - 1)

! CCHU: IMPORTANT CHANGE HERE:
! ONLY FOR PFF-DART (the inverse part has already been calculated in the increment):
if (obs_prior_var > 0.0_r8) then
!   reg_coef = obs_state_cov/obs_prior_var
   reg_coef = obs_state_cov
else
   reg_coef = 0.0_r8
endif

! If correl_out is present, need correl for adaptive inflation
! Also needed for file correction below.

! WARNING: we have had several different numerical problems in this
! section, especially with users running in single precision floating point.
! Be very cautious if changing any code in this section, taking into
! account underflow and overflow for 32 bit floats.

if(present(correl_out) .or. sampling_error_correction) then
   if (obs_state_cov == 0.0_r8 .or. obs_prior_var <= 0.0_r8) then
      correl = 0.0_r8
   else
      state_var = sum((state - state_mean)**2) / (ens_size - 1)
      if (state_var <= 0.0_r8) then
         correl = 0.0_r8
      else
         intermed = sqrt(obs_prior_var) * sqrt(state_var)
         if (intermed <= 0.0_r8) then
            correl = 0.0_r8
         else
            correl = obs_state_cov / intermed
         endif
      endif
   endif
   if(correl >  1.0_r8) correl =  1.0_r8
   if(correl < -1.0_r8) correl = -1.0_r8
endif
if(present(correl_out)) correl_out = correl


! Get the expected actual correlation and the regression weight reduction factor
if(sampling_error_correction) then
   call get_correction_from_table(correl, mean_factor, exp_true_correl, ens_size)
   ! Watch out for division by zero; if correl is really small regression is safely 0
   if(abs(correl) > 0.001_r8) then
      reg_coef = reg_coef * (exp_true_correl / correl) * mean_factor
   else
      reg_coef = 0.0_r8
   endif
   correl = exp_true_correl
endif



! Then compute the increment as product of reg_coef and observation space increment
state_inc = reg_coef * obs_inc
!
! FIXME: craig schwartz has a degenerate case involving externally computed
! forward operators in which the obs prior variance is in fact exactly 0.
! adding this test allowed him to continue to  use spread restoration
! without numerical problems.  we don't know if this is sufficient;
! for now we'll leave the original code but it needs to be revisited.
!
! Spread restoration algorithm option.
!if(spread_restoration .and. obs_prior_var > 0.0_r8) then
!

! Spread restoration algorithm option.
if(spread_restoration) then
   ! Don't use this to reduce spread at present (should revisit this line)
   net_a = min(net_a_in, 1.0_r8)

   ! Default restoration increment is 0.0
   restoration_inc = 0.0_r8

   ! Compute the factor by which to inflate
   ! These come from correl_error.f90 in system_simulation and the files ens??_pairs and
   ! ens_pairs_0.5 in work under system_simulation. Assume a linear reduction from 1
   ! as a function of the net_a. Assume that the slope of this reduction is a function of
   ! the reciprocal of the ensemble_size (slope = 0.80 / ens_size). These are empirical
   ! for now. See also README in spread_restoration_paper documentation.
   !!!factor = 1.0_r8 / (1.0_r8 + (net_a - 1.0_r8) * (0.8_r8 / ens_size)) - 1.0_r8
   factor = 1.0_r8 / (1.0_r8 + (net_a - 1.0_r8) / (-2.4711_r8 + 1.6386_r8 * ens_size)) - 1.0_r8
   !!!factor = 1.0_r8 / (1.0_r8 + (net_a**2 - 1.0_r8) * (-0.0111_r8 + .8585_r8 / ens_size)) - 1.0_r8

   ! Variance restoration
   state_mean = sum(state) / ens_size
   restoration_inc = factor * (state - state_mean)
   state_inc = state_inc + restoration_inc
endif

!! NOTE: if requested to be returned, correl_out is set further up in the
!! code, before the sampling error correction, if enabled, is applied.
!! this means it's returning a different larger value than the correl
!! being returned here.  it's used by the adaptive inflation and so the
!! inflation will see a slightly different correlation value.  it isn't
!! clear that this is a bad thing; it means the inflation might be a bit
!! larger than it would otherwise.  before we move any code this would
!! need to be studied to see what the real impact would be.

end subroutine update_from_obs_inc


!------------------------------------------------------------------------

subroutine get_correction_from_table(scorrel, mean_factor, expected_true_correl, ens_size)

real(r8),  intent(in) :: scorrel
real(r8), intent(out) :: mean_factor, expected_true_correl
integer,  intent(in)  :: ens_size

! Uses interpolation to get correction factor into the table

integer             :: low_indx, high_indx
real(r8)            :: correl, fract, low_correl, low_exp_correl, low_alpha
real(r8)            :: high_correl, high_exp_correl, high_alpha

logical, save :: first_time = .true.

if (first_time) then
   call read_sampling_error_correction(ens_size, exp_true_correl, alpha)
   first_time = .false.
endif

! Interpolate to get values of expected correlation and mean_factor
if(scorrel < -1.0_r8) then
   correl = -1.0_r8
   mean_factor = 1.0_r8
else if(scorrel > 1.0_r8) then
   correl = 1.0_r8
   mean_factor = 1.0_r8
else if(scorrel <= -0.995_r8) then
   fract = (scorrel + 1.0_r8) / 0.005_r8
   correl = (exp_true_correl(1) + 1.0_r8) * fract - 1.0_r8
   mean_factor = (alpha(1) - 1.0_r8) * fract + 1.0_r8
else if(scorrel >= 0.995_r8) then
   fract = (scorrel - 0.995_r8) / 0.005_r8
   correl = (1.0_r8 - exp_true_correl(sec_table_size)) * fract + exp_true_correl(sec_table_size)
   mean_factor = (1.0_r8 - alpha(sec_table_size)) * fract + alpha(sec_table_size)
else
   ! given the ifs above, the floor() computation below for low_indx
   ! should always result in a value in the range 1 to 199.  but if this
   ! code is compiled with r8=r4 (single precision reals) it turns out
   ! to be possible to get values a few bits below 0 which results in
   ! a very large negative integer.  the limit tests below ensure the
   ! index stays in a legal range.
   low_indx = floor((scorrel + 0.995_r8) / 0.01_r8 + 1.0_r8)
   if (low_indx <   1) low_indx =   1
   if (low_indx > 199) low_indx = 199
   low_correl = -0.995_r8 + (low_indx - 1) * 0.01_r8
   low_exp_correl = exp_true_correl(low_indx)
   low_alpha = alpha(low_indx)
   high_indx = low_indx + 1
   high_correl = low_correl + 0.01_r8
   high_exp_correl = exp_true_correl(high_indx)
   high_alpha = alpha(high_indx)
   fract = (scorrel - low_correl) / (high_correl - low_correl)
   correl = (high_exp_correl - low_exp_correl) * fract + low_exp_correl
   mean_factor = (high_alpha - low_alpha) * fract + low_alpha
endif

expected_true_correl = correl

! Don't want Monte Carlo interpolation problems to put us outside of a
! ratio between 0 and 1 for expected_true_correl / sample_correl
! If they have different signs, expected should just be 0
if(expected_true_correl * scorrel <= 0.0_r8) then
   expected_true_correl = 0.0_r8
else if(abs(expected_true_correl) > abs(scorrel)) then
   ! If same sign, expected should not be bigger in absolute value
   expected_true_correl = scorrel
endif

end subroutine get_correction_from_table



subroutine obs_increment_boxcar(ens, ens_size, obs, obs_var, obs_inc, rel_weight)
!------------------------------------------------------------------------
!
! An observation space update that uses a set of boxcar kernels plus two
! half-gaussians on the wings to represent the prior distribution. If N is
! the ensemble size, 1/(N+1) of the mass is placed between each ensemble
! member. This is reminiscent of the ranked historgram approach for
! evaluating ensembles. The prior distribution on the wings is
! represented by a half gaussian with mean being the outermost ensemble
! member (left or right) and variance being somewhat arbitrarily chosen
! as half the total ensemble sample variance. A particle
! filter like algorithm is then used for the update. The weight associated
! with each prior ensemble member is computed by evaluating the likelihood.
! For the interior, the domain for each boxcar is divided in half and each
! half is associated with the nearest ensemble member. The updated mass in
! each half box is the product of the prior mass and the ensemble weight.
! In the wings, the observation likelihood gaussian is convolved with the
! prior gaussian to get an updated weighted gaussian that is assumed to
! represent the posterior outside of the outermost ensemble members. The
! updated ensemble members are chosen so that 1/(N+1) of the updated
! mass is between each member and also on the left and right wings. This
! algorithm is able to deal well with outliers, bimodality and other
! non-gaussian behavior in observation space. It could also be modified to
! deal with non-gaussian likelihoods in the future.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)
real(r8), intent(out) :: rel_weight(ens_size)

integer  :: i, e_ind(ens_size), lowest_box, j
real(r8) :: sx, prior_mean, prior_var, prior_var_d2
real(r8) :: var_ratio, new_var, new_sd, umass, left_weight, right_weight
real(r8) :: mass(2*ens_size), weight(ens_size), cumul_mass(0:2*ens_size)
real(r8) :: new_mean_left, new_mean_right, prod_weight_left, prod_weight_right
real(r8) :: new_ens(ens_size), mass_sum, const_term
real(r8) :: x(1:2*ens_size - 1), sort_inc(ens_size)

! The factor a is not defined for this filter for now (could it be???)

! The relative weights could be used for a multi-dimensional particle-type
! update using update_ens_from_weights. There are algorithmic challenges
! with outliers so this is not currently a supported option. For now,
! rel_weight is simply set to 0 and is unused elsewhere.
rel_weight = 0.0_r8

! Do an index sort of the ensemble members; Need sorted ensemble
call index_sort(ens, e_ind, ens_size)

! Prior distribution is boxcar in the central bins with 1/(n+1) density
! in each intermediate bin. BUT, distribution on the wings is a normal with
! 1/(n + 1) of the mass on each side.

! Begin by computing a weight for each of the prior ensemble membersA
! This is just evaluating the gaussian likelihood
const_term = 1.0_r8 / (sqrt(2.0_r8 * PI) * sqrt(obs_var))
do i = 1, ens_size
   weight(i) = const_term * exp(-1.0_r8 * (ens(i) - obs)**2 / (2.0_r8 * obs_var))
end do

! Compute the points that bound all the updated mass boxes; start with ensemble
do i = 1, ens_size
   x(2*i - 1) = ens(e_ind(i))
end do
! Compute the mid-point interior boundaries; these are halfway between ensembles
do i = 2, 2*ens_size - 2, 2
   x(i) = (x(i - 1) + x(i + 1)) / 2.0_r8
end do

! Compute the s.d. of the ensemble for getting the gaussian wings
sx         = sum(ens)
prior_mean = sx / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)

! Need to normalize the wings so they have 1/(ens_size + 1) mass outside
! Since 1/2 of a normal is outside, need to multiply by 2 / (ens_size + 1)

! Need some sort of width for the boundary kernel, try 1/2 the VAR for now
prior_var_d2 = prior_var / 2.0_r8

! Compute the product of the obs error gaussian with the prior gaussian (EAKF)
! Left wing first
var_ratio = obs_var / (prior_var_d2 + obs_var)
new_var = var_ratio * prior_var_d2
new_sd = sqrt(new_var)
new_mean_left  = var_ratio * (ens(e_ind(1))  + prior_var_d2*obs / obs_var)
new_mean_right  = var_ratio * (ens(e_ind(ens_size))  + prior_var_d2*obs / obs_var)
! REMEMBER, this product has an associated weight which must be taken into account
! See Anderson and Anderson for this weight term (or tutorial kernel filter)
prod_weight_left =  2.71828_r8 ** (-0.5_r8 * (ens(e_ind(1))**2 / prior_var_d2 + &
      obs**2 / obs_var - new_mean_left**2 / new_var)) / sqrt(2.0_r8 * PI)

prod_weight_right =  2.71828_r8 ** (-0.5_r8 * (ens(e_ind(ens_size))**2 / prior_var_d2 + &
      obs**2 / obs_var - new_mean_right**2 / new_var)) / sqrt(2.0_r8 * PI)

! Split into 2*ens_size domains; mass in each is computed
! Start by computing mass in the outermost (gaussian) regions
mass(1) = norm_cdf(ens(e_ind(1)), new_mean_left, new_sd) * &
   prod_weight_left * (2.0_r8 / (ens_size + 1.0_r8))
mass(2*ens_size) = (1.0_r8 - norm_cdf(ens(e_ind(ens_size)), new_mean_right, &
   new_sd)) * prod_weight_right * (2.0_r8 / (ens_size + 1.0_r8))

! Compute mass in the inner half boxes that have ensemble point on the left
do i = 2, 2*ens_size - 2, 2
   mass(i) = (1.0_r8 / (2.0_r8 * (ens_size + 1.0_r8))) * weight(e_ind(i/2))
end do

! Now right inner half boxes
do i = 3, 2*ens_size - 1, 2
   mass(i) = (1.0_r8 / (2.0_r8 * (ens_size + 1.0_r8))) * weight(e_ind(i/2 + 1))
end do

! Now normalize the mass in the different bins
mass_sum = sum(mass)
mass = mass / mass_sum

! Find cumulative mass at each box boundary and middle boundary
cumul_mass(0) = 0.0_r8
do i = 1, 2*ens_size
   cumul_mass(i) = cumul_mass(i - 1) + mass(i)
end do

! Get resampled ensemble, Need 1/(ens_size + 1) between each
umass = 1.0_r8 / (ens_size + 1.0_r8)

! Begin search at bottom of lowest box, but then update for efficiency
lowest_box = 1

! Find each new ensemble members location
do i = 1, ens_size
   ! If it's in the inner or outer range have to use normal
   if(umass < cumul_mass(1)) then
      ! In the first normal box
      left_weight = (1.0_r8 / mass_sum) * prod_weight_left * (2.0_r8 / (ens_size + 1.0_r8))
      call weighted_norm_inv(left_weight, new_mean_left, new_sd, umass, new_ens(i))
   else if(umass > cumul_mass(2*ens_size - 1)) then
      ! In the last normal box; Come in from the outside
      right_weight = (1.0_r8 / mass_sum) * prod_weight_right * (2.0_r8 / (ens_size + 1.0_r8))
      call weighted_norm_inv(right_weight, new_mean_right, new_sd, 1.0_r8 - umass, new_ens(i))
      new_ens(i) = new_mean_right + (new_mean_right - new_ens(i))
   else
      ! In one of the inner uniform boxes.
      FIND_BOX:do j = lowest_box, 2 * ens_size - 2
         ! Find the box that this mass is in
         if(umass >= cumul_mass(j) .and. umass <= cumul_mass(j + 1)) then
            new_ens(i) = x(j) + ((umass - cumul_mass(j)) / (cumul_mass(j+1) - cumul_mass(j))) * &
               (x(j + 1) - x(j))
            ! Don't need to search lower boxes again
            lowest_box = j
            exit FIND_BOX
         end if
      end do FIND_BOX
   endif
   ! Want equally partitioned mass in update with exception that outermost boxes have half
   umass = umass + 1.0_r8 / (ens_size + 1.0_r8)
end do

! Can now compute sorted increments
do i = 1, ens_size
   sort_inc(i) = new_ens(i) - ens(e_ind(i))
end do

! Now, need to convert to increments for unsorted
do i = 1, ens_size
   obs_inc(e_ind(i)) = sort_inc(i)
end do

end subroutine obs_increment_boxcar



subroutine obs_increment_rank_histogram(ens, ens_size, prior_var, &
   obs, obs_var, obs_inc)
!------------------------------------------------------------------------
!
! Revised 14 November 2008
!
! Does observation space update by approximating the prior distribution by
! a rank histogram. Prior and posterior are assumed to have 1/(n+1) probability
! mass between each ensemble member. The tails are assumed to be gaussian with
! a variance equal to sample variance of the entire ensemble and a mean
! selected so that 1/(n+1) of the mass is in each tail.
!
! The likelihood between the extreme ensemble members is approximated by
! quadrature. Two options are available and controlled by the namelist entry
! rectangular_quadrature. If this namelist is true than the likelihood between
! a pair of ensemble members is assumed to be uniform with the average of
! the likelihood computed at the two ensemble members. If it is false then
! the likelihood between two ensemble members is approximated by a line
! connecting the values of the likelihood computed at each of the ensemble
! members (trapezoidal quadrature).
!
! Two options are available for approximating the likelihood on the tails.
! If gaussian_likelihood_tails is true that the likelihood is assumed to
! be N(obs, obs_var) on the tails. If this is false, then the likelihood
! on the tails is taken to be uniform (to infinity) with the value at the
! outermost ensemble members.
!
! A product of the approximate prior and approximate posterior is taken
! and new ensemble members are located so that 1/(n+1) of the mass is between
! each member and on the tails.

! This code is still under development. Please contact Jeff Anderson at
! jla@ucar.edu if you are interested in trying it.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), prior_var, obs, obs_var
real(r8), intent(out) :: obs_inc(ens_size)

integer  :: i, e_ind(ens_size), lowest_box, j
real(r8) :: prior_sd, var_ratio, umass, left_amp, right_amp
real(r8) :: left_sd, left_var, right_sd, right_var, left_mean, right_mean
real(r8) :: mass(ens_size + 1), like(ens_size), cumul_mass(0:ens_size + 1)
real(r8) :: nmass(ens_size + 1)
real(r8) :: new_mean_left, new_mean_right, prod_weight_left, prod_weight_right
real(r8) :: new_var_left, new_var_right, new_sd_left, new_sd_right
real(r8) :: new_ens(ens_size), mass_sum
real(r8) :: x(ens_size)
real(r8) :: like_dense(2:ens_size), height(2:ens_size)
real(r8) :: dist_for_unit_sd
real(r8) :: a, b, c, hright, hleft, r1, r2, adj_r1, adj_r2

! Do an index sort of the ensemble members; Will want to do this very efficiently
call index_sort(ens, e_ind, ens_size)

do i = 1, ens_size
   ! The boundaries of the interior bins are just the sorted ensemble members
   x(i) = ens(e_ind(i))
   ! Compute likelihood for each ensemble member; just evaluate the gaussian
   ! No need to compute the constant term since relative likelihood is what matters
   like(i) = exp(-1.0_r8 * (x(i) - obs)**2 / (2.0_r8 * obs_var))
end do

! Prior distribution is boxcar in the central bins with 1/(n+1) density
! in each intermediate bin. BUT, distribution on the tails is a normal with
! 1/(n + 1) of the mass on each side.

! Can now compute the mean likelihood density in each interior bin
do i = 2, ens_size
   like_dense(i) = ((like(i - 1) + like(i)) / 2.0_r8)
end do

! Compute the s.d. of the ensemble for getting the gaussian tails
prior_sd = sqrt(prior_var)

! For unit normal, find distance from mean to where cdf is 1/(n+1)
! Lots of this can be done once in first call and then saved
call weighted_norm_inv(1.0_r8, 0.0_r8, 1.0_r8, &
   1.0_r8 / (ens_size + 1.0_r8), dist_for_unit_sd)
dist_for_unit_sd = -1.0_r8 * dist_for_unit_sd

! Have variance of tails just be sample prior variance
! Mean is adjusted so that 1/(n+1) is outside
left_mean = x(1) + dist_for_unit_sd * prior_sd
left_var = prior_var
left_sd = prior_sd
! Same for right tail
right_mean = x(ens_size) - dist_for_unit_sd * prior_sd
right_var = prior_var
right_sd = prior_sd

if(gaussian_likelihood_tails) then
   !*************** Block to do Gaussian-Gaussian on tail **************
   ! Compute the product of the obs likelihood gaussian with the priors
   ! Left tail gaussian first
   var_ratio = obs_var / (left_var + obs_var)
   new_var_left = var_ratio * left_var
   new_sd_left = sqrt(new_var_left)
   new_mean_left  = var_ratio * (left_mean  + left_var*obs / obs_var)
   ! REMEMBER, this product has an associated weight which must be taken into account
   ! See Anderson and Anderson for this weight term (or tutorial kernel filter)
   ! NOTE: The constant term has been left off the likelihood so we don't have
   ! to divide by sqrt(2 PI) in this expression
   prod_weight_left =  exp(-0.5_r8 * (left_mean**2 / left_var + &
         obs**2 / obs_var - new_mean_left**2 / new_var_left)) / &
         sqrt(left_var + obs_var)
   ! Determine how much mass is in the updated tails by computing gaussian cdf
   mass(1) = norm_cdf(x(1), new_mean_left, new_sd_left) * prod_weight_left

   ! Same for the right tail
   var_ratio = obs_var / (right_var + obs_var)
   new_var_right = var_ratio * right_var
   new_sd_right = sqrt(new_var_right)
   new_mean_right  = var_ratio * (right_mean  + right_var*obs / obs_var)
   ! NOTE: The constant term has been left off the likelihood so we don't have
   ! to divide by sqrt(2 PI) in this expression
   prod_weight_right =  exp(-0.5_r8 * (right_mean**2 / right_var + &
         obs**2 / obs_var - new_mean_right**2 / new_var_right)) / &
         sqrt(right_var + obs_var)
   ! Determine how much mass is in the updated tails by computing gaussian cdf
   mass(ens_size + 1) = (1.0_r8 - norm_cdf(x(ens_size), new_mean_right, &
      new_sd_right)) * prod_weight_right
   !************ End Block to do Gaussian-Gaussian on tail **************
else
   !*************** Block to do flat tail for likelihood ****************
   ! Flat tails: THIS REMOVES ASSUMPTIONS ABOUT LIKELIHOOD AND CUTS COST
   new_var_left = left_var
   new_sd_left = left_sd
   new_mean_left = left_mean
   prod_weight_left = like(1)
   mass(1) = like(1) / (ens_size + 1.0_r8)

   ! Same for right tail
   new_var_right = right_var
   new_sd_right = right_sd
   new_mean_right = right_mean
   prod_weight_right = like(ens_size)
   mass(ens_size + 1) = like(ens_size) / (ens_size + 1.0_r8)
   !*************** End block to do flat tail for likelihood ****************
endif

! The mass in each interior box is the height times the width
! The height of the likelihood is like_dense
! For the prior, mass is 1/(n+1),   and mass = height x width so...
! The height of the prior is 1 / ((n+1) width);   multiplying by width leaves 1/(n+1)

! In prior, have 1/(n+1) mass in each bin, multiply by mean likelihood density
! to get approximate mass in updated bin
do i = 2, ens_size
   mass(i) = like_dense(i) / (ens_size + 1.0_r8)
   ! Height of prior in this bin is mass/width; Only needed for trapezoidal
   ! If two ensemble members are the same, set height to -1 as flag
   if(x(i) == x(i - 1)) then
      height(i) = -1.0_r8
   else
      height(i) = 1.0_r8 / ((ens_size + 1.0_r8) * (x(i) - x(i-1)))
   endif
end do

! Now normalize the mass in the different bins to get a pdf
mass_sum = sum(mass)
nmass = mass / mass_sum

! Get the weight for the final normalized tail gaussians
! This is the same as left_amp=(ens_size + 1)*nmass(1)
left_amp = prod_weight_left / mass_sum
! This is the same as right_amp=(ens_size + 1)*nmass(ens_size + 1)
right_amp = prod_weight_right / mass_sum

! Find cumulative mass at each box boundary and middle boundary
cumul_mass(0) = 0.0_r8
do i = 1, ens_size + 1
   cumul_mass(i) = cumul_mass(i - 1) + nmass(i)
end do

! Begin intenal box search at bottom of lowest box, update for efficiency
lowest_box = 1

! Find each new ensemble members location
do i = 1, ens_size
   ! Each update ensemble member has 1/(n+1) mass before it
   umass = (1.0_r8 * i) / (ens_size + 1.0_r8)

   ! If it is in the inner or outer range have to use normal
   if(umass < cumul_mass(1)) then
      ! It's in the left tail
      ! Get position of x in weighted gaussian where the cdf has value umass
      call weighted_norm_inv(left_amp, new_mean_left, new_sd_left, &
         umass, new_ens(i))
   else if(umass > cumul_mass(ens_size)) then
      ! It's in the right tail
      ! Get position of x in weighted gaussian where the cdf has value umass
      call weighted_norm_inv(right_amp, new_mean_right, new_sd_right, &
         1.0_r8 - umass, new_ens(i))
      ! Coming in from the right, use symmetry after pretending its on left
      new_ens(i) = new_mean_right + (new_mean_right - new_ens(i))
   else
      ! In one of the inner uniform boxes.
      FIND_BOX:do j = lowest_box, ens_size - 1
         ! Find the box that this mass is in
         if(umass >= cumul_mass(j) .and. umass <= cumul_mass(j + 1)) then

            if(rectangular_quadrature) then
               !********* Block for rectangular quadrature *******************
               ! Linearly interpolate in mass
               new_ens(i) = x(j) + ((umass - cumul_mass(j)) / &
                  (cumul_mass(j+1) - cumul_mass(j))) * (x(j + 1) - x(j))
               !********* End block for rectangular quadrature *******************

            else

               !********* Block for trapezoidal interpolation *******************
               ! Assume that mass has linear profile, quadratic interpolation
               ! If two ensemble members are the same, just keep that value
               if(height(j + 1) < 0) then
                  new_ens(i) = x(j)
               else
                  ! Height on left side and right side
                  hleft = height(j + 1) * like(j) / mass_sum
                  hright = height(j + 1) * like(j + 1) / mass_sum
                  ! Will solve a quadratic for desired x-x(j)
                  ! a is 0.5(hright - hleft) / (x(j+1) - x(j))
                  a = 0.5_r8 * (hright - hleft) / (x(j+1) - x(j))
                  ! b is hleft
                  b = hleft
                  ! c is cumul_mass(j) - umass
                  c = cumul_mass(j) - umass
                  ! Use stable quadratic solver
                  call solve_quadratic(a, b, c, r1, r2)
                  adj_r1 = r1 + x(j)
                  adj_r2 = r2 + x(j)
                  if(adj_r1 >= x(j) .and. adj_r1 <= x(j+1)) then
                     new_ens(i) = adj_r1
                  elseif (adj_r2 >= x(j) .and. adj_r2 <= x(j+1)) then
                     new_ens(i) = adj_r2
                  else
                     msgstring = 'Did not get a satisfactory quadratic root'
                     call error_handler(E_ERR, 'obs_increment_rank_histogram', msgstring, &
                        source)
                  endif
               endif
               !********* End block for quadratic interpolation *******************

            endif

            ! Don't need to search lower boxes again
            lowest_box = j
            exit FIND_BOX
         end if
      end do FIND_BOX
   endif
end do

! Convert to increments for unsorted
do i = 1, ens_size
   obs_inc(e_ind(i)) = new_ens(i) - x(i)
end do

end subroutine obs_increment_rank_histogram




subroutine update_ens_from_weights(ens, ens_size, rel_weight, ens_inc)
!------------------------------------------------------------------------
! Given relative weights for an ensemble, compute increments for the
! ensemble members. Assumes that prior distributon is equal uniform mass
! between each ensemble member. On the edges, have a normal with the
! sample mean and s.d. BUT normalized by a factor alpha so that only
! 1/(2*ens_size) of the total mass lies on each flank.

integer,  intent(in)  :: ens_size
real(r8), intent(in)  :: ens(ens_size), rel_weight(ens_size)
real(r8), intent(out) :: ens_inc(ens_size)

integer  :: i, j, lowest_box
integer  :: e_ind(ens_size)
real(r8) :: x(1:2*ens_size - 1), cumul_mass(1:2*ens_size - 1), new_ens(ens_size)
real(r8) :: sort_inc(ens_size), updated_mass(2 * ens_size)
real(r8) :: sx, prior_mean, prior_var, prior_sd, mass
real(r8) :: total_mass_left, total_mass_right, alpha(2)

! Initialize assim_tools_module if needed
if (.not. module_initialized) call assim_tools_init()

call error_handler(E_ERR,'update_ens_from_weight','Routine needs testing.', &
           source, text2='Talk to Jeff before using.')

! Do an index sort of the ensemble members
call index_sort(ens, e_ind, ens_size)

! Have half boxes between all ensembles in the interior
! Total number of mass boxes is 2*ens_size

! Compute the points that bound all the updated mass boxes; start with ensemble
do i = 1, ens_size
   x(2*i - 1) = ens(e_ind(i))
end do
! Compute the mid-point interior boundaries; these are halfway between ensembles
do i = 2, 2*ens_size - 2, 2
   x(i) = (x(i - 1) + x(i + 1)) / 2.0_r8
end do

! Compute the mean and s.d. of the prior ensemble to handle wings
sx         = sum(ens)
prior_mean = sx / ens_size
prior_var  = sum((ens - prior_mean)**2) / (ens_size - 1)
prior_sd = sqrt(prior_var)

! Need to normalize the wings so they have 1/(2*ens_size) mass outside
! Use cdf to find out how much mass is left of 1st member, right of last
total_mass_left = norm_cdf(ens(e_ind(1)), prior_mean, prior_sd)
total_mass_right = 1.0_r8 - norm_cdf(ens(e_ind(ens_size)), prior_mean, prior_sd)

! Find the mass in each division given the initial equal partition and the weights
updated_mass(1) = rel_weight(e_ind(1)) / (2.0_r8 * ens_size)
updated_mass(2 * ens_size) = rel_weight(e_ind(ens_size)) / (2.0_r8 * ens_size)
do i = 2, 2*ens_size - 2, 2
   updated_mass(i) = rel_weight(e_ind(i / 2)) / (2.0_r8 * ens_size)
end do
do i = 3, 2*ens_size - 1, 2
   updated_mass(i) = rel_weight(e_ind((i+1) / 2)) / (2.0_r8 * ens_size)
end do

! Normalize the mass; (COULD IT EVER BE 0 necessitating error check?)
updated_mass = updated_mass / sum(updated_mass)

! Find a normalization factor to get tail mass right
if(total_mass_left > 0.0_r8) then
   alpha(1) = updated_mass(1) / total_mass_left
else
   alpha(1) = 0.0_r8
endif
if(total_mass_right > 0.0_r8) then
   alpha(2) = updated_mass(2 * ens_size) / total_mass_right
else
   alpha(2) = 0.0_r8
endif

! Find cumulative mass at each box boundary and middle boundary
cumul_mass(1) = updated_mass(1)
do i = 2, 2*ens_size - 1
   cumul_mass(i) = cumul_mass(i - 1) + updated_mass(i)
end do

! Get resampled position an inefficient way
! Need 1/ens_size between each EXCEPT for outers which get half of this
mass = 1.0_r8 / (2.0_r8 * ens_size)

do i = 1, ens_size
   ! If it's in the inner or outer range have to use normal
   if(mass < cumul_mass(1)) then
      ! In the first normal box
      call weighted_norm_inv(alpha(1), prior_mean, prior_sd, mass, new_ens(i))
   else if(mass > cumul_mass(2*ens_size - 1)) then
      ! In the last normal box; Come in from the outside
      call weighted_norm_inv(alpha(2), prior_mean, prior_sd, 1.0_r8 - mass, new_ens(i))
      new_ens(i) = prior_mean + (prior_mean - new_ens(i))
   else
      ! In one of the inner uniform boxes. Make this much more efficient search?
      lowest_box = 1
      FIND_BOX:do j = lowest_box, 2 * ens_size - 2
         ! Find the box that this mass is in
         if(mass >= cumul_mass(j) .and. mass <= cumul_mass(j + 1)) then
            new_ens(i) = x(j) + ((mass - cumul_mass(j)) / (cumul_mass(j+1) - cumul_mass(j))) * &
               (x(j + 1) - x(j))
            ! Don't need to search lower boxes again
            lowest_box = j
            exit FIND_BOX
         end if
      end do FIND_BOX
   endif
   ! Want equally partitioned mass in update with exception that outermost boxes have half
   mass = mass + 1.0_r8 / ens_size
end do

! Can now compute sorted increments
do i = 1, ens_size
   sort_inc(i) = new_ens(i) - ens(e_ind(i))
end do

! Now, need to convert to increments for unsorted
do i = 1, ens_size
   ens_inc(e_ind(i)) = sort_inc(i)
end do

end subroutine update_ens_from_weights
!---------------------------------------------------------------

subroutine obs_updates_ens(ens_size, num_groups, ens, ens_loc, ens_kind, &
   obs_prior, obs_inc, obs_prior_mean, obs_prior_var, obs_loc, obs_type, obs_time,    &
   net_a, grp_size, grp_beg, grp_end, reg_factor_obs_index,         &
   reg_factor_ens_index, final_factor, correl, correl_needed, inflate_only, &
   state_prior, one_inner_to_state_inc)

integer,             intent(in)  :: ens_size
integer,             intent(in)  :: num_groups

! CCHU: important change below:
! Note: ens = state_handle%copies (the ensemble member for a "current" state)
! which is not used at all in PFF-DART
!real(r8),            intent(inout)  :: ens(ens_size) ! current state variable
real(r8),            intent(in)  :: ens(ens_size)

type(location_type), intent(in)  :: ens_loc
integer,             intent(in)  :: ens_kind
real(r8),            intent(in)  :: obs_prior(ens_size)
real(r8),            intent(in)  :: obs_inc(ens_size)
real(r8),            intent(in)  :: obs_prior_mean(num_groups)
real(r8),            intent(in)  :: obs_prior_var(num_groups)
type(location_type), intent(in)  :: obs_loc
integer,             intent(in)  :: obs_type
type(time_type),     intent(in)  :: obs_time
real(r8),            intent(in)  :: net_a(num_groups)
integer,             intent(in)  :: grp_size
integer,             intent(in)  :: grp_beg(num_groups)
integer,             intent(in)  :: grp_end(num_groups)
integer,             intent(in)  :: reg_factor_obs_index
integer(i8),         intent(in)  :: reg_factor_ens_index
real(r8),            intent(inout) :: final_factor
real(r8),            intent(out) :: correl(num_groups)
logical,             intent(in)  :: correl_needed
logical,             intent(in)  :: inflate_only

real(r8) :: reg_coef(num_groups), increment(ens_size)
real(r8) :: reg_factor
integer  :: group, grp_bot, grp_top

! CCHU: two new argument just for PFF-DART:
real(r8), intent(in),  optional:: state_prior(ens_size) ! prior state variable
real(r8), intent(out), optional:: one_inner_to_state_inc(ens_size)

! Loop through groups to update the state variable ensemble members
do group = 1, num_groups
   grp_bot = grp_beg(group); grp_top = grp_end(group)
   ! Do update of state, correl only needed for varying ss inflate
   if(correl_needed) then
     ! call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
     !    obs_prior_var(group), obs_inc(grp_bot:grp_top), ens(grp_bot:grp_top), grp_size, &
     !    increment(grp_bot:grp_top), reg_coef(group), net_a(group), correl(group))
      call update_from_obs_inc(obs_prior(grp_bot:grp_top),obs_prior_mean(group), &
         obs_prior_var(group), obs_inc(grp_bot:grp_top), state_prior, grp_size, &
         increment(grp_bot:grp_top), reg_coef(group), net_a(group),correl(group))

   else
     ! call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
     !    obs_prior_var(group), obs_inc(grp_bot:grp_top), ens(grp_bot:grp_top), grp_size, &
     !    increment(grp_bot:grp_top), reg_coef(group), net_a(group))
      call update_from_obs_inc(obs_prior(grp_bot:grp_top), obs_prior_mean(group), &
         obs_prior_var(group), obs_inc(grp_bot:grp_top), state_prior, grp_size, &
         increment(grp_bot:grp_top), reg_coef(group), net_a(group))
   endif
end do

if(num_groups > 1) then
   reg_factor = comp_reg_factor(num_groups, reg_coef, obs_time, &
      reg_factor_obs_index, reg_factor_ens_index)
   final_factor = min(final_factor, reg_factor)
endif

! Get the updated ensemble
!if(.not. inflate_only) ens = ens + final_factor * increment

! CCWU: important change here!
if (.not. inflate_only) one_inner_to_state_inc = final_factor * increment

end subroutine obs_updates_ens

!-------------------------------------------------------------

function cov_and_impact_factors(base_obs_loc, base_obs_type, state_loc, state_kind, &
dist, cutoff_rev)

! Computes the cov_factor and multiplies by obs_impact_factor if selected

real(r8) :: cov_and_impact_factors
type(location_type), intent(in) :: base_obs_loc
integer, intent(in) :: base_obs_type
type(location_type), intent(in) :: state_loc
integer, intent(in) :: state_kind
real(r8), intent(in) :: dist
real(r8), intent(in) :: cutoff_rev

real(r8) :: impact_factor, cov_factor

! Get external impact factors, cycle if impact of this ob on this state is zero
if (adjust_obs_impact) then
   ! Get the impact factor from the table if requested
   impact_factor = obs_impact_table(base_obs_type, state_kind)
   if(impact_factor <= 0.0_r8) then
      ! Avoid the cost of computing cov_factor if impact is 0
      cov_and_impact_factors = 0.0_r8
      return
   endif
else
   impact_factor = 1.0_r8
endif

! Compute the covariance factor
cov_factor = comp_cov_factor(dist, cutoff_rev, &
   base_obs_loc, base_obs_type, state_loc, state_kind)

! Combine the impact_factor and the cov_factor
cov_and_impact_factors = cov_factor * impact_factor

end function cov_and_impact_factors


!------------------------------------------------------------------------

function norm_cdf(x_in, mean, sd)

! Approximate cumulative distribution function for normal
! with mean and sd evaluated at point x_in
! Only works for x>= 0.

real(r8)             :: norm_cdf
real(r8), intent(in) :: x_in, mean, sd

real(digits12) :: x, p, b1, b2, b3, b4, b5, t, density, nx

! Convert to a standard normal
nx = (x_in - mean) / sd

x = abs(nx)


! Use formula from Abramowitz and Stegun to approximate
p = 0.2316419_digits12
b1 = 0.319381530_digits12
b2 = -0.356563782_digits12
b3 = 1.781477937_digits12
b4 = -1.821255978_digits12
b5 = 1.330274429_digits12

t = 1.0_digits12 / (1.0_digits12 + p * x)

density = (1.0_digits12 / sqrt(2.0_digits12 * PI)) * exp(-x*x / 2.0_digits12)

norm_cdf = 1.0_digits12 - density * &
   ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t

if(nx < 0.0_digits12) norm_cdf = 1.0_digits12 - norm_cdf

!write(*, *) 'cdf is ', norm_cdf

end function norm_cdf


!------------------------------------------------------------------------

subroutine weighted_norm_inv(alpha, mean, sd, p, x)

! Find the value of x for which the cdf of a N(mean, sd) multiplied times
! alpha has value p.

real(r8), intent(in)  :: alpha, mean, sd, p
real(r8), intent(out) :: x

real(r8) :: np

! Can search in a standard normal, then multiply by sd at end and add mean
! Divide p by alpha to get the right place for weighted normal
np = p / alpha

! Find spot in standard normal
call norm_inv(np, x)

! Add in the mean and normalize by sd
x = mean + x * sd

end subroutine weighted_norm_inv


!------------------------------------------------------------------------

subroutine norm_inv(p, x)

real(r8), intent(in)  :: p
real(r8), intent(out) :: x

! normal inverse
! translate from http://home.online.no/~pjacklam/notes/invnorm
! a routine written by john herrero

real(r8) :: p_low,p_high
real(r8) :: a1,a2,a3,a4,a5,a6
real(r8) :: b1,b2,b3,b4,b5
real(r8) :: c1,c2,c3,c4,c5,c6
real(r8) :: d1,d2,d3,d4
real(r8) :: q,r
a1 = -39.69683028665376_digits12
a2 =  220.9460984245205_digits12
a3 = -275.9285104469687_digits12
a4 =  138.357751867269_digits12
a5 = -30.66479806614716_digits12
a6 =  2.506628277459239_digits12
b1 = -54.4760987982241_digits12
b2 =  161.5858368580409_digits12
b3 = -155.6989798598866_digits12
b4 =  66.80131188771972_digits12
b5 = -13.28068155288572_digits12
c1 = -0.007784894002430293_digits12
c2 = -0.3223964580411365_digits12
c3 = -2.400758277161838_digits12
c4 = -2.549732539343734_digits12
c5 =  4.374664141464968_digits12
c6 =  2.938163982698783_digits12
d1 =  0.007784695709041462_digits12
d2 =  0.3224671290700398_digits12
d3 =  2.445134137142996_digits12
d4 =  3.754408661907416_digits12
p_low  = 0.02425_digits12
p_high = 1_digits12 - p_low
! Split into an inner and two outer regions which have separate fits
if(p < p_low) then
   q = sqrt(-2.0_digits12 * log(p))
   x = (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
      ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_digits12)
else if(p > p_high) then
   q = sqrt(-2.0_digits12 * log(1.0_digits12 - p))
   x = -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) / &
      ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0_digits12)
else
   q = p - 0.5_digits12
   r = q*q
   x = (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6)*q / &
      (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0_digits12)
endif

end subroutine norm_inv

!------------------------------------------------------------------------

subroutine set_assim_tools_trace(execution_level, timestamp_level)
 integer, intent(in) :: execution_level
 integer, intent(in) :: timestamp_level

! set module local vars from the calling code to indicate how much
! output we should generate from this code.  execution level is
! intended to make it easier to figure out where in the code a crash
! is happening; timestamp level is intended to help with gross levels
! of overall performance profiling.  eventually, a level of 1 will
! print out only basic info; level 2 will be more detailed.
! (right now, only > 0 prints anything and it doesn't matter how
! large the value is.)

! Initialize assim_tools_module if needed
if (.not. module_initialized) call assim_tools_init()

print_trace_details = execution_level
print_timestamps    = timestamp_level

end subroutine set_assim_tools_trace

!--------------------------------------------------------------------

function revised_distance(orig_dist, newcount, oldcount, base, cutfloor)
 real(r8),            intent(in) :: orig_dist
 integer,             intent(in) :: newcount, oldcount
 type(location_type), intent(in) :: base
 real(r8),            intent(in) :: cutfloor

 real(r8)                        :: revised_distance

! take the ratio of the old and new counts, and revise the
! original cutoff distance to match.

! for now, only allow the code to do a 2d area adaption.
! to experiment with other schemes, set this local variable
! to .false. at the top of the file and recompile.

if (only_area_adapt) then

   revised_distance = orig_dist * sqrt(real(newcount, r8) / oldcount)

   ! allow user to set a minimum cutoff, so even if there are very dense
   ! observations the cutoff distance won't go below this floor.
   if (revised_distance < cutfloor) revised_distance = cutfloor
   return

endif

! alternatives for different dimensionalities and schemes

! Change the cutoff radius to get the appropriate number
if (LocationDims == 1) then
   ! linear (be careful of cyclic domains; if > domain, this is
   ! not going to be right)
   revised_distance = orig_dist * real(newcount, r8) / oldcount

else if (LocationDims == 2) then
   ! do an area scaling
   revised_distance = orig_dist * sqrt(real(newcount, r8) / oldcount)

else if (LocationDims == 3) then
   ! do either a volume or area scaling (depending on whether we are
   ! localizing in the vertical or not.)   if surface obs, assume a hemisphere
   ! and shrink more.

   if (vertical_localization_on()) then
      ! cube root for volume
      revised_distance = orig_dist * ((real(newcount, r8) / oldcount) &
                                      ** 0.33333333333333333333_r8)

      ! Cut the adaptive localization threshold in half again for 'surface' obs
      if (is_vertical(base, "SURFACE")) then
         revised_distance = revised_distance * (0.5_r8 ** 0.33333333333333333333_r8)
      endif
   else
      ! do an area scaling, even if 3d obs
      revised_distance = orig_dist * sqrt(real(newcount, r8) / oldcount)

      ! original code was:
      !cutoff_rev =  sqrt((2.0_r8*cutoff)**2 * adaptive_localization_threshold / &
      !   total_num_close_obs) / 2.0_r8

      ! original comment
      ! Need to get thinning out of assim_tools and into something about locations
   endif
else
   call error_handler(E_ERR, 'revised_distance', 'unknown locations dimension, not 1, 2 or 3', &
      source)
endif

! allow user to set a minimum cutoff, so even if there are very dense
! observations the cutoff distance won't go below this floor.
if (revised_distance < cutfloor) revised_distance = cutfloor

end function revised_distance

!--------------------------------------------------------------------

function count_close(num_close, index_list, my_types, dist, maxdist)
 integer, intent(in)  :: num_close, index_list(:), my_types(:)
 real(r8), intent(in) :: dist(:), maxdist
 integer :: count_close

! return the total number of items from the index_list which
! are types which are going to be assimilated, and within distance.
! this excludes items on the eval list only, not listed, or
! items too far away.   this routine does a global communication
! so if any MPI tasks make this call, all must.

integer :: k, thistype, local_count

local_count = 0
do k=1, num_close

   ! only accept items closer than limit
   if (dist(k) > maxdist) cycle

   ! include identity obs, plus types on assim list.
   ! you have to do the if tests separately because fortran allows
   ! both parts of an if(a .or. b) test to be eval'd at the same time.
   ! you'd be using a negative index if it was an identity obs.
   thistype = my_types(index_list(k))
   if (thistype < 0) then
      local_count = local_count + 1
   else if (assimilate_this_type_of_obs(thistype)) then
      local_count = local_count + 1
   endif
end do

! broadcast sums from all tasks to compute new total
call sum_across_tasks(local_count, count_close)

end function count_close

!----------------------------------------------------------------------
! Revise the cutoff for this observation if adaptive localization is required
! Output diagnostics for localization if requested

subroutine adaptive_localization_and_diags(cutoff_orig, cutoff_rev, adaptive_localization_threshold, &
   adaptive_cutoff_floor, num_close_obs, close_obs_ind, close_obs_dist, my_obs_type, &
   base_obs_index, base_obs_loc, obs_def, out_unit)

real(r8),            intent(in)  :: cutoff_orig
real(r8),            intent(out) :: cutoff_rev
integer,             intent(in)  :: adaptive_localization_threshold
real(r8),            intent(in)  :: adaptive_cutoff_floor
integer,             intent(in)  :: num_close_obs
integer,             intent(in)  :: close_obs_ind(:)
real(r8),            intent(in)  :: close_obs_dist(:)
integer,             intent(in)  :: my_obs_type(:)
integer,             intent(in)  :: base_obs_index
type(location_type), intent(in)  :: base_obs_loc
type(obs_def_type),  intent(in)  :: obs_def
integer,             intent(in)  :: out_unit

integer :: total_num_close_obs, rev_num_close_obs, secs, days
type(time_type) :: this_obs_time
character(len = 200) :: base_loc_text   ! longer than longest location formatting possible

! Default is that cutoff is not revised
cutoff_rev = cutoff_orig

! For adaptive localization, need number of other obs close to the chosen observation
if(adaptive_localization_threshold > 0) then
   ! this does a cross-task sum, so all tasks must make this call.
   total_num_close_obs = count_close(num_close_obs, close_obs_ind, my_obs_type, &
                                     close_obs_dist, cutoff_rev*2.0_r8)

   ! Want expected number of close observations to be reduced to some threshold;
   ! accomplish this by cutting the size of the cutoff distance.
   if(total_num_close_obs > adaptive_localization_threshold) then
      cutoff_rev = revised_distance(cutoff_rev*2.0_r8, adaptive_localization_threshold, &
                                    total_num_close_obs, base_obs_loc, &
                                    adaptive_cutoff_floor*2.0_r8) / 2.0_r8
   endif
endif

if ( output_localization_diagnostics ) then
   ! Warning, this can be costly and generate large output
   ! This is referred to as revised in case adaptive localization was done
   rev_num_close_obs = count_close(num_close_obs, close_obs_ind, my_obs_type, &
                                     close_obs_dist, cutoff_rev*2.0_r8)

   ! Output diagnostic information about the number of close obs
   if (my_task_id() == 0) then
      this_obs_time = get_obs_def_time(obs_def)
      call get_time(this_obs_time,secs,days)
      call write_location(-1, base_obs_loc, charstring=base_loc_text)

      ! If adaptive localization did something, output info about what it did
      ! Probably would be more consistent to just output for all observations
      if(adaptive_localization_threshold > 0 .and. &
         total_num_close_obs > adaptive_localization_threshold) then
         write(out_unit,'(i12,1x,i5,1x,i8,1x,A,2(f14.5,1x,i12))') base_obs_index, &
            secs, days, trim(base_loc_text), cutoff_orig, total_num_close_obs, cutoff_rev, &
            rev_num_close_obs
      else
         write(out_unit,'(i12,1x,i5,1x,i8,1x,A,f14.5,1x,i12)') base_obs_index, &
            secs, days, trim(base_loc_text), cutoff_rev, rev_num_close_obs
      endif
   endif
endif

end subroutine adaptive_localization_and_diags

!----------------------------------------------------------------------
!> gets the location of of all my observations
subroutine get_my_obs_loc(obs_ens_handle, obs_seq, keys, my_obs_loc, my_obs_kind, my_obs_type, my_obs_time)

type(ensemble_type),      intent(in)  :: obs_ens_handle
type(obs_sequence_type),  intent(in)  :: obs_seq
integer,                  intent(in)  :: keys(:)
type(location_type),      intent(out) :: my_obs_loc(:)
integer,                  intent(out) :: my_obs_type(:), my_obs_kind(:)
type(time_type),          intent(out) :: my_obs_time

type(obs_type) :: observation
type(obs_def_type)   :: obs_def
integer :: this_obs_key
integer i
type(location_type) :: dummyloc

Get_Obs_Locations: do i = 1, obs_ens_handle%my_num_vars

   this_obs_key = keys(obs_ens_handle%my_vars(i)) ! if keys becomes a local array, this will need changing
   call get_obs_from_key(obs_seq, this_obs_key, observation)
   call get_obs_def(observation, obs_def)
   my_obs_loc(i)  = get_obs_def_location(obs_def)
   my_obs_type(i) = get_obs_def_type_of_obs(obs_def)
   if (my_obs_type(i) > 0) then
         my_obs_kind(i) = get_quantity_for_type_of_obs(my_obs_type(i))
   else
      call get_state_meta_data(-1 * int(my_obs_type(i),i8), dummyloc, my_obs_kind(i))
   endif
end do Get_Obs_Locations

! Need the time for regression diagnostics potentially; get from first observation
my_obs_time = get_obs_def_time(obs_def)

end subroutine get_my_obs_loc

!--------------------------------------------------------------------
!> Get close obs from cache if appropriate. Cache new get_close_obs info
!> if requested.

subroutine get_close_obs_cached(gc_obs, base_obs_loc, base_obs_type, &
   my_obs_loc, my_obs_kind, my_obs_type, num_close_obs, close_obs_ind, close_obs_dist,  &
   ens_handle, last_base_obs_loc, last_num_close_obs, num_close_obs_cached,               &
   num_close_obs_calls_made)

type(get_close_type),          intent(in)  :: gc_obs
type(location_type),           intent(inout) :: base_obs_loc, my_obs_loc(:)
integer,                       intent(in)  :: base_obs_type, my_obs_kind(:), my_obs_type(:)
integer,                       intent(out) :: num_close_obs
integer,                       intent(inout) :: close_obs_ind(:)
real(r8),                      intent(inout) :: close_obs_dist(:)
type(ensemble_type),           intent(in)  :: ens_handle
type(location_type), intent(inout) :: last_base_obs_loc
integer, intent(inout) :: last_num_close_obs
integer, intent(inout) :: num_close_obs_cached, num_close_obs_calls_made

! This logic could be arranged to make code less redundant
if (.not. close_obs_caching) then
   call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                      my_obs_loc, my_obs_kind, my_obs_type, &
                      num_close_obs, close_obs_ind, close_obs_dist, ens_handle)
else
   if (base_obs_loc == last_base_obs_loc) then
      num_close_obs     = last_num_close_obs
      num_close_obs_cached = num_close_obs_cached + 1
   else
      call get_close_obs(gc_obs, base_obs_loc, base_obs_type, &
                         my_obs_loc, my_obs_kind, my_obs_type, &
                         num_close_obs, close_obs_ind, close_obs_dist, ens_handle)

      last_base_obs_loc      = base_obs_loc
      last_num_close_obs     = num_close_obs
      num_close_obs_calls_made = num_close_obs_calls_made +1
   endif
endif

end subroutine get_close_obs_cached

!--------------------------------------------------------------------
!> Get close state from cache if appropriate. Cache new get_close_state info
!> if requested.

subroutine get_close_state_cached(gc_state, base_obs_loc, base_obs_type, &
   my_state_loc, my_state_kind, my_state_indx, num_close_states, close_state_ind, close_state_dist,  &
   ens_handle, last_base_states_loc, last_num_close_states, num_close_states_cached,               &
   num_close_states_calls_made)

type(get_close_type),          intent(in)    :: gc_state
type(location_type),           intent(inout) :: base_obs_loc, my_state_loc(:)
integer,                       intent(in)    :: base_obs_type, my_state_kind(:)
integer(i8),                   intent(in)    :: my_state_indx(:)
integer,                       intent(out)   :: num_close_states
integer,                       intent(inout) :: close_state_ind(:)
real(r8),                      intent(inout) :: close_state_dist(:)
type(ensemble_type),           intent(in)    :: ens_handle
type(location_type), intent(inout) :: last_base_states_loc
integer, intent(inout) :: last_num_close_states
integer, intent(inout) :: num_close_states_cached, num_close_states_calls_made

! This logic could be arranged to make code less redundant
if (.not. close_obs_caching) then
   call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                      my_state_loc, my_state_kind, my_state_indx, &
                      num_close_states, close_state_ind, close_state_dist, ens_handle)
else
   if (base_obs_loc == last_base_states_loc) then
      num_close_states     = last_num_close_states
      num_close_states_cached = num_close_states_cached + 1
   else
      call get_close_state(gc_state, base_obs_loc, base_obs_type, &
                         my_state_loc, my_state_kind, my_state_indx, &
                         num_close_states, close_state_ind, close_state_dist, ens_handle)

      last_base_states_loc      = base_obs_loc
      last_num_close_states     = num_close_states
      num_close_states_calls_made = num_close_states_calls_made +1
   endif
endif

end subroutine get_close_state_cached

!--------------------------------------------------------------------
!> log what the user has selected via the namelist choices

subroutine log_namelist_selections(num_special_cutoff, cache_override)

integer, intent(in) :: num_special_cutoff
logical, intent(in) :: cache_override

integer :: i

select case (filter_kind)
 case (1)
   msgstring = 'Ensemble Adjustment Kalman Filter (EAKF)'
 case (2)
   msgstring = 'Ensemble Kalman Filter (ENKF)'
 case (3)
   msgstring = 'Kernel filter'
 case (4)
   msgstring = 'observation space particle filter'
 case (5)
   msgstring = 'random draw from posterior'
 case (6)
   msgstring = 'deterministic draw from posterior with fixed kurtosis'
 case (7)
   msgstring = 'Boxcar'
 case (8)
   msgstring = 'Rank Histogram Filter'
 case (9)
   msgstring = 'Particle Flow Filter'
 case default
   call error_handler(E_ERR, 'assim_tools_init:', 'illegal filter_kind value, valid values are 1-8', &
                      source)
end select
call error_handler(E_MSG, 'assim_tools_init:', 'Selected filter type is '//trim(msgstring))

if (adjust_obs_impact) then
   call allocate_impact_table(obs_impact_table)
   call read_impact_table(obs_impact_filename, obs_impact_table, allow_any_impact_values, "allow_any_impact_values")
   call error_handler(E_MSG, 'assim_tools_init:', &
                      'Using observation impact table from file "'//trim(obs_impact_filename)//'"')
endif

write(msgstring,  '(A,F18.6)') 'The cutoff namelist value is ', cutoff
write(msgstring2, '(A)') 'cutoff is the localization half-width parameter,'
write(msgstring3, '(A,F18.6)') 'so the effective localization radius is ', cutoff*2.0_r8
call error_handler(E_MSG,'assim_tools_init:', msgstring, text2=msgstring2, text3=msgstring3)

if (has_special_cutoffs) then
   call error_handler(E_MSG, '', '')
   call error_handler(E_MSG,'assim_tools_init:','Observations with special localization treatment:')
   call error_handler(E_MSG,'assim_tools_init:','(type name, specified cutoff distance, effective localization radius)')

   do i = 1, num_special_cutoff
      write(msgstring, '(A32,F18.6,F18.6)') special_localization_obs_types(i), &
            special_localization_cutoffs(i), special_localization_cutoffs(i)*2.0_r8
      call error_handler(E_MSG,'assim_tools_init:', msgstring)
   end do
   call error_handler(E_MSG,'assim_tools_init:','all other observation types will use the default cutoff distance')
   call error_handler(E_MSG, '', '')
endif

if (cache_override) then
   call error_handler(E_MSG,'assim_tools_init:','Disabling the close obs caching because specialized localization')
   call error_handler(E_MSG,'assim_tools_init:','distances are enabled. ')
endif

if(adaptive_localization_threshold > 0) then
   write(msgstring, '(A,I10,A)') 'Using adaptive localization, threshold ', &
                                  adaptive_localization_threshold, ' obs'
   call error_handler(E_MSG,'assim_tools_init:', msgstring)
   if(adaptive_cutoff_floor > 0.0_r8) then
      write(msgstring, '(A,F18.6)') 'Minimum cutoff will not go below ', &
                                     adaptive_cutoff_floor
      call error_handler(E_MSG,'assim_tools_init:', 'Using adaptive localization cutoff floor.', &
                         text2=msgstring)
   endif
endif

if(output_localization_diagnostics) then
   call error_handler(E_MSG,'assim_tools_init:', 'Writing localization diagnostics to file:')
   call error_handler(E_MSG,'assim_tools_init:', trim(localization_diagnostics_file))
endif

if(sampling_error_correction) then
   call error_handler(E_MSG,'assim_tools_init:', 'Using Sampling Error Correction')
endif

if (task_count() > 1) then
    if(distribute_mean) then
       msgstring  = 'Distributing one copy of the ensemble mean across all tasks'
       msgstring2 = 'uses less memory per task but may run slower if doing vertical '
    else
       msgstring  = 'Replicating a copy of the ensemble mean on every task'
       msgstring2 = 'uses more memory per task but may run faster if doing vertical '
    endif
    call error_handler(E_MSG,'assim_tools_init:', msgstring, text2=msgstring2, &
                       text3='coordinate conversion; controlled by namelist item "distribute_mean"')
endif

if (has_vertical_choice()) then
   if (.not. vertical_localization_on()) then
      msgstring = 'Not doing vertical localization, no vertical coordinate conversion required'
      call error_handler(E_MSG,'assim_tools_init:', msgstring)
   else
      msgstring = 'Doing vertical localization, vertical coordinate conversion may be required'
      if (convert_all_state_verticals_first) then
         msgstring2 = 'Converting all state vector verticals to localization coordinate first.'
      else
         msgstring2 = 'Converting all state vector verticals only as needed.'
      endif
      if (convert_all_obs_verticals_first) then
         msgstring3 = 'Converting all observation verticals to localization coordinate first.'
      else
         msgstring3 = 'Converting all observation verticals only as needed.'
      endif
      call error_handler(E_MSG,'assim_tools_init:', msgstring, text2=msgstring2, text3=msgstring3)
   endif
endif

end subroutine log_namelist_selections

!===========================================================
! TEST FUNCTIONS BELOW THIS POINT
!-----------------------------------------------------------
!> test get_state_meta_data
!> Write out the resutls of get_state_meta_data for each task
!> They should be the same as the Trunk version
subroutine test_get_state_meta_data(locations, num_vars)

type(location_type), intent(in) :: locations(:)
integer,             intent(in) :: num_vars

character*20  :: task_str !< string to hold the task number
character*129 :: file_meta !< output file name
character(len=128) :: locinfo
integer :: i

write(task_str, '(i10)') my_task_id()
file_meta = TRIM('test_get_state_meta_data' // TRIM(ADJUSTL(task_str)))

open(15, file=file_meta, status = 'unknown')

do i = 1, num_vars
   call write_location(-1, locations(i), charstring=locinfo)
   write(15,*) trim(locinfo)
enddo

close(15)


end subroutine test_get_state_meta_data

!--------------------------------------------------------
!> dump out the copies array for the state ens handle
subroutine test_state_copies(state_ens_handle, information)

type(ensemble_type), intent(in) :: state_ens_handle
character(len=*),        intent(in) :: information

character*20  :: task_str !< string to hold the task number
character*129 :: file_copies !< output file name
integer :: i

write(task_str, '(i10)') state_ens_handle%my_pe
file_copies = TRIM('statecopies_'  // TRIM(ADJUSTL(information)) // '.' // TRIM(ADJUSTL(task_str)))
open(15, file=file_copies, status ='unknown')

do i = 1, state_ens_handle%num_copies - state_ens_handle%num_extras
   write(15, *) state_ens_handle%copies(i,:)
enddo

close(15)

end subroutine test_state_copies

!--------------------------------------------------------
!> dump out the distances calculated in get_close_obs
subroutine test_close_obs_dist(distances, num_close, ob)

real(r8), intent(in) :: distances(:) !< array of distances calculated in get_close
integer,  intent(in) :: num_close !< number of close obs
integer,  intent(in) :: ob

character*20  :: task_str !< string to hold the task number
character*20  :: ob_str !< string to hold ob number
character*129 :: file_dist !< output file name
integer :: i

write(task_str, '(i10)') my_task_id()
write(ob_str, '(i20)') ob
file_dist = TRIM('distances'   // TRIM(ADJUSTL(task_str)) // '.' // TRIM(ADJUSTL(ob_str)))
open(15, file=file_dist, status ='unknown')

write(15, *) num_close

do i = 1, num_close
   write(15, *) distances(i)
enddo

close(15)

end subroutine test_close_obs_dist

!> @}

!========================================================================
! end module assim_tools_mod
!========================================================================

end module assim_tools_mod

