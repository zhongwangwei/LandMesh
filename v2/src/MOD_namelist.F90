!-- @brief Namelist processing module
!-- @details This module contains subroutines for reading namelist files
!-- and processing mode4 mesh creation
module MOD_namelist
    implicit none
    
    public :: read_nl

contains

!-- @brief Reads namelist settings for grid generation and refinement.
subroutine read_nl(nlfile)
    use consts_coms 
    use refine_vars 
    use MOD_mask_processing, only : Mask_make
    implicit none

    character(*), intent(in) :: nlfile     
    integer :: i, pos, iostat             
    logical :: fexists                    
    character(pathlen) :: path, fprefix, filename, lndname

    namelist /mkgrd/ nl    
    namelist /mkrefine/ rl 
    
    inquire(file = nlfile, exist = fexists) 
    write(io6, *) nlfile
    if (.not. fexists) then
        write(*, *) "The namelist file " // trim(nlfile) // " is missing."
        stop "Stopping model run."
    endif
    open(iunit, status = 'OLD', file = nlfile) 
    
    REWIND(iunit)
    read(iunit, nml = mkgrd)   
    close(iunit)
    write(*, nml = mkgrd)      
    write(io6, *) ""

    expnme               = nl%expnme           
    nxp                  = nl%nxp              
    GXR                  = nl%GXR              
    base_dir             = nl%base_dir         
    source_dir           = nl%source_dir       
    mesh_type            = nl%mesh_type        
    mode_grid            = nl%mode_grid        
    mode_file            = nl%mode_file        
    mode_file_description= nl%mode_file_description
    refine               = nl%refine           
    lcs                  = nl%lcs              
    openmp               = nl%openmp           
    mask_sea_ratio       = nl%mask_sea_ratio   
    mask_restart         = nl%mask_restart     
    mask_domain_type     = nl%mask_domain_type 
    mask_domain_fprefix  = nl%mask_domain_fprefix
    mask_patch_on        = nl%mask_patch_on    
    mask_patch_type      = nl%mask_patch_type  
    mask_patch_fprefix   = nl%mask_patch_fprefix 
    file_dir             = trim(base_dir) // trim(expnme) // '/'

    if (GXR < 0) stop "ERROR! GXR must >= 0" 

    if (mask_restart) then
        if (mask_patch_on) CALL Mask_make('mask_patch', mask_patch_type, mask_patch_fprefix) 
        return
    end if

    CALL execute_command_line('rm -rf '//trim(file_dir)) 
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"contain/")   
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"gridfile/")  
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"patchtype/") 
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"result/")    
    CALL execute_command_line('mkdir -p '//trim(file_dir)//"tmpfile/")   
    CALL execute_command_line('rm *_filelist.txt') 
    
    CALL Mask_make('mask_domain', mask_domain_type, mask_domain_fprefix)
    
    if (mask_patch_on) CALL Mask_make('mask_patch', mask_patch_type, mask_patch_fprefix)

    if (refine) then
        CALL execute_command_line('mkdir -p '//trim(file_dir)//"threshold/") 
        open(iunit, status = 'OLD', file = nlfile) 
        REWIND(iunit)
        read(iunit, nml = mkrefine) 
        close(iunit)
        write(*, nml = mkrefine)    

        weak_concav_eliminate = rl%weak_concav_eliminate
        Istransition          = rl%Istransition         
        max_sa_iter           = rl%max_sa_iter          
        halo                  = rl%halo                 
        max_transition_row    = rl%max_transition_row   
        if (Istransition == .false. .and. mode_grid /= 'tri') STOP "ERROR! not Istransition can only use in the tri" 

        refine_spc            = rl%refine_spc           
        refine_cal            = rl%refine_cal           
        if (refine_spc) max_iter_spc = rl%max_iter_spc  
        if (refine_cal) max_iter_cal = rl%max_iter_cal  

        if ((refine_spc .eqv. .TRUE.) .and. (refine_cal .eqv. .TRUE.)) then
            refine_setting = 'mixed'    
        else if ((refine_spc .eqv. .TRUE.)  .and. (refine_cal .eqv. .FALSE.)) then
            refine_setting = 'specified'
        else if ((refine_spc .eqv. .FALSE.) .and. (refine_cal .eqv. .TRUE.)) then
            refine_setting = 'calculate'
        else
            stop "ERROR! MUst one of TRUE in the refine_spc and refine_cal when refine is TRUE"
        end if
        write(io6, *) "refine_setting = ", refine_setting

        if (refine_setting == 'specified' .or. refine_setting == 'mixed') then
            mask_refine_spc_type       = RL%mask_refine_spc_type   
            mask_refine_spc_fprefix    = RL%mask_refine_spc_fprefix
            CALL Mask_make('mask_refine', mask_refine_spc_type, mask_refine_spc_fprefix) 
            if (mask_refine_ndm(max_iter_spc) == 0) then 
                write(io6, *) "max_iter_spc = ", max_iter_spc
                write(io6, *) "mask_refine_ndm(max_iter_spc) = ", mask_refine_ndm(max_iter_spc)
                stop "ERROR! mask_refine_ndm(max_iter_spc) must larger then one, please modify max_iter_spc"
            end if
        end if

        if (refine_setting == 'calculate' .or. refine_setting == 'mixed') then
            if ((mesh_type == 'landmesh') .or. (mesh_type == 'earthmesh')) then
                refine_num_landtypes      = rl%refine_num_landtypes   
                refine_area_mainland      = rl%refine_area_mainland   
                
                refine_onelayer_Lnd( 1)   = rl%refine_lai_m       
                refine_onelayer_Lnd( 2)   = rl%refine_lai_s       
                refine_onelayer_Lnd( 3)   = rl%refine_slope_m     
                refine_onelayer_Lnd( 4)   = rl%refine_slope_s     
                
                refine_twolayer_Lnd( 1)   = rl%refine_k_s_m       
                refine_twolayer_Lnd( 2)   = rl%refine_k_s_s       
                refine_twolayer_Lnd( 3)   = rl%refine_k_solids_m  
                refine_twolayer_Lnd( 4)   = rl%refine_k_solids_s  
                refine_twolayer_Lnd( 5)   = rl%refine_tkdry_m     
                refine_twolayer_Lnd( 6)   = rl%refine_tkdry_s     
                refine_twolayer_Lnd( 7)   = rl%refine_tksatf_m    
                refine_twolayer_Lnd( 8)   = rl%refine_tksatf_s    
                refine_twolayer_Lnd( 9)   = rl%refine_tksatu_m    
                refine_twolayer_Lnd(10)   = rl%refine_tksatu_s    

                th_num_landtypes          = rl%th_num_landtypes
                th_area_mainland          = rl%th_area_mainland
                th_onelayer_Lnd( 1)       = rl%th_lai_m
                th_onelayer_Lnd( 2)       = rl%th_lai_s
                th_onelayer_Lnd( 3)       = rl%th_slope_m
                th_onelayer_Lnd( 4)       = rl%th_slope_s
                th_twolayer_Lnd( 1, 1:2)  = rl%th_k_s_m        
                th_twolayer_Lnd( 2, 1:2)  = rl%th_k_s_s
                th_twolayer_Lnd( 3, 1:2)  = rl%th_k_solids_m
                th_twolayer_Lnd( 4, 1:2)  = rl%th_k_solids_s
                th_twolayer_Lnd( 5, 1:2)  = rl%th_tkdry_m
                th_twolayer_Lnd( 6, 1:2)  = rl%th_tkdry_s
                th_twolayer_Lnd( 7, 1:2)  = rl%th_tksatf_m
                th_twolayer_Lnd( 8, 1:2)  = rl%th_tksatf_s
                th_twolayer_Lnd( 9, 1:2)  = rl%th_tksatu_m
                th_twolayer_Lnd(10, 1:2)  = rl%th_tksatu_s
            end if

            if ((mesh_type == 'oceanmesh') .or. (mesh_type == 'earthmesh')) then
                refine_sea_ratio          = rl%refine_sea_ratio       
                refine_Rossby_radius      = rl%refine_Rossby_radius   
                
                refine_onelayer_Ocn(1)    = rl%refine_sst_m         
                refine_onelayer_Ocn(2)    = rl%refine_sst_s         
                refine_onelayer_Ocn(3)    = rl%refine_ssh_m         
                refine_onelayer_Ocn(4)    = rl%refine_ssh_s         
                refine_onelayer_Ocn(5)    = rl%refine_eke_m         
                refine_onelayer_Ocn(6)    = rl%refine_eke_s         
                refine_onelayer_Ocn(7)    = rl%refine_sea_slope_m   
                refine_onelayer_Ocn(8)    = rl%refine_sea_slope_s   

                th_sea_ratio              = rl%th_sea_ratio
                th_Rossby_radius          = rl%th_Rossby_radius
                th_onelayer_Ocn(1)        = rl%th_sst_m
                th_onelayer_Ocn(2)        = rl%th_sst_s
                th_onelayer_Ocn(3)        = rl%th_ssh_m
                th_onelayer_Ocn(4)        = rl%th_ssh_s
                th_onelayer_Ocn(5)        = rl%th_eke_m
                th_onelayer_Ocn(6)        = rl%th_eke_s
                th_onelayer_Ocn(7)        = rl%th_sea_slope_m
                th_onelayer_Ocn(8)        = rl%th_sea_slope_s
            end if

            if (mesh_type == 'earthmesh') then
                refine_onelayer_Earth( 1) = rl%refine_typhoon_m   
                refine_onelayer_Earth( 2) = rl%refine_typhoon_s   
                th_onelayer_Earth( 1)     = rl%th_typhoon_m
                th_onelayer_Earth( 2)     = rl%th_typhoon_s
            end if

            if (refine_setting == 'calculate' .or. refine_setting == 'mixed') then
                if (mesh_type == 'landmesh') then
                    if ((refine_num_landtypes .eqv. .false.) .and. &
                        (refine_area_mainland .eqv. .false.) .and. &
                        (all(refine_onelayer_Lnd  .eqv. .false.)).and. &
                        (all(refine_twolayer_Lnd  .eqv. .false.))) then
                        stop "Error! MUst one of TRUE in the refine_num_landtypes or &
                                refine_area_mainland or refine_onelayer_Lnd or refine_twolayer_Lnd &
                                when refine is TRUE and meshtype = landmesh"
                    end if

                else if (mesh_type == 'oceanmesh') then
                    if ((refine_sea_ratio .eqv. .false.) .and. &
                        (refine_Rossby_radius .eqv. .false.) .and. &
                        (all(refine_onelayer_Ocn .eqv. .false.))) then
                        stop "ERROR! MUst one of TRUE in the refine_sea_ratio or refine_onelayer_Ocn when refine is TRUE and meshtype = oceanmesh"
                    end if

                else if (mesh_type == 'earthmesh') then
                    if ((refine_num_landtypes .eqv. .false.) .and. &
                        (refine_area_mainland .eqv. .false.) .and. &
                        (refine_sea_ratio     .eqv. .false.) .and. &
                        (refine_Rossby_radius .eqv. .false.) .and. &
                        (all(refine_onelayer_Lnd  .eqv. .false.)) .and. &
                        (all(refine_twolayer_Lnd  .eqv. .false.)) .and. &
                        (all(refine_onelayer_Ocn  .eqv. .false.)) .and. &
                        (all(refine_onelayer_Earth .eqv. .false.))) then
                        write(io6, *) "refine_num_landtypes = ", refine_num_landtypes
                        write(io6, *) "refine_area_mainland = ", refine_area_mainland
                        write(io6, *) "refine_sea_ratio = ", refine_sea_ratio
                        write(io6, *) "refine_Rossby_radius = ", refine_Rossby_radius
                        stop "Error! MUst one of TRUE in the refine_sea_ratio or refine_Rossby_radius or &
                                refine_num_landtypes or &
                                refine_area_mainland or refine_onelayer_Lnd or refine_twolayer_Lnd or &
                                refine_onelayer_Ocn or refine_onelayer_Earth &
                                when refine is TRUE and meshtype = earthmesh"
                    end if
                end if
            end if
            mask_refine_cal_type       = RL%mask_refine_cal_type   
            mask_refine_cal_fprefix    = RL%mask_refine_cal_fprefix
            CALL Mask_make('mask_refine', mask_refine_cal_type, mask_refine_cal_fprefix) 
        end if

        ! Validation checks for consistency between refinement flags and threshold values
        if ((refine_onelayer_Lnd( 1) .eqv. .true.) .and. (th_onelayer_Lnd( 1) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Lnd( 1) and   th_onelayer_Lnd( 1) "
        if ((refine_onelayer_Lnd( 2) .eqv. .true.) .and. (th_onelayer_Lnd( 2) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Lnd( 2) and   th_onelayer_Lnd( 2) "
        if ((refine_onelayer_Lnd( 3) .eqv. .true.) .and. (th_onelayer_Lnd( 3) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Lnd( 3) and   th_onelayer_Lnd( 3) "
        if ((refine_onelayer_Lnd( 4) .eqv. .true.) .and. (th_onelayer_Lnd( 4) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Lnd( 4) and   th_onelayer_Lnd( 4) "
        
        if ((refine_twolayer_Lnd( 1) .eqv. .true.) .and. any(th_twolayer_Lnd( 1, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 1)  and      th_twolayer_Lnd( 1, 1:2) "
        if ((refine_twolayer_Lnd( 2) .eqv. .true.) .and. any(th_twolayer_Lnd( 2, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 2)  and      th_twolayer_Lnd( 2, 1:2) "
        if ((refine_twolayer_Lnd( 3) .eqv. .true.) .and. any(th_twolayer_Lnd( 3, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 3)  and      th_twolayer_Lnd( 3, 1:2) "
        if ((refine_twolayer_Lnd( 4) .eqv. .true.) .and. any(th_twolayer_Lnd( 4, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 4)  and      th_twolayer_Lnd( 4, 1:2) "
        if ((refine_twolayer_Lnd( 5) .eqv. .true.) .and. any(th_twolayer_Lnd( 5, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 5)  and      th_twolayer_Lnd( 5, 1:2) "
        if ((refine_twolayer_Lnd( 6) .eqv. .true.) .and. any(th_twolayer_Lnd( 6, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 6)  and      th_twolayer_Lnd( 6, 1:2) "
        if ((refine_twolayer_Lnd( 7) .eqv. .true.) .and. any(th_twolayer_Lnd( 7, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 7)  and      th_twolayer_Lnd( 7, 1:2) "
        if ((refine_twolayer_Lnd( 8) .eqv. .true.) .and. any(th_twolayer_Lnd( 8, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 8)  and      th_twolayer_Lnd( 8, 1:2) "
        if ((refine_twolayer_Lnd( 9) .eqv. .true.) .and. any(th_twolayer_Lnd( 9, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd( 9)  and      th_twolayer_Lnd( 9, 1:2) "
        if ((refine_twolayer_Lnd(10) .eqv. .true.) .and. any(th_twolayer_Lnd(10, 1:2) == 999.))  stop "stop for &
        mismatch between refine_twolayer_Lnd(10)  and      th_twolayer_Lnd(10, 1:2) "

        if ((refine_onelayer_Ocn( 1) .eqv. .true.) .and. (th_onelayer_Ocn( 1) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 1) and   th_onelayer_Ocn( 1) "
        if ((refine_onelayer_Ocn( 2) .eqv. .true.) .and. (th_onelayer_Ocn( 2) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 2) and   th_onelayer_Ocn( 2) "
        if ((refine_onelayer_Ocn( 3) .eqv. .true.) .and. (th_onelayer_Ocn( 3) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 3) and   th_onelayer_Ocn( 3) "
        if ((refine_onelayer_Ocn( 4) .eqv. .true.) .and. (th_onelayer_Ocn( 4) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 4) and   th_onelayer_Ocn( 4) "
        if ((refine_onelayer_Ocn( 5) .eqv. .true.) .and. (th_onelayer_Ocn( 5) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 5) and   th_onelayer_Ocn( 5) "
        if ((refine_onelayer_Ocn( 6) .eqv. .true.) .and. (th_onelayer_Ocn( 6) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 6) and   th_onelayer_Ocn( 6) "
        if ((refine_onelayer_Ocn( 7) .eqv. .true.) .and. (th_onelayer_Ocn( 7) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 7) and   th_onelayer_Ocn( 7) "
        if ((refine_onelayer_Ocn( 8) .eqv. .true.) .and. (th_onelayer_Ocn( 8) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Ocn( 8) and   th_onelayer_Ocn( 8) "

        if ((refine_onelayer_Earth( 1) .eqv. .true.) .and. (th_onelayer_Earth( 1) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Earth( 1) and   th_onelayer_Earth( 1) "
        if ((refine_onelayer_Earth( 2) .eqv. .true.) .and. (th_onelayer_Earth( 2) == 999.) ) stop "stop for &
            mismatch between refine_onelayer_Earth( 2) and   th_onelayer_Earth( 2) "

    end if

end subroutine read_nl

end module MOD_namelist 
