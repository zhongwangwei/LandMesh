MODULE MOD_Threshold_Read
   
    use consts_coms, only: r8, source_dir, nlons_source, nlats_source
    use refine_vars, only: refine_onelayer, refine_twolayer
    use netcdf
    implicit none
    character(LEN = 10), dimension(:), public :: RL_onelayer(2) = (/"lai", "slope_avg"/) ! Add by Rui Zhang
    character(LEN = 10), dimension(:), public :: RL_twolayer(5) = (/"k_s", "k_solids", "tkdry", "tksatf", "tksatu"/) ! Add by Rui Zhang
    type :: var_data2d
        real(r8), allocatable :: var2d(:, :)
    end type
    type(var_data2d), allocatable, public :: input2d(:)

    type :: var_data3d
        real(r8), allocatable :: var3d(:, :, :)
    end type
    type(var_data3d), allocatable, public :: input3d(:)

    !interface data_read
    !    module procedure data_read_onelayer
    !    module procedure data_read_twolayer
    !end interface data_read

    contains

    SUBROUTINE Threshold_Read()
        ! only use for refine_onelayer and refine_twolayer but not refine_num_landtypes and refine_area_mainland
        IMPLICIT NONE
        integer :: num_dataset
        character(len = 256) :: lndname
        character(len = 20) :: varname_select
        integer :: i
        real(r8), allocatable :: input2d_temp(:, :), input3d_temp(:, :, :)

        
        if (all(refine_onelayer .eqv. .false.) .and. &
            all(refine_twolayer .eqv. .false.)) then
            print*, "all false in refine_onelayer and refine_twolayer and return!"
            return
        end if
        print*, "true exist in refine_onelayer or refine_twolayer and go on!"

        ! refine_onelayer
        if (any(refine_onelayer .eqv. .true.)) then
            num_dataset = 0 ! 确定需要的阈值文件个数
            allocate(input2d(size(refine_onelayer)/2)) !%var2d(nlons_source, nlats_source) ! 因为最多只有两个一层数据
            allocate(input2d_temp(nlons_source, nlats_source)); input2d_temp = 0.

            ! 还需要分配inputdata的数组大小, 目前还没有数据裁剪，后面会加上
            do i = 1, size(refine_onelayer)/2, 1 ! 这个7在未来可以更加智能化
                if ((refine_onelayer(2*i-1) .eqv. .true.) .or. &
                    (refine_onelayer(2*i)   .eqv. .true.)) then! 说明该数据集需要读入
                    num_dataset = num_dataset + 1
                    varname_select = RL_onelayer(i)
                    lndname = trim(source_dir) // trim(varname_select) //'.nc' ! slope 应该为 slope_avg.nc
                    print*,lndname
                    allocate(input2d(i)%var2d(nlons_source, nlats_source))
                    CALL data_read_onelayer(lndname, i, varname_select, input2d_temp)
                    input2d(i)%var2d = input2d_temp
                end if
            end do
            deallocate(input2d_temp)
            print*, "onelayer num_dataset = ", num_dataset
        end if

        ! refine_twolayer
        if (any(refine_twolayer .eqv. .true.)) then
            num_dataset = 0 ! 确定需要的阈值文件个数
            allocate(input3d(size(refine_twolayer)/2)) ! 因为最多只有FIVE个一层数据
            allocate(input3d_temp(2, nlons_source, nlats_source)); input3d_temp = 0.

            ! 还需要分配inputdata的数组大小, 目前还没有数据裁剪，后面会加上
            do i = 1, size(refine_twolayer)/2, 1 ! 这个7在未来可以更加智能化
                if ((refine_twolayer(2*i-1) .eqv. .true.) .or. &
                    (refine_twolayer(2*i)   .eqv. .true.)) then! 说明该数据集需要读入
                    num_dataset = num_dataset + 1
                    varname_select = RL_twolayer(i)
                    lndname = trim(source_dir) // trim(varname_select) //'.nc' ! slope 应该为 slope_avg.nc
                    print*,lndname
                    allocate(input3d(i)%var3d(2, nlons_source, nlats_source))
                    CALL data_read_twolayer(lndname, i, varname_select, input3d_temp)
                    input3d(i)%var3d = input3d_temp
                end if
            end do
            deallocate(input3d_temp)
            print*, "twolayer num_dataset = ", num_dataset
        end if

    END SUBROUTINE Threshold_Read

    SUBROUTINE data_read_onelayer(lndname, i, varname_select, input2d_temp)

        IMPLICIT NONE
        character(len = 256), intent(in) :: lndname
        integer, intent(in) :: i
        character(len = 20), intent(in) :: varname_select
        integer :: ncid, varid
        real(r8), dimension(nlons_source, nlats_source), intent(out) :: input2d_temp

        CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid)) ! 文件打开
        CALL CHECK(NF90_INQ_VARID(ncid, trim(varname_select), varid))
        CALL CHECK(NF90_GET_VAR(ncid, varid, input2d_temp))
        CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件
        print*, varname_select, minval(input2d_temp), maxval(input2d_temp)

    END SUBROUTINE data_read_onelayer

    SUBROUTINE data_read_twolayer(lndname, i, varname_select, input3d_temp)

        IMPLICIT NONE
        character(len = 256), intent(in) :: lndname
        integer, intent(in) :: i
        character(len = 20), intent(in) :: varname_select
        character(len = 20)  :: varname_new ! 用于存放需要读取数据的数据集名字
        integer :: k, ncid, varid(2)
        real(r8), dimension(2, nlons_source, nlats_source), intent(out) :: input3d_temp

        CALL CHECK(NF90_OPEN(trim(lndname), nf90_nowrite, ncid)) ! 文件打开
        do k = 1, 2, 1
            if (k == 1) then! ["k_s", "k_solids", "tkdry", "tksatf", "tksatu"] ! 双层信息
                varname_new =  trim(varname_select)//"_l1"
            else
                varname_new =  trim(varname_select)//"_l2"
            end if
            CALL CHECK(NF90_INQ_VARID(ncid,varname_new,varid(k))) 
            CALL CHECK(NF90_GET_VAR(ncid, varid(k), input3d_temp(k, :, :)))
            print*, varname_new, minval(input3d_temp(k, :, :)), maxval(input3d_temp(k, :, :))
        end do
        CALL CHECK(NF90_CLOSE(ncid))! 7. NF90_CLOSE关闭文件

    END SUBROUTINE data_read_twolayer

END MODULE MOD_Threshold_Read
