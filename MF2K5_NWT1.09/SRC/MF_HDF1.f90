   Module HDF_Vars
        use HDF5
        
        !LOGICAL,dimension(3) ::L_HDF  !This is a variable to record if the HDF file format is defined in the MODFLOW name file
        !integer, dimension(3)::File_UN
        
        CHARACTER(LEN=8), PARAMETER :: groupname1  = "ModelDim" ! Group name
        CHARACTER(LEN=8), PARAMETER :: groupname2 = "Stepinfo"  ! Group name
        CHARACTER(LEN=4), PARAMETER :: groupname3 ="Time"  ! Group nam
        CHARACTER(LEN=5), PARAMETER :: groupname4 = "Array" ! Group name
        !CHARACTER(LEN=20) ::groupTname1 = "Array/head/T" ! Group name for timely head
        !CHARACTER(LEN=20) ::groupTname2 = "Array/DD/T" ! Group name for timely drawdown data
        CHARACTER(LEN=30) ::groupTname = "Array/CON/T" ! Group name for timely concentration data
        
        
        character(len=10)::datasetT  ! dataset name for the timely dataset
        CHARACTER(LEN=10)datasetL !dataset name for the array data for each layer
        character(LEN=30):: dataset1 !dataset name for the single data
        
        
        
        INTEGER           :: chunk0   = 64  !Chunk size for 2D DATA
        INTEGER           :: chunk1   = 64 !Chunk size for 2D DATA
        INTEGER           :: chunk2   = 1 !Chunk size for 2D DATA
        
        INTEGER :: hdferr  
        LOGICAL  :: avail
        !INTEGER(HID_T),dimension(3)  :: file_id=(/0,0,0/)  !//id to store the HDF file ID
        INteger::FileCount=0 ! Counting the number of HDF file in the name file
        INTEGER(HID_T)  :: space, dset, dcpl ! Handles
        INTEGER, DIMENSION(1:1) :: cd_values
        INTEGER :: filter_id
        INTEGER :: filter_info_both
        INTEGER(HSIZE_T), DIMENSION(1:2) :: dims , chunk 
        INTEGER(HSIZE_T), DIMENSION(1) :: dims_time = (/2/), chunk_time =(/1/)
        INTEGER(HSIZE_T), DIMENSION(1) :: dims_step = (/3/), chunk_step =(/1/)
        integer(HSIZE_T), dimension(10)::int_dim
        INTEGER(HID_T)  :: space_time, dset_time , dcpl_time,space_step, dset_step , dcpl_step ! Handles
        INTEGER(HSIZE_T), DIMENSION(1:2) :: maxdims
    
        INTEGER(SIZE_T) :: nelmts
        INTEGER :: flags, filter_info
        real, DIMENSION(:, :, :),allocatable :: wdata, & ! Write buffer 
                                               rdata    ! Read buffer
                                          
         integer,dimension(3)::datastep
         real,dimension(2)::datatime
         logical IsHDF
         
         !integer::ICOUNT=0
         
   Type:: HDF_headDD
        !CHARACTER(LEN=8), PARAMETER :: groupname1  = "ModelDim" ! Group name
        !CHARACTER(LEN=8), PARAMETER :: groupname2 = "Stepinfo"  ! Group name
        !CHARACTER(LEN=4), PARAMETER :: groupname3 ="Time"  ! Group nam
        !CHARACTER(LEN=5), PARAMETER :: groupname4 = "Array" ! Group name
        INTEGER(HID_T):: group1_id, group2_id, group3_id, group4_id, groupT_id
        INTEGER(HID_T)::fileid=0, I_TCOUNT(17)=0, File_UN=0
   End type HDF_HeadDD  
   
   type (HDF_headDD)::hdf_file(20)
           
   end module HDF_Vars
    
  
        
        
    
    Subroutine HDF1AR()
    !--Variable definition and allocation
    !--write the model dimension information to the HDF file 
    
    use GLOBAL,only:    NCOL,NROW,NLAY,HNEW,NPER
    use HDF_Vars
    
    dims=(/NCOL,NROW/)
    
    chunk0=64
    CHUNK1=64
    CHUNK2=1
    
    chunk=(/chunk0,chunk1/) !, chunk2/)
    
    END SUBROUTINE 
        
    !////////////////////////////////
    subroutine HDF1WriteG(fileid)
    !/////////////////////////////////
    ! write general information, such as model dimensions
    
    use GLOBAL,only:    NCOL,NROW,NLAY,HNEW,NPER
    use HDF_Vars
    
    integer(HID_T) fileid
    !write the model dimension variable

    int_dim(1)=1
    CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
    call h5screate_simple_f(1,int_dim,space,hdferr)
    dataset1=groupname1//"/"//"NROW"
        
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
    CALL h5dcreate_f(fileid, dataset1,  H5T_NATIVE_INTEGER, space, dset, hdferr, dcpl)
    CALL h5dwrite_f(dset, H5T_NATIVE_INTEGER, NROW, int_dim, hdferr) 
    call h5sclose_f(space,hdferr)
  
    CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
    dataset1=groupname1//"/"//"NCOL"
    CALL h5dcreate_f(fileid, dataset1,  H5T_NATIVE_INTEGER, space, dset, hdferr, dcpl)
    CALL h5dwrite_f(dset, H5T_NATIVE_INTEGER, NCOL, int_dim, hdferr) 
    call h5sclose_f(space,hdferr)
        
    CALL h5screate_f(H5S_SCALAR_F, space, hdferr)
    dataset1=groupname1//"/"//"NLAY"
    CALL h5dcreate_f(fileid, dataset1,  H5T_NATIVE_INTEGER, space, dset, hdferr, dcpl)
    CALL h5dwrite_f(dset, H5T_NATIVE_INTEGER, NLAY, int_dim, hdferr) 
    call h5sclose_f(space,hdferr)
    call h5dclose_f(dset,hdferr)
    call h5pclose_f(dcpl,hdferr)

    END SUBROUTINE
    
    Subroutine HDF1PF(filename,IFLEN, findex)
    !Prepare the HDF files, and create the groups
       
        use HDF_Vars
        
        character(len=iflen)::filename
        !integer(HID_T) fileid
        integer::IFLEN
        integer::findex
        avail=.false.
        call h5open_f(hdferr)
        
        !print *, filename, len(filename)

       
        CALL h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, avail, hdferr)
    
        IF (.NOT.avail) THEN
        WRITE(*,'("gzip filter not available.",/)')
        STOP
        ENDIF
        !CALL h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, hdferr)
        CALL h5zget_filter_info_f(H5Z_FILTER_SHUFFLE_F, filter_info, hdferr)
        filter_info_both=IOR(H5Z_FILTER_ENCODE_ENABLED_F,H5Z_FILTER_DECODE_ENABLED_F)
        IF (filter_info .NE. filter_info_both) THEN
            WRITE(*,'("gzip filter not available for encoding and decoding.",/)')
            STOP
        ENDIF
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf_file(findex)%fileid, hdferr)
        CALL h5gcreate_f(hdf_file(findex)%fileid, groupname1, hdf_file(findex)%group1_id, hdferr)
        CALL h5gcreate_f(hdf_file(findex)%fileid, groupname2, hdf_file(findex)%group2_id, hdferr)
        CALL h5gcreate_f(hdf_file(findex)%fileid, groupname3, hdf_file(findex)%group3_id, hdferr)
        CALL h5gcreate_f(hdf_file(findex)%fileid, groupname4, hdf_file(findex)%group4_id, hdferr)
        
     
    end Subroutine
    
    
    Subroutine HDF1WriteHead_DD(idUnit, KSTP, KPER,varType)
    !Write the head/drawdown data for one time step
    
    use HDF_vars
    USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,BUFF,IBOUND,IOUT
    Use GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBHDSV,CHEDFM,IOFLG
    

    integer::idUnit, KSTP, KPER
    integer(HID_T) fileid
    INTEGER(HID_T):: group1_id, group2_id, group3_id, group4_id, groupT_id
    
    integer::i, j
    character*10:: tempchar
    integer::I_TCOUNT
    Integer:: varType ! a variable to indiate what data are writing; 1 -- head data, 2 -- drawdown data
    
    !searching the fileID for output
    do i=1,FileCount
        if(idUnit.eq.HDF_FILE(i)%File_UN)then
            fileid=hdf_file(i)%fileid
            !fileID is found, start to export the data once a time
            if(varType.eq.1)then
                Print *, "write head data to file unit", idUnit
            elseif(varType.eq.2)then    
                Print *, "write drawdown data to file unit", idUnit
            endif
                
            group1_id=hdf_file(i)%group1_id
            group2_id=hdf_file(i)%group2_id
            group3_id=hdf_file(i)%group3_id
            group4_id=hdf_file(i)%group4_id
            groupT_id=hdf_file(i)%groupT_id
            
            
            
              
  ! Write the data to the dataset.
            maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
            CALL h5screate_simple_f(2, dims, space, hdferr, maxdims)
            CALL h5screate_simple_f(1, dims_time, space_time, hdferr)
            CALL h5screate_simple_f(1, dims_step, space_step, hdferr)

  !
  ! Create the dataset creation property list, add the gzip
  ! compression filter and set the chunk size.
  !
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_time, hdferr)
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_step, hdferr)
  
            if(avail)then
                CALL h5pset_chunk_f(dcpl, 2, chunk, hdferr)
                CALL h5pset_deflate_f(dcpl, 4, hdferr)

                CALL h5pset_chunk_f(dcpl_time, 1, chunk_time, hdferr)                
                CALL h5pset_deflate_f(dcpl_time, 4, hdferr)

                CALL h5pset_chunk_f(dcpl_step, 1, chunk_step, hdferr)  
                CALL h5pset_deflate_f(dcpl_step, 4, hdferr)

            endif    
            
            hdf_file(i)%I_TCOUNT(1)=hdf_file(i)%I_TCOUNT(1)+1
            
            I_TCOUNT=hdf_file(i)%I_TCOUNT(1)
  
  
            write(tempchar,'(I5.1)')I_TCOUNT
            groupTname=trim(groupname4)//"/"//"T"
            groupTname=trim(groupTname)//trim(adjustl(tempchar))
            CALL h5gcreate_f(fileid, trim(groupTname), groupT_id, hdferr)
    
            datastep(1)=0
            datastep(2)=KSTP
            datastep(3)=KPER
    
            datatime(1)=TOTIM
            datatime(2)=Pertim
    
    
            do K=1,NLAY
                write(tempchar,'(I5.1)')K
                dataset1="DS"//adjustl(tempchar)
        
                CALL h5dcreate_f(groupT_id, dataset1,  H5T_IEEE_F32BE, space, dset, hdferr, dcpl)
                CALL h5dwrite_f(dset, H5T_NATIVE_REAL, buff(:,:,k), dims, hdferr) 
            enddo
    
            !call h5gclose_f(groupT_id,hdferr)
    
    
            write(tempchar,'(I5.1)')I_TCOUNT
            dataset1="T"
            dataset1=trim(dataset1)//adjustl(tempchar)
    
    
            CALL h5dcreate_f(group2_id, dataset1,  H5T_NATIVE_INTEGER, space_step, dset_step, hdferr, dcpl_step)
            CALL h5dwrite_f(dset_step, H5T_NATIVE_INTEGER, DATASTEP, dims_step, hdferr)  
    
            CALL h5dcreate_f(group3_id, dataset1,  H5T_NATIVE_REAL, space_time, dset_time, hdferr, dcpl_time)
            CALL h5dwrite_f(dset_time, H5T_NATIVE_REAL, DATATIME, dims_time, hdferr)  
            
            CALL h5pclose_f(dcpl , hdferr)
            CALL h5dclose_f(dset , hdferr)
            CALL h5sclose_f(space, hdferr)
    
            CALL h5pclose_f(dcpl_time , hdferr)
            CALL h5dclose_f(dset_time , hdferr)
            CALL h5sclose_f(space_time, hdferr)
            CALL h5pclose_f(dcpl_step , hdferr)
            CALL h5dclose_f(dset_step , hdferr)
            CALL h5sclose_f(space_Step, hdferr)
            
            RETURN
            
        endif
    enddo
    !can not find file ID, report error
    Write(*,*)" can not find the fileid of the HDF file for output"
    stop
    
    End subroutine
    
    Subroutine HDF1WriteCBB(idUnit, KSTP, KPER,varType)
    !Write the CBB data for one time step
    use HDF_vars
    USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IXSEC,HNEW,BUFF,IBOUND,IOUT
    Use GWFBASMODULE,ONLY:PERTIM,TOTIM,IHEDFM,IHEDUN,LBHDSV,CHEDFM,IOFLG
    

    integer::idUnit, KSTP, KPER
    integer(HID_T) fileid
    INTEGER(HID_T):: group1_id, group2_id, group3_id, group4_id, groupT_id
    
    integer::i, j
    character*10:: tempchar
    integer::I_TCOUNT
    Integer:: varType ! a variable to indiate what data are writing; 1 -- head data, 2 -- drawdown data
    character*16::TEXT
    
    
    !searching the fileID for output
    do i=1,FileCount
        if(idUnit.eq.HDF_FILE(i)%File_UN)then
            fileid=hdf_file(i)%fileid
            !fileID is found, start to export the data once a time
            
            select case(varType)
                case (1)
                    TEXT='         STORAGE'
                    Print *, "write STORAGE cell-by-cell flow data to file unit", idUnit
                Case(2)
                    TEXT='   CONSTANT HEAD'
                    Print *, "write CONSTANT HEAD cell-by-cell flow data to file unit", idUnit
                Case(3)
                    TEXT='FLOW RIGHT FACE '
                    Print *, "write FLOW RIGHT FACE cell-by-cell flow data to file unit", idUnit
                Case(4)
                    TEXT='FLOW FRONT FACE '
                    Print *, "FLOW FRONT FACE flow data to file unit", idUnit
                Case(5)                    
                    TEXT='FLOW LOWER FACE '
                    Print *, "write FLOW LOWER FACE cell-by-cell flow data to file unit", idUnit
                Case(6)
                    TEXT='           WELLS'
                    Print *, "write WELLS cell-by-cell flow data to file unit", idUnit
                Case(7)
                    TEXT='          DRAINS'
                    Print *, "write DRAIN cell-by-cell flow data to file unit", idUnit
                Case(8)
                    TEXT='        RECHARGE'
                    Print *, "write RECHARGE cell-by-cell flow data to file unit", idUnit
                Case(9)
                    TEXT='     ET SEGMENTS'
                    Print *, "write ET cell-by-cell flow data to file unit", idUnit
                Case(10)
                    TEXT='  STREAM LEAKAGE'
                    Print *, "write STREAM LEAKAGE cell-by-cell flow data to file unit", idUnit
                Case(11)
                    TEXT='STREAM FLOW OUT '
                    Print *, "write STREAM FLOW OUT cell-by-cell flow data to file unit", idUnit
                Case(12)
                    TEXT=' HEAD DEP BOUNDS'
                    Print *, "write GHB cell-by-cell flow data to file unit", idUnit
                Case(13)
                    TEXT='INTERBED STORAGE'
                    Print *, "write INTERBED STORAGE cell-by-cell flow data to file unit", idUnit
                Case(14)
                    TEXT='   LAKE  SEEPAGE'                    
                    Print *, "write LAKE  SEEPAGE cell-by-cell flow data to file unit", idUnit
                Case(15)
                    TEXT='     RIPARIAN ET'
                    Print *, "write RIPARIAN ET cell-by-cell flow data to file unit", idUnit
                Case(16)
                    TEXT='   RIVER LEAKAGE'
                    Print *, "write RIVER LEAKAGE cell-by-cell flow data to file unit", idUnit
                Case(17)
                    TEXT='              ET'
                    Print *, "write ET LEAKAGE cell-by-cell flow data to file unit", idUnit
                    
            end select
            
                
            group1_id=hdf_file(i)%group1_id
            group2_id=hdf_file(i)%group2_id
            group3_id=hdf_file(i)%group3_id
            group4_id=hdf_file(i)%group4_id
            groupT_id=hdf_file(i)%groupT_id
            
            
            
              
  ! Write the data to the dataset.
            maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
            CALL h5screate_simple_f(2, dims, space, hdferr, maxdims)
            CALL h5screate_simple_f(1, dims_time, space_time, hdferr)
            CALL h5screate_simple_f(1, dims_step, space_step, hdferr)

  !
  ! Create the dataset creation property list, add the gzip
  ! compression filter and set the chunk size.
  !
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, hdferr)
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_time, hdferr)
            CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_step, hdferr)
  
            if(avail)then
                CALL h5pset_chunk_f(dcpl, 2, chunk, hdferr)
                CALL h5pset_deflate_f(dcpl, 4, hdferr)

                CALL h5pset_chunk_f(dcpl_time, 1, chunk_time, hdferr)                
                CALL h5pset_deflate_f(dcpl_time, 4, hdferr)

                CALL h5pset_chunk_f(dcpl_step, 1, chunk_step, hdferr)  
                CALL h5pset_deflate_f(dcpl_step, 4, hdferr)

            endif    
            
            hdf_file(i)%I_TCOUNT(varType)=hdf_file(i)%I_TCOUNT(varType)+1
            
            I_TCOUNT=hdf_file(i)%I_TCOUNT(varType)
  
  
            write(tempchar,'(I5.1)')I_TCOUNT
            groupTname=trim(groupname4)//"/"//TEXT//"_T"
            groupTname=trim(groupTname)//trim(adjustl(tempchar))
            CALL h5gcreate_f(fileid, trim(groupTname), groupT_id, hdferr)
    
            datastep(1)=0
            datastep(2)=KSTP
            datastep(3)=KPER
    
            datatime(1)=TOTIM
            datatime(2)=Pertim
    
    
            do K=1,NLAY
                write(tempchar,'(I5.1)')K
                dataset1="DS"//adjustl(tempchar)
        
                CALL h5dcreate_f(groupT_id, dataset1,  H5T_IEEE_F32BE, space, dset, hdferr, dcpl)
                CALL h5dwrite_f(dset, H5T_NATIVE_REAL, buff(:,:,k), dims, hdferr) 
            enddo
    
            !call h5gclose_f(groupT_id,hdferr)
    
    
            write(tempchar,'(I5.1)')I_TCOUNT
            dataset1=TEXT//"_T"
            dataset1=trim(dataset1)//adjustl(tempchar)
    
    
            CALL h5dcreate_f(group2_id, dataset1,  H5T_NATIVE_INTEGER, space_step, dset_step, hdferr, dcpl_step)
            CALL h5dwrite_f(dset_step, H5T_NATIVE_INTEGER, DATASTEP, dims_step, hdferr)  
    
            CALL h5dcreate_f(group3_id, dataset1,  H5T_NATIVE_REAL, space_time, dset_time, hdferr, dcpl_time)
            CALL h5dwrite_f(dset_time, H5T_NATIVE_REAL, DATATIME, dims_time, hdferr)  
            
            CALL h5pclose_f(dcpl , hdferr)
            CALL h5dclose_f(dset , hdferr)
            CALL h5sclose_f(space, hdferr)
    
            CALL h5pclose_f(dcpl_time , hdferr)
            CALL h5dclose_f(dset_time , hdferr)
            CALL h5sclose_f(space_time, hdferr)
            CALL h5pclose_f(dcpl_step , hdferr)
            CALL h5dclose_f(dset_step , hdferr)
            CALL h5sclose_f(space_Step, hdferr)
            
            RETURN
            
        endif
    enddo
    !can not find file ID, report error
    Write(*,*)" can not find the fileid of the HDF file for output"
    stop
    
    End subroutine
    
    
    Subroutine HDF1CheckGZip()
    !This subroutine is to check if the GZIp library is available
    
    
    End subroutine
    
    Subroutine HDF1_CheckifHDF(iFileUN)
    !This function is to check if the file format is HDf format.
    
    use HDF_Vars
    Use GWFBASMODULE,ONLY:IHEDFM,IHEDUN
    
    integer::iFileUN
    integer::i
    IsHDF=.false.
    
    
    do i=1,FileCount
        !print *, HDF_File(i)%File_UN, iFileUN
        if(HDF_File(i)%File_UN.EQ.iFileUN)then
            IsHDF=.true.
            exit 
        else
            IsHDF=.false.
            continue
        endif
    enddo
    end subroutine
    