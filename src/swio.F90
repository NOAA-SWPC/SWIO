module SWIO

  !-----------------------------------------------------------------------------
  ! IPE Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS           => SetServices,          &
    model_label_DataInitialize => label_DataInitialize, &
    model_label_Advance        => label_Advance,        &
    model_label_Finalize       => label_Finalize,       &
    model_label_CheckImport    => label_CheckImport

  use COMIO
  use mpi

  implicit none

  type type_InternalStateStruct
    integer                 :: logLevel
    integer                 :: fieldCount
    character(ESMF_MAXSTR)  :: file_prefix
    character(ESMF_MAXSTR)  :: file_suffix
    type(ESMF_Config)       :: config
    class(COMIO_T), pointer :: io
  end type

  type type_InternalState
    type(type_InternalStateStruct), pointer :: wrap
  end type

  integer, parameter :: nlon = 90
  integer, parameter :: nlat = 45
  integer, parameter :: nlev = 25

  private

  public :: SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! Provide InitializeP0 to switch from default IPDv00 to IPDv01
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv02p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, &
      specLabel=model_label_DataInitialize, specRoutine=InitializeData, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    ! -- advance method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    ! -- do not check time stamp of imported fields
    call ESMF_MethodRemove(gcomp, model_label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
      specRoutine=NUOPC_NoOp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -- finalize method
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
     specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetServices
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    ! local variables
    integer                    :: stat
    character(len=5)           :: value
    character(len=ESMF_MAXSTR) :: msgString
    type(type_InternalState)   :: is
    
    rc = ESMF_SUCCESS

    ! Allocate memory for the internal state and set it in the component
    allocate(is % wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of internal state memory failed.", &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridCompSetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! Get component's attributes
    ! - Verbosity
    call ESMF_AttributeGet(gcomp, name="Verbosity", value=value, &
      defaultValue="max", convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    ! convert value to logLevel
    is % wrap % logLevel = ESMF_UtilString2Int(value, &
      specialStringList=(/"min","max"/), specialValueList=(/0,255/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    write(msgString,'("IPE: logLevel = ",i0)') is % wrap % logLevel
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! Load component's configuration
    is % wrap % config = ESMF_ConfigCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ConfigLoadFile(is % wrap % config, "swio.conf", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv02p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP0

  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables
    logical :: complete
    integer :: item
!   integer :: columnCount
    character(len=ESMF_MAXSTR) :: importFieldName
    type(type_InternalState) :: is

    ! begin
    rc = ESMF_SUCCESS

    ! get internal state to access component's configuration
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! get import fields from Config
    ! -- read in fields in table
    call ESMF_ConfigFindLabel(is % wrap % config, "import_fields::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    complete = .false.
    item = 0
    do while (.not.complete)
      item = item + 1
      ! -- get input field name
      call ESMF_ConfigNextLine(is % wrap % config, tableEnd=complete, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      if (complete) exit
      call ESMF_ConfigGetAttribute(is % wrap % config, importFieldName, rc=rc )
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      ! -- advertise field
      call NUOPC_Advertise(importState, StandardName=importFieldName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
    end do

  end subroutine InitializeAdvertise
  
  !-----------------------------------------------------------------------------


  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field) :: field
    type(ESMF_Grid)  :: grid

    integer          :: item
    integer          :: stat
    logical          :: isConnected
    character(ESMF_MAXSTR), pointer :: importFieldNames(:)
    type(type_InternalState) :: is

    rc = ESMF_SUCCESS

    ! remove unconnected fields from component's import state
    nullify(importFieldNames)
    call NUOPC_GetStateMemberLists(importState, &
      StandardNameList=importFieldNames, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (.not.associated(importFieldNames)) return

    do item = 1, size(importFieldNames)
      isConnected = NUOPC_IsConnected(importState, &
        fieldName=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (.not.isConnected) then
        call ESMF_StateRemove(importState, (/ importFieldNames(item) /), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end if
    end do

    deallocate(importFieldNames, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Allocation of internal state memory failed.", &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    nullify(importFieldNames)
    
    ! get internal state to access component's configuration
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! retrieve all connected import fields
    call NUOPC_GetStateMemberLists(importState, &
      StandardNameList=importFieldNames, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (associated(importFieldNames)) then
      is % wrap % fieldCount = size(importFieldNames)
    else
      return
    end if

    call IOGridCreateReg(grid, nlon, nlat, nlev, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! realize connected Fields in the importState
    do item = 1, is % wrap % fieldCount
      field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, &
        name=trim(importFieldNames(item)), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_Realize(importState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end do

    deallocate(importFieldNames)

  end subroutine InitializeRealize
  
  !-----------------------------------------------------------------------------

  subroutine InitializeData(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_VM)            :: vm
    type(type_InternalState) :: is
    integer :: comm, iofmt
    character(len=ESMF_MAXSTR) :: svalue
    
    
    rc = ESMF_SUCCESS
    
    ! get internal state to access component's configuration
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (is % wrap % fieldCount > 0) then
      ! get I/O format
      call ESMF_ConfigGetAttribute(is % wrap % config, svalue, &
        label="io_format:", default="none", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      print *,'SWIO: io_format: ' // trim(svalue)
      ! convert text to COMIO I/O selector
      select case (trim(svalue))
        case ("hdf5", "HDF5")
          iofmt = COMIO_FMT_HDF5
          is % wrap % file_suffix = "hd5"
        case ("pnetcdf", "parallel-netcdf", "PNETCDF")
          iofmt = COMIO_FMT_PNETCDF
          is % wrap % file_suffix = "nc"
        case default
          call ESMF_LogSetError(rcToCheck=ESMF_RC_VAL_OUTOFRANGE, &
            msg="Invalid I/O format", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
      end select
      ! get file info
      call ESMF_ConfigGetAttribute(is % wrap % config, is % wrap % file_prefix, &
        label="output_file_prefix:", default="data", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out

      ! get local VM
      ! -- get VM for current gridded component
      call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      ! -- get VM for current gridded component
      call ESMF_VMGet(vm, mpiCommunicator=comm, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
  
      is % wrap % io => COMIO_T(fmt=iofmt, comm=comm, info=MPI_INFO_NULL)
      if (ESMF_LogFoundError(rcToCheck=is % wrap % io % err % rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
    end if

    ! -> set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data
    call NUOPC_CompAttributeSet(gcomp, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
        
  end subroutine InitializeData

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_State)   :: importState
    type(ESMF_Field)   :: field
!---
    logical :: isTimeValid
    integer :: item
    integer :: mm, dd
    integer, dimension(3) :: mstart, mcount
    integer(ESMF_KIND_I4) :: yy, h, m, s
    character(ESMF_MAXSTR) :: fileName
    character(ESMF_MAXSTR), pointer :: importFieldNames(:)
    type(type_InternalState) :: is
    type(ESMF_Time) :: timeStamp
    real(ESMF_KIND_R8), pointer :: fptr(:,:,:)

    real(8), dimension(:,:,:), allocatable :: buf

    ! begin
    rc = ESMF_SUCCESS

    ! get internal state to access component's configuration
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (is % wrap % fieldCount == 0) return

    ! query the component for its importState
    call ESMF_GridCompGet(gcomp, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! retrieve and write imported fields
    nullify(importFieldNames)
    call NUOPC_GetStateMemberLists(importState, &
      StandardNameList=importFieldNames, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (.not.associated(importFieldNames)) return

    do item = 1, size(importFieldNames)
      call ESMF_StateGet(importState, trim(importFieldNames(item)), field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      print *,'SWIO: field: '//trim(importFieldNames(item)) ; flush 6
      if (item == 1) then
        call NUOPC_GetTimestamp(field, isValid=isTimeValid, &
          time=timeStamp, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        print *,'SWIO: field: '//trim(importFieldNames(item)), isTimeValid ; flush 6
        if (.not.isTimeValid) exit
        call ESMF_TimeGet(timeStamp, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        write(fileName, '(a,"_",i4.4,2i2.2,"_",3i2.2,".",a)') &
          trim(is % wrap % file_prefix), yy, mm, dd, h, m, s, trim(is % wrap % file_suffix)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        print *,'SWIO: output file: '//trim(fileName) ; flush 6
        call ESMF_FieldGetBounds(field, computationalLBound=mstart, computationalUBound=mcount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
!       mcount = mcount - mstart + 1
        call is % wrap % io % open(fileName, "c")
        if (ESMF_LogFoundError(rcToCheck=is % wrap % io % err % rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
!       call is % wrap % io % domain((/ nlon, nlat, nlev /), mstart, mcount)
        call is % wrap % io % domain((/ nlon, nlat, nlev /), mstart, mcount - mstart + 1)
        if (ESMF_LogFoundError(rcToCheck=is % wrap % io % err % rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
      nullify(fptr)
      call ESMF_FieldGet(field, farrayPtr=fptr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      allocate(buf(mstart(1):mcount(1),mstart(2):mcount(2),mstart(3):mcount(3)))
      buf = fptr(mstart(1):mcount(1),mstart(2):mcount(2),mstart(3):mcount(3))
      call is % wrap % io % write(trim(importFieldNames(item)), buf)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      deallocate(buf)
    end do
    call is % wrap % io % close()
    if (ESMF_LogFoundError(rcToCheck=is % wrap % io % err % rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    deallocate(importFieldNames)

  end subroutine ModelAdvance
 
  !-----------------------------------------------------------------------------

  subroutine Finalize(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! -- local variables
    integer                  :: stat
    type(type_InternalState) :: is
    

    rc = ESMF_SUCCESS

    ! -- get internal state
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -- destroy config
    call ESMF_ConfigDestroy(is % wrap % config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -- shotdown I/O layer
    call is % wrap % io % shutdown()
    if (ESMF_LogFoundError(rcToCheck=is % wrap % io % err % rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -- finally, deallocate internal state
    deallocate(is % wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Allocation of internal state memory failed.", &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine Finalize

  !------------------------------------------------------

  subroutine IOGridCreateReg(grid, nlon, nlat, nlev, rc)

    type(ESMF_Grid) :: grid
    integer, intent(in) :: nlon, nlat, nlev
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: i, n, localrc
    integer, dimension(1) :: lbnd, ubnd
    real(ESMF_KIND_R8), dimension(:), pointer :: coordPtr
    real(ESMF_KIND_R8), dimension(3) :: coordbase, coordIncr


    if (present(rc)) rc = ESMF_SUCCESS

    coordBase = (/ 0._ESMF_KIND_R8, -90._ESMF_KIND_R8, 90._ESMF_KIND_R8 /)
    coordIncr = (/ 360._ESMF_KIND_R8/(nlon-1), 180._ESMF_KIND_R8/(nlat-1), &
                   800._ESMF_KIND_R8/(nlev-1) /)

    ! -- create 2D grid
    grid = ESMF_GridCreate1PeriDim( &
             maxIndex = (/ nlon, nlat, nlev /), &
             coordSys  = ESMF_COORDSYS_SPH_DEG, &
             coordDep1 = (/ 1 /), &
             coordDep2 = (/ 2 /), &
             coordDep3 = (/ 3 /), &
             indexflag = ESMF_INDEX_GLOBAL, &
             rc = localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    ! -- add coordinates
    call ESMF_GridAddCoord(grid, &
           staggerloc = ESMF_STAGGERLOC_CENTER, rc = localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    do n = 1, size(coordIncr)
      nullify(coordPtr)
      call ESMF_GridGetCoord(grid, coordDim = n, localDE = 0, &
             staggerloc = ESMF_STAGGERLOC_CENTER, &
             computationalLBound = lbnd, computationalUBound = ubnd, &
             farrayPtr = coordPtr, rc = localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return
      do i = lbnd(1), ubnd(1)
        coordPtr(i) = (i-1)*coordIncr(n) + coordBase(n)
      end do
      if (n == 3) then
        do i = lbnd(1), ubnd(1)
!         coordPtr(i) = 1._ESMF_KIND_R8 + coordPtr(i)/6371.2_ESMF_KIND_R8
          coordPtr(i) = 1000._ESMF_KIND_R8 * coordPtr(i)
        end do
      end if
    end do

    nullify(coordPtr)

  end subroutine IOGridCreateReg

end module SWIO
